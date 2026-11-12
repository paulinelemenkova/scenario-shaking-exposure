import csv

import numpy as np
import rasterio
from openquake.hazardlib.contexts import ContextMaker
from openquake.hazardlib.geo import Point
from openquake.hazardlib.geo.surface import PlanarSurface
from openquake.hazardlib.gsim.akkar_2014 import AkkarEtAlRjb2014
from openquake.hazardlib.gsim.boore_2014 import BooreEtAl2014
from openquake.hazardlib.gsim.cauzzi_2014 import CauzziEtAl2014
from openquake.hazardlib.gsim.kuehn_2020 import KuehnEtAl2020SSlab
from openquake.hazardlib.gsim.parker_2020 import ParkerEtAl2020SSlab
from openquake.hazardlib.site import Site, SiteCollection
from openquake.hazardlib.source.rupture import BaseRupture
from openquake.hazardlib.tom import PoissonTOM
from pyproj import Geod
from rasterio.transform import from_origin

geod = Geod(ellps="WGS84")
D = 1 / 120
CITIES = {
    "kabul": (69.181, 34.519),
    "kandahar": (65.710, 31.620),
    "herat": (62.200, 34.348),
    "mazar": (67.110, 36.709),
    "jalalabad": (70.450, 34.430),
    "kunduz": (68.860, 36.728),
}
HX, HY = 0.30, 0.24

vsrc = rasterio.open("data/vs30.tif")
vsa = vsrc.read(1).astype(float)


def z1pt0_cy14(v: np.ndarray) -> np.ndarray:

    return np.exp(-7.15 / 4 * np.log((v**4 + 571**4) / (1360**4 + 571**4)))


def z2pt5_cb14(v: np.ndarray) -> np.ndarray:

    return np.exp(7.089 - 1.144 * np.log(v))


def surface(
    clon: float, clat: float, cz: float, strike: float, dip: float, L: float, Wd: float
) -> PlanarSurface:

    hw = Wd / 2 * np.cos(np.radians(dip)) * 1e3
    dz = Wd / 2 * np.sin(np.radians(dip))
    tlo, tla, _ = geod.fwd(clon, clat, strike - 90, hw)
    blo, bla, _ = geod.fwd(clon, clat, strike + 90, hw)
    pts = []
    for lo, la, z in ((tlo, tla, cz - dz), (blo, bla, cz + dz)):
        a = geod.fwd(lo, la, strike + 180, L / 2 * 1e3)
        b = geod.fwd(lo, la, strike, L / 2 * 1e3)
        pts.append((Point(a[0], a[1], max(z, 0.0)), Point(b[0], b[1], max(z, 0.0))))
    (tl, tr), (bl, br) = pts
    return PlanarSurface.from_corner_points(tl, tr, br, bl)


CRUSTAL = [(BooreEtAl2014(), 0.4), (CauzziEtAl2014(), 0.3), (AkkarEtAlRjb2014(), 0.3)]
SLAB = [(ParkerEtAl2020SSlab(), 0.5), (KuehnEtAl2020SSlab(), 0.5)]
KEYS = ("centre_lon", "centre_lat", "centre_depth_km", "strike", "dip", "L_km", "W_km")
with open("ruptures.csv") as f:
    SCEN = list(csv.DictReader(f))

summary = []
for city, (clo, cla) in CITIES.items():
    lons = np.arange(round((clo - HX) / D) * D, clo + HX + D / 2, D)
    lats = np.arange(round((cla + HY) / D) * D, cla - HY - D / 2, -D)
    LON, LAT = np.meshgrid(lons, lats)
    r, c = rasterio.transform.rowcol(vsrc.transform, LON.ravel(), LAT.ravel())
    vs = vsa[np.clip(r, 0, vsrc.height - 1), np.clip(c, 0, vsrc.width - 1)]
    vs = np.where(vs <= 0, 760.0, vs)
    sites = SiteCollection(
        [
            Site(Point(lo, la), vs30=v, vs30measured=False, z1pt0=z1, z2pt5=z2)
            for lo, la, v, z1, z2 in zip(
                LON.ravel(),
                LAT.ravel(),
                vs,
                z1pt0_cy14(vs),
                z2pt5_cb14(vs),
                strict=True,
            )
        ]
    )
    stack = np.zeros((len(SCEN), LON.size))
    for k, row in enumerate(SCEN):
        slab = row["scenario"] == "Hindu Kush"
        gsims = SLAB if slab else CRUSTAL
        trt = "Subduction IntraSlab" if slab else "Active Shallow Crust"
        hypo = Point(
            float(row["centre_lon"]),
            float(row["centre_lat"]),
            float(row["centre_depth_km"]),
        )
        rup = BaseRupture(
            float(row["Mw"]),
            float(row["rake"]),
            trt,
            hypo,
            surface(*(float(row[q]) for q in KEYS)),
        )
        rup.tom = PoissonTOM(1.0)
        rup.rup_id = 0
        cm = ContextMaker(
            trt,
            [g for g, _ in gsims],
            {
                "imtls": {"PGA": [0.1]},
                "maximum_distance": {"default": [(1, 2000), (10, 2000)]},
                "truncation_level": 0,
            },
        )
        ctx = next(iter(cm.get_ctxs([rup], sites)))
        mean, _, _, _ = cm.get_mean_stds([ctx])
        w = np.array([wt for _, wt in gsims])[:, None, None]
        g = np.zeros(LON.size)
        g[ctx.sids] = np.exp((w * mean).sum(0)[0])
        stack[k] = g
    env = stack.max(0).reshape(LON.shape)
    gov = stack.argmax(0).reshape(LON.shape)
    tr = from_origin(lons[0] - D / 2, lats[0] + D / 2, D, D)
    for nm, arr, dt in (
        (f"{city}_env.tif", env, "float32"),
        (f"{city}_gov.tif", gov, "int16"),
    ):
        with rasterio.open(
            f"data/{nm}",
            "w",
            driver="GTiff",
            height=LON.shape[0],
            width=LON.shape[1],
            count=1,
            dtype=dt,
            crs="EPSG:4326",
            transform=tr,
            compress="deflate",
        ) as o:
            o.write(arr.astype(dt), 1)

    i = np.argmin(np.abs(lats - cla))
    j = np.argmin(np.abs(lons - clo))
    k0 = int(gov[i, j])
    summary.append(
        [
            city,
            SCEN[k0]["scenario"],
            SCEN[k0]["Mw"],
            f"{env[i, j]:.3f}",
            f"{np.percentile(env, 5):.3f}",
            f"{np.percentile(env, 95):.3f}",
            f"{env.max():.3f}",
            int(vs.reshape(LON.shape)[i, j]),
        ]
    )
    print(summary[-1])
with open("data/city_pga.csv", "w", newline="") as f:
    wcsv = csv.writer(f)
    wcsv.writerow(
        [
            "city",
            "governing_scenario",
            "Mw",
            "pga_centre_g",
            "pga_p05_g",
            "pga_p95_g",
            "pga_max_g",
            "vs30_centre",
        ]
    )
    wcsv.writerows(summary)
