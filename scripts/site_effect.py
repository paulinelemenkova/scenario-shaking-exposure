import csv
import os

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
HX, HY = 0.30, 0.24
CITIES = {
    "kabul": (69.181, 34.519),
    "herat": (62.200, 34.348),
    "mazar": (67.110, 36.709),
    "kunduz": (68.860, 36.728),
    "jalalabad": (70.450, 34.430),
    "kandahar": (65.710, 31.620),
}
RANGE = {1: (180.0, 360.0), 2: (270.0, 560.0), 3: (360.0, 1e4)}
CRUSTAL = [(BooreEtAl2014(), 0.4), (CauzziEtAl2014(), 0.3), (AkkarEtAlRjb2014(), 0.3)]
SLAB = [(ParkerEtAl2020SSlab(), 0.5), (KuehnEtAl2020SSlab(), 0.5)]
KEYS = ("centre_lon", "centre_lat", "centre_depth_km", "strike", "dip", "L_km", "W_km")
os.makedirs("out", exist_ok=True)


def z1(v: np.ndarray) -> np.ndarray:
    return np.exp(-7.15 / 4 * np.log((v**4 + 571**4) / (1360**4 + 571**4)))


def z25(v: np.ndarray) -> np.ndarray:
    return np.exp(7.089 - 1.144 * np.log(v))


def surface(clon, clat, cz, strike, dip, L, Wd) -> PlanarSurface:
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


def median_pga(
    row: dict, lon: np.ndarray, lat: np.ndarray, vs: np.ndarray
) -> np.ndarray:
    slab = row["scenario"] == "Hindu Kush"
    gs = SLAB if slab else CRUSTAL
    trt = "Subduction IntraSlab" if slab else "Active Shallow Crust"
    p = [float(row[k]) for k in KEYS]
    rup = BaseRupture(
        float(row["Mw"]), float(row["rake"]), trt, Point(p[0], p[1], p[2]), surface(*p)
    )
    rup.tom = PoissonTOM(1.0)
    rup.rup_id = 0
    cm = ContextMaker(
        trt,
        [g for g, _ in gs],
        {
            "imtls": {"PGA": [0.1]},
            "maximum_distance": {"default": [(1, 3000), (10, 3000)]},
            "truncation_level": 0,
        },
    )
    sites = SiteCollection(
        [
            Site(Point(a, b), vs30=v, vs30measured=False, z1pt0=c, z2pt5=d)
            for a, b, v, c, d in zip(lon, lat, vs, z1(vs), z25(vs), strict=True)
        ]
    )
    ctx = next(iter(cm.get_ctxs([rup], sites)))
    mean, _, _, _ = cm.get_mean_stds([ctx])
    w = np.array([x for _, x in gs])[:, None]
    out = np.zeros(lon.size)
    out[ctx.sids] = np.exp((w * mean[:, 0]).sum(0))
    return out


CLS_DIR = os.environ.get("CLS_DIR", "cls")

with open("ruptures.csv") as f:
    SCEN = {r["scenario"]: r for r in csv.DictReader(f)}
with open(os.environ.get("CITY_PGA", "city_pga.csv")) as f:
    GOV = {r["city"]: r["governing_scenario"] for r in csv.DictReader(f)}

vsrc = rasterio.open(os.environ.get("VS30", "vs30.tif"))
vsa = vsrc.read(1).astype(float)
rows = []
for city, (clo, cla) in CITIES.items():
    lons = np.arange(round((clo - HX) / D) * D, clo + HX + D / 2, D)
    lats = np.arange(round((cla + HY) / D) * D, cla - HY - D / 2, -D)
    LON, LAT = np.meshgrid(lons, lats)
    tr = from_origin(lons[0] - D / 2, lats[0] + D / 2, D, D)
    r, c = rasterio.transform.rowcol(vsrc.transform, LON.ravel(), LAT.ravel())
    vs = vsa[np.clip(r, 0, vsrc.height - 1), np.clip(c, 0, vsrc.width - 1)]
    vs = np.where(vs <= 0, 760.0, vs)
    with rasterio.open(f"{CLS_DIR}/{city}_cls.tif") as cs:
        cg = cs.read(1).ravel()
    va = vs.copy()
    for k, (lo, hi) in RANGE.items():
        m = cg == k
        va[m] = np.clip(vs[m], lo, hi)
    row = SCEN[GOV[city]]
    p0 = median_pga(row, LON.ravel(), LAT.ravel(), vs)
    p1 = median_pga(row, LON.ravel(), LAT.ravel(), va)
    ratio = (p1 / p0).reshape(LON.shape)
    with rasterio.open(
        f"out/{city}_ratio.tif",
        "w",
        driver="GTiff",
        height=LON.shape[0],
        width=LON.shape[1],
        count=1,
        dtype="float32",
        crs="EPSG:4326",
        transform=tr,
        compress="deflate",
    ) as o:
        o.write(ratio.astype("float32"), 1)
    i = np.argmin(np.abs(lats - cla))
    j = np.argmin(np.abs(lons - clo))
    changed = va != vs
    rr = ratio.ravel()
    rows.append(
        {
            "city": city,
            "scenario": GOV[city],
            "nodes": LON.size,
            "pct_changed": f"{100 * changed.mean():.1f}",
            "pct_alluvium_capped": f"{100 * ((cg == 1) & (vs > 360)).mean():.1f}",
            "pct_bedrock_raised": f"{100 * ((cg == 3) & (vs < 360)).mean():.1f}",
            "vs30_centre": f"{vs.reshape(LON.shape)[i, j]:.0f}",
            "vs30_adj_centre": f"{va.reshape(LON.shape)[i, j]:.0f}",
            "ratio_centre": f"{ratio[i, j]:.3f}",
            "ratio_min": f"{rr.min():.3f}",
            "ratio_max": f"{rr.max():.3f}",
            "ratio_mean_changed": f"{rr[changed].mean():.3f}"
            if changed.any()
            else "1.000",
            "pga_centre_global_g": f"{p0.reshape(LON.shape)[i, j]:.3f}",
            "pga_centre_adj_g": f"{p1.reshape(LON.shape)[i, j]:.3f}",
        }
    )
    print(rows[-1])
with open("out/site_effect.csv", "w", newline="") as f:
    wr = csv.DictWriter(f, fieldnames=list(rows[0]))
    wr.writeheader()
    wr.writerows(rows)
