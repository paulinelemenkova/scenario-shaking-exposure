import csv
import importlib.metadata as md

import numpy as np
import rasterio
from openquake.hazardlib import imt as oqimt
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
W, E, S, N, D = 60.0, 75.2, 29.0, 39.2, 0.05
lons = np.arange(W, E + D / 2, D)
lats = np.arange(N, S - D / 2, -D)
LON, LAT = np.meshgrid(lons, lats)

with rasterio.open("data/vs30.tif") as src:
    rows, cols = rasterio.transform.rowcol(src.transform, LON.ravel(), LAT.ravel())
    vs = src.read(1).astype(float)[
        np.clip(rows, 0, src.height - 1), np.clip(cols, 0, src.width - 1)
    ]
vs = np.where(vs <= 0, 760.0, vs)


def z1pt0_cy14(v: np.ndarray) -> np.ndarray:

    return np.exp(-7.15 / 4 * np.log((v**4 + 571**4) / (1360**4 + 571**4)))


def z2pt5_cb14(v: np.ndarray) -> np.ndarray:

    return np.exp(7.089 - 1.144 * np.log(v))


sites = SiteCollection(
    [
        Site(Point(lo, la), vs30=v, vs30measured=False, z1pt0=z1, z2pt5=z2)
        for lo, la, v, z1, z2 in zip(
            LON.ravel(), LAT.ravel(), vs, z1pt0_cy14(vs), z2pt5_cb14(vs), strict=True
        )
    ]
)

CRUSTAL = [(BooreEtAl2014(), 0.4), (CauzziEtAl2014(), 0.3), (AkkarEtAlRjb2014(), 0.3)]
SLAB = [(ParkerEtAl2020SSlab(), 0.5), (KuehnEtAl2020SSlab(), 0.5)]


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


tag = {
    "Chaman": "chaman",
    "Herat": "herat",
    "Hindu Kush": "hindukush",
    "Kholm-type": "kholm",
}
summary = []
with open("ruptures.csv") as f:
    for r in csv.DictReader(f):
        name = r["scenario"]
        slab = name == "Hindu Kush"
        gsims = SLAB if slab else CRUSTAL
        trt = "Subduction IntraSlab" if slab else "Active Shallow Crust"
        surf = surface(
            *(
                float(r[k])
                for k in (
                    "centre_lon",
                    "centre_lat",
                    "centre_depth_km",
                    "strike",
                    "dip",
                    "L_km",
                    "W_km",
                )
            )
        )
        hypo = Point(
            float(r["centre_lon"]), float(r["centre_lat"]), float(r["centre_depth_km"])
        )
        rup = BaseRupture(float(r["Mw"]), float(r["rake"]), trt, hypo, surf)
        rup.tom = PoissonTOM(1.0)
        rup.rup_id = 0
        cm = ContextMaker(
            trt,
            [g for g, _ in gsims],
            {"imtls": {"PGA": [0.1], "PGV": [1.0]}, "truncation_level": 0},
        )
        ctx = next(iter(cm.get_ctxs([rup], sites)))
        mean, sig, _, _ = cm.get_mean_stds([ctx])
        w = np.array([wt for _, wt in gsims])[:, None, None]
        lnmed = (w * mean).sum(0)
        sid = ctx.sids
        out = {}
        for m, key in enumerate(["pga", "pgv"]):
            g = np.full(LON.size, np.nan)
            g[sid] = np.exp(lnmed[m])
            out[key] = g.reshape(LON.shape)
            with rasterio.open(
                f"data/{tag[name]}_{key}.tif",
                "w",
                driver="GTiff",
                height=LON.shape[0],
                width=LON.shape[1],
                count=1,
                dtype="float32",
                crs="EPSG:4326",
                transform=from_origin(W - D / 2, N + D / 2, D, D),
            ) as o:
                o.write(out[key].astype("float32"), 1)
        summary.append(
            (
                name,
                np.nanmax(out["pga"]),
                np.nanmax(out["pgv"]),
                float(np.mean(sig[:, 0])),
            )
        )
        print(
            f"{name:12s} PGAmax {np.nanmax(out['pga']):.3f} g  PGVmax {np.nanmax(out['pgv']):.1f} cm/s  "
            f"sigma(lnPGA) {np.mean(sig[:, 0]):.2f}"
        )
_ = oqimt
with open("data/gmm_versions.txt", "w") as f:
    f.write(f"openquake.engine {md.version('openquake.engine')}\n")

with (
    open("ruptures.csv") as f,
    open("data/rup.gmt", "w") as fr,
    open("data/top.gmt", "w") as ft,
):
    for k, r in enumerate(csv.DictReader(f)):
        s = surface(
            *(
                float(r[c])
                for c in (
                    "centre_lon",
                    "centre_lat",
                    "centre_depth_km",
                    "strike",
                    "dip",
                    "L_km",
                    "W_km",
                )
            )
        )
        lo, la = s.corner_lons, s.corner_lats
        ring = [
            (lo[0], la[0]),
            (lo[1], la[1]),
            (lo[3], la[3]),
            (lo[2], la[2]),
            (lo[0], la[0]),
        ]
        fr.write(f"> -Z{k}\n" + "".join(f"{x:.4f} {y:.4f}\n" for x, y in ring))
        ft.write(f"> -Z{k}\n{lo[0]:.4f} {la[0]:.4f}\n{lo[1]:.4f} {la[1]:.4f}\n")
