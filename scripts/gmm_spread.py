import csv

import numpy as np
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

geod = Geod(ellps="WGS84")
VS = 300.0
CITIES = {
    "Kabul": (69.181, 34.519),
    "Kandahar": (65.710, 31.620),
    "Herat": (62.200, 34.348),
    "Mazar-i-Sharif": (67.110, 36.709),
    "Jalalabad": (70.450, 34.430),
    "Kunduz": (68.860, 36.728),
}
CRUSTAL = [
    ("BSSA14", BooreEtAl2014(), 0.4),
    ("CZ15", CauzziEtAl2014(), 0.3),
    ("ASB14", AkkarEtAlRjb2014(), 0.3),
]
SLAB = [("PEA20", ParkerEtAl2020SSlab(), 0.5), ("KBCG20", KuehnEtAl2020SSlab(), 0.5)]
KEYS = ("centre_lon", "centre_lat", "centre_depth_km", "strike", "dip", "L_km", "W_km")


def z1(v: float) -> float:
    return float(np.exp(-7.15 / 4 * np.log((v**4 + 571**4) / (1360**4 + 571**4))))


def z25(v: float) -> float:
    return float(np.exp(7.089 - 1.144 * np.log(v)))


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


def run(
    cm: ContextMaker, rup: BaseRupture, sites: SiteCollection
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    ctx = next(iter(cm.get_ctxs([rup], sites)))
    mean, sig, _, _ = cm.get_mean_stds([ctx])
    return mean[:, 0], sig[:, 0], ctx.rrup


offsets = np.geomspace(0.2, 700.0, 46)
rows, crow = [], []
with open("ruptures.csv") as f:
    for r in csv.DictReader(f):
        name = r["scenario"]
        slab = name == "Hindu Kush"
        gs = SLAB if slab else CRUSTAL
        trt = "Subduction IntraSlab" if slab else "Active Shallow Crust"
        p = [float(r[k]) for k in KEYS]
        rup = BaseRupture(
            float(r["Mw"]), float(r["rake"]), trt, Point(p[0], p[1], p[2]), surface(*p)
        )
        rup.tom = PoissonTOM(1.0)
        rup.rup_id = 0
        cm = ContextMaker(
            trt,
            [g for _, g, _ in gs],
            {
                "imtls": {"PGA": [0.1]},
                "maximum_distance": {"default": [(1, 3000), (10, 3000)]},
                "truncation_level": 0,
            },
        )
        pts = [geod.fwd(p[0], p[1], p[3] - 90, d * 1e3)[:2] for d in offsets]
        sites = SiteCollection(
            [
                Site(
                    Point(lo, la),
                    vs30=VS,
                    vs30measured=False,
                    z1pt0=z1(VS),
                    z2pt5=z25(VS),
                )
                for lo, la in pts
            ]
        )
        mu, sg, rr = run(cm, rup, sites)
        w = np.array([x for _, _, x in gs])[:, None]
        lnmed = (w * mu).sum(0)
        sigw = np.sqrt((w * sg**2).sum(0))
        order = np.argsort(rr)
        for i in order:
            row = {
                "scenario": name,
                "rrup_km": f"{rr[i]:.2f}",
                "median_g": f"{np.exp(lnmed[i]):.4g}",
                "p16_g": f"{np.exp(lnmed[i] - sigw[i]):.4g}",
                "p84_g": f"{np.exp(lnmed[i] + sigw[i]):.4g}",
                "sigma_ln": f"{sigw[i]:.3f}",
            }
            for k, (lab, _, _) in enumerate(gs):
                row[lab] = f"{np.exp(mu[k, i]):.4g}"
            rows.append(row)
        cs = SiteCollection(
            [
                Site(
                    Point(lo, la),
                    vs30=VS,
                    vs30measured=False,
                    z1pt0=z1(VS),
                    z2pt5=z25(VS),
                )
                for lo, la in CITIES.values()
            ]
        )
        _, _, rrc = run(cm, rup, cs)
        for cname, d in zip(CITIES, rrc, strict=True):
            crow.append({"scenario": name, "city": cname, "rrup_km": f"{d:.1f}"})

        for d0 in (10.0, 30.0, 100.0, 300.0):
            if rr.min() <= d0 <= rr.max():
                vals = [
                    np.interp(np.log(d0), np.log(rr[order]), mu[k][order])
                    for k in range(len(gs))
                ]
                print(
                    f"{name:11s} R={d0:5.0f} km: max/min GMM = {np.exp(max(vals) - min(vals)):.2f}, "
                    f"sigma_ln = {np.interp(np.log(d0), np.log(rr[order]), sigw[order]):.2f}"
                )
cols = [
    "scenario",
    "rrup_km",
    "median_g",
    "p16_g",
    "p84_g",
    "sigma_ln",
    "BSSA14",
    "CZ15",
    "ASB14",
    "PEA20",
    "KBCG20",
]
with open("data/gmm_curves.csv", "w", newline="") as f:
    wr = csv.DictWriter(f, fieldnames=cols, restval="")
    wr.writeheader()
    wr.writerows(rows)
with open("data/city_rrup.csv", "w", newline="") as f:
    wr = csv.DictWriter(f, fieldnames=["scenario", "city", "rrup_km"])
    wr.writeheader()
    wr.writerows(crow)
