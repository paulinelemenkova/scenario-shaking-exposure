from __future__ import annotations

import csv
import json
from itertools import pairwise
from typing import Any

import numpy as np
import numpy.typing as npt
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
from pyproj import Geod, Transformer
from rasterio.transform import rowcol
from shapely.geometry import Point as ShPoint
from shapely.geometry import shape
from shapely.ops import transform

FArr = npt.NDArray[np.float64]
geod = Geod(ellps="WGS84")
CITIES = {
    "Kabul": (69.181, 34.519),
    "Kandahar": (65.710, 31.620),
    "Herat": (62.200, 34.348),
    "Mazar-i-Sharif": (67.110, 36.709),
    "Jalalabad": (70.450, 34.430),
    "Kunduz": (68.860, 36.728),
}
CRUSTAL = [(BooreEtAl2014(), 0.4), (CauzziEtAl2014(), 0.3), (AkkarEtAlRjb2014(), 0.3)]
SLAB = [(ParkerEtAl2020SSlab(), 0.5), (KuehnEtAl2020SSlab(), 0.5)]
KEYS = ("centre_lon", "centre_lat", "centre_depth_km", "strike", "dip", "L_km", "W_km")


def z1pt0_cy14(v: FArr) -> FArr:

    return np.exp(-7.15 / 4 * np.log((v**4 + 571**4) / (1360**4 + 571**4)))


def z2pt5_cb14(v: FArr) -> FArr:

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


def contour_value(geojson: str, lon: float, lat: float) -> tuple[float, float, float]:

    aeqd = Transformer.from_crs(
        "EPSG:4326",
        f"+proj=aeqd +lat_0={lat} +lon_0={lon} +ellps=WGS84",
        always_xy=True,
    )
    with open(geojson) as f:
        feats = json.load(f)["features"]
    origin = ShPoint(0.0, 0.0)
    dist: dict[float, float] = {}
    for ft in feats:
        lev = float(ft["properties"]["ACC_VAL"])
        geom = transform(aeqd.transform, shape(ft["geometry"]))
        dist[lev] = min(geom.distance(origin), dist.get(lev, np.inf))
    levels = sorted(v for v in dist if v > 0)
    best: tuple[float, float, float] | None = None
    score = np.inf
    for lo, hi in pairwise(levels):
        s = dist[lo] + dist[hi]
        if s < score:
            score, best = s, (lo, hi, lo + (hi - lo) * dist[lo] / (dist[lo] + dist[hi]))
    assert best is not None
    lo, hi, val = best
    return val, lo, hi


def scenario_at(lon: float, lat: float, vs30: float) -> list[tuple[str, float, float]]:

    v = np.array([vs30])
    sites = SiteCollection(
        [
            Site(
                Point(lon, lat),
                vs30=vs30,
                vs30measured=False,
                z1pt0=float(z1pt0_cy14(v)[0]),
                z2pt5=float(z2pt5_cb14(v)[0]),
            )
        ]
    )
    out = []
    with open("ruptures.csv") as f:
        for row in csv.DictReader(f):
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
                    "truncation_level": 0,
                    "maximum_distance": {"default": [(1, 2000), (10, 2000)]},
                },
            )
            ctx = next(iter(cm.get_ctxs([rup], sites)))
            mean, sig, _, _ = cm.get_mean_stds([ctx])
            w = np.array([wt for _, wt in gsims])
            out.append(
                (
                    row["scenario"],
                    float(np.exp((w * mean[:, 0, 0]).sum())),
                    float((w * sig[:, 0, 0]).sum()),
                )
            )
    return out


def main() -> None:
    gem = rasterio.open("gem.tif")
    garr = gem.read(1).astype(np.float64)
    with rasterio.open("data/vs30.tif") as vsrc:
        vsa = vsrc.read(1).astype(np.float64)
        vtr = vsrc.transform
    rows_p: list[list[Any]] = []
    rows_s: list[list[Any]] = []
    for city, (lon, lat) in CITIES.items():
        fr, fc = ~gem.transform * (lon, lat)
        r0, c0 = int(np.floor(fr - 0.5)), int(np.floor(fc - 0.5))
        tr_, tc_ = fr - 0.5 - r0, fc - 0.5 - c0
        q = garr[r0 : r0 + 2, c0 : c0 + 2]
        g475 = (
            q[0, 0] * (1 - tr_) * (1 - tc_)
            + q[0, 1] * (1 - tr_) * tc_
            + q[1, 0] * tr_ * (1 - tc_)
            + q[1, 1] * tr_ * tc_
        )
        u475, l475, h475 = contour_value("afpga050_arc.json", lon, lat)
        u2475, l2475, h2475 = contour_value("afpga250_arc.json", lon, lat)
        rows_p.append(
            [
                city,
                lon,
                lat,
                f"{100 * g475:.2f}",
                f"{u475:.2f}",
                l475,
                h475,
                f"{u2475:.2f}",
                l2475,
                h2475,
            ]
        )
        print(
            f"{city:15s} GEM475 {100 * g475:5.1f}  USGS475 {u475:5.1f} [{l475:g}-{h475:g}]  USGS2475 {u2475:5.1f} [{l2475:g}-{h2475:g}] % g"
        )
        vr, vc = rowcol(vtr, lon, lat)
        vs_site = (
            float(vsa[vr, vc])
            if np.isfinite(vsa[vr, vc]) and vsa[vr, vc] > 0
            else 760.0
        )
        for vs, tag in ((760.0, "rock"), (vs_site, "site")):
            for scen, med, sg in scenario_at(lon, lat, vs):
                rows_s.append(
                    [city, scen, tag, f"{vs:.0f}", f"{100 * med:.3f}", f"{sg:.3f}"]
                )
                print(
                    f"   {scen:11s} {tag:4s} Vs30 {vs:4.0f}  median {100 * med:6.2f} % g  sigma {sg:.2f}"
                )
    with open("data/psha_cities.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "city",
                "lon",
                "lat",
                "gem475_pctg",
                "usgs475_pctg",
                "usgs475_lo",
                "usgs475_hi",
                "usgs2475_pctg",
                "usgs2475_lo",
                "usgs2475_hi",
            ]
        )
        w.writerows(rows_p)
    with open("data/scenario_cities.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["city", "scenario", "site", "vs30", "pga_median_pctg", "sigma_ln"])
        w.writerows(rows_s)


if __name__ == "__main__":
    main()
