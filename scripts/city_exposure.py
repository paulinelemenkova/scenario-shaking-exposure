from __future__ import annotations

import csv
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
from pyproj import Geod
from scipy import ndimage

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
LIM = [2.76, 6.2, 11.5, 21.5, 40.1, 74.7]
LAB = ["<IV", "IV", "V", "VI", "VII", "VIII", "IX+"]
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


def urban_centres(pop: FArr, lat: FArr) -> tuple[npt.NDArray[np.int32], FArr]:

    area = (6371.0088 * np.pi / 180 / 120) ** 2 * np.cos(np.radians(lat))[:, None]
    dense = pop / area >= 1500.0
    lab, _ = ndimage.label(dense)
    return lab.astype(np.int32), area


def main() -> None:
    with rasterio.open("data/pop.tif") as src:
        pop = np.nan_to_num(src.read(1).astype(np.float64)).clip(min=0)
        tr = src.transform
    rows_lat = tr.f + tr.e * (np.arange(pop.shape[0]) + 0.5)
    lab, _ = urban_centres(pop, rows_lat)
    with rasterio.open("data/vs30.tif") as vsrc:
        vsa = vsrc.read(1).astype(np.float64)
        vtr, vh, vw = vsrc.transform, vsrc.height, vsrc.width
    with open("ruptures.csv") as f:
        scen = list(csv.DictReader(f))

    out: list[list[Any]] = []
    cells: list[list[Any]] = []
    for city, (clo, cla) in CITIES.items():
        r0, c0 = rasterio.transform.rowcol(tr, clo, cla)
        k = int(lab[r0, c0])
        if k == 0:
            rr, cc = np.nonzero(lab)
            i = int(np.argmin((rr - r0) ** 2 + (cc - c0) ** 2))
            k = int(lab[rr[i], cc[i]])
        sel = ndimage.binary_fill_holes(lab == k)
        rr, cc = np.nonzero(sel)
        p = pop[rr, cc]
        lon = tr.c + tr.a * (cc + 0.5)
        lat = tr.f + tr.e * (rr + 0.5)
        vr, vc = rasterio.transform.rowcol(vtr, lon, lat)
        vs = vsa[np.clip(vr, 0, vh - 1), np.clip(vc, 0, vw - 1)]
        vs = np.where(~np.isfinite(vs) | (vs <= 0), 760.0, vs)
        print(
            f"{city:15s} cells {len(p):5d}  people {p.sum():11.0f}  centre found {lab[r0, c0] > 0}"
        )
        assert p.sum() >= 50_000
        sites = SiteCollection(
            [
                Site(Point(a, b), vs30=v, vs30measured=False, z1pt0=z1, z2pt5=z2)
                for a, b, v, z1, z2 in zip(
                    lon, lat, vs, z1pt0_cy14(vs), z2pt5_cb14(vs), strict=True
                )
            ]
        )
        for row in scen:
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
            w = np.array([wt for _, wt in gsims])[:, None]
            lnm = np.full(len(p), -np.inf)
            sg = np.zeros(len(p))
            lnm[ctx.sids] = (w * mean[:, 0, :]).sum(0)
            sg[ctx.sids] = (w * sig[:, 0, :]).sum(0)
            for lev, z in (("p16", -1.0), ("p50", 0.0), ("p84", 1.0)):
                pga = 100.0 * np.exp(lnm + z * sg)
                cls = np.digitize(pga, LIM)
                e = [float(p[cls == j].sum()) for j in range(len(LAB))]
                out.append(
                    [
                        city,
                        row["scenario"],
                        lev,
                        f"{p.sum():.0f}",
                        *[f"{v:.0f}" for v in e],
                    ]
                )
            pm = 100.0 * np.exp(lnm)
            print(
                f"   {row['scenario']:11s} PGA median {np.median(pm):5.1f} % g  (range {pm.min():.1f}-{pm.max():.1f})"
            )
            cells += [
                [
                    city,
                    row["scenario"],
                    f"{a:.5f}",
                    f"{b:.5f}",
                    f"{q:.0f}",
                    f"{v:.0f}",
                    f"{g:.2f}",
                ]
                for a, b, q, v, g in zip(lon, lat, p, vs, pm, strict=True)
            ]
    with open("data/city_exposure.csv", "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["city", "scenario", "level", "population", *LAB])
        wr.writerows(out)
    with open("data/urban_centres.csv", "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(
            ["city", "scenario", "lon", "lat", "people", "vs30", "pga_median_pctg"]
        )
        wr.writerows(cells)


if __name__ == "__main__":
    main()
