import csv
import json
import os
import urllib.request
import xml.etree.ElementTree as ET

import numpy as np
import rasterio
from openquake.hazardlib.contexts import ContextMaker
from openquake.hazardlib.geo import Point
from openquake.hazardlib.geo.surface import PlanarSurface
from openquake.hazardlib.gsim.akkar_2014 import AkkarEtAlRjb2014
from openquake.hazardlib.gsim.boore_2014 import BooreEtAl2014
from openquake.hazardlib.gsim.cauzzi_2014 import CauzziEtAl2014
from openquake.hazardlib.site import Site, SiteCollection
from openquake.hazardlib.source.rupture import BaseRupture
from openquake.hazardlib.tom import PoissonTOM
from pyproj import Geod
from rasterio.transform import from_origin

geod = Geod(ellps="WGS84")
GSIMS = [(BooreEtAl2014(), 0.4), (CauzziEtAl2014(), 0.3), (AkkarEtAlRjb2014(), 0.3)]
TRT = "Active Shallow Crust"
EVENTS = {
    "us7000hj3u": ("2022 Paktika", None),
    "us6000ldpg": ("2023 Herat", (253.13, 26.19, 74.93)),
    "us7000qsvj": ("2025 Kunar", None),
    "us6000rl31": ("2025 Kholm", (241.3, 70.21, 37.28)),
}
STEP = 3
HALF = 1.5
os.makedirs("out", exist_ok=True)


def get(url: str) -> bytes:
    with urllib.request.urlopen(url, timeout=120) as r:
        return r.read()


def z1(v: np.ndarray) -> np.ndarray:
    return np.exp(-7.15 / 4 * np.log((v**4 + 571**4) / (1360**4 + 571**4)))


def z25(v: np.ndarray) -> np.ndarray:
    return np.exp(7.089 - 1.144 * np.log(v))


def wc_dims(mag: float) -> tuple[float, float]:

    return 10 ** (-2.44 + 0.59 * mag), 10 ** (-1.01 + 0.32 * mag)


def plane(lon, lat, z, strike, dip, L, W):

    if z > 1:
        W = min(W, (z - 0.5) / np.sin(np.radians(dip)) * 2)
    hw = W / 2 * np.cos(np.radians(dip)) * 1e3
    dz = W / 2 * np.sin(np.radians(dip))
    tlo, tla, _ = geod.fwd(lon, lat, strike - 90, hw)
    blo, bla, _ = geod.fwd(lon, lat, strike + 90, hw)
    pts = []
    for lo, la, zz in ((tlo, tla, z - dz), (blo, bla, z + dz)):
        a = geod.fwd(lo, la, strike + 180, L / 2 * 1e3)
        b = geod.fwd(lo, la, strike, L / 2 * 1e3)
        pts.append((Point(a[0], a[1], max(zz, 0.0)), Point(b[0], b[1], max(zz, 0.0))))
    (tl, tr), (bl, br) = pts
    return PlanarSurface.from_corner_points(tl, tr, br, bl), (L, W)


def quad(coords):

    t1, t2, b2, b1 = (np.array(c[:3], dtype=float) for c in coords[:4])
    _, _, ltop = geod.inv(t1[0], t1[1], t2[0], t2[1])
    tm = (t1 + t2) / 2
    bm = (b1 + b2) / 2
    ddir, _, h = geod.inv(tm[0], tm[1], bm[0], bm[1])
    strike = ddir - 90.0
    dz = bm[2] - tm[2]
    dip = float(np.degrees(np.arctan2(dz, h / 1e3)))
    W = float(np.hypot(h / 1e3, dz))
    c = (tm + bm) / 2
    return plane(c[0], c[1], c[2], strike % 360, dip, ltop / 1e3, W)


stats, nodes = [], []
for ev, (label, np_) in EVENTS.items():
    det = json.loads(
        get(
            f"https://earthquake.usgs.gov/fdsnws/event/1/query?eventid={ev}&format=geojson"
        )
    )
    lon, lat, dep = det["geometry"]["coordinates"]
    mag = det["properties"]["mag"]
    sm = det["properties"]["products"]["shakemap"][0]
    c = sm["contents"]
    rup_js = json.loads(get(c["download/rupture.json"]["url"]))
    g = rup_js["features"][0]["geometry"]
    mt = det["properties"]["products"]["moment-tensor"][0]["properties"]
    rake = float(mt["nodal-plane-1-rake"]) if np_ is None else np_[2]
    if g["type"] == "MultiPolygon":
        ring = g["coordinates"][0][0]
        surf, (L, W) = quad(ring)
        rsrc = f"ShakeMap finite rupture, L={L:.1f} km, W={W:.1f} km"
        poly = ring
    else:
        assert np_ is not None
        surf, (L, W) = plane(lon, lat, dep, np_[0], np_[1], *wc_dims(mag))
        rsrc = f"W-phase plane {np_[0]:.0f}/{np_[1]:.0f}/{np_[2]:.0f}, L={L:.1f} km, W={W:.1f} km"
        cl, ca = surf.corner_lons, surf.corner_lats
        poly = [
            [cl[0], ca[0]],
            [cl[1], ca[1]],
            [cl[3], ca[3]],
            [cl[2], ca[2]],
            [cl[0], ca[0]],
        ]
    with open(f"out/{ev}_rupture.gmt", "w") as f:
        f.write(">\n" + "".join(f"{p[0]:.4f} {p[1]:.4f}\n" for p in poly))

    root = ET.fromstring(get(c["download/grid.xml"]["url"]))
    ns = {"s": root.tag.split("}")[0].strip("{")}
    spec_el = root.find("s:grid_specification", ns)
    data_el = root.find("s:grid_data", ns)
    assert spec_el is not None and data_el is not None and data_el.text is not None
    spec = spec_el.attrib
    nlon, nlat = int(spec["nlon"]), int(spec["nlat"])
    dx, dy = float(spec["nominal_lon_spacing"]), float(spec["nominal_lat_spacing"])
    fields = [f.attrib["name"] for f in root.findall("s:grid_field", ns)]
    data = np.array(data_el.text.split(), dtype=float).reshape(-1, len(fields))
    grid = {k: data[:, i].reshape(nlat, nlon) for i, k in enumerate(fields)}
    lo0, la0 = float(spec["lon_min"]), float(spec["lat_max"])

    LON, LAT = grid["LON"][::STEP, ::STEP], grid["LAT"][::STEP, ::STEP]
    keep_c = (LON[0] >= lon - HALF) & (LON[0] <= lon + HALF)
    keep_r = (LAT[:, 0] >= lat - HALF) & (LAT[:, 0] <= lat + HALF)
    sub = {k: v[::STEP, ::STEP][keep_r][:, keep_c] for k, v in grid.items()}
    LON, LAT = sub["LON"], sub["LAT"]
    pga_sm = sub["PGA"] / 100.0
    vs = sub["SVEL"]
    sites = SiteCollection(
        [
            Site(Point(a, b), vs30=v, vs30measured=False, z1pt0=zz1, z2pt5=zz2)
            for a, b, v, zz1, zz2 in zip(
                LON.ravel(),
                LAT.ravel(),
                vs.ravel(),
                z1(vs.ravel()),
                z25(vs.ravel()),
                strict=True,
            )
        ]
    )
    rup = BaseRupture(mag, rake, TRT, Point(lon, lat, dep), surf)
    rup.tom = PoissonTOM(1.0)
    rup.rup_id = 0
    cm = ContextMaker(
        TRT,
        [gg for gg, _ in GSIMS],
        {
            "imtls": {"PGA": [0.1]},
            "maximum_distance": {"default": [(1, 3000), (10, 3000)]},
            "truncation_level": 0,
        },
    )
    ctx = next(iter(cm.get_ctxs([rup], sites)))
    mean, sig, _, _ = cm.get_mean_stds([ctx])
    w = np.array([x for _, x in GSIMS])[:, None]
    model = np.full(LON.size, np.nan)
    model[ctx.sids] = np.exp((w * mean[:, 0]).sum(0))
    rjb = np.full(LON.size, np.nan)
    rjb[ctx.sids] = ctx.rjb
    model = model.reshape(LON.shape)
    rjb = rjb.reshape(LON.shape)
    res = np.log(pga_sm / model)
    tr = from_origin(
        LON[0, 0] - STEP * dx / 2, LAT[0, 0] + STEP * dy / 2, STEP * dx, STEP * dy
    )
    for nm, arr in (("sm", pga_sm), ("model", model), ("res", res)):
        with rasterio.open(
            f"out/{ev}_{nm}.tif",
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
            o.write(arr.astype("float32"), 1)

    st = json.loads(get(c["download/stationlist.json"]["url"]))
    with open(f"out/{ev}_stations.txt", "w") as f:
        for s in st["features"]:
            x, y = s["geometry"]["coordinates"][:2]
            typ = (
                "seismic"
                if s["properties"].get("station_type") == "seismic"
                else "macroseismic"
            )
            f.write(f"{x:.4f} {y:.4f} {typ}\n")
    nseis = sum(
        1 for s in st["features"] if s["properties"].get("station_type") == "seismic"
    )
    m = np.isfinite(res) & (rjb <= 200) & (pga_sm >= 0.005)
    r = res[m]
    near = m & (rjb <= 50)
    stats.append(
        {
            "event": ev,
            "label": label,
            "mag": mag,
            "depth_km": dep,
            "rupture": rsrc,
            "shakemap_version": sm["properties"].get("version"),
            "n_seismic": nseis,
            "n_macroseismic": len(st["features"]) - nseis,
            "n_nodes": int(m.sum()),
            "bias_ln": f"{r.mean():.3f}",
            "sd_ln": f"{r.std(ddof=1):.3f}",
            "bias_ln_rjb50": f"{res[near].mean():.3f}",
            "rmse_ln": f"{np.sqrt((r**2).mean()):.3f}",
            "pga_sm_max_g": f"{np.nanmax(pga_sm):.3f}",
            "pga_model_max_g": f"{np.nanmax(model):.3f}",
            "w_lon": f"{LON.min():.3f}",
            "e_lon": f"{LON.max():.3f}",
            "s_lat": f"{LAT.min():.3f}",
            "n_lat": f"{LAT.max():.3f}",
            "ep_lon": lon,
            "ep_lat": lat,
        }
    )
    print(stats[-1])
    for a, b, rr, x, y, p1, p2 in zip(
        LON[m], LAT[m], rjb[m], res[m], vs[m], pga_sm[m], model[m], strict=True
    ):
        nodes.append(
            [
                ev,
                f"{a:.4f}",
                f"{b:.4f}",
                f"{rr:.2f}",
                f"{x:.4f}",
                f"{y:.0f}",
                f"{p1:.5f}",
                f"{p2:.5f}",
            ]
        )
with open("out/validation_stats.csv", "w", newline="") as f:
    wr = csv.DictWriter(f, fieldnames=list(stats[0]))
    wr.writeheader()
    wr.writerows(stats)
with open("out/validation_nodes.csv", "w", newline="") as f:
    wn = csv.writer(f)
    wn.writerow(
        ["event", "lon", "lat", "rjb_km", "res_ln", "vs30", "pga_sm_g", "pga_model_g"]
    )
    wn.writerows(nodes)
