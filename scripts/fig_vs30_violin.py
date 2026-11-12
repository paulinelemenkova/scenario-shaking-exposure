from __future__ import annotations

import csv
import os
import re
import subprocess
import zipfile
from pathlib import Path
from typing import Any, cast

import geopandas as gpd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import rasterio
from matplotlib.axes import Axes
from matplotlib.font_manager import FontProperties, findfont
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
from rasterio.features import rasterize
from rasterio.transform import Affine, from_origin
from scipy import stats

FArr = npt.NDArray[np.float64]
OUT = Path(os.environ.get("OUT", "."))
WORK = Path(os.environ.get("WORK", "work_violin"))
GDRIVE_ID = os.environ.get("GDRIVE_ID", "1TLTNVIpkC4p7Mff4VkBo8mUVRfA3Yh6v")
VS30_URL = "https://apps.usgs.gov/shakemap_geodata/vs30/global_vs30.grd"
RADIUS_KM = 20.0

CLASSES = ["alluvium", "fan", "bedrock"]
LABELS = {"alluvium": "Alluvium", "fan": "Alluvial fans", "bedrock": "Bedrock"}
SHORT = {"alluvium": "Alluvium", "fan": "Fans", "bedrock": "Bedrock"}

FILL = {"alluvium": "#FAE178", "fan": "#EB963C", "bedrock": "#A59687"}
DARK = {"alluvium": "#8C7414", "fan": "#8A4A0B", "bedrock": "#4A4038"}
MARK = {"alluvium": "o", "fan": "s", "bedrock": "^"}
CITIES = [
    ("Kabul", 69.181, 34.519),
    ("Kandahar", 65.710, 31.620),
    ("Herat", 62.200, 34.348),
    ("Mazar-i-Sharif", 67.110, 36.709),
    ("Jalalabad", 70.450, 34.430),
    ("Kunduz", 68.860, 36.728),
]
NEHRP = [(180.0, "E | D"), (360.0, "D | C"), (760.0, "C | B")]


def set_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Nimbus Sans", "Helvetica"],
            "font.size": 8.5,
            "axes.labelsize": 9.5,
            "axes.linewidth": 0.6,
            "axes.labelpad": 2,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "xtick.minor.width": 0.4,
            "ytick.minor.width": 0.4,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "lines.linewidth": 0.8,
            "pdf.fonttype": 42,
            "mathtext.fontset": "custom",
            "mathtext.rm": "Nimbus Sans",
            "mathtext.it": "Nimbus Sans:italic",
            "mathtext.bf": "Nimbus Sans:bold",
        }
    )
    fp = FontProperties(family=plt.rcParams["font.sans-serif"])
    assert "dejavu" not in findfont(fp).lower(), "DejaVu fallback: Nimbus Sans missing"


def load_vs30() -> tuple[FArr, Affine]:

    WORK.mkdir(exist_ok=True)
    local = Path("global_vs30.grd")
    d = 1.0 / 120.0
    if local.exists():
        with rasterio.open(local) as src:
            win = src.window(59.5 - d / 2, 28.5 - d / 2, 75.7 + d / 2, 39.7 + d / 2)
            arr = src.read(1, window=win).astype(np.float64)
            return arr, src.window_transform(win)
    raw = WORK / "raw.tif"
    env = dict(os.environ, GDAL_SKIP="netCDF")
    cmd = ["gdal_translate", "-q", "-srcwin", "28740", "10140", "1945", "1345"]
    subprocess.run([*cmd, f"/vsicurl/{VS30_URL}", str(raw)], check=True, env=env)
    with rasterio.open(raw) as src:
        arr = src.read(1).astype(np.float64)[::-1]
    return arr, from_origin(59.5 - d / 2, 39.7 + d / 2, d, d)


def classify(sym: str | None) -> str:
    s = (sym or "").strip()
    if s.lower().startswith("w"):
        return "water"
    if re.match(r"^Q\d*(ac|t)$", s):
        return "fan"
    if re.match(r"^Q\d*(a|lap|lac|sm|e|loe)$", s):
        return "alluvium"
    return "bedrock"


def load_classes(shape: tuple[int, int], transform: Affine) -> npt.NDArray[np.int8]:

    zpath = WORK / "geol.zip"
    url = f"https://drive.usercontent.google.com/download?id={GDRIVE_ID}&export=download&confirm=t"
    subprocess.run(
        ["curl", "-s", "-f", "-m", "600", "-L", "-o", str(zpath), url], check=True
    )
    with zipfile.ZipFile(zpath) as z:
        members = [m for m in z.namelist() if "/Geologic age/" in m]
        z.extractall(WORK, members)
    shp = next(WORK.glob("*/Geologic age/naf_disk/shapes/geology/glg*.shp"))
    g = gpd.read_file(shp)
    g["CLS"] = g["UNIT_SYMBO"].map(classify)
    grid = np.zeros(shape, dtype=np.int8)
    for k, c in enumerate(CLASSES, start=1):
        geoms = [(geom, 1) for geom in g.loc[g["CLS"] == c, "geometry"]]
        m = rasterize(
            geoms, out_shape=shape, transform=transform, fill=0, dtype="uint8"
        )
        grid[m == 1] = k
    return grid


def node_coords(shape: tuple[int, int], transform: Affine) -> tuple[FArr, FArr]:
    rows, cols = np.indices(shape)
    lon = transform.c + (cols + 0.5) * transform.a
    lat = transform.f + (rows + 0.5) * transform.e
    return lon.astype(np.float64), lat.astype(np.float64)


def haversine_km(lon: FArr, lat: FArr, lon0: float, lat0: float) -> FArr:
    p1, p2 = np.radians(lat), np.radians(lat0)
    dl = np.radians(lon - lon0)
    a = np.sin((p2 - p1) / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return 2 * 6371.0 * np.arcsin(np.sqrt(a))


def violin(ax: Axes, x: float, v: FArr, c: str, width: float) -> None:

    if v.size < 20:
        ax.plot(
            np.full(v.size, x), v, MARK[c], ms=2.5, mfc=FILL[c], mec=DARK[c], mew=0.3
        )
        return
    parts = ax.violinplot(v, positions=[x], widths=width, showextrema=False)
    bodies = cast("list[Any]", parts["bodies"])
    for body in bodies:
        body.set_facecolor(FILL[c])
        body.set_edgecolor(DARK[c])
        body.set_linewidth(0.6)
        body.set_alpha(1.0)
    q1, med, q3 = np.percentile(v, [25, 50, 75])
    ax.plot([x, x], [q1, q3], color=DARK[c], lw=1.0, solid_capstyle="butt", zorder=3)
    ax.plot(x, med, MARK[c], ms=4, mfc="white", mec=DARK[c], mew=0.8, zorder=4)


def style_axis(ax: Axes) -> None:
    ax.yaxis.set_major_locator(MultipleLocator(100))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.tick_params(which="both", top=False, right=True, direction="in")
    ax.tick_params(which="major", length=4, pad=2)
    ax.tick_params(which="minor", length=2)
    ax.tick_params(axis="x", which="both", length=0)
    ax.grid(axis="y", which="major", lw=0.3, color="0.85")
    ax.set_axisbelow(True)
    for yv, _ in NEHRP:
        ax.axhline(yv, color="0.35", lw=0.5, ls=(0, (4, 2)), zorder=1)


def main() -> None:
    set_style()
    vs30, tr = load_vs30()
    cls = load_classes(vs30.shape, tr)
    lon, lat = node_coords(vs30.shape, tr)
    ok = np.isfinite(vs30)

    samples: dict[tuple[str, str], FArr] = {}
    for k, c in enumerate(CLASSES, start=1):
        samples[("Domain", c)] = vs30[ok & (cls == k)]
    for name, lo, la in CITIES:
        near = haversine_km(lon, lat, lo, la) <= RADIUS_KM
        for k, c in enumerate(CLASSES, start=1):
            samples[(name, c)] = vs30[ok & near & (cls == k)]

    rows = []
    for (area, c), v in samples.items():
        if v.size:
            q1, med, q3 = np.percentile(v, [25, 50, 75])
            rows.append(
                [
                    area,
                    c,
                    v.size,
                    round(med),
                    round(q1),
                    round(q3),
                    round(100 * float(np.mean(v < 360)), 1),
                    round(100 * float(np.mean((v >= 360) & (v < 760))), 1),
                    round(100 * float(np.mean(v >= 760)), 1),
                ]
            )
        else:
            rows.append([area, c, 0, "", "", "", "", "", ""])
    with open(OUT / "vs30_stats.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            ["area", "class", "n", "median", "q1", "q3", "pct_D", "pct_C", "pct_B"]
        )
        w.writerows(rows)
        nat = [samples[("Domain", c)] for c in CLASSES]
        h, p = stats.kruskal(*nat)
        w.writerow([])
        w.writerow(["kruskal_wallis_H", round(float(h), 1), "p", f"{p:.3g}"])
        for i in range(3):
            for j in range(i + 1, 3):
                u, pu = stats.mannwhitneyu(nat[i], nat[j])
                a12 = u / (nat[i].size * nat[j].size)
                w.writerow(
                    [
                        f"{CLASSES[i]}_vs_{CLASSES[j]}",
                        "U",
                        float(u),
                        "p",
                        f"{pu:.3g}",
                        "A12",
                        round(float(a12), 3),
                    ]
                )
    for r in rows:
        print(r)
    print("Kruskal-Wallis H =", round(float(h), 1), "p =", p)

    fig = plt.figure(figsize=(17 / 2.54, 8.6 / 2.54))
    gs = fig.add_gridspec(
        1,
        2,
        width_ratios=[1.45, 3.6],
        wspace=0.035,
        left=0.075,
        right=0.955,
        bottom=0.14,
        top=0.86,
    )
    ax_a = fig.add_subplot(gs[0])
    ax_b = fig.add_subplot(gs[1], sharey=ax_a)
    ylo, yhi = 150.0, 1060.0

    for j, c in enumerate(CLASSES):
        violin(ax_a, float(j), samples[("Domain", c)], c, 0.8)
    ax_a.set_xlim(-0.6, 2.6)
    ax_a.set_xticks(
        range(3),
        [f"{SHORT[c]}\n({samples[('Domain', c)].size:,})" for c in CLASSES],
        fontsize=8,
    )
    ax_a.set_ylabel("$V_{S30}$ (m/s)")
    style_axis(ax_a)

    step, off = 3.6, (-1.0, 0.0, 1.0)
    for i, (name, _, _) in enumerate(CITIES):
        for j, c in enumerate(CLASSES):
            v = samples[(name, c)]
            x = i * step + off[j]
            if v.size:
                violin(ax_b, x, v, c, 0.85)
            ax_b.text(
                x,
                915,
                f"{v.size}",
                ha="center",
                va="bottom",
                fontsize=8,
                color=DARK[c],
                rotation=90,
            )
        if i:
            ax_b.axvline(i * step - step / 2, color="0.75", lw=0.5, zorder=1)
    ax_b.set_xlim(-step / 2 - 0.9, (len(CITIES) - 0.5) * step)
    ax_b.set_xticks(
        [i * step for i in range(len(CITIES))],
        [c[0].replace("Mazar-i-Sharif", "Mazar-\ni-Sharif") for c in CITIES],
    )
    style_axis(ax_b)
    ax_b.tick_params(axis="y", labelleft=False)
    ax_a.set_ylim(ylo, yhi)

    for yv, lab in NEHRP:
        ax_b.text(
            1.008,
            yv,
            lab,
            transform=ax_b.get_yaxis_transform(),
            ha="left",
            va="center",
            fontsize=8,
            color="0.25",
        )
    ax_b.text(
        1.008,
        270,
        "D",
        transform=ax_b.get_yaxis_transform(),
        ha="left",
        va="center",
        fontsize=8,
        fontweight="bold",
        color="0.25",
    )
    ax_b.text(
        1.008,
        560,
        "C",
        transform=ax_b.get_yaxis_transform(),
        ha="left",
        va="center",
        fontsize=8,
        fontweight="bold",
        color="0.25",
    )
    ax_b.text(
        1.008,
        870,
        "B",
        transform=ax_b.get_yaxis_transform(),
        ha="left",
        va="center",
        fontsize=8,
        fontweight="bold",
        color="0.25",
    )

    for ax, tag in ((ax_a, "(a)"), (ax_b, "(b)")):
        ax.text(
            0.015,
            0.985,
            tag,
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=10.5,
            fontweight="bold",
        )

    handles: list[Patch | Line2D] = [
        Patch(facecolor=FILL[c], edgecolor=DARK[c], lw=0.6, label=LABELS[c])
        for c in CLASSES
    ]
    handles.append(
        Line2D(
            [],
            [],
            marker="o",
            ls="-",
            color="0.25",
            mfc="white",
            ms=4,
            lw=1.0,
            label="Median and interquartile range",
        )
    )
    handles.append(
        Line2D(
            [],
            [],
            color="0.35",
            lw=0.5,
            ls=(0, (4, 2)),
            label="NEHRP site-class boundary",
        )
    )
    fig.legend(
        handles=handles,
        loc="upper center",
        ncol=5,
        frameon=False,
        fontsize=9,
        bbox_to_anchor=(0.5, 1.0),
        handlelength=1.6,
        columnspacing=1.2,
    )

    for ext in ("pdf", "png"):
        fig.savefig(
            OUT / f"fig_vs30_violin.{ext}",
            dpi=600,
            bbox_inches="tight",
            pad_inches=0.02,
        )


if __name__ == "__main__":
    main()
