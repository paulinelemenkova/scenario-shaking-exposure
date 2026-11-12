from __future__ import annotations

import csv
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import itertools

import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
from matplotlib.axes import Axes
from matplotlib.font_manager import FontProperties, findfont
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import LogLocator, NullFormatter

FArr = npt.NDArray[np.float64]
IN = Path(os.environ.get("IN", "."))
OUT = Path(os.environ.get("OUT", "."))

SCEN = [
    ("Chaman", r"a) Chaman $M_\mathrm{w}$ 7.0, strike-slip"),
    ("Herat", r"b) Herat $M_\mathrm{w}$ 7.5, strike-slip"),
    ("Hindu Kush", r"c) Hindu Kush $M_\mathrm{w}$ 7.5, intraslab, 231 km"),
    ("Kholm-type", r"d) Kholm-type $M_\mathrm{w}$ 6.5, oblique reverse"),
]

GMM = {
    "BSSA14": ("Boore et al. (2014)", "#0072B2", (0, (5, 1.5)), 400.0),
    "CZ15": ("Cauzzi et al. (2015)", "#D55E00", (0, (1.2, 1.2)), 150.0),
    "ASB14": ("Akkar et al. (2014)", "#009E73", (0, (6, 1.5, 1.2, 1.5)), 200.0),
    "PEA20": ("Parker et al. (2022)", "#CC79A7", (0, (5, 1.5)), 1000.0),
    "KBCG20": ("Kuehn et al. (2023)", "#E69F00", (0, (1.2, 1.2)), 1000.0),
}
CRUSTAL = ["BSSA14", "CZ15", "ASB14"]
SLAB = ["PEA20", "KBCG20"]
WEIGHT = {"BSSA14": 0.4, "CZ15": 0.3, "ASB14": 0.3, "PEA20": 0.5, "KBCG20": 0.5}

MMI_B = [0.0005, 0.003, 0.0276, 0.062, 0.115, 0.215, 0.401, 0.747, 1.39]
MMI_L = ["I", "II-III", "IV", "V", "VI", "VII", "VIII", "IX"]
CITY_SHORT = {"Mazar-i-Sharif": "Mazar"}
XLIM = (0.7, 900.0)
YLIM = (0.0015, 2.5)


def set_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Nimbus Sans", "Helvetica"],
            "font.size": 8,
            "axes.labelsize": 9,
            "axes.linewidth": 0.6,
            "axes.labelpad": 2,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "xtick.minor.width": 0.4,
            "ytick.minor.width": 0.4,
            "xtick.major.size": 3.5,
            "ytick.major.size": 3.5,
            "xtick.minor.size": 2.0,
            "ytick.minor.size": 2.0,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "lines.linewidth": 0.8,
            "legend.fontsize": 7,
            "pdf.fonttype": 42,
            "mathtext.fontset": "custom",
            "mathtext.rm": "Nimbus Sans",
            "mathtext.it": "Nimbus Sans:italic",
            "mathtext.bf": "Nimbus Sans:bold",
        }
    )
    fp = FontProperties(family=plt.rcParams["font.sans-serif"])
    assert "dejavu" not in findfont(fp).lower(), "DejaVu fallback: Nimbus Sans missing"


def load_curves() -> dict[str, dict[str, FArr]]:
    out: dict[str, dict[str, list[float]]] = {}
    with open(IN / "gmm_curves.csv") as f:
        for r in csv.DictReader(f):
            d = out.setdefault(r["scenario"], {})
            for k, v in r.items():
                if k != "scenario" and v:
                    d.setdefault(k, []).append(float(v))
    res: dict[str, dict[str, FArr]] = {}
    for s, d in out.items():
        a = {k: np.asarray(v) for k, v in d.items()}

        keep = np.r_[True, np.diff(a["rrup_km"]) > 0.05]
        res[s] = {k: v[keep] for k, v in a.items()}
    return res


def load_cities() -> dict[str, list[tuple[str, float]]]:
    out: dict[str, list[tuple[str, float]]] = {}
    with open(IN / "city_rrup.csv") as f:
        for r in csv.DictReader(f):
            out.setdefault(r["scenario"], []).append((r["city"], float(r["rrup_km"])))
    return out


def mmi_axis(ax: Axes, show_label: bool) -> None:

    for b in MMI_B[2:-1]:
        ax.axhline(b, color="0.55", lw=0.4, ls=(0, (2, 2)), zorder=1)
    tw = ax.twinx()
    tw.set_yscale("log")
    tw.set_ylim(*YLIM)
    mids = [np.sqrt(a * b) for a, b in itertools.pairwise(MMI_B)]
    keep = [
        (m, lab) for m, lab in zip(mids, MMI_L, strict=True) if YLIM[0] < m < YLIM[1]
    ]
    tw.set_yticks([m for m, _ in keep])
    tw.set_yticklabels([lab for _, lab in keep], fontsize=7, color="0.3")
    tw.yaxis.set_minor_locator(LogLocator(numticks=1))
    tw.yaxis.set_minor_formatter(NullFormatter())
    tw.tick_params(axis="y", length=0, pad=2)
    if show_label:
        tw.set_ylabel("Instrumental intensity, MMI", color="0.3", labelpad=3)


def city_marks(
    ax: Axes, r: FArr, med: FArr, p84: FArr, cities: list[tuple[str, float]]
) -> None:

    vis = sorted(
        (d, c) for c, d in cities if XLIM[0] < d < XLIM[1] and r.min() <= d <= r.max()
    )
    pts = [
        (d, float(np.exp(np.interp(np.log(d), np.log(r), np.log(med)))), c)
        for d, c in vis
    ]
    groups: list[list[tuple[float, float, str]]] = []
    for p in pts:
        if groups and np.log10(p[0]) - np.log10(groups[-1][-1][0]) < 0.25:
            groups[-1].append(p)
        else:
            groups.append([p])
    lead = {"arrowstyle": "-", "lw": 0.5, "color": "0.25", "shrinkA": 1, "shrinkB": 2.5}
    for g in groups:
        for d, y, _ in g:
            ax.plot(d, y, "o", ms=3.6, mfc="white", mec="black", mew=0.7, zorder=6)
        if len(g) == 1:
            d, y, c = g[0]
            near_edge = np.log10(d) > np.log10(XLIM[1]) - 0.35
            ax.annotate(
                CITY_SHORT.get(c, c),
                xy=(d, y),
                xytext=(-4 if near_edge else 0, 14),
                textcoords="offset points",
                ha="right" if near_edge else "center",
                va="bottom",
                fontsize=6.8,
                arrowprops=lead,
                zorder=7,
            )
            continue
        xcol = g[0][0] / 10**0.07
        span = (r >= xcol / 3) & (r <= xcol)
        top = float(p84[span].max()) if span.any() else 0.0
        ylab = max(top * 1.35, max(y for _, y, _ in g) * 2.0)
        for d, y, c in g:
            ax.annotate(
                CITY_SHORT.get(c, c),
                xy=(d, y),
                xytext=(xcol, ylab),
                textcoords="data",
                ha="right",
                va="center",
                fontsize=6.8,
                arrowprops={
                    **lead,
                    "relpos": (1.0, 0.5),
                    "connectionstyle": "angle,angleA=0,angleB=90,rad=0",
                },
                zorder=7,
            )
            ylab *= 1.55


def panel(
    ax: Axes,
    name: str,
    title: str,
    cv: dict[str, FArr],
    cities: list[tuple[str, float]],
) -> None:
    r = cv["rrup_km"]
    keys = SLAB if name == "Hindu Kush" else CRUSTAL
    ax.fill_between(r, cv["p16_g"], cv["p84_g"], color="0.82", lw=0, zorder=2)
    for k in keys:
        _lab, col, dash, rmax = GMM[k]
        inside = r <= rmax
        ax.plot(r[inside], cv[k][inside], color=col, ls=dash, lw=0.8, zorder=3)
        if (~inside).any():
            j = max(int(inside.sum()) - 1, 0)
            ax.plot(r[j:], cv[k][j:], color=col, ls=dash, lw=0.8, alpha=0.35, zorder=3)
    ax.plot(r, cv["median_g"], color="black", lw=1.0, zorder=4)
    city_marks(ax, r, cv["median_g"], cv["p84_g"], cities)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(*XLIM)
    ax.set_ylim(*YLIM)
    ax.grid(True, which="major", color="0.88", lw=0.3)
    ax.set_axisbelow(True)
    ax.xaxis.set_minor_locator(
        LogLocator(subs=[2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0], numticks=12)
    )
    ax.yaxis.set_minor_locator(
        LogLocator(subs=[2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0], numticks=12)
    )
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.tick_params(which="both", top=True)
    ax.text(
        0.975,
        0.965,
        title,
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=8,
        fontweight="bold",
        bbox={"fc": "white", "ec": "none", "pad": 1.2},
        zorder=8,
    )


def legend(ax: Axes, keys: list[str]) -> None:
    h: list[Line2D | Patch] = [Line2D([], [], color="black", lw=1.0)]
    lab = ["Weighted median"]
    h.append(Patch(fc="0.82", ec="none"))
    lab.append(r"16th$-$84th percentile")
    for k in keys:
        h.append(Line2D([], [], color=GMM[k][1], ls=GMM[k][2], lw=0.8))
        lab.append(f"{GMM[k][0]}, w = {WEIGHT[k]:.1f}")
    h.append(
        Line2D([], [], ls="none", marker="o", ms=3.6, mfc="white", mec="black", mew=0.7)
    )
    lab.append("City (rupture distance)")
    ax.legend(
        h,
        lab,
        loc="lower left",
        frameon=True,
        framealpha=0.95,
        edgecolor="0.6",
        fancybox=False,
        handlelength=2.6,
        borderpad=0.4,
        labelspacing=0.25,
    ).set_zorder(9)


def main() -> None:
    set_style()
    cv = load_curves()
    cities = load_cities()
    fig, axs = plt.subplots(
        2, 2, figsize=(17.4 / 2.54, 13.0 / 2.54), sharex=True, sharey=True
    )
    for ax, (name, title) in zip(axs.flat, SCEN, strict=True):
        panel(ax, name, title, cv[name], cities[name])
    for i, ax in enumerate(axs.flat):
        mmi_axis(ax, show_label=(i % 2 == 1))
    legend(axs[0, 0], CRUSTAL)
    legend(axs[1, 0], SLAB)
    for ax in axs[1, :]:
        ax.set_xlabel("Rupture distance, $R_\\mathrm{rup}$ (km)")
    for ax in axs[:, 0]:
        ax.set_ylabel("Median PGA (g), $V_\\mathrm{S30}$ = 300 m/s")
    fig.subplots_adjust(
        left=0.085, right=0.925, bottom=0.085, top=0.985, wspace=0.17, hspace=0.07
    )
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"fig_gmm_spread.{ext}", dpi=600 if ext == "png" else None)

    fig.canvas.draw()
    for ax in axs.flat:
        boxes = [t.get_window_extent() for t in ax.texts]
        if ax.get_legend() is not None:
            boxes.append(ax.get_legend().get_window_extent())
        n = 0
        for ln in ax.lines:
            if ln.get_transform() != ax.transData or len(ln.get_xdata()) < 3:
                continue
            xy = ax.transData.transform(np.column_stack(ln.get_data()).astype(float))
            for b in boxes:
                n += int(
                    (
                        (xy[:, 0] > b.x0)
                        & (xy[:, 0] < b.x1)
                        & (xy[:, 1] > b.y0)
                        & (xy[:, 1] < b.y1)
                    ).sum()
                )
        print(
            ax.texts[0].get_text()[:12] if ax.texts else "",
            "points under text/legend:",
            n,
        )


if __name__ == "__main__":
    main()
