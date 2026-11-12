from __future__ import annotations

import csv
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.font_manager import FontProperties, findfont
from matplotlib.patches import Patch
from matplotlib.ticker import AutoMinorLocator, MultipleLocator

IN = Path(os.environ.get("IN", "."))
OUT = Path(os.environ.get("OUT", "."))
CLS = ["IV", "V", "VI", "VII", "VIII", "IX+"]

COL = ["#1AD3D1", "#6DFE62", "#D9E436", "#FE8F28", "#D13005", "#7A0403"]
NAME = {
    "Chaman": "Chaman M$_\\mathrm{w}$ 7.0",
    "Herat": "Herat M$_\\mathrm{w}$ 7.5",
    "Hindu Kush": "Hindu Kush M$_\\mathrm{w}$ 7.5",
    "Kholm-type": "Kholm-type M$_\\mathrm{w}$ 6.5",
}
LEVELS = [("p16", "16th"), ("p50", "median"), ("p84", "84th")]
HGT = {"p16": 0.2, "p50": 0.3, "p84": 0.2}
OFF = {"p16": -0.29, "p50": 0.0, "p84": 0.29}


def set_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Nimbus Sans", "Helvetica"],
            "font.size": 7.5,
            "axes.labelsize": 8.5,
            "axes.linewidth": 0.6,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "xtick.minor.width": 0.4,
            "xtick.major.size": 3.5,
            "xtick.minor.size": 2.0,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "legend.fontsize": 7.5,
            "pdf.fonttype": 42,
            "mathtext.fontset": "custom",
            "mathtext.rm": "Nimbus Sans",
            "mathtext.it": "Nimbus Sans:italic",
        }
    )
    fp = FontProperties(family=plt.rcParams["font.sans-serif"])
    assert "dejavu" not in findfont(fp).lower(), "DejaVu fallback: Nimbus Sans missing"


def millions(v: float) -> str:

    if v < 500:
        return "0"
    return f"{v / 1e6:.2f} M" if v >= 1e6 else f"{v / 1e3:.0f} k"


def main() -> None:
    set_style()
    with open(IN / "city_exposure.csv") as f:
        rows = list(csv.DictReader(f))
    data: dict[tuple[str, str], dict[str, list[float]]] = {}
    popn: dict[str, float] = {}
    for r in rows:
        key = (r["city"], r["scenario"])
        data.setdefault(key, {})[r["level"]] = [float(r[c]) for c in CLS]
        popn[r["city"]] = float(r["population"])
    order = sorted(popn, key=lambda c: -popn[c])
    groups = [(c, s) for c in order for s in NAME if sum(data[(c, s)]["p84"]) > 0]
    silent = [c for c in order if not any(g[0] == c for g in groups)]

    fig, ax = plt.subplots(figsize=(17.4 / 2.54, 11.2 / 2.54))
    yy = np.arange(len(groups))[::-1].astype(float)
    for y, (c, s) in zip(yy, groups, strict=True):
        for lev, _ in LEVELS:
            share = 100.0 * np.array(data[(c, s)][lev]) / popn[c]
            left = 0.0
            yb = y + OFF[lev]
            ax.barh(
                yb,
                100.0,
                height=HGT[lev],
                left=0,
                color="white",
                edgecolor="0.55",
                lw=0.3,
                zorder=2,
            )
            for k, w in enumerate(share):
                if w > 0:
                    ax.barh(
                        yb,
                        w,
                        height=HGT[lev],
                        left=left,
                        color=COL[k],
                        edgecolor="none",
                        zorder=3,
                    )
                left += w
            ax.barh(
                yb,
                100.0,
                height=HGT[lev],
                left=0,
                color="none",
                edgecolor="0.25",
                lw=0.3,
                zorder=4,
            )
        ge6 = sum(data[(c, s)]["p50"][2:])
        lo6 = sum(data[(c, s)]["p16"][2:])
        hi6 = sum(data[(c, s)]["p84"][2:])
        ax.text(
            102.5,
            y,
            f"{millions(ge6)}\n({millions(lo6)}$-${millions(hi6)})",
            va="center",
            ha="left",
            fontsize=7,
            linespacing=1.15,
        )

    ax.set_yticks(yy)
    ax.set_yticklabels(
        [f"{c} ({millions(popn[c])})\n{NAME[s]}" for c, s in groups], linespacing=1.2
    )
    ax.tick_params(axis="y", length=0, pad=4)

    for i in range(1, len(groups)):
        if groups[i][0] != groups[i - 1][0]:
            ax.axhline(yy[i] + 0.5, color="0.6", lw=0.4, zorder=1)
    ax.set_xlim(0, 100)
    ax.set_ylim(yy[-1] - 0.55, yy[0] + 0.55)
    ax.xaxis.set_major_locator(MultipleLocator(20))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.tick_params(axis="x", which="both", top=True)
    ax.grid(True, axis="x", which="major", color="0.88", lw=0.3)
    ax.set_axisbelow(True)
    ax.set_xlabel("Share of the urban-centre population (%)")
    ax.text(
        102.5,
        yy[0] + 0.55,
        "People at MMI $\\geq$ VI\nmedian (16th$-$84th)",
        va="bottom",
        ha="left",
        fontsize=7,
        linespacing=1.15,
        clip_on=False,
    )
    handles = [
        Patch(facecolor=COL[k], edgecolor="none", label=CLS[k]) for k in range(len(CLS))
    ]
    handles.append(Patch(facecolor="white", edgecolor="0.55", lw=0.3, label="<IV"))
    ax.legend(
        handles=handles,
        loc="lower left",
        bbox_to_anchor=(0.06, 1.0),
        ncol=len(handles),
        frameon=False,
        handlelength=1.6,
        columnspacing=1.2,
        borderaxespad=0.0,
    )
    ax.text(
        0.0,
        1.0 + 0.012,
        "MMI",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=7.5,
    )
    note = "Bars of each pair: 84th percentile (top), median (thick), 16th percentile (bottom)."
    if silent:
        note += f"\n{', '.join(silent)}: below MMI IV in all scenarios, also at the 84th percentile."
    if note:
        ax.text(
            0.0,
            -0.125,
            note,
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=7,
            color="0.2",
        )
    fig.subplots_adjust(left=0.215, right=0.835, bottom=0.17, top=0.93)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"fig_exposure_bars.{ext}", dpi=600 if ext == "png" else None)


if __name__ == "__main__":
    main()
