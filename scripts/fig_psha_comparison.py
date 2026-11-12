from __future__ import annotations

import csv
import os
from itertools import pairwise
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.font_manager import FontProperties, findfont
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, LogLocator, NullFormatter

IN = Path(os.environ.get("IN", "."))
OUT = Path(os.environ.get("OUT", "."))
CITIES = ["Kabul", "Herat", "Kandahar", "Mazar-i-Sharif", "Jalalabad", "Kunduz"]
SCEN = [
    ("Chaman", "Chaman M$_\\mathrm{w}$ 7.0", "#D55E00", "o"),
    ("Herat", "Herat M$_\\mathrm{w}$ 7.5", "#0072B2", "s"),
    ("Hindu Kush", "Hindu Kush M$_\\mathrm{w}$ 7.5", "#009E73", "^"),
    ("Kholm-type", "Kholm-type M$_\\mathrm{w}$ 6.5", "#CC79A7", "D"),
]
MMI = [
    (2.76, "IV"),
    (6.2, "V"),
    (11.5, "VI"),
    (21.5, "VII"),
    (40.1, "VIII"),
    (74.7, "IX"),
]


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
            "ytick.minor.width": 0.4,
            "xtick.major.size": 3.5,
            "ytick.major.size": 3.5,
            "ytick.minor.size": 2.0,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "legend.fontsize": 7,
            "pdf.fonttype": 42,
            "mathtext.fontset": "custom",
            "mathtext.rm": "Nimbus Sans",
            "mathtext.it": "Nimbus Sans:italic",
        }
    )
    fp = FontProperties(family=plt.rcParams["font.sans-serif"])
    assert "dejavu" not in findfont(fp).lower(), "DejaVu fallback: Nimbus Sans missing"


def main() -> None:
    set_style()
    with open(IN / "psha_cities.csv") as f:
        psha = {r["city"]: r for r in csv.DictReader(f)}
    with open(IN / "scenario_cities.csv") as f:
        sc = {
            (r["city"], r["scenario"]): (
                float(r["pga_median_pctg"]),
                float(r["sigma_ln"]),
            )
            for r in csv.DictReader(f)
            if r["site"] == "rock"
        }
    fig, ax = plt.subplots(figsize=(17.4 / 2.54, 8.6 / 2.54))
    hw = 0.42
    for i, c in enumerate(CITIES):
        p = psha[c]
        ax.hlines(
            float(p["gem475_pctg"]), i - hw, i + hw, color="black", lw=1.0, zorder=4
        )
        ax.hlines(
            float(p["usgs475_pctg"]),
            i - hw,
            i + hw,
            color="0.35",
            lw=0.8,
            ls=(0, (4, 2)),
            zorder=4,
        )
        ax.hlines(
            float(p["usgs2475_pctg"]),
            i - hw,
            i + hw,
            color="0.35",
            lw=0.8,
            ls=(0, (1, 1.2)),
            zorder=4,
        )
        for k, (s, _, col, mk) in enumerate(SCEN):
            med, sg = sc[(c, s)]
            x = i + (k - 1.5) * 0.17
            ax.vlines(
                x, med * np.exp(-sg), med * np.exp(sg), color=col, lw=0.8, zorder=5
            )
            ax.plot(
                x,
                med,
                marker=mk,
                ms=4.6,
                mfc=col,
                mec="black",
                mew=0.4,
                ls="none",
                zorder=6,
            )
        if i:
            ax.axvline(i - 0.5, color="0.8", lw=0.4, zorder=1)

    ax.set_yscale("log")
    ax.set_ylim(0.02, 160)
    for lim, _ in MMI:
        ax.axhline(lim, color="0.88", lw=0.4, zorder=0)
    ax2 = ax.twinx()
    ax2.set_yscale("log")
    ax2.set_ylim(ax.get_ylim())
    mids = [np.sqrt(a * b) for (a, _), (b, _) in pairwise(MMI)]
    ax2.yaxis.set_major_locator(FixedLocator([*mids, 74.7 * 1.45]))
    ax2.set_yticklabels([lab for _, lab in MMI])
    ax2.yaxis.set_minor_locator(FixedLocator([]))
    ax2.tick_params(axis="y", length=0, labelsize=7, colors="0.3")
    ax2.set_ylabel("MMI (ShakeMap PGA limits)", color="0.3")
    ax.set_xlim(-0.5, len(CITIES) - 0.5)
    ax.set_xticks(range(len(CITIES)))
    ax.set_xticklabels(CITIES)
    ax.tick_params(axis="x", length=0)
    ax.yaxis.set_major_locator(LogLocator(base=10, numticks=10))
    ax.yaxis.set_minor_locator(
        LogLocator(base=10, subs=[2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0], numticks=10)
    )
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_yticks([0.1, 1, 10, 100])
    ax.set_yticklabels(["0.1", "1", "10", "100"])
    ax.tick_params(axis="y", which="both", left=True)
    ax.set_ylabel("PGA on rock, $V_\\mathrm{S30}$ = 760 m s$^{-1}$ (% $g$)")
    h_s = [
        Line2D(
            [],
            [],
            color=col,
            marker=mk,
            ms=4.6,
            mfc=col,
            mec="black",
            mew=0.4,
            lw=0.8,
            label=lab,
        )
        for _, lab, col, mk in SCEN
    ]
    h_p = [
        Line2D([], [], color="black", lw=1.0, label="GEM 2023, 475 yr"),
        Line2D([], [], color="0.35", lw=0.8, ls=(0, (4, 2)), label="USGS 2007, 475 yr"),
        Line2D(
            [], [], color="0.35", lw=0.8, ls=(0, (1, 1.2)), label="USGS 2007, 2475 yr"
        ),
    ]
    blank = Line2D([], [], color="none", label=" ")
    order = [h_s[0], h_p[0], h_s[1], h_p[1], h_s[2], h_p[2], h_s[3], blank]
    ax.legend(
        handles=order,
        loc="lower left",
        bbox_to_anchor=(0.0, 1.01),
        ncol=4,
        frameon=False,
        borderaxespad=0.0,
        handlelength=2.4,
        columnspacing=1.4,
        labelspacing=0.35,
    )
    fig.subplots_adjust(left=0.085, right=0.915, bottom=0.09, top=0.86)
    for ext in ("pdf", "png"):
        fig.savefig(
            OUT / f"fig_psha_comparison.{ext}", dpi=600 if ext == "png" else None
        )


if __name__ == "__main__":
    main()
