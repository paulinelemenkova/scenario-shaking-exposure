from __future__ import annotations

import csv
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
from matplotlib.axes import Axes
from matplotlib.font_manager import FontProperties, findfont
from matplotlib.lines import Line2D
from matplotlib.ticker import (
    AutoMinorLocator,
    LogLocator,
    NullFormatter,
    ScalarFormatter,
)

FArr = npt.NDArray[np.float64]
IN = Path(os.environ.get("IN", "."))
OUT = Path(os.environ.get("OUT", "."))
EVS = [
    ("us7000hj3u", "2022 Paktika", "#0072B2", "o"),
    ("us6000ldpg", "2023 Herat", "#009E73", "s"),
    ("us7000qsvj", "2025 Kunar", "#D55E00", "^"),
    ("us6000rl31", "2025 Kholm", "#CC79A7", "D"),
]
YLIM = (-0.45, 1.45)


def set_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Nimbus Sans", "Helvetica"],
            "font.size": 8,
            "axes.labelsize": 8.5,
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


def read_csv(name: str) -> list[dict[str, str]]:
    with open(IN / name) as f:
        return list(csv.DictReader(f))


def tag(ax: Axes, text: str) -> None:
    ax.text(
        0.025,
        0.975,
        text,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8.5,
        fontweight="bold",
        bbox={"fc": "white", "ec": "none", "pad": 1.0},
        zorder=8,
    )


def residual_panel(
    ax: Axes, key: str, pts: list[dict[str, str]], binned: list[dict[str, str]]
) -> None:
    ax.axhline(0.0, color="0.25", lw=0.6, zorder=2)
    for k, (ev, _lab, col, mk) in enumerate(EVS):
        x = np.array([float(p[key]) for p in pts if p["event"] == ev])
        y = np.array([float(p["res_ln"]) for p in pts if p["event"] == ev])
        ax.scatter(
            x, y, s=4, marker=mk, facecolor=col, edgecolor="none", alpha=0.28, zorder=3
        )
        b = [r for r in binned if r["event"] == ev and r["x"] == key]
        xc = np.array([float(r["centre"]) for r in b])
        m = np.array([float(r["mean"]) for r in b])
        s = np.array([float(r["sd"]) for r in b])
        off = 1.0 + 0.025 * (k - 1.5) if key == "rjb_km" else 4.0 * (k - 1.5)
        xo = xc * off if key == "rjb_km" else xc + off
        ax.errorbar(
            xo, m, yerr=s, fmt="none", ecolor=col, elinewidth=0.6, capsize=1.5, zorder=4
        )
        ax.plot(
            xo,
            m,
            color=col,
            lw=0.9,
            marker=mk,
            ms=3.2,
            mfc=col,
            mec="black",
            mew=0.4,
            zorder=5,
        )
    ax.set_ylim(*YLIM)
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.grid(True, which="major", color="0.88", lw=0.3)
    ax.set_axisbelow(True)
    ax.tick_params(which="both", top=True, right=True)


def main() -> None:
    set_style()
    pts = read_csv("points.csv")
    binned = read_csv("binned.csv")
    imp = read_csv("importance.csv")
    score = dict(kv.split("=") for kv in (IN / "rf_score.txt").read_text().split())
    fig, axs = plt.subplots(
        1,
        3,
        figsize=(17.4 / 2.54, 6.4 / 2.54),
        gridspec_kw={"width_ratios": [1.0, 1.0, 0.9]},
    )
    a, b, c = axs
    residual_panel(a, "rjb_km", pts, binned)
    a.set_xscale("log")
    a.set_xlim(1.5, 230)
    a.xaxis.set_major_formatter(ScalarFormatter())
    a.xaxis.set_minor_locator(
        LogLocator(subs=[2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0], numticks=10)
    )
    a.xaxis.set_minor_formatter(NullFormatter())
    a.set_xlabel(r"Joyner$-$Boore distance, $R_\mathrm{JB}$ (km)")
    a.set_ylabel(r"$\ln(\mathrm{PGA_{ShakeMap}}\,/\,\mathrm{PGA_{model}})$")
    tag(a, "a)")
    residual_panel(b, "vs30", pts, binned)
    b.set_xlim(150, 950)
    b.xaxis.set_minor_locator(AutoMinorLocator(2))
    b.set_xlabel(r"$V_\mathrm{S30}$ (m s$^{-1}$)")
    b.tick_params(labelleft=False)
    tag(b, "b)")
    hs = [
        Line2D(
            [], [], color=col, lw=0.9, marker=mk, ms=3.2, mfc=col, mec="black", mew=0.4
        )
        for _ev, _lab, col, mk in EVS
    ]
    b.legend(
        hs,
        [lab for _ev, lab, _c, _m in EVS],
        loc="upper right",
        ncol=2,
        frameon=True,
        framealpha=0.95,
        edgecolor="0.6",
        fancybox=False,
        handlelength=1.8,
        columnspacing=0.8,
        borderpad=0.35,
        labelspacing=0.25,
    )

    names = [r["feature"] for r in imp][::-1]
    mean = np.array([float(r["mean"]) for r in imp])[::-1]
    sd = np.array([float(r["sd"]) for r in imp])[::-1]
    yy = np.arange(len(names))
    c.barh(
        yy,
        mean,
        xerr=sd,
        height=0.6,
        color="#56B4E9",
        edgecolor="black",
        linewidth=0.4,
        error_kw={"elinewidth": 0.6, "capsize": 1.8, "ecolor": "black"},
        zorder=3,
    )
    labels = {
        "log10 R_JB": r"$R_\mathrm{JB}$",
        "Vs30": r"$V_\mathrm{S30}$",
        "Distance to observation": "Distance to\nobservation",
        "Azimuth to strike": "Azimuth\nto strike",
    }
    c.set_yticks(yy)
    c.set_yticklabels([labels.get(n, n) for n in names])
    c.set_xlim(0, max(mean + sd) * 1.12)
    c.set_ylim(-0.6, len(names) + 0.75)
    c.xaxis.set_minor_locator(AutoMinorLocator(2))
    c.tick_params(axis="y", which="both", length=0)
    c.tick_params(axis="x", which="both", top=True)
    c.grid(True, axis="x", which="major", color="0.88", lw=0.3)
    c.set_axisbelow(True)
    c.set_xlabel(r"Importance ($\Delta R^2$)")
    tag(c, "c)")
    c.text(
        0.97,
        0.975,
        f"RF, within-event residual\n$n$ = {int(score['n']):,}; CV $R^2$ = "
        f"{float(score['r2_cv']):.2f}",
        transform=c.transAxes,
        ha="right",
        va="top",
        fontsize=6.3,
        bbox={"fc": "white", "ec": "0.6", "lw": 0.4, "pad": 1.6},
        zorder=8,
    )
    fig.subplots_adjust(left=0.075, right=0.99, bottom=0.19, top=0.97, wspace=0.08)
    pos = c.get_position()
    c.set_position((pos.x0 + 0.075, pos.y0, pos.width - 0.085, pos.height))
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"fig_residuals.{ext}", dpi=600 if ext == "png" else None)


if __name__ == "__main__":
    main()
