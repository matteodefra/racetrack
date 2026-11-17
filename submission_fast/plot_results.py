#!/usr/bin/env python3
"""
Plots of output/sweeps_merged.csv (made by run_sweeps_parallel.sh):
DP (= ALGOWIN = Bekos et al., 1 thread) against DP-parallel (T threads).

Main plots (mean time per instance, log scale):
    plot_time_vs_n(df, L)          L fixed, n on the x axis            (sweep "n")
    plot_time_vs_L(df, n)          n fixed, L on the x axis            (sweep "L")
Speed-up plots (speed-up = t(DP) / t(DP-parallel), same instances):
    plot_speedup_vs_n(df)          against n, one line per L, with the ideal value
    plot_speedup_vs_L(df, n)       against L (n fixed), with the ideal value
    plot_thread_scaling(scaling)   speed-up and efficiency against the number of
                                   threads (sweep "T"), with the ideal line

Usage:
    python3 plot_results.py output/sweeps_merged.csv [--out plots] [--format pdf]
It makes every plot that the data allows. The functions can also be imported.
"""
import argparse
import os

import matplotlib.pyplot as plt
import pandas as pd

DP_COLOR, PAR_COLOR = "#2a78d6", "#eb6834"          # blue, orange
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]
MARKERS = "os^Dv*"
IDEAL = dict(color="#5f5e5a", linewidth=1, linestyle=":")

plt.rcParams.update({
    "font.size": 11, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": "#e4e3df", "grid.linewidth": 0.8,
    "axes.edgecolor": "#8a8984", "axes.labelcolor": "#2b2b29",
    "xtick.color": "#5f5e5a", "ytick.color": "#5f5e5a", "legend.frameon": False,
})


# --------------------------------------------------------------------- data
def load_results(path):
    """Read the CSV, keep the rows with status ok, add dp_s, dp_parallel_s, speedup.
    Returns (main, scaling): main = sweeps n and L, one row per (n, L) (a point in two
    sweeps is averaged); scaling = sweep T, one row per number of threads."""
    df = pd.read_csv(path)
    bad = df[df["status"] != "ok"]
    if len(bad):
        print(f"ignored {len(bad)} rows (status != ok):")
        print(bad[["sweep", "status", "n", "L", "threads"]].to_string(index=False))
    df = df[df["status"] == "ok"].copy()
    df["dp_s"] = df["dp_us"] / 1e6
    df["dp_parallel_s"] = df["dp_parallel_us"] / 1e6
    df["speedup"] = df["dp_us"] / df["dp_parallel_us"]
    num = [c for c in df.select_dtypes("number").columns if c not in ("n", "L", "threads")]
    main = df[df["sweep"].isin(["n", "L"])]
    main = main.groupby(["n", "L", "threads"], as_index=False)[num].mean()
    main["speedup"] = main["dp_us"] / main["dp_parallel_us"]
    scaling = df[df["sweep"] == "T"].sort_values("threads")
    return main, scaling


def _new_ax(ax):
    return ax or plt.subplots(figsize=(6.5, 4.2))[1]


def _threads(df):
    return int(df["threads"].iloc[0])


def _time_lines(ax, sub, x):
    sub = sub.sort_values(x)
    ax.plot(sub[x], sub["dp_s"], color=DP_COLOR, marker="o", markersize=6,
            linewidth=2, label="DP (1 thread)")
    ax.plot(sub[x], sub["dp_parallel_s"], color=PAR_COLOR, marker="s", linestyle="--",
            markersize=6, linewidth=2, label=f"DP parallel ({_threads(sub)} threads)")
    ax.set_yscale("log")
    ax.set_ylabel("mean time per instance (s)")
    ax.legend(loc="best")


def _ideal(ax, value, text):
    ax.axhline(value, **IDEAL)
    ax.text(0.99, value, text, transform=ax.get_yaxis_transform(),
            ha="right", va="bottom", fontsize=8, color="#5f5e5a")


def _top_smax(ax, sub):
    top = ax.secondary_xaxis("top")
    top.set_xticks(sub["L"].tolist(), [str(int(s)) for s in sub["smax"]])
    top.set_xlabel("smax", fontsize=9)
    top.tick_params(labelsize=8)


def _sweeps_in_n(df):
    return sorted(L for L, g in df.groupby("L") if g["n"].nunique() >= 3)


# --------------------------------------------------------------------- main plots
def plot_time_vs_n(df, L, ax=None):
    """L fixed, n on the x axis."""
    sub = df[df["L"] == L]
    if sub.empty:
        raise ValueError(f"no data for L={L}")
    ax = _new_ax(ax)
    _time_lines(ax, sub, "n")
    ax.set_xlabel("n (number of cities)")
    ax.set_title(f"Time against n  (L = {L}, smax = {int(sub['smax'].iloc[0])})", loc="left")
    return ax


def plot_time_vs_L(df, n, ax=None):
    """n fixed, L on the x axis (smax = floor(sqrt(3L)) on a top axis)."""
    sub = df[df["n"] == n].sort_values("L")
    if sub.empty:
        raise ValueError(f"no data for n={n}")
    ax = _new_ax(ax)
    _time_lines(ax, sub, "L")
    ax.set_xlabel("L (side of the box)")
    ax.set_title(f"Time against L  (n = {n}, smax = ⌊√(3L)⌋)", loc="left")
    _top_smax(ax, sub)
    return ax


# --------------------------------------------------------------------- speed-up plots
def plot_speedup_vs_n(df, Ls=None, ax=None):
    """t(DP) / t(DP parallel) against n, one line per L."""
    Ls = Ls or _sweeps_in_n(df)
    ax = _new_ax(ax)
    for i, L in enumerate(Ls):
        sub = df[df["L"] == L].sort_values("n")
        ax.plot(sub["n"], sub["speedup"], color=SERIES[i % 6], marker=MARKERS[i % 6],
                markersize=6, linewidth=2, label=f"L = {L}")
    thr = _threads(df)
    _ideal(ax, thr, f"ideal ({thr} threads)")
    ax.set_ylim(0, thr * 1.1)
    ax.set_xlabel("n (number of cities)")
    ax.set_ylabel("speed-up  t(DP) / t(DP parallel)")
    ax.set_title(f"Speed-up of the parallel DP ({thr} threads) against n", loc="left")
    ax.legend(loc="best")
    return ax


def plot_speedup_vs_L(df, n, ax=None):
    """t(DP) / t(DP parallel) against L, n fixed."""
    sub = df[df["n"] == n].sort_values("L")
    ax = _new_ax(ax)
    ax.plot(sub["L"], sub["speedup"], color=PAR_COLOR, marker="s", markersize=6, linewidth=2)
    thr = _threads(sub)
    _ideal(ax, thr, f"ideal ({thr} threads)")
    ax.set_ylim(0, thr * 1.1)
    ax.set_xlabel("L (side of the box)")
    ax.set_ylabel("speed-up  t(DP) / t(DP parallel)")
    ax.set_title(f"Speed-up of the parallel DP ({thr} threads) against L  (n = {n})", loc="left")
    _top_smax(ax, sub)
    return ax


def plot_thread_scaling(scaling, axes=None):
    """Sweep T: speed-up and efficiency (speed-up / threads) against the number of
    threads, with the ideal values. Two panels."""
    if axes is None:
        _, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    s = scaling.sort_values("threads")
    T, sp = s["threads"], s["speedup"]
    n, L = int(s["n"].iloc[0]), int(s["L"].iloc[0])

    ax = axes[0]
    ax.plot(T, T, label="ideal", **IDEAL)
    ax.plot(T, sp, color=PAR_COLOR, marker="s", linestyle="--", markersize=6,
            linewidth=2, label="DP parallel")
    ax.set_xscale("log", base=2); ax.set_yscale("log", base=2)
    ax.set_xticks(T.tolist(), [str(int(t)) for t in T])
    ax.set_xlabel("threads"); ax.set_ylabel("speed-up  t(DP) / t(DP parallel)")
    ax.set_title(f"Scaling in threads  (n = {n}, L = {L})", loc="left")
    ax.legend(loc="upper left")

    ax = axes[1]
    ax.plot(T, sp / T, color=PAR_COLOR, marker="s", linestyle="--", markersize=6, linewidth=2)
    ax.axhline(1, **IDEAL)
    ax.set_xscale("log", base=2)
    ax.set_xticks(T.tolist(), [str(int(t)) for t in T])
    ax.set_ylim(0, 1.1)
    ax.set_xlabel("threads"); ax.set_ylabel("efficiency  (speed-up / threads)")
    ax.set_title("Parallel efficiency", loc="left")
    return axes


# --------------------------------------------------------------------- all plots
def main():
    p = argparse.ArgumentParser()
    p.add_argument("csv")
    p.add_argument("--out", default="plots")
    p.add_argument("--format", default="pdf", help="pdf, png, svg ...")
    a = p.parse_args()
    os.makedirs(a.out, exist_ok=True)
    df, scaling = load_results(a.csv)

    def save(obj, name):
        fig = obj[0].figure if hasattr(obj, "__len__") else obj.figure
        fig.tight_layout()
        path = os.path.join(a.out, f"{name}.{a.format}")
        fig.savefig(path, dpi=200)
        plt.close(fig)
        print("wrote", path)

    for L in _sweeps_in_n(df):                       # sweeps in n
        save(plot_time_vs_n(df, L), f"time_vs_n_L{L}")
    if _sweeps_in_n(df):
        save(plot_speedup_vs_n(df), "speedup_vs_n")
    for n, g in df.groupby("n"):                     # sweeps in L
        if g["L"].nunique() >= 3:
            save(plot_time_vs_L(df, n), f"time_vs_L_n{n}")
            save(plot_speedup_vs_L(df, n), f"speedup_vs_L_n{n}")
    if len(scaling) >= 2:                            # sweep in threads
        save(plot_thread_scaling(scaling), "thread_scaling")


if __name__ == "__main__":
    main()
