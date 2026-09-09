"""Shared publication style + palette for the sign-reversal paper figures."""
import os, matplotlib as mpl, matplotlib.pyplot as plt

PAL = dict(
    REST="#2c6fbb", STIM="#e07b39",
    LOST="#c1272d",     # up in Rest -> down in Stim48  (edge lost on stimulation)
    GAINED="#2a9d8f",   # down in Rest -> up in Stim48  (edge gained on stimulation)
    GREY="#9aa0a6", CONC="#c7ccd1", CORE="#6a51a3",
    INK="#22262b", HAIR="#dfe2e6", HL="#f2c14e",
)

_SCALE = [1.0]   # current text scale, so panel_tag matches setup()


def setup(scale=1.0):
    """scale < 1 shrinks every text element, for dense multi-panel figures that must
    fit a fixed page size (e.g. the A4 composite Figure 1)."""
    _SCALE[0] = scale

    def s(v):
        return round(v * scale, 2)
    mpl.rcParams.update({
        "figure.dpi": 120, "savefig.dpi": 300,
        "pdf.fonttype": 42, "ps.fonttype": 42,
        "font.family": "DejaVu Sans", "font.size": s(9),
        "axes.titlesize": s(10), "axes.titleweight": "bold", "axes.labelsize": s(9),
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.edgecolor": "#565b61", "axes.linewidth": 0.9,
        "xtick.labelsize": s(8), "ytick.labelsize": s(8), "xtick.color": "#565b61",
        "ytick.color": "#565b61", "legend.fontsize": s(7.5), "legend.frameon": False,
        "axes.grid": False, "figure.facecolor": "white", "axes.facecolor": "white",
    })


# A4 portrait canvas in inches. Figure 1 is saved at exactly this size, so it must NOT be
# written with bbox_inches="tight" (that crops to content and silently discards the page size).
A4_PORTRAIT = (8.27, 11.69)

OUT = ("/Users/sabih/CoDEG_Tcell/paper/figures_final")
os.makedirs(OUT, exist_ok=True)

def panel_tag(ax, s, dx=-0.04, dy=1.06):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=12 * _SCALE[0], fontweight="bold",
            va="top", ha="right", color="#22262b")

def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"), bbox_inches="tight",
                    facecolor="white")
    print("  saved", name + ".pdf/.png")
