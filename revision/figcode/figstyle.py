"""Shared figure style and palettes."""
import matplotlib as mpl

RC = {
    "font.family":      "sans-serif",
    "font.weight":      "bold",
    "axes.labelweight": "bold",
    "axes.titleweight": "bold",
    "font.size":        9.5,
    "axes.linewidth":   1.2,
    "axes.labelsize":   10.5,
    "xtick.labelsize":  9.0,
    "ytick.labelsize":  9.0,
    "legend.fontsize":  9,
    "pdf.fonttype":     42,
    "ps.fonttype":      42,
    "svg.fonttype":     "none",
}

PAL = {
    "Conservation":    "#2ca02c",
    "Destabilization": "#ff7f0e",
    "Ablation":        "#d62728",
    "Emergence":       "#1f77b4",
    "Baseline":        "#7f7f7f",
}

TOPO_ORDER = ["N-Terminal Loss", "Internal Deletion / Splice", "C-Terminal Loss", "Whole-chain"]
TOPO_SHORT = {
    "N-Terminal Loss":            "N-term",
    "Internal Deletion / Splice": "Internal",
    "C-Terminal Loss":            "C-term",
    "Whole-chain":                "Both ends",
}
TOPO_COLOR = {
    "N-Terminal Loss":            "#1f77b4",
    "Internal Deletion / Splice": "#7f7f7f",
    "C-Terminal Loss":            "#d62728",
    "Whole-chain":                "#9467bd",
}

METRIC_COLOR = {"ipTM": "#8c6bb1", "ipSAE_min": "#238b45"}


def apply():
    mpl.rcParams.update(RC)


def save(fig, basename, outdir, dpi=400):
    import os
    os.makedirs(outdir, exist_ok=True)
    for ext in ("png", "pdf", "svg"):
        fig.savefig(os.path.join(outdir, "%s.%s" % (basename, ext)),
                    dpi=dpi, bbox_inches="tight")
