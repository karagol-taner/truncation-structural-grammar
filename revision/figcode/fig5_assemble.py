"""Figure 5. Four predicted homodimers, canonical beside isoform, coloured by pLDDT.

Inputs:  auxiliary/fig5_manifest.json, auxiliary/fig5_structures/*.cif,
         and the panel renders written by fig5_render_3dmol.py
Outputs: figures/Figure5_interface_structures.{png,pdf,svg}
         auxiliary/fig5_plddt.json
"""
import json, os, sys, warnings
warnings.filterwarnings("ignore")
from PIL import Image
import numpy as np
from Bio.PDB import MMCIFParser
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.patches as mpatches

ROOT = os.environ.get("AFDB_REPO_DIR", os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")))
MANIFEST = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "auxiliary", "fig5_manifest.json")
PANELS = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, "auxiliary", "fig5_panels")
OUTDIR = sys.argv[3] if len(sys.argv) > 3 else os.path.join(ROOT, "figures")

pairs = json.load(open(MANIFEST))
SHORT = {"TYW5": "conservation", "PDCD6": "ablation",
         "EPHB4": "emergence", "GFER": "metric shift"}

CIF_KEY = {"canonical": "canonical_cif", "isoform": "isoform_cif"}

def content_box(im, pad=26):
    mask = im.convert("L").point(lambda v: 0 if v > 247 else 255)
    box = mask.getbbox()
    if box is None:
        return (0, 0, im.size[0], im.size[1])
    x0, y0, x1, y1 = box
    return (max(0, x0 - pad), max(0, y0 - pad),
            min(im.size[0], x1 + pad), min(im.size[1], y1 + pad))


def row_images(pair):
    ims = [Image.open(os.path.join(PANELS, "%s_%s.png" % (pair["gene"], role)))
           for role in ("canonical", "isoform")]
    boxes = [content_box(im) for im in ims]
    shared = (min(b[0] for b in boxes), min(b[1] for b in boxes),
              max(b[2] for b in boxes), max(b[3] for b in boxes))
    return [im.crop(shared) for im in ims]


pl = {}
for p in pairs:
    for role in ("canonical", "isoform"):
        cif = os.path.join(ROOT, p[CIF_KEY[role]])
        st = MMCIFParser(QUIET=True).get_structure("x", cif)
        pl[(p["gene"], role)] = float(np.mean([a.get_bfactor() for a in st.get_atoms()]))

plt.rcParams.update({"font.family": "sans-serif", "font.weight": "bold",
                     "axes.titleweight": "bold"})
rows = [row_images(p) for p in pairs]
ratios = [r[0].size[1] / r[0].size[0] for r in rows]
panel_w = 9.6 / 2
height = panel_w * sum(ratios) + 0.42 * len(rows) + 0.5
fig, axes = plt.subplots(len(pairs), 2, figsize=(9.6, height), dpi=200,
                         gridspec_kw={"height_ratios": ratios})
for i, p in enumerate(pairs):
    for j, role in enumerate(("canonical", "isoform")):
        ax = axes[i][j]
        ax.imshow(rows[i][j])
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values(): s.set_visible(False)
        ipt = p["canon_ipTM"] if role == "canonical" else p["iso_ipTM"]
        ips = p["canon_ipSAE"] if role == "canonical" else p["iso_ipSAE"]
        ax.set_title("%s   ipTM %.2f   ipSAE %.3f   pLDDT %.1f"
                     % (role, ipt, ips, pl[(p["gene"], role)]), fontsize=8.5, pad=3)
    axes[i][0].text(-0.075, 0.5, "%s\n%s" % (p["gene"], SHORT[p["gene"]]),
                    transform=axes[i][0].transAxes, rotation=90, va="center",
                    ha="center", fontsize=9.5, fontweight="bold")

h = [mpatches.Patch(color="#0053d6", label="pLDDT > 90"),
     mpatches.Patch(color="#65cbf3", label="70 to 90"),
     mpatches.Patch(color="#ffdb13", label="50 to 70"),
     mpatches.Patch(color="#ff7d45", label="< 50"),
     mpatches.Patch(color="#cccccc", label="removed by truncation, shown faded")]
fig.legend(handles=h, loc="lower center", ncol=5, frameon=False,
           bbox_to_anchor=(0.5, 0.004), fontsize=8.2)
fig.tight_layout(rect=[0.055, 0.03, 1, 1])
os.makedirs(OUTDIR, exist_ok=True)
for ext in ("png", "pdf", "svg"):
    fig.savefig(os.path.join(OUTDIR, "Figure5_interface_structures." + ext),
                dpi=400, bbox_inches="tight")

json.dump({"%s_%s" % (g.lower(), r): round(v, 1) for (g, r), v in sorted(pl.items())},
          open(os.path.join(ROOT, "auxiliary", "fig5_plddt.json"), "w"), indent=1)
print("wrote Figure5_interface_structures to %s, %d rows" % (OUTDIR, len(pairs)))
for k in sorted(pl):
    print("    %-8s %-10s mean pLDDT %.1f" % (k[0], k[1], pl[k]))
