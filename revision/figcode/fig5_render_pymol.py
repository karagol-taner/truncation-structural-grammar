"""Render each predicted homodimer coloured by pLDDT, with the removed segment faded.

An alternative to fig5_render_3dmol.py writing the same panels. Both read the
models and the shared scale written by fig5_prepare_models.py.

Inputs:  auxiliary/fig5_manifest.json, auxiliary/fig5_oriented/*.pdb,
         auxiliary/fig5_oriented/frames.json
Outputs: auxiliary/fig5_panels/<gene>_<role>.png

Usage:
    python fig5_render_pymol.py
"""
import json, os

import pymol2

ROOT = os.environ.get("AFDB_REPO_DIR", os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")))
WIDTH, HEIGHT = 2200, 1900
BANDS = (("plddt_very_low", 0xFF7D45, None),
         ("plddt_low", 0xFFDB13, 50),
         ("plddt_high", 0x65CBF3, 70),
         ("plddt_very_high", 0x0053D6, 90))
REMOVED = ("removed", 0xCCCCCC)


def rgb(value):
    return [(value >> 16 & 255) / 255.0, (value >> 8 & 255) / 255.0, (value & 255) / 255.0]


def spans(window):
    out = []
    for part in (window or "").split(";"):
        if "-" in part:
            a, b = part.split("-")
            out.append((int(a), int(b)))
    return out


def frame(cmd, radius):
    cmd.delete("frame")
    for axis in range(3):
        for sign in (1, -1):
            pos = [0.0, 0.0, 0.0]
            pos[axis] = sign * radius
            cmd.pseudoatom("frame", pos=pos)
    cmd.zoom("frame", buffer=0.0)
    cmd.delete("frame")


def style(cmd, obj, deleted):
    cmd.hide("everything", obj)
    cmd.show("cartoon", obj)
    cmd.color(BANDS[0][0], obj)
    for name, _, floor in BANDS[1:]:
        cmd.color(name, "%s and b > %d" % (obj, floor))
    if deleted:
        sel = "%s and (%s)" % (obj, " or ".join("resi %d-%d" % s for s in deleted))
        cmd.color(REMOVED[0], sel)
        cmd.set("cartoon_transparency", 0.67, sel)


def render(cmd):
    for name, value, _ in BANDS:
        cmd.set_color(name, rgb(value))
    cmd.set_color(REMOVED[0], rgb(REMOVED[1]))
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.set("antialias", 2)
    cmd.set("ray_shadows", 0)
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("cartoon_transparency", 0.0)

    oriented = os.path.join(ROOT, "auxiliary", "fig5_oriented")
    panels = os.path.join(ROOT, "auxiliary", "fig5_panels")
    os.makedirs(panels, exist_ok=True)
    frames = json.load(open(os.path.join(oriented, "frames.json")))

    for pair in json.load(open(os.path.join(ROOT, "auxiliary", "fig5_manifest.json"))):
        gene = pair["gene"]
        for role in ("canonical", "isoform"):
            model = os.path.join(oriented, "%s_%s.pdb" % (gene, role))
            obj = "%s_%s" % (gene, role)
            cmd.delete("all")
            cmd.load(model, obj)
            cmd.remove("solvent or hydrogens")
            style(cmd, obj, spans(pair["deletion_window"]) if role == "canonical" else None)
            frame(cmd, frames[gene]["common_radius"])
            out = os.path.join(panels, "%s.png" % obj)
            cmd.png(out, width=WIDTH, height=HEIGHT, dpi=300, ray=1)
            print("  %s" % os.path.basename(out))


def main():
    with pymol2.PyMOL() as session:
        render(session.cmd)


if __name__ == "__main__":
    main()
