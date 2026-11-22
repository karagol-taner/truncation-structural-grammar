"""Place each canonical and its truncated isoform in one frame at one scale.

The isoform is superposed on the canonical through one chain, so a change in the
quaternary arrangement stays visible instead of being averaged away.

Inputs:  auxiliary/fig5_manifest.json, auxiliary/fig5_structures/*.cif
Outputs: auxiliary/fig5_oriented/<gene>_<role>.pdb
         auxiliary/fig5_oriented/frames.json
"""
import json, os, sys, warnings
warnings.filterwarnings("ignore")
import numpy as np
from Bio.PDB import MMCIFParser, PDBIO, Superimposer

ROOT = os.environ.get("AFDB_REPO_DIR", os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")))
MANIFEST = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "auxiliary", "fig5_manifest.json")
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, "auxiliary", "fig5_oriented")


def spans(w):
    out = []
    for part in (w or "").split(";"):
        if "-" in part:
            a, b = part.split("-")
            out.append((int(a), int(b)))
    return out


def ca_by_chain(structure):
    return {ch.id: [r for r in ch if "CA" in r] for ch in structure[0]}


def principal_rotation(coords):
    centred = coords - coords.mean(axis=0)
    _, _, vt = np.linalg.svd(centred, full_matrices=False)
    rot = vt
    if np.linalg.det(rot) < 0:
        rot[2] *= -1
    return rot


def main():
    os.makedirs(OUT, exist_ok=True)
    frames = {}
    for pair in json.load(open(MANIFEST)):
        gene = pair["gene"]
        canon = MMCIFParser(QUIET=True).get_structure("c", os.path.join(ROOT, pair["canonical_cif"]))
        iso = MMCIFParser(QUIET=True).get_structure("i", os.path.join(ROOT, pair["isoform_cif"]))
        deleted = set()
        for a, b in spans(pair["deletion_window"]):
            deleted.update(range(a, b + 1))

        cca, ica = ca_by_chain(canon), ca_by_chain(iso)
        order = sorted(cca)
        best = None
        for ikey in order:
            kept = [r for r in cca[order[0]] if r.id[1] not in deleted]
            if len(kept) != len(ica[ikey]):
                continue
            sup = Superimposer()
            sup.set_atoms([r["CA"] for r in kept], [r["CA"] for r in ica[ikey]])
            if best is None or sup.rms < best[0]:
                best = (sup.rms, sup.rotran)
        if best is None:
            raise SystemExit("no residue correspondence for %s" % gene)
        rms, rotran = best
        for atom in iso.get_atoms():
            atom.transform(*rotran)

        canon_ca = np.array([r["CA"].coord for ch in cca.values() for r in ch])
        iso_ca = np.array([r["CA"].coord for ch in ca_by_chain(iso).values() for r in ch])
        origin = canon_ca.mean(axis=0)
        rot = principal_rotation(np.vstack([canon_ca, iso_ca]))
        for struct in (canon, iso):
            for atom in struct.get_atoms():
                atom.set_coord(rot.dot(atom.coord - origin))

        ordered_radius, full_radius = {}, {}
        for role, struct in (("canonical", canon), ("isoform", iso)):
            ca = np.array([a.coord for a in struct.get_atoms() if a.get_id() == "CA"])
            bf = np.array([a.get_bfactor() for a in struct.get_atoms() if a.get_id() == "CA"])
            dist = np.linalg.norm(ca, axis=1)
            ordered = dist[bf >= 50]
            if len(ordered) < 0.2 * len(ca):
                ordered = dist
            ordered_radius[role] = float(ordered.max())
            full_radius[role] = float(dist.max())
            writer = PDBIO()
            writer.set_structure(struct)
            writer.save(os.path.join(OUT, "%s_%s.pdb" % (gene, role)))
        common = max(ordered_radius.values())
        frames[gene] = {"common_radius": round(common, 2),
                        "zoom": {r: round(full_radius[r] / common, 4) for r in full_radius},
                        "ordered_radius": {r: round(v, 1) for r, v in ordered_radius.items()},
                        "superposition_rmsd": round(float(rms), 2)}
        print("  %-6s chain superposition RMSD %.2f A   ordered radius canonical %.1f  isoform %.1f"
              % (gene, rms, ordered_radius["canonical"], ordered_radius["isoform"]))
    json.dump(frames, open(os.path.join(OUT, "frames.json"), "w"), indent=1)
    print("wrote oriented models to %s" % OUT)


if __name__ == "__main__":
    main()
