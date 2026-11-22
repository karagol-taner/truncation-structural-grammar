"""Render each predicted homodimer coloured by pLDDT, with the removed segment faded.

Both members of a pair are drawn in the frame and at the scale written by
fig5_prepare_models.py, so the isoform can be compared with its canonical.

Inputs:  auxiliary/fig5_manifest.json, auxiliary/fig5_oriented/*.pdb,
         auxiliary/fig5_oriented/frames.json
Outputs: auxiliary/fig5_panels/<gene>_<role>.png

Usage:
    python fig5_render_3dmol.py
"""
import json, os, sys, asyncio

ROOT = os.environ.get("AFDB_REPO_DIR", os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")))
MOL3D = os.environ.get(
    "MOL3D_JS", "https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.4.2/3Dmol-min.js")
WIDTH, HEIGHT = 1100, 950


def spans(w):
    o = []
    for p in (w or "").split(";"):
        if "-" in p:
            s, e = p.split("-"); o.append((int(s), int(e)))
    return o


def html(model, dele, zoom, out):
    txt = open(model).read()
    delset = set()
    for s, e in dele: delset.update(range(s, e + 1))
    dl = "[" + ",".join(str(i) for i in sorted(delset)) + "]"
    js = """
    let v = $3Dmol.createViewer(document.getElementById('v'), {backgroundColor:'white'});
    v.addModel(DATA, 'pdb');
    const plddt = function(atom){
      const b = atom.b;
      if (b >= 90) return '0x0053d6';
      if (b >= 70) return '0x65cbf3';
      if (b >= 50) return '0xffdb13';
      return '0xff7d45';
    };
    v.setViewStyle({style:'outline', color:'black', width:0.03});
    v.setStyle({}, {cartoon:{colorfunc: plddt, thickness:0.3, arrows:true}});
    const DEL = %s;
    if (DEL.length) {
      v.setStyle({resi: DEL}, {cartoon:{color:'0xcccccc', opacity:0.33, thickness:0.14}});
    }
    v.zoomTo();
    v.zoom(%f);
    v.render();
    window.__done = true;
    """ % (dl, zoom)
    page = ("<html><head><meta charset='utf-8'><script src='" + MOL3D + "'></script></head>"
            "<body style='margin:0'><div id='v' style='width:%dpx;height:%dpx;position:relative'></div>"
            % (WIDTH, HEIGHT)
            + "<script>const DATA=" + json.dumps(txt) + ";" + js + "</script></body></html>")
    open(out, "w").write(page)


async def shot(htmlfile, png):
    from playwright.async_api import async_playwright
    async with async_playwright() as p:
        br = await p.chromium.launch()
        pg = await br.new_page(viewport={"width": WIDTH, "height": HEIGHT},
                               device_scale_factor=2)
        await pg.goto("file://" + os.path.abspath(htmlfile))
        await pg.wait_for_function("window.__done === true", timeout=90000)
        await pg.wait_for_timeout(1400)
        await pg.locator("#v").screenshot(path=png)
        await br.close()


def main():
    pairs = json.load(open(os.path.join(ROOT, "auxiliary", "fig5_manifest.json")))
    oriented = os.path.join(ROOT, "auxiliary", "fig5_oriented")
    frames = json.load(open(os.path.join(oriented, "frames.json")))
    panels = os.path.join(ROOT, "auxiliary", "fig5_panels")
    os.makedirs(panels, exist_ok=True)
    for p in pairs:
        for role in ("canonical", "isoform"):
            window = p["deletion_window"] if role == "canonical" else ""
            png = os.path.join(panels, "%s_%s.png" % (p["gene"], role))
            h = png.replace(".png", ".html")
            html(os.path.join(oriented, "%s_%s.pdb" % (p["gene"], role)),
                 spans(window), frames[p["gene"]]["zoom"][role], h)
            asyncio.run(shot(h, png))
            print("  %s at scale %.3f" % (os.path.basename(png),
                                          frames[p["gene"]]["zoom"][role]))


if __name__ == "__main__":
    main()
