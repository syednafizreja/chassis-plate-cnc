#!/usr/bin/env python3
"""
Independent G-code simulator for chassis_plate_acrylic.nc.
Reads the .nc file only (not the generator), cuts a height-map of the acrylic blank
with flat end mills, then checks the result against the drawing geometry.

    python3 simulate_nc.py chassis_plate_acrylic.nc preview.png
"""
import math
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- drawing (measured from the PDF/PNG) ---------------------------------------------------
W, H, C, THICK = 160.0, 130.0, 15.0, 4.0
OCT = [(C, 0), (W - C, 0), (W, C), (W, H - C), (W - C, H), (C, H), (0, H - C), (0, C)]
FEATS = {
    "slot": [(78, 40), (82, 40), (82, 90), (78, 90)],
    "rect_L1": [(5, 55), (17, 55), (17, 61), (5, 61)],
    "rect_L2": [(5, 70), (17, 70), (17, 76), (5, 76)],
    "rect_R1": [(143, 55), (155, 55), (155, 61), (143, 61)],
    "rect_R2": [(143, 70), (155, 70), (155, 76), (143, 76)],
    "tri_top": [(60, 110), (80, 95), (100, 110)],
    "tri_bot": [(60, 20), (100, 20), (80, 35)],
}
TOOL_DIA = {1: 6.0, 2: 3.0}          # T01 / T02 as stated in the program header
BLANK = (-10.0, 170.0, -10.0, 140.0)  # recommended blank 180 x 150, part centred
RES = 0.05                            # grid size, mm
STEP = 0.2                            # stamp spacing along a move, mm


# ---- parser --------------------------------------------------------------------------------
def parse(path):
    segs, tool = [], None
    pos = [None, None, None]
    mode, f = None, None
    rapids_low = []                    # rapid XY moves made below Z_SAFE-ish
    min_z = 0.0
    for raw in open(path):
        line = re.sub(r"\(.*?\)", "", raw).strip()
        if not line or line == "%" or re.match(r"O\d+", line):
            continue
        words = re.findall(r"([A-Z])\s*(-?\d*\.?\d+)", line)
        gs = [int(float(v)) for k, v in words if k == "G"]
        if 28 in gs:                   # reference return: position unknown afterwards
            for k, v in words:
                if k == "Z":
                    pos[2] = None
                if k in "XY":
                    pos[0] = pos[1] = None
            continue
        new = list(pos)
        for k, v in words:
            v = float(v)
            if k == "G" and int(v) in (0, 1):
                mode = int(v)
            elif k == "T":
                tool = int(v)
            elif k == "F":
                f = v
            elif k in "XYZ":
                new["XYZ".index(k)] = v
        if 43 in gs:                   # tool length offset line sets Z only
            pos = new
            continue
        moved = new != pos
        if moved and mode in (0, 1) and None not in pos and None not in new:
            segs.append(dict(mode=mode, p0=tuple(pos), p1=tuple(new), f=f, tool=tool))
            if mode == 0 and (pos[0] != new[0] or pos[1] != new[1]) and min(pos[2], new[2]) < 1.0:
                rapids_low.append((tuple(pos), tuple(new)))
            min_z = min(min_z, new[2])
        pos = new
    return segs, rapids_low, min_z


# ---- height-map ----------------------------------------------------------------------------
def simulate(segs):
    x0, x1, y0, y1 = BLANK
    xs = np.arange(x0 - 5, x1 + 5, RES, dtype=np.float32)
    ys = np.arange(y0 - 5, y1 + 5, RES, dtype=np.float32)
    h = np.full((len(xs), len(ys)), -99.0, dtype=np.float32)       # no material
    ix = (xs >= x0) & (xs <= x1)
    iy = (ys >= y0) & (ys <= y1)
    h[np.ix_(ix, iy)] = 0.0                                        # top of blank = Z0
    max_doc = {}

    def stamp(cx, cy, z, r, tool):
        i0 = max(int((cx - r - xs[0]) / RES) - 1, 0)
        i1 = min(int((cx + r - xs[0]) / RES) + 2, len(xs))
        j0 = max(int((cy - r - ys[0]) / RES) - 1, 0)
        j1 = min(int((cy + r - ys[0]) / RES) + 2, len(ys))
        sub = h[i0:i1, j0:j1]
        mask = (xs[i0:i1, None] - cx) ** 2 + (ys[None, j0:j1] - cy) ** 2 <= r * r
        vals = (sub - z)[mask & (sub > -50)]
        k = int(0.1 / (RES * RES))                      # 0.1 mm2 of footprint
        doc = float(np.partition(vals, -k)[-k]) if vals.size >= k else 0.0
        if doc > max_doc.get(tool, 0.0):
            max_doc[tool] = doc
        np.minimum(sub, np.float32(z), out=sub, where=mask)

    cut_len = {}
    for s in segs:
        if s["mode"] != 1:
            continue
        r = TOOL_DIA[s["tool"]] / 2
        (xa, ya, za), (xb, yb, zb) = s["p0"], s["p1"]
        d = max(math.hypot(xb - xa, yb - ya), abs(zb - za))
        cut_len[s["tool"]] = cut_len.get(s["tool"], 0.0) + math.dist(s["p0"], s["p1"])
        n = max(1, math.ceil(d / STEP))
        for k in range(n + 1):
            t = k / n
            stamp(xa + (xb - xa) * t, ya + (yb - ya) * t, za + (zb - za) * t, r, s["tool"])
    return xs, ys, h, max_doc, cut_len


# ---- convex polygon helpers ---------------------------------------------------------------
def ccw(p):
    a = sum(p[i][0] * p[(i + 1) % len(p)][1] - p[(i + 1) % len(p)][0] * p[i][1] for i in range(len(p)))
    return p if a > 0 else p[::-1]


def offset_poly(pts, d):
    pts = ccw(pts)
    n = len(pts)
    lines = []
    for i in range(n):
        p, q = pts[i], pts[(i + 1) % n]
        dx, dy = q[0] - p[0], q[1] - p[1]
        L = math.hypot(dx, dy)
        lines.append(((p[0] + dy / L * d, p[1] - dx / L * d), (dx, dy)))
    out = []
    for i in range(n):
        (p1, d1), (p2, d2) = lines[i - 1], lines[i]
        den = d1[0] * d2[1] - d1[1] * d2[0]
        t = ((p2[0] - p1[0]) * d2[1] - (p2[1] - p1[1]) * d2[0]) / den
        out.append((p1[0] + t * d1[0], p1[1] + t * d1[1]))
    return out


def inside(pts, xs, ys):
    pts = ccw(pts)
    m = np.ones((len(xs), len(ys)), dtype=bool)
    for i in range(len(pts)):
        a, b = pts[i], pts[(i + 1) % len(pts)]
        cross = (b[0] - a[0]) * (ys[None, :] - a[1]) - (b[1] - a[1]) * (xs[:, None] - a[0])
        m &= cross >= 0
    return m


def fillet_area(poly, r):
    """Area a round tool of radius r cannot remove in the corners of a convex polygon."""
    p = ccw(poly)
    tot = 0.0
    for i in range(len(p)):
        a, b, c = p[i - 1], p[i], p[(i + 1) % len(p)]
        v1 = (a[0] - b[0], a[1] - b[1])
        v2 = (c[0] - b[0], c[1] - b[1])
        th = math.acos((v1[0] * v2[0] + v1[1] * v2[1]) / (math.hypot(*v1) * math.hypot(*v2)))
        tot += r * r * (1 / math.tan(th / 2) - (math.pi - th) / 2)
    return tot


def main():
    nc = sys.argv[1]
    png = sys.argv[2] if len(sys.argv) > 2 else "preview.png"
    segs, rapids_low, min_z = parse(nc)
    xs, ys, h, max_doc, cut_len = simulate(segs)

    through = h <= -THICK + 0.01
    on_blank = h > -50

    d = 0.10                                          # pixel-noise band around ideal edges
    keep = inside(offset_poly(OCT, -d), xs, ys)
    for f in FEATS.values():
        keep &= ~inside(offset_poly(f, +d), xs, ys)
    gouge = keep & (h < -0.01)                        # ANY removal on the finished plate
    overcut = keep & through

    groove = inside(offset_poly(OCT, +2.9), xs, ys) & ~inside(offset_poly(OCT, +d), xs, ys)
    groove_missed = groove & ~through

    per_feat, sim_tot, th_tot, extra_tot = {}, 0.0, 0.0, 0.0
    for name, f in FEATS.items():
        xa = [p[0] for p in f]; ya = [p[1] for p in f]
        i0 = int((min(xa) - 2 - xs[0]) / RES); i1 = int((max(xa) + 2 - xs[0]) / RES)
        j0 = int((min(ya) - 2 - ys[0]) / RES); j1 = int((max(ya) + 2 - ys[0]) / RES)
        sx, sy = xs[i0:i1], ys[j0:j1]
        want = inside(offset_poly(f, -d), sx, sy)            # must be removed
        # cells a D3 round cutter can physically reach inside this feature
        c = offset_poly(f, -1.5)
        reach = inside(c, sx, sy)
        P = np.stack(np.meshgrid(sx, sy, indexing="ij"), -1)
        dmin = np.full(reach.shape, 1e9)
        for k in range(len(c)):
            a_, b_ = np.array(c[k]), np.array(c[(k + 1) % len(c)])
            ab = b_ - a_
            t = np.clip(((P - a_) @ ab) / (ab @ ab), 0, 1)
            dmin = np.minimum(dmin, np.linalg.norm(P - (a_ + t[..., None] * ab), axis=-1))
        reach |= dmin <= 1.5 + 0.02
        sim_missed = want & ~through[i0:i1, j0:j1]
        theory = want & ~reach
        extra = sim_missed & ~theory
        per_feat[name] = tuple(float(m.sum()) * RES * RES for m in (sim_missed, theory, extra))
        sim_tot += per_feat[name][0]; th_tot += per_feat[name][1]; extra_tot += per_feat[name][2]

    print("=== G-code checks ===")
    print(f"cutting moves             : {sum(1 for s in segs if s['mode']==1)}")
    print(f"lowest Z reached          : {min_z:.3f}  (stock {THICK} + 0.2 overcut = -4.2)")
    print(f"rapid XY moves below Z1   : {len(rapids_low)}  (must be 0)")
    for t in sorted(max_doc):
        print(f"T{t:02d} max depth of cut     : {max_doc[t]:.2f} mm   (path {cut_len[t]:.0f} mm)")
    print("=== final shape vs drawing ===")
    print(f"finished plate surface gouged anywhere : {float(gouge.sum())*RES*RES:.3f} mm2  (must be 0)")
    print(f"plate material cut through wrongly      : {float(overcut.sum())*RES*RES:.3f} mm2  (must be 0)")
    print(f"scrap not cleared, 0-2.9 mm outside edge: {float(groove_missed.sum())*RES*RES:.3f} mm2  (must be 0)")
    print(f"left inside features (cutter corner radii): sim {sim_tot:.2f} mm2 vs theory {th_tot:.2f} mm2")
    print(f"left inside features beyond what a D3 round cutter can reach: {extra_tot:.3f} mm2  (must be ~0)")
    for k, (m, t, e) in per_feat.items():
        print(f"   {k:8s} sim {m:5.2f}  theory {t:5.2f}  extra {e:5.3f}")

    # ---- picture ---------------------------------------------------------------------------
    fig, ax = plt.subplots(1, 2, figsize=(15, 7.2))
    ext = (xs[0], xs[-1], ys[0], ys[-1])
    from matplotlib.colors import ListedColormap
    cls = np.where(~on_blank, 0, np.where(through, 2, 1))
    cmap = ListedColormap(["#3a3a3a", "#8ab6d6", "#ffffff"])
    ax[1].imshow(cls.T, origin="lower", extent=ext, cmap=cmap, vmin=0, vmax=2, interpolation="nearest")
    ox = [p[0] for p in OCT] + [OCT[0][0]]
    oy = [p[1] for p in OCT] + [OCT[0][1]]
    for a in ax:
        a.plot(ox, oy, "k--", lw=0.8)
        a.set_aspect("equal")
        a.set_xlim(-14, 174)
        a.set_ylim(-14, 144)
        a.set_xlabel("X (mm)")
        a.set_ylabel("Y (mm)")
    cols = {1: "#d95f02", 2: "#1b6ca8"}
    for s in segs:
        if s["mode"] == 1 and s["p1"][2] <= -2.0 and s["p0"][2] <= -2.0:
            ax[0].plot([s["p0"][0], s["p1"][0]], [s["p0"][1], s["p1"][1]], color=cols[s["tool"]], lw=0.9)
    ax[0].plot([], [], color=cols[2], label="T02  D3 - internal features")
    ax[0].plot([], [], color=cols[1], label="T01  D6 - outline (tool-centre path)")
    h_, l_ = ax[0].get_legend_handles_labels()
    fig.legend(h_, l_, loc="lower center", ncol=2, fontsize=9, frameon=False)
    ax[0].set_title("Tool-centre paths (cutting moves below Z-2)")
    ax[1].set_title("Simulated result: white = cut through, blue = acrylic left, dashed = drawing")
    fig.suptitle("Chassis plate 160 x 130 - acrylic - simulation of chassis_plate_acrylic.nc", fontsize=12)
    fig.tight_layout(rect=[0, 0.04, 1, 0.95])
    fig.savefig(png, dpi=130)
    print("saved", png)


if __name__ == "__main__":
    main()
