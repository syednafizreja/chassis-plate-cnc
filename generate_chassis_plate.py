#!/usr/bin/env python3
"""
Chassis / mounting plate - CNC mill program generator (ACRYLIC)
Southwest Petroleum University - Engineering Training Vehicle (two-shaft, four-wheel)

Writes a Fanuc-style ISO G-code program (chassis_plate_acrylic.nc).
All geometry is measured from the plate drawing (160 x 130 mm).
Edit the PART / PROCESS blocks below and re-run to regenerate the program.

    python3 generate_chassis_plate.py [output.nc]
"""
import math
import sys
from bisect import bisect_right
from pathlib import Path

# =====================================================================
# PART GEOMETRY  (mm, origin X0 Y0 = lower-left corner of FINISHED plate)
# =====================================================================
PLATE_W, PLATE_H = 160.0, 130.0
THICK = 4.0            # acrylic sheet thickness (drawing says 3-5 mm; 4 mm assumed)
CHAMFER = 15.0         # four corner cuts, 15 x 15 as drawn (drawing label says "R5")

SLOT = (78.0, 40.0, 82.0, 90.0)            # centre alignment slot 4 x 50 (x0,y0,x1,y1)
RECTS = [                                  # side mounting windows 12 x 6
    (5.0, 55.0, 17.0, 61.0),
    (5.0, 70.0, 17.0, 76.0),
    (143.0, 70.0, 155.0, 76.0),
    (143.0, 55.0, 155.0, 61.0),
]
TRIS = [                                   # lightening cut-outs, 40 wide x 15 high
    [(60.0, 110.0), (100.0, 110.0), (80.0, 95.0)],   # top, pointing down
    [(60.0, 20.0), (100.0, 20.0), (80.0, 35.0)],     # bottom, pointing up
]

# =====================================================================
# PROCESS  (acrylic - keep chip load up, never dwell, no flood coolant)
# =====================================================================
OVERCUT = 0.2          # break through into the spoilboard
DOC_MAX = 2.1          # max depth per pass  -> 2 passes for 4 mm + 0.2
RAMP_DEG = 4.0         # ramp-in angle (no straight plunging)
Z_SAFE = 5.0           # rapid height
Z_APP = 1.0            # feed-down start height
PLUNGE_F = 300         # feed from Z_APP down to Z0 (in air)

TOOLS = {
    # H-number = T-number.  Carbide single-flute (O-flute) cutters for plastics.
    "T02": dict(n=2, dia=3.0, rpm=8000, feed=350, ramp=200, step=1.5,
                name="D3.0 SINGLE FLUTE O-FLUTE CARBIDE"),
    "T01": dict(n=1, dia=6.0, rpm=8000, feed=600, ramp=250, step=3.0,
                name="D6.0 SINGLE FLUTE O-FLUTE CARBIDE"),
}
START_X_OUTLINE = 120.0   # outline ramp-in starts on the bottom edge, in the scrap groove


# =====================================================================
# geometry helpers
# =====================================================================
def area2(p):
    return sum(p[i][0] * p[(i + 1) % len(p)][1] - p[(i + 1) % len(p)][0] * p[i][1]
               for i in range(len(p)))


def ccw(p):
    return list(p) if area2(p) > 0 else list(p)[::-1]


def rect_poly(r):
    x0, y0, x1, y1 = r
    return [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]


def offset_poly(pts, d):
    """Mitred offset of a CONVEX CCW polygon. d>0 outward, d<0 inward (exact for convex shapes)."""
    n = len(pts)
    lines = []
    for i in range(n):
        p, q = pts[i], pts[(i + 1) % n]
        dx, dy = q[0] - p[0], q[1] - p[1]
        L = math.hypot(dx, dy)
        nx, ny = dy / L, -dx / L                      # outward normal of a CCW polygon
        lines.append(((p[0] + nx * d, p[1] + ny * d), (dx, dy)))
    out = []
    for i in range(n):
        (p1, d1), (p2, d2) = lines[i - 1], lines[i]
        den = d1[0] * d2[1] - d1[1] * d2[0]
        t = ((p2[0] - p1[0]) * d2[1] - (p2[1] - p1[1]) * d2[0]) / den
        out.append((p1[0] + t * d1[0], p1[1] + t * d1[1]))
    return out


def inradius(tri):
    a = abs(area2(tri)) / 2.0
    per = sum(math.dist(tri[i], tri[(i + 1) % 3]) for i in range(3))
    return 2 * a / per


class Loop:
    """Closed polyline that can be walked by arc length (wraps around)."""

    def __init__(self, pts):
        self.p = [(float(x), float(y)) for x, y in pts]
        self.n = len(self.p)
        self.cum = [0.0]
        for i in range(self.n):
            self.cum.append(self.cum[-1] + math.dist(self.p[i], self.p[(i + 1) % self.n]))
        self.L = self.cum[-1]

    def at(self, s):
        s2 = s % self.L
        i = min(bisect_right(self.cum, s2) - 1, self.n - 1)
        seg = self.cum[i + 1] - self.cum[i]
        t = (s2 - self.cum[i]) / seg if seg > 0 else 0.0
        a, b = self.p[i], self.p[(i + 1) % self.n]
        return a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t

    def path(self, s0, s1):
        """[(x, y, s)] from s0 to s1 including every vertex in between."""
        pts = [(*self.at(s0), s0)]
        for lap in range(int(s0 // self.L), int(s1 // self.L) + 1):
            for i in range(self.n):
                sv = self.cum[i] + lap * self.L
                if s0 + 1e-9 < sv < s1 - 1e-9:
                    pts.append((*self.p[i], sv))
        pts.append((*self.at(s1), s1))
        return pts

    def nearest(self, xy):
        return min(range(self.n), key=lambda i: math.dist(self.p[i], xy))


# =====================================================================
# G-code writer
# =====================================================================
def f3(v):
    v = round(v, 3)
    if abs(v) < 0.0005:
        v = 0.0
    return f"{v:.3f}".rstrip("0")        # always keeps the decimal point: 80.  79.5  -2.1


class NC:
    def __init__(self):
        self.out = []
        self.pos = [None, None, None]
        self.mode = None
        self.f = None
        self.time = 0.0          # cutting time, minutes
        self.cut_len = 0.0

    def raw(self, s):
        self.out.append(s)

    def cmt(self, s):
        self.out.append(f"({s})")

    def _move(self, g, x, y, z, f=None):
        parts = []
        new = list(self.pos)
        for i, (name, v) in enumerate((("X", x), ("Y", y), ("Z", z))):
            if v is None:
                continue
            v = round(v, 3)
            if self.pos[i] is None or abs(v - self.pos[i]) >= 0.0005:
                parts.append(f"{name}{f3(v)}")
            new[i] = v
        if not parts:
            return
        if g == "G01" and f is not None and f != self.f:
            parts.append(f"F{int(f)}")
            self.f = f
        if g == "G01":
            d = math.dist([p if p is not None else n for p, n in zip(self.pos, new)], new)
            self.cut_len += d
            self.time += d / (f or self.f)
        head = g + " " if g != self.mode else ""
        self.mode = g
        self.out.append(head + " ".join(parts))
        self.pos = new

    def rapid(self, x=None, y=None, z=None):
        self._move("G00", x, y, z)

    def feed(self, x=None, y=None, z=None, f=None):
        self._move("G01", x, y, z, f)


# =====================================================================
# toolpath building blocks
# =====================================================================
def make_levels():
    total = THICK + OVERCUT
    n = math.ceil(round(total / DOC_MAX, 6))
    return [-total * (k + 1) / n for k in range(n)]


def start_tool(nc, key, x, y, first=False):
    t = TOOLS[key]
    if not first:
        nc.rapid(z=Z_SAFE)
        nc.raw("M05")
        nc.raw("G91 G28 Z0.")
        nc.raw("G90")
    nc.cmt(f"{key} = {t['name']}")
    nc.raw(f"T{t['n']:02d} M06")
    nc.raw(f"G00 G90 G54 X{f3(x)} Y{f3(y)} S{t['rpm']} M03")
    nc.raw(f"G43 H{t['n']:02d} Z{f3(Z_SAFE)}")
    nc.mode = "G00"
    nc.pos = [round(x, 3), round(y, 3), Z_SAFE]
    nc.f = None


def approach(nc, tool, x, y):
    nc.rapid(z=Z_SAFE)
    nc.rapid(x=x, y=y)
    nc.rapid(z=Z_APP)
    nc.feed(z=0.0, f=PLUNGE_F)


def lap(nc, loop, s0, z, feed):
    for x, y, _ in loop.path(s0, s0 + loop.L)[1:]:
        nc.feed(x, y, z, feed)


def ramp_lap(nc, loop, s0, z_from, z_to, tool):
    """Ramp down along the loop from z_from to z_to, then one full flat lap at z_to."""
    D = (z_from - z_to) / math.tan(math.radians(RAMP_DEG))
    s1 = s0 + D
    for x, y, s in loop.path(s0, s1)[1:]:
        nc.feed(x, y, z_from + (z_to - z_from) * (s - s0) / D, tool["ramp"])
    lap(nc, loop, s1, z_to, tool["feed"])
    return s1


def ring_feature(nc, tool, rings, entry, levels, center=None):
    """
    rings : closed tool-centre loops, rings[0] = finishing ring (outermost), last = innermost
    entry : ring used for the ramp-in.  Order = entry -> inner rings -> back out to rings[0],
            so the finishing ring is always cut last and the wall stays attached to the plate.
    center: optional ((x,y),(x,y)) clean-up pass along the centre line at full depth.
    """
    order = list(range(entry, len(rings))) + list(range(entry - 1, -1, -1))
    z_prev = 0.0
    for k, z in enumerate(levels):
        for j, ri in enumerate(order):
            loop = rings[ri]
            i0 = loop.nearest((nc.pos[0], nc.pos[1]))
            s0 = loop.cum[i0]
            x0, y0 = loop.p[i0]
            if j == 0:
                if k == 0:
                    approach(nc, tool, x0, y0)
                else:
                    nc.feed(x0, y0, z_prev, tool["feed"])
                ramp_lap(nc, loop, s0, z_prev, z, tool)
            else:
                nc.feed(x0, y0, z, tool["feed"])
                lap(nc, loop, s0, z, tool["feed"])
        z_prev = z
    if center:
        (xa, ya), (xb, yb) = center
        nc.feed(xa, ya, levels[-1], tool["feed"])
        nc.feed(xb, yb, levels[-1], tool["feed"])
    nc.rapid(z=Z_SAFE)


def feature_slot(tool):
    r = tool["dia"] / 2
    ring = Loop(offset_poly(ccw(rect_poly(SLOT)), -r))
    return dict(name="CENTRE SLOT 4 X 50", rings=[ring], entry=0, center=None)


def feature_rect(tool, rc):
    r = tool["dia"] / 2
    x0, y0, x1, y1 = rc
    ring = Loop(offset_poly(ccw(rect_poly(rc)), -r))
    ym = (y0 + y1) / 2
    return dict(name=f"WINDOW 12 X 6 AT X{f3(x0)} Y{f3(y0)}", rings=[ring], entry=0,
                center=((x0 + r, ym), (x1 - r, ym)))


def feature_tri(tool, tri, label):
    r, s = tool["dia"] / 2, tool["step"]
    tri = ccw(tri)
    rho = inradius(tri)
    ds, d = [], r
    while d < rho:
        ds.append(d)
        d += s
    rings = [Loop(offset_poly(tri, -d)) for d in ds]
    entry = max(i for i, lp in enumerate(rings) if lp.L >= 25.0)
    return dict(name=f"{label} TRIANGLE CUT-OUT", rings=rings, entry=entry, center=None)


def outline_loop(r):
    """Tool-centre path r outside the plate, CLOCKWISE (climb milling, M03), starting on the bottom edge."""
    W, H, C = PLATE_W, PLATE_H, CHAMFER
    octa = [(C, 0), (W - C, 0), (W, C), (W, H - C), (W - C, H), (C, H), (0, H - C), (0, C)]
    o = offset_poly(octa, +r)
    start = (START_X_OUTLINE, -r)
    return Loop([start, o[0], o[7], o[6], o[5], o[4], o[3], o[2], o[1]])


# =====================================================================
# program
# =====================================================================
def build():
    nc = NC()
    levels = make_levels()
    t3, t6 = TOOLS["T02"], TOOLS["T01"]

    feats = [
        feature_tri(t3, TRIS[0], "TOP"),
        feature_slot(t3),
        feature_tri(t3, TRIS[1], "BOTTOM"),
        feature_rect(t3, RECTS[0]),
        feature_rect(t3, RECTS[1]),
        feature_rect(t3, RECTS[2]),
        feature_rect(t3, RECTS[3]),
    ]

    # ---- OP1: D3 internal features -------------------------------------------------
    first = feats[0]["rings"][feats[0]["entry"]].p[0]
    nc.cmt("OP1 - INTERNAL FEATURES, PART STAYS ATTACHED TO BLANK")
    start_tool(nc, "T02", first[0], first[1], first=True)
    for ft in feats:
        nc.cmt(ft["name"])
        ring_feature(nc, t3, ft["rings"], ft["entry"], levels, ft["center"])

    # ---- OP2: D6 outline LAST -------------------------------------------------------
    ol = outline_loop(t6["dia"] / 2)
    nc.cmt("OP2 - OUTLINE LAST, CLIMB CUT CLOCKWISE")
    start_tool(nc, "T01", *ol.p[0])
    nc.cmt("OUTLINE 160 X 130 WITH 4 CORNER CUTS 15 X 15")
    approach(nc, t6, *ol.p[0])
    s, zp = 0.0, 0.0
    for z in levels:
        s = ramp_lap(nc, ol, s, zp, z, t6)
        zp = z
    nc.rapid(z=Z_SAFE)

    # ---- end ------------------------------------------------------------------------
    nc.raw("M05")
    nc.raw("G91 G28 Z0.")
    nc.raw("G28 X0. Y0.")
    nc.raw("G90")
    nc.raw("M30")
    return nc, levels


def header(nc, levels):
    t1, t2 = TOOLS["T01"], TOOLS["T02"]
    lv = " ".join(f"Z{f3(z)}" for z in levels)
    return [
        "%",
        "O0002 (CHASSIS PLATE 160 X 130 - ACRYLIC)",
        "(SOUTHWEST PETROLEUM UNIV - TWO-SHAFT FOUR-WHEEL VEHICLE)",
        f"(MATERIAL: CAST ACRYLIC {THICK:g} MM THICK - CHECK YOUR SHEET)",
        "(BLANK: MIN 180 X 150, PART CENTRED ON IT, ON SPOILBOARD)",
        "(HOLD BLANK WITH TAPE OR CLAMPS OUTSIDE 8 MM OF PLATE EDGE)",
        "(G54: X0 Y0 = LOWER LEFT CORNER OF FINISHED PLATE)",
        "(G54: Z0 = TOP FACE OF ACRYLIC)",
        f"(T01 H01 = {t1['name']}  S{t1['rpm']} F{t1['feed']})",
        f"(T02 H02 = {t2['name']}  S{t2['rpm']} F{t2['feed']})",
        f"(DEPTH PASSES: {lv}  = THICKNESS + {OVERCUT:g} INTO SPOILBOARD)",
        "(PATHS ARE TOOL-CENTRE PATHS, ALREADY OFFSET: NO G41 G42)",
        "(AIR BLAST OR MIST ONLY, NO FLOOD COOLANT, NO DWELLING)",
        f"(EST CUTTING TIME {nc.time:.1f} MIN, FEEDS ARE STARTING VALUES)",
        "(PROVE OUT: DRY RUN 20 MM ABOVE PART, THEN SINGLE BLOCK)",
        "G21 G17 G40 G49 G80 G90 G94",
        "G54",
    ]


def main():
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).with_name("chassis_plate_acrylic.nc")
    nc, levels = build()
    text = "\n".join(header(nc, levels) + nc.out + ["%"]) + "\n"
    out.write_text(text)
    print(f"wrote {out}  ({len(nc.out)} blocks, cut length {nc.cut_len:.0f} mm, ~{nc.time:.1f} min cutting)")


if __name__ == "__main__":
    main()
