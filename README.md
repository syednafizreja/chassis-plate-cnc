# Chassis Plate CNC Program (Acrylic)

Fanuc-style G-code, a Python generator and a material-removal simulator for the **160 × 130 mm chassis / mounting plate** of a two-shaft, four-wheel training vehicle.

Part of the *Engineering Training Vehicle Design & Racing Project* at Southwest Petroleum University.

<p align="center">
  <img src="docs/code_screenshot.png" width="46%" alt="Syntax-highlighted excerpt of chassis_plate_acrylic.nc">
  <img src="docs/preview.png" width="52%" alt="Simulated toolpaths and resulting plate">
</p>

## What is in this repo

| File | Purpose |
|---|---|
| `chassis_plate_acrylic.nc` | The CNC mill program (Fanuc-style ISO G-code, 289 lines) |
| `generate_chassis_plate.py` | Generates the program from the part geometry and cutting parameters (standard library only) |
| `simulate_nc.py` | Reads the `.nc` file, simulates material removal on a height map and checks the result against the drawing |
| `docs/` | Code screenshot and simulation preview |

## The part

Measured from the technical drawing. Origin `X0 Y0` is the lower-left corner of the finished plate.

| Feature | Geometry |
|---|---|
| Outline | 160 × 130 mm, four corner cuts of 15 × 15 mm |
| Centre alignment slot | 4 × 50 mm, X78–82, Y40–90 |
| Side windows (4×) | 12 × 6 mm at X5–17 and X143–155, Y55–61 and Y70–76 |
| Lightening cut-outs (2×) | Triangles 40 wide × 15 high: (60,110) (100,110) (80,95) and (60,20) (100,20) (80,35) |
| Material | Cast acrylic, 4 mm assumed (drawing allows 3–5 mm) |

> The drawing labels the corners "R5", but they are drawn as straight 15 × 15 mm chamfers. The program follows the drawn geometry. Change `CHAMFER` in the generator if your design differs.

## Machining strategy

| Op | Tool | Cuts |
|---|---|---|
| 1 | T02, Ø3 mm single-flute carbide | Two triangles, centre slot, four side windows |
| 2 | T01, Ø6 mm single-flute carbide | Outline, cut last so the part stays attached while internal features are machined |

- **Tool-centre paths.** Every path is already offset by the cutter radius, so no `G41`/`G42` and no radius wear offset is needed.
- **Two depth passes** of 2.1 mm each, to Z-4.2 (4 mm sheet + 0.2 mm into the spoilboard).
- **Ramp entry at 4°** on every feature. No straight plunging into acrylic.
- **Climb milling.** Outline runs clockwise, internal features counter-clockwise (spindle `M03`).
- **Triangles are cleared inside-out** in 50 % stepover rings, so no loose slug is left behind. The finishing ring is always last.
- **Decimal points on every value**, because many Fanuc controls read `X78` as 0.078 mm.
- Acrylic settings: air blast or mist only, no flood coolant, no dwelling in the cut.

### Cutting data (starting values, prove out on scrap)

| Tool | Spindle | Feed | Ramp feed |
|---|---|---|---|
| T02 Ø3 | 8000 rpm | 350 mm/min | 200 mm/min |
| T01 Ø6 | 8000 rpm | 600 mm/min | 250 mm/min |

## Setup

- **Work offset G54:** `X0 Y0` at the lower-left corner of the finished plate, `Z0` at the top face of the acrylic.
- **Blank:** at least 180 × 150 mm, with the plate centred on it, on a spoilboard.
- **Holding:** tape or clamps in the frame only. The Ø6 cutter reaches 6 mm outside the plate edge, so keep clamps at least 8 mm away.
- **Tools:** T01 / H01 = Ø6, T02 / H02 = Ø3.

## Usage

```bash
pip install -r requirements.txt          # only the simulator needs numpy and matplotlib

python3 generate_chassis_plate.py        # writes chassis_plate_acrylic.nc
python3 simulate_nc.py chassis_plate_acrylic.nc docs/preview.png
```

Part size, thickness, feature positions, feeds, speeds, depth per pass and ramp angle are constants at the top of `generate_chassis_plate.py`.

## Verification

The simulator reads only the `.nc` file, cuts a 0.05 mm height map of the blank with flat end mills, and compares the result with the drawing geometry.

| Check | Result |
|---|---|
| Finished plate surface gouged anywhere | 0.000 mm² |
| Plate material cut through that should stay | 0.000 mm² |
| Scrap within 2.9 mm of the outline left uncut | 0.000 mm² |
| Rapid XY moves below Z1 | 0 |
| Deepest Z reached | -4.200 mm |
| Max depth of cut, Ø3 and Ø6 | 2.10 mm |
| Left in feature corners (Ø3 corner radius) | 18.2 mm², theory 17.6 mm² |
| Estimated cutting time | about 7.5 min |

**This is a simulation, not a machine proof-out.** Run a dry run 20 mm above the part, then single block, before cutting. Feeds and speeds depend on your machine, cutter and sheet.

## Adapting to another controller

The program is Fanuc-style (`G43 Hxx`, `T## M06`, `G91 G28 Z0.`). Check tool change, home return and dwell conventions for Siemens, Mach3, LinuxCNC and similar controllers.

## License

Add a license of your choice (GitHub: *Add file → Create new file → name it `LICENSE` → Choose a license template*).
