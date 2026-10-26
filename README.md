# scene-mechanics

Signed distances for a mug and a spoon on a table. Given the two bodies' poses, `SignedDistance` returns the scene's **gaps**:
the signed distance between two pieces of geometry.

- positive: the distance when they are apart;
- 0: on contact;
- negative: minus the penetration depth when they overlap.

A gap is continuous with slope 1 across contact. The geometry is the MuJoCo model in [`models/mug_spoon.xml`](models/mug_spoon.xml):

- the table is a plane;
- the mug is a cylinder base, 24 box wall segments and 12 capsule handle segments;
- the spoon is an ellipsoid bowl and a capsule shaft.

## Build and run the demo

You need a C++20 compiler and MuJoCo 3.x (headers and shared library). The pip package provides both:

```bash
pip install mujoco
```

```bash
./build.sh
```

```bash
./build/demo
```

`build.sh` finds MuJoCo through the Python package. To use another install, set `MUJOCO_DIR` to a directory holding
`include/mujoco/mujoco.h` and `libmujoco`:

```bash
MUJOCO_DIR=/path/to/mujoco ./build.sh
```

## How to call it

The library is header-only: include [`signed_distance.hpp`](include/scene_mechanics/signed_distance.hpp) (it includes
[`geom3d.hpp`](include/scene_mechanics/geom3d.hpp)) and link MuJoCo.

```cpp
#include "scene_mechanics/signed_distance.hpp"
using namespace ms3;

SignedDistance sd("models/mug_spoon.xml");   // optional 2nd argument: FAR, the cap on pair gaps (default 0.01 m)

// a pose: mug (px py pz qw qx qy qz), spoon (px py pz qw qx qy qz); metres, unit quaternions (w first)
Pose x{0, 0, 0, 1, 0, 0, 0,   0.02, 0, 0.031, 1, 0, 0, 0};   // mug upright on the table, spoon standing in it

// 1. every gap of the scene at x: g[i] is the gap of feature sd.F[i] (its name in sd.F[i].name)
std::vector<double> g;
sd.gaps(x, {}, g);                    // 413 values

// 2. a subset of features, by index into sd.F
sd.gaps(x, {0, 7, 120}, g);

// 3. one pair of geoms, or one point against a geom, after placing the bodies
sd.set_pose(x);
double bowl_wall = sd.pair_gap(sd.geom("bowl"), sd.geom("wall_0"));   // 0.007 here
const double p[3] = {0.0, 0.0, 0.02};                                 // a world point
double into_base = sd.point_gap(sd.geom("base"), p, 0.0);             // signed distance of p to the base cylinder
double one = sd.gap(sd.F[5]);                                         // one feature at the placed pose
```

| call | returns |
|---|---|
| `SignedDistance(xml, far = 0.01)` | loads the model and builds the feature list `F` |
| `set_pose(x)` | places both bodies (MuJoCo forward kinematics) |
| `gaps(x, idx, out, cut = 1e300)` | gaps of the features in `idx` (all if empty) at `x`, into `out` |
| `gap(feature, cut = 1e300)` | gap of one feature at the pose last placed |
| `pair_gap(g1, g2, cut = 1e300)` | signed distance between geoms `g1` and `g2` (ids from `geom(name)`), capped at `FAR` |
| `point_gap(g, p, r)` | signed distance from world point `p` to geom `g`, less radius `r` (a sphere of radius `r` at `p`) |
| `geom(name)` | geom id of a named geom |

Notes:

- **The cap.** Pair gaps larger than `FAR` are returned as `FAR`, because only gaps near contact matter for contacts. Pass a larger `far` to the constructor for true distances (the demo uses 0.1 m). Point gaps are not capped.
- **`cut`.** When a cheap separating-axis lower bound on a pair gap already exceeds `cut`, that bound is returned instead of the exact value. Callers that only need gaps below a threshold pass it.
- **Threads.** One `SignedDistance` holds one `mjData` and is not thread-safe. Make one per thread.

## The features

There are 413 gaps. The pair gaps alone see every contact, including the curved ones. The point gaps are needed because a single
number per pair cannot express a face lying flush: a flush face is several point gaps at 0.

| kind | features | count |
|---|---|---|
| pair | table to every mug and spoon geom; every mug geom to both spoon geoms | 113 |
| point | the 8 corners of every wall box, and 8 points of the base's bottom rim, against the table | 200 |
| point | both ends of every handle capsule against the table | 24 |
| point | both ends of the spoon shaft (spheres of the shaft's radius) against the table and every mug geom | 76 |

Each pair is computed by one of these methods:

| pair | method |
|---|---|
| point against plane, box, cylinder or capsule | closed-form signed distance |
| plane and cylinder (table and base) | closed form: the height of the lowest rim point |
| ellipsoid and plane | support function of the ellipsoid |
| ellipsoid and box, cylinder or capsule | point-to-ellipsoid distance (safeguarded Newton on the Lagrange equation), minimised over the other body's faces, edges, axis and rim circles |
| capsule and cylinder | the segment's distance to the cylinder's ends, rims and lateral surface, less the radius |
| any other pair (box, capsule, plane combinations) | MuJoCo's `mj_geomDistance` (analytic for these types) |

The exact routines exist because `mj_geomDistance` is unreliable on curved pairs. In section 4 of the demo it returns 0 for an
ellipsoid–box pair that is 9.6 mm apart.

## What the demo shows

Each section compares the function with an answer known independently:

1. The feature list.
2. **Mug at height h.** The pair gap `table|base` and the point gap `base.rim0|table` both equal h, for h = 1 mm, 1 µm, 0 and −0.5 mm.
3. **Spoon slid into the wall.** The spoon stands on its bowl inside the mug and slides into the wall. `bowl|wall_0` follows 0.027 − x, `shaft|wall_0` follows 0.038 − x and `bowl|base` stays 0. The match is to 1e-18 m, through contact and into 3 mm of overlap.
4. **Random spoon poses.** Each exact pair routine is compared with a brute-force minimum over the moving surface (a fine grid, refined four times). Where the bodies are apart, the two agree to 5e-17 m for every pair type. MuJoCo's distance is off by 9.6 mm in one ellipsoid–box case.
5. **Lowering a tilted spoon.** A tilted spoon is lowered into the mug by stepping down by its smallest gap each time. Because a gap is a true distance, this never overlaps and stops at first touch: 7 steps, and the bowl touches the base with gap 3e-18 m.

The whole demo runs in under a second.

## Known limitation

The ellipsoid–cylinder and capsule–cylinder routines are built for bodies that are apart or overlap slightly. They miss a deep
overlap through the base disc, which is only 8 mm thick.

- **Ellipsoid–cylinder:** the bowl's centre lies inside the disc's slab.
- **Capsule–cylinder:** the shaft pierces the disc near its axis.

In both cases the routine returns the positive distance to the rim instead of a negative value. In section 4 of the demo this
happens in 13 of 60 capsule–cylinder cases and 4 of 60 ellipsoid–cylinder cases. These are reported in the "overlapping, missed"
column. MuJoCo's collision pass flags all 17 of these poses as penetrating.
