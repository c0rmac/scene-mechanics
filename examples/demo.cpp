// What the signed distance function returns, on configurations whose answer is known.
//
//   ./build/demo [models/mug_spoon.xml]
//
// 1. the scene's features (the gaps it measures);
// 2. the mug on the table, lifted and pushed in: pair gap and point gap against the closed form;
// 3. the spoon standing in the mug, slid into the wall: exact ellipsoid-box gap, MuJoCo capsule-box gap, both against the
//    closed form, through contact into overlap;
// 4. random spoon poses near the mug: every exact pair routine against a brute-force minimum over the moving surface, and
//    MuJoCo's own distance on the same pairs; this also shows the routines' known failure (a deep overlap through the base);
// 5. the spoon lowered into the mug by stepping down by its smallest gap until it touches (sphere tracing), and the gaps at
//    the touching pose.
#include "scene_mechanics/signed_distance.hpp"
#include <algorithm>
#include <cstdio>
#include <map>
#include <random>
#include <string>
#include <vector>

using namespace ms3;

// world point of a point given in a geom's frame
static void world(const Frame &f, const double *l, double *p) { for (int k = 0; k < 3; ++k) { p[k] = f.c[k]; for (int j = 0; j < 3; ++j) p[k] += f.R[3 * k + j] * l[j]; } }

// brute force for an ellipsoid geom ge against geom g: the smallest exact point distance to g over a grid on the ellipsoid's
// surface, refined four times around the best grid point (equals the gap when the two are apart: the closest point of a convex
// body to a disjoint convex body is on its surface)
static double brute_ellipsoid(const SignedDistance &sd, int ge, int g) {
    const double *a = sd.m->geom_size + 3 * ge; const Frame E = sd.frame(ge);
    auto f = [&](double th, double ph) { const double l[3] = {a[0] * std::sin(th) * std::cos(ph), a[1] * std::sin(th) * std::sin(ph), a[2] * std::cos(th)}; double p[3]; world(E, l, p); return sd.point_gap(g, p, 0); };
    double best = 1e300, tb = 0, pb = 0, dt = M_PI / 200, dp = 2 * M_PI / 400;
    for (int i = 0; i <= 200; ++i) for (int j = 0; j < 400; ++j) { const double v = f(i * dt, j * dp); if (v < best) best = v, tb = i * dt, pb = j * dp; }
    for (int level = 0; level < 4; ++level) {
        const double t0 = tb, p0 = pb, wt = 2 * dt, wp = 2 * dp; dt = 2 * wt / 100; dp = 2 * wp / 100;
        for (int i = 0; i <= 100; ++i) for (int j = 0; j <= 100; ++j) { const double th = t0 - wt + i * dt, ph = p0 - wp + j * dp, v = f(th, ph); if (v < best) best = v, tb = th, pb = ph; }
    }
    return best;
}

// brute force for a capsule geom gc against geom g: the smallest exact distance from its axis segment to g, less its radius
static double brute_capsule(const SignedDistance &sd, int gc, int g) {
    const double r = sd.m->geom_size[3 * gc], h = sd.m->geom_size[3 * gc + 1]; const Frame C = sd.frame(gc);
    auto f = [&](double s) { const double l[3] = {0, 0, s}; double p[3]; world(C, l, p); return sd.point_gap(g, p, r); };
    double best = 1e300, sb = 0, ds = 2 * h / 2000;
    for (int i = 0; i <= 2000; ++i) { const double s = -h + i * ds, v = f(s); if (v < best) best = v, sb = s; }
    for (int level = 0; level < 4; ++level) {
        const double s0 = sb, w = 2 * ds; ds = 2 * w / 1000;
        for (int i = 0; i <= 1000; ++i) { const double s = std::max(-h, std::min(h, s0 - w + i * ds)), v = f(s); if (v < best) best = v, sb = s; }
    }
    return best;
}

int main(int argc, char **argv) {
    const std::string xml = argc > 1 ? argv[1] : "models/mug_spoon.xml";
    SignedDistance sd(xml, 0.10);   // pair gaps exact up to 10 cm
    const int base = sd.geom("base"), wall0 = sd.geom("wall_0"), bowl = sd.geom("bowl"), shaft = sd.geom("shaft");
    auto feature = [&](const std::string &name) { for (size_t i = 0; i < sd.F.size(); ++i) if (sd.F[i].name == name) return static_cast<int>(i); return -1; };

    // 1. features
    int npair = 0, npoint = 0; for (const Feature &f : sd.F) (f.type == Feature::GEOM_PAIR ? npair : npoint)++;
    std::printf("== 1. The scene has %zu gaps: %d between pairs of geoms, %d between a point and a geom\n", sd.F.size(), npair, npoint);
    std::printf("   for example  %s,  %s,  %s\n\n", sd.F[0].name.c_str(), sd.F[npair].name.c_str(), sd.F.back().name.c_str());

    // 2. the mug upright on the table: its base's bottom is at the mug's origin, so every base gap equals the mug's height
    std::printf("== 2. Mug upright at height h (spoon far away): table|base (pair) and base.rim0|table (point) should both equal h\n");
    std::printf("   %10s  %14s  %14s\n", "h (m)", "table|base", "base.rim0|table");
    const int rim0 = feature("base.rim0|table");
    for (double h : {1e-3, 1e-6, 0.0, -5e-4}) {
        const Pose x{0, 0, h, 1, 0, 0, 0, 0.3, 0, 0.1, 1, 0, 0, 0}; sd.set_pose(x);
        std::printf("   %+10.1e  %+14.6e  %+14.6e\n", h, sd.pair_gap(sd.table, base), sd.gap(sd.F[rim0]));
    }
    std::printf("\n");

    // 3. the spoon upright on its bowl inside the mug, slid along +x towards wall_0 (inner face at x = 0.041). Closed forms: the bowl
    //    (semi-axis 0.014 along x) meets the wall at x = 0.027, the shaft (radius 0.003) at x = 0.038, and the bowl's bottom sits on
    //    the base's top (z = 0.008) when the spoon's origin is at z = 0.031
    std::printf("== 3. Spoon standing in the mug at (x, 0, 0.031), slid into the wall\n");
    std::printf("   %8s  %14s  %14s  %14s  %14s  %14s\n", "x (m)", "bowl|wall_0", "0.027 - x", "shaft|wall_0", "0.038 - x", "bowl|base");
    for (double X : {0.0, 0.01, 0.02, 0.026, 0.027, 0.0275, 0.03}) {
        const Pose x{0, 0, 0, 1, 0, 0, 0, X, 0, 0.031, 1, 0, 0, 0}; sd.set_pose(x);
        std::printf("   %8.4f  %+14.6e  %+14.6e  %+14.6e  %+14.6e  %+14.6e\n", X, sd.pair_gap(bowl, wall0), 0.027 - X, sd.pair_gap(shaft, wall0), 0.038 - X, sd.pair_gap(bowl, base));
    }
    std::printf("   (negative = the depth of overlap; the gap keeps unit slope through contact)\n\n");

    // 4. random spoon poses near the upright mug
    std::printf("== 4. Random spoon poses near the mug: pairs whose exact gap is 0.1 mm to 1 cm, checked by brute force\n");
    struct Tally { int n = 0, apart = 0, missed = 0; double max_err = 0; int mj_bad = 0; double mj_max = 0; };
    std::map<std::string, Tally> tally;
    const std::map<int, std::string> kind{{mjGEOM_PLANE, "ellipsoid-plane"}, {mjGEOM_BOX, "ellipsoid-box"}, {mjGEOM_CYLINDER, "ellipsoid-cylinder"}, {mjGEOM_CAPSULE, "ellipsoid-capsule"}};
    std::vector<int> others{sd.table}; for (int g = 0; g < sd.m->ngeom; ++g) if (sd.m->geom_bodyid[g] == sd.body_mug) others.push_back(g);
    std::mt19937_64 rg(7); std::uniform_real_distribution<double> ux(-0.09, 0.11), uy(-0.09, 0.09), uz(-0.01, 0.13); std::normal_distribution<double> nq;
    const int CAP = 60;
    auto record = [&](const std::string &k, double exact, double brute, double mj) {
        Tally &t = tally[k]; ++t.n;
        if (brute < 0) { ++t.missed; return; }   // the bodies overlap, though the exact routine reports them apart
        ++t.apart; t.max_err = std::max(t.max_err, std::fabs(exact - brute));
        if (std::fabs(mj - brute) > 1e-6) ++t.mj_bad; t.mj_max = std::max(t.mj_max, std::fabs(mj - brute));
    };
    for (int it = 0; it < 20000; ++it) {
        double q[4], qn = 0; for (double &v : q) { v = nq(rg); qn += v * v; } qn = std::sqrt(qn);
        const Pose x{0, 0, 0, 1, 0, 0, 0, ux(rg), uy(rg), uz(rg), q[0] / qn, q[1] / qn, q[2] / qn, q[3] / qn}; sd.set_pose(x);
        for (int g : others) {
            const std::string &k = kind.at(sd.m->geom_type[g]); if (tally[k].n >= CAP) continue;
            const double e = sd.pair_gap(bowl, g); if (e < 1e-4 || e > 0.01) continue;
            record(k, e, brute_ellipsoid(sd, bowl, g), mj_geomDistance(sd.m, sd.d, bowl, g, sd.FAR, nullptr));
        }
        if (tally["capsule-cylinder"].n < CAP) {
            const double e = sd.pair_gap(shaft, base);
            if (e >= 1e-4 && e <= 0.01) record("capsule-cylinder", e, brute_capsule(sd, shaft, base), mj_geomDistance(sd.m, sd.d, shaft, base, sd.FAR, nullptr));
        }
        bool full = true; for (const char *k : {"ellipsoid-plane", "ellipsoid-box", "ellipsoid-cylinder", "ellipsoid-capsule", "capsule-cylinder"}) full = full && tally[k].n >= CAP;
        if (full) break;
    }
    std::printf("   %-20s %6s  %6s  %20s  %28s  %22s\n", "pair", "cases", "apart", "max |exact - brute|", "MuJoCo off by > 1 um (max)", "overlapping, missed");
    for (const auto &[k, t] : tally) std::printf("   %-20s %6d  %6d  %18.2e m  %8d of %-3d (%.2e m)  %22d\n", k.c_str(), t.n, t.apart, t.max_err, t.mj_bad, t.apart, t.mj_max, t.missed);
    std::printf("   apart: brute force > 0, so the exact gap must equal it. overlapping, missed: brute force < 0 (the spoon passes through the\n"
                "   8 mm base disc) but the routine reports a positive gap; see the README's known limitation\n\n");

    // 5. lower a tilted spoon into the mug: the smallest gap involving the spoon is a safe step (no point moves farther than the
    //    translation, so no gap drops by more than it); stepping down by it converges to the first touch without ever overlapping
    std::vector<int> spoon_feats;
    for (size_t i = 0; i < sd.F.size(); ++i) { const Feature &f = sd.F[i];
        const bool on_spoon = f.type == Feature::GEOM_PAIR ? (sd.m->geom_bodyid[f.g1] == sd.body_spoon || sd.m->geom_bodyid[f.g2] == sd.body_spoon) : f.body == sd.body_spoon;
        if (on_spoon) spoon_feats.push_back(static_cast<int>(i)); }
    const double th = 25 * M_PI / 180;   // spoon tilted 25 degrees about the y axis
    Pose x{0, 0, 0, 1, 0, 0, 0, -0.01, 0, 0.15, std::cos(th / 2), 0, std::sin(th / 2), 0};
    std::vector<double> g; int steps = 0; double mn = 1;
    for (; steps < 500; ++steps) { sd.gaps(x, spoon_feats, g); mn = *std::min_element(g.begin(), g.end()); if (mn < 1e-12) break; x[9] -= mn; }
    std::printf("== 5. Tilted spoon lowered into the mug by its smallest gap: %d steps, spoon origin now at z = %.9f m, smallest gap %.1e m\n", steps, x[9], mn);
    std::vector<int> order(g.size()); for (size_t k = 0; k < order.size(); ++k) order[k] = static_cast<int>(k);
    std::sort(order.begin(), order.end(), [&](int a, int b) { return g[a] < g[b]; });
    std::printf("   the six smallest of its %zu gaps:\n", g.size());
    for (int k = 0; k < 6; ++k) std::printf("     %+.6e m  %s\n", g[order[k]], sd.F[spoon_feats[order[k]]].name.c_str());
    return 0;
}
