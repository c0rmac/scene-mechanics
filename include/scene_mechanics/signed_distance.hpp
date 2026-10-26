// Signed distances between a mug, a spoon and a table (the gap functions of the scene-construction project, 26 Oct 2026).
//
// A pose x = (p_M, q_M, p_S, q_S): each body's position and unit quaternion (w, x, y, z). A gap is the signed distance between
// two pieces of geometry: the distance when they are apart, 0 on contact, minus the penetration depth when they overlap.
// It is continuous with unit slope across contact. The scene has two kinds of gap (its features):
//  - geom-pair gaps: the table to every mug and spoon geom, every mug geom to both spoon geoms. They see any contact, curved
//    ones included, but a single number per pair cannot express a face lying flush;
//  - point gaps, for face and line contacts: the 8 corners of every wall box and 8 points of the base's bottom rim against
//    the table; both ends of every handle capsule against the table; both ends of the spoon shaft (spheres of the shaft's
//    radius) against the table, every wall box, the base and every handle capsule. A flush face is several point gaps at 0.
// Pair gaps come from exact formulas (geom3d.hpp) where MuJoCo's convex solver is unreliable (the ellipsoid against anything,
// the capsule against the cylinder, the cylinder against the plane) and from MuJoCo's analytic distance elsewhere.
#pragma once
#include "geom3d.hpp"
#include <mujoco/mujoco.h>
#include <array>
#include <stdexcept>
#include <string>
#include <vector>

namespace ms3 {

using Pose = std::array<double, 14>;   // mug (px py pz qw qx qy qz), spoon (px py pz qw qx qy qz)

struct Feature {
    enum Type { GEOM_PAIR, POINT } type;
    int g1 = -1, g2 = -1;          // GEOM_PAIR
    int body = -1;                 // POINT: body carrying the point
    double local[3] = {0, 0, 0};   // POINT: point in the body frame
    double radius = 0;             // POINT: sphere radius (capsule end), 0 for a corner
    int target = -1;               // POINT: geom measured against (another body, or the table)
    std::string name;
    int geom = -1;                 // POINT: geom the point belongs to
};

inline void quat_rot(const double *q, const double *v, double *o) {   // rotate v by unit quaternion q
    const double w = q[0], x = q[1], y = q[2], z = q[3];
    const double R[9] = {1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y), 2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x), 2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)};
    for (int i = 0; i < 3; ++i) o[i] = R[3 * i] * v[0] + R[3 * i + 1] * v[1] + R[3 * i + 2] * v[2];
}

struct SignedDistance {
    mjModel *m = nullptr;
    mjData *d = nullptr;
    int body_mug = -1, body_spoon = -1, table = -1;
    std::vector<Feature> F;
    double FAR = 0.01;   // pair gaps beyond FAR are reported as FAR (only contacts matter; point gaps are not capped)

    explicit SignedDistance(const std::string &xml, double far = 0.01) : FAR(far) {
        char err[1000] = "";
        m = mj_loadXML(xml.c_str(), nullptr, err, sizeof err);
        if (!m) throw std::runtime_error(std::string("mj_loadXML: ") + err);
        d = mj_makeData(m);
        body_mug = mj_name2id(m, mjOBJ_BODY, "mug"); body_spoon = mj_name2id(m, mjOBJ_BODY, "spoon");
        table = mj_name2id(m, mjOBJ_GEOM, "table");
        std::vector<int> mug, spoon;
        for (int g = 0; g < m->ngeom; ++g) { if (m->geom_bodyid[g] == body_mug) mug.push_back(g); else if (m->geom_bodyid[g] == body_spoon) spoon.push_back(g); }
        auto gname = [&](int g) { return std::string(mj_id2name(m, mjOBJ_GEOM, g)); };
        for (int g : mug) F.push_back({Feature::GEOM_PAIR, table, g, -1, {0, 0, 0}, 0, -1, "table|" + gname(g)});
        for (int g : spoon) F.push_back({Feature::GEOM_PAIR, table, g, -1, {0, 0, 0}, 0, -1, "table|" + gname(g)});
        for (int g : mug) for (int s : spoon) F.push_back({Feature::GEOM_PAIR, g, s, -1, {0, 0, 0}, 0, -1, gname(g) + "|" + gname(s)});
        // a point given in a geom's frame, expressed in its body's frame
        auto body_point = [&](int g, const double *pl, double *out) { double r[3]; quat_rot(m->geom_quat + 4 * g, pl, r); for (int k = 0; k < 3; ++k) out[k] = m->geom_pos[3 * g + k] + r[k]; };
        auto add_point = [&](int g, const double *pl, double radius, int target, const std::string &nm) {
            Feature f{Feature::POINT, -1, -1, m->geom_bodyid[g], {0, 0, 0}, radius, target, nm, g}; body_point(g, pl, f.local); F.push_back(f);
        };
        for (int g : mug) {
            const double *sz = m->geom_size + 3 * g;
            if (m->geom_type[g] == mjGEOM_BOX) {
                for (int c = 0; c < 8; ++c) { const double pl[3] = {(c & 1 ? 1 : -1) * sz[0], (c & 2 ? 1 : -1) * sz[1], (c & 4 ? 1 : -1) * sz[2]}; add_point(g, pl, 0, table, gname(g) + ".corner" + std::to_string(c) + "|table"); }
            } else if (m->geom_type[g] == mjGEOM_CYLINDER) {
                for (int a = 0; a < 8; ++a) { const double th = 2 * M_PI * a / 8, pl[3] = {sz[0] * std::cos(th), sz[0] * std::sin(th), -sz[1]}; add_point(g, pl, 0, table, gname(g) + ".rim" + std::to_string(a) + "|table"); }
            } else if (m->geom_type[g] == mjGEOM_CAPSULE) {
                for (int e = 0; e < 2; ++e) { const double pl[3] = {0, 0, (e ? 1 : -1) * sz[1]}; add_point(g, pl, sz[0], table, gname(g) + ".end" + std::to_string(e) + "|table"); }
            }
        }
        for (int s : spoon) {
            if (m->geom_type[s] != mjGEOM_CAPSULE) continue;
            const double *sz = m->geom_size + 3 * s;
            for (int e = 0; e < 2; ++e) {
                const double pl[3] = {0, 0, (e ? 1 : -1) * sz[1]};
                add_point(s, pl, sz[0], table, gname(s) + ".end" + std::to_string(e) + "|table");
                for (int g : mug) add_point(s, pl, sz[0], g, gname(s) + ".end" + std::to_string(e) + "|" + gname(g));
            }
        }
    }
    ~SignedDistance() { if (d) mj_deleteData(d); if (m) mj_deleteModel(m); }
    SignedDistance(const SignedDistance &) = delete;
    SignedDistance &operator=(const SignedDistance &) = delete;

    int geom(const char *name) const { return mj_name2id(m, mjOBJ_GEOM, name); }

    // place both bodies at pose x (geom frames in d->geom_xpos, d->geom_xmat)
    void set_pose(const Pose &x) const { for (int i = 0; i < 14; ++i) d->qpos[i] = x[i]; mj_kinematics(m, d); }

    // signed distance from world point p to geom g (negative inside), less radius r
    double point_gap(int g, const double *p, double r) const {
        const double *c = d->geom_xpos + 3 * g, *R = d->geom_xmat + 9 * g, *sz = m->geom_size + 3 * g;
        double q[3]; for (int k = 0; k < 3; ++k) { q[k] = 0; for (int j = 0; j < 3; ++j) q[k] += R[3 * j + k] * (p[j] - c[j]); }
        switch (m->geom_type[g]) {
        case mjGEOM_PLANE: return q[2] - r;
        case mjGEOM_BOX: {
            double out = 0, in = -1e300; for (int k = 0; k < 3; ++k) { const double dk = std::fabs(q[k]) - sz[k]; out += std::max(dk, 0.0) * std::max(dk, 0.0); in = std::max(in, dk); }
            return (in > 0 ? std::sqrt(out) : in) - r;
        }
        case mjGEOM_CYLINDER: {
            const double dr = std::sqrt(q[0] * q[0] + q[1] * q[1]) - sz[0], dz = std::fabs(q[2]) - sz[1];
            return std::min(std::max(dr, dz), 0.0) + std::sqrt(std::max(dr, 0.0) * std::max(dr, 0.0) + std::max(dz, 0.0) * std::max(dz, 0.0)) - r;
        }
        case mjGEOM_CAPSULE: {
            const double t = std::max(-sz[1], std::min(sz[1], q[2]));
            return std::sqrt(q[0] * q[0] + q[1] * q[1] + (q[2] - t) * (q[2] - t)) - sz[0] - r;
        }
        default: return 1e300;
        }
    }

    Frame frame(int g) const { Frame f; for (int k = 0; k < 3; ++k) f.c[k] = d->geom_xpos[3 * g + k]; for (int k = 0; k < 9; ++k) f.R[k] = d->geom_xmat[9 * g + k]; return f; }
    double bound_radius(int g) const {
        const double *sz = m->geom_size + 3 * g;
        switch (m->geom_type[g]) { case mjGEOM_ELLIPSOID: return std::max(sz[0], std::max(sz[1], sz[2])); case mjGEOM_BOX: return std::sqrt(sz[0] * sz[0] + sz[1] * sz[1] + sz[2] * sz[2]);
            case mjGEOM_CYLINDER: return std::sqrt(sz[0] * sz[0] + sz[1] * sz[1]); case mjGEOM_CAPSULE: return sz[0] + sz[1]; default: return 1e300; }
    }
    // gap between two geoms: exact formulas (geom3d.hpp) where MuJoCo's convex solver fails, MuJoCo's analytic distance elsewhere
    // cut: when a cheap lower bound (separating axes) exceeds cut, the bound is returned instead of the exact gap; callers that
    // only need the gaps below a threshold (the active set, penetration) pass it
    double pair_gap(int g1, int g2, double cut = 1e300) const {
        int t1 = m->geom_type[g1], t2 = m->geom_type[g2];
        if (t2 == mjGEOM_ELLIPSOID || (t1 == mjGEOM_CYLINDER && t2 == mjGEOM_CAPSULE)) { std::swap(g1, g2); std::swap(t1, t2); }   // g1: the ellipsoid, or the capsule against a cylinder
        if ((t1 == mjGEOM_PLANE && t2 == mjGEOM_CYLINDER) || (t1 == mjGEOM_CYLINDER && t2 == mjGEOM_PLANE)) {   // MuJoCo's routine loses
            const int gp = t1 == mjGEOM_PLANE ? g1 : g2, gc = t1 == mjGEOM_PLANE ? g2 : g1;                       // flush discs: the lowest
            const double *Rp = d->geom_xmat + 9 * gp, *Rc = d->geom_xmat + 9 * gc, *cp = d->geom_xpos + 3 * gp, *cc = d->geom_xpos + 3 * gc;   // rim point
            const double n[3] = {Rp[2], Rp[5], Rp[8]}, a[3] = {Rc[2], Rc[5], Rc[8]}, R = m->geom_size[3 * gc], h = m->geom_size[3 * gc + 1];
            const double na[3] = {n[1] * a[2] - n[2] * a[1], n[2] * a[0] - n[0] * a[2], n[0] * a[1] - n[1] * a[0]};
            const double sn = std::sqrt(na[0] * na[0] + na[1] * na[1] + na[2] * na[2]), na_ = n[0] * a[0] + n[1] * a[1] + n[2] * a[2];
            const double hc = n[0] * (cc[0] - cp[0]) + n[1] * (cc[1] - cp[1]) + n[2] * (cc[2] - cp[2]);
            return std::min(FAR, hc - h * std::fabs(na_) - R * sn);   // lower end centre at hc - h |n.a|, rim reaching R sin(angle) below it
        }
        const bool custom = t1 == mjGEOM_ELLIPSOID || (t1 == mjGEOM_CAPSULE && t2 == mjGEOM_CYLINDER);
        if (!custom) return mj_geomDistance(m, d, g1, g2, FAR, nullptr);
        if (t2 != mjGEOM_PLANE) {
            double dc = 0; for (int k = 0; k < 3; ++k) dc += (d->geom_xpos[3 * g1 + k] - d->geom_xpos[3 * g2 + k]) * (d->geom_xpos[3 * g1 + k] - d->geom_xpos[3 * g2 + k]);
            if (std::sqrt(dc) - bound_radius(g1) - bound_radius(g2) > FAR) return FAR;
        }
        const Frame F1 = frame(g1), F2 = frame(g2); const double *s1 = m->geom_size + 3 * g1, *s2 = m->geom_size + 3 * g2;
        double v;
        if (t1 == mjGEOM_ELLIPSOID) {
            if (t2 == mjGEOM_PLANE) { const double nz[3] = {F2.R[2], F2.R[5], F2.R[8]}; double off = 0; for (int k = 0; k < 3; ++k) off += nz[k] * (F1.c[k] - F2.c[k]); double mneg[3] = {-nz[0], -nz[1], -nz[2]}; v = off - ell_support(s1, F1, mneg, nullptr); }
            else if (t2 == mjGEOM_BOX) {
                const Frame E = relative(F2, F1); double lb = -1e300;
                for (int k = 0; k < 3; ++k) { double mk[3] = {0, 0, 0}; mk[k] = 1; lb = std::max(lb, std::fabs(E.c[k]) - s2[k] - ell_support(s1, E, mk, nullptr)); }   // separating box axes
                if (lb > cut) return std::min(lb, FAR);
                v = gap_ellipsoid_box(s1, E, s2);
            } else if (t2 == mjGEOM_CYLINDER) {
                const Frame E = relative(F2, F1); const double mz[3] = {0, 0, 1}, rho = std::sqrt(E.c[0] * E.c[0] + E.c[1] * E.c[1]);
                double lb = std::fabs(E.c[2]) - s2[1] - ell_support(s1, E, mz, nullptr);   // along the axis, and radially
                if (rho > 0) { const double mr[3] = {E.c[0] / rho, E.c[1] / rho, 0}; lb = std::max(lb, rho - s2[0] - ell_support(s1, E, mr, nullptr)); }
                if (lb > cut) return std::min(lb, FAR);
                v = gap_ellipsoid_cylinder(s1, E, s2[0], s2[1]);
            } else if (t2 == mjGEOM_CAPSULE) {
                double P[3], Q[3], A[3], B[3]; for (int k = 0; k < 3; ++k) { P[k] = F2.c[k] - s2[1] * F2.R[3 * k + 2]; Q[k] = F2.c[k] + s2[1] * F2.R[3 * k + 2]; }
                const double lb = point_segment(F1.c, P, Q) - std::max(s1[0], std::max(s1[1], s1[2])) - s2[0];
                if (lb > cut) return std::min(lb, FAR);
                to_local(F1, P, A); to_local(F1, Q, B); v = seg_min_ellipsoid(s1, A, B) - s2[0];
            }
            else return mj_geomDistance(m, d, g1, g2, FAR, nullptr);
        } else { double P[3], Q[3], A[3], B[3]; for (int k = 0; k < 3; ++k) { P[k] = F1.c[k] - s1[1] * F1.R[3 * k + 2]; Q[k] = F1.c[k] + s1[1] * F1.R[3 * k + 2]; } to_local(F2, P, A); to_local(F2, Q, B); v = gap_capsule_cylinder(A, B, s1[0], s2[0], s2[1]); }
        return std::min(v, FAR);
    }

    // gap of one feature at the pose last placed by set_pose
    double gap(const Feature &f, double cut = 1e300) const {
        if (f.type == Feature::GEOM_PAIR) return pair_gap(f.g1, f.g2, cut);
        double p[3]; const double *xp = d->xpos + 3 * f.body, *xm = d->xmat + 9 * f.body;
        for (int k = 0; k < 3; ++k) p[k] = xp[k] + xm[3 * k] * f.local[0] + xm[3 * k + 1] * f.local[1] + xm[3 * k + 2] * f.local[2];
        return point_gap(f.target, p, f.radius);
    }

    // gaps of the features in idx (all if empty) at pose x
    void gaps(const Pose &x, const std::vector<int> &idx, std::vector<double> &out, double cut = 1e300) const {
        set_pose(x);
        const int n = idx.empty() ? static_cast<int>(F.size()) : static_cast<int>(idx.size());
        out.resize(n);
        for (int a = 0; a < n; ++a) out[a] = gap(F[idx.empty() ? a : idx[a]], cut);
    }
};

}  // namespace ms3
