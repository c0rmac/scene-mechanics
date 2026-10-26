// Exact signed gaps for the geom pairs where MuJoCo's convex solver is unreliable (21 Oct 2026): ellipsoid against box,
// cylinder and capsule, and capsule against cylinder. On all of Majid's poses mj_geomDistance (native solver) returned 0 for
// 1,749 pairs whose closest points were centimetres apart, and libccd disagreed by up to 2.6 cm; the analytic pairs (plane, box
// and capsule combinations) had no such failure.
// A gap g(P, Q) is the distance when the bodies are apart, 0 on contact, and minus the depth of the deepest point of one inside
// the other when they overlap slightly: continuous, with unit slope across contact, which is what the strata machinery needs.
// Everything rests on the signed distance from a point to an ellipsoid (closest point by safeguarded Newton on the
// Lagrange-multiplier equation), minimised over the features of the other body: faces through the ellipsoid's support point,
// edges and the axis as convex one-dimensional problems, rim circles by sampling and refinement.
#pragma once
#include <algorithm>
#include <cmath>
#include <array>
#include <algorithm>

namespace ms3 {

// GEOM_FAST (23 Oct): the ellipsoid distance's Newton solve starts from the previous solution in a chain of nearby queries (along a
// rim circle or a segment), and the angle / segment-parameter refinements stop at 1e-9 rad / 1e-11 instead of 1e-14 / 1e-15. The
// root is the same (unique on its branch, the safeguarded Newton converges from any start in the bracket); the minimum values move
// by far less than 1e-12 m. Off by default (tools keep the original path); pcbo3d switches it on (env GEOM_FAST, default 1).
inline bool GEOM_FAST = false;

struct Frame { double c[3]; double R[9]; };   // world = c + R local, R row-major (columns: the local axes in world)

inline void to_local(const Frame &f, const double *y, double *z) { for (int k = 0; k < 3; ++k) { z[k] = 0; for (int j = 0; j < 3; ++j) z[k] += f.R[3 * j + k] * (y[j] - f.c[j]); } }
inline void dir_local(const Frame &f, const double *v, double *z) { for (int k = 0; k < 3; ++k) { z[k] = 0; for (int j = 0; j < 3; ++j) z[k] += f.R[3 * j + k] * v[j]; } }
// frame of b expressed in the frame of a
inline Frame relative(const Frame &a, const Frame &b) { Frame r; to_local(a, b.c, r.c); for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) { double v = 0; for (int k = 0; k < 3; ++k) v += a.R[3 * k + i] * b.R[3 * k + j]; r.R[3 * i + j] = v; } return r; }

// signed distance from y (ellipsoid frame) to the ellipsoid with semi-axes a; n = outward unit normal at the closest point
inline double sdf_ellipsoid(const double *a, const double *y, double *n = nullptr, double *tws = nullptr) {   // tws: warm start in/out (GEOM_FAST)
    double s = 0; for (int i = 0; i < 3; ++i) s += (y[i] / a[i]) * (y[i] / a[i]);
    int im = 0; for (int i = 1; i < 3; ++i) if (a[i] < a[im]) im = i;
    const double amin2 = a[im] * a[im];
    auto F = [&](double t) { double v = -1; for (int i = 0; i < 3; ++i) { const double q = a[i] * y[i] / (t + a[i] * a[i]); v += q * q; } return v; };
    auto dF = [&](double t) { double v = 0; for (int i = 0; i < 3; ++i) { const double u = t + a[i] * a[i]; v -= 2 * (a[i] * y[i]) * (a[i] * y[i]) / (u * u * u); } return v; };
    double t;
    if (s == 1.0) { t = 0; if (tws) *tws = NAN; }
    else {
        double lo, hi;
        if (s > 1) { lo = 0; double ny = 0; for (int i = 0; i < 3; ++i) ny += y[i] * y[i]; hi = std::max(a[0], std::max(a[1], a[2])) * std::sqrt(ny); }
        else { hi = 0; lo = -amin2 * (1 - 1e-12); if (F(lo) <= 0) lo = -amin2; }
        if (s < 1 && F(-amin2 * (1 - 1e-12)) <= 0) {   // deep inside, closest point off the root branch (Eberly's degenerate case)
            double x[3], r = 0; for (int i = 0; i < 3; ++i) if (i != im) { x[i] = a[i] * a[i] * y[i] / (a[i] * a[i] - amin2); r += (x[i] / a[i]) * (x[i] / a[i]); }
            x[im] = a[im] * std::sqrt(std::max(0.0, 1 - r)) * (y[im] < 0 ? -1 : 1);
            double d2 = 0, g[3], gn = 0; for (int i = 0; i < 3; ++i) { d2 += (y[i] - x[i]) * (y[i] - x[i]); g[i] = x[i] / (a[i] * a[i]); gn += g[i] * g[i]; }
            if (n) for (int i = 0; i < 3; ++i) n[i] = g[i] / std::sqrt(gn);
            if (tws) *tws = NAN;
            return -std::sqrt(d2);
        }
        t = s > 1 ? 0.5 * hi : 0.5 * lo;
        if (GEOM_FAST && tws && std::isfinite(*tws) && *tws > lo && *tws < hi) t = *tws;   // warm start on the same branch
        for (int it = 0; it < 200; ++it) {   // safeguarded Newton on F(t) = 0, F decreasing on (-amin2, inf)
            const double f = F(t);
            if (f > 0) lo = t; else hi = t;
            if (f == 0 || hi - lo <= 1e-17 * std::max(1.0, std::fabs(t)) * amin2) break;
            double tn = t - f / dF(t);
            if (!(tn > lo && tn < hi)) tn = 0.5 * (lo + hi);
            if (tn == t) break;
            t = tn;
        }
        if (tws) *tws = t;
    }
    double d2 = 0, g[3], gn = 0;
    for (int i = 0; i < 3; ++i) { const double u = t + a[i] * a[i], dx = y[i] * t / u; d2 += dx * dx; g[i] = y[i] / u; gn += g[i] * g[i]; }   // y - x = y t / (t + a^2); normal along x / a^2 = y / (t + a^2)
    if (n) { gn = std::sqrt(gn); for (int i = 0; i < 3; ++i) n[i] = gn > 0 ? g[i] / gn : (i == im ? 1.0 : 0.0); }
    return (s < 1 ? -1 : 1) * std::sqrt(d2);
}

// minimum over t in [0, 1] of the ellipsoid's signed distance at A + t (B - A) (A, B in the ellipsoid frame): convex in t, so
// its derivative n(x(t)) . (B - A) is increasing; its root by regula falsi (Illinois) safeguarded with bisection
inline double seg_min_ellipsoid(const double *a, const double *A, const double *B, double *tstar = nullptr) {
    double d[3] = {B[0] - A[0], B[1] - A[1], B[2] - A[2]};
    double tw = NAN;
    auto at = [&](double t, double *dv) { double y[3] = {A[0] + t * d[0], A[1] + t * d[1], A[2] + t * d[2]}, n[3]; const double v = sdf_ellipsoid(a, y, n, GEOM_FAST ? &tw : nullptr); *dv = n[0] * d[0] + n[1] * d[1] + n[2] * d[2]; return v; };
    double d0, d1; const double f0 = at(0, &d0), f1 = at(1, &d1);
    if (d0 >= 0) { if (tstar) *tstar = 0; return f0; }
    if (d1 <= 0) { if (tstar) *tstar = 1; return f1; }
    double lo = 0, hi = 1, dlo = d0, dhi = d1, best = std::min(f0, f1), tb = f0 < f1 ? 0 : 1; int side = 0;
    for (int it = 0; it < 80 && hi - lo > (GEOM_FAST ? 1e-11 : 1e-15); ++it) {
        double m = lo - dlo * (hi - lo) / (dhi - dlo);
        if (!(m > lo && m < hi) || it % 8 == 7) m = 0.5 * (lo + hi);   // every eighth step a bisection guarantees progress
        double dm; const double fm = at(m, &dm);
        if (fm < best) best = fm, tb = m;
        if (dm == 0) break;
        if (dm > 0) { hi = m; dhi = dm; if (side == 1) dlo *= 0.5; side = 1; } else { lo = m; dlo = dm; if (side == -1) dhi *= 0.5; side = -1; }
    }
    if (tstar) *tstar = tb;
    return best;
}

// minimum over the circle (radius R, height z, axis = local z) of a function f(point): sampled, then refined by golden section
template <class Fn> inline double circle_min(double R, double z, Fn f) {
    const int N = 96; double best = 1e300; int kb = 0;
    for (int k = 0; k < N; ++k) { const double th = 2 * M_PI * k / N, p[3] = {R * std::cos(th), R * std::sin(th), z}; const double v = f(p); if (v < best) best = v, kb = k; }
    double lo = 2 * M_PI * (kb - 1) / N, hi = 2 * M_PI * (kb + 1) / N; const double gr = 0.5 * (std::sqrt(5.0) - 1);
    auto g = [&](double th) { const double p[3] = {R * std::cos(th), R * std::sin(th), z}; return f(p); };
    double x1 = hi - gr * (hi - lo), x2 = lo + gr * (hi - lo), f1 = g(x1), f2 = g(x2);
    for (int it = 0; it < 80 && hi - lo > (GEOM_FAST ? 1e-9 : 1e-14); ++it) { if (f1 < f2) { hi = x2; x2 = x1; f2 = f1; x1 = hi - gr * (hi - lo); f1 = g(x1); } else { lo = x1; x1 = x2; f1 = f2; x2 = lo + gr * (hi - lo); f2 = g(x2); } }
    return std::min(best, std::min(f1, f2));
}

// GEOM_FAST form for a function with a cheap lower bound lb(p) <= f(p): the 96 samples are visited in increasing order of their bound
// and evaluated only while the bound is below the best value found (exact: a skipped sample cannot be the smallest), each evaluation
// warm-starting from the previous; then the same golden-section refinement around the best sample
template <class Fn, class Lb> inline double circle_min_lb(double R, double z, Fn f, Lb lb) {
    const int N = 96; std::array<std::pair<double, int>, N> order;
    for (int k = 0; k < N; ++k) { const double th = 2 * M_PI * k / N, p[3] = {R * std::cos(th), R * std::sin(th), z}; order[k] = {lb(p), k}; }
    std::sort(order.begin(), order.end());
    double best = 1e300; int kb = 0;
    for (const auto &[b, k] : order) { if (b >= best) break; const double th = 2 * M_PI * k / N, p[3] = {R * std::cos(th), R * std::sin(th), z}; const double v = f(p); if (v < best) best = v, kb = k; }
    double lo = 2 * M_PI * (kb - 1) / N, hi = 2 * M_PI * (kb + 1) / N; const double gr = 0.5 * (std::sqrt(5.0) - 1);
    auto g = [&](double th) { const double p[3] = {R * std::cos(th), R * std::sin(th), z}; return f(p); };
    double x1 = hi - gr * (hi - lo), x2 = lo + gr * (hi - lo), f1 = g(x1), f2 = g(x2);
    for (int it = 0; it < 80 && hi - lo > 1e-9; ++it) { if (f1 < f2) { hi = x2; x2 = x1; f2 = f1; x1 = hi - gr * (hi - lo); f1 = g(x1); } else { lo = x1; x1 = x2; f1 = f2; x2 = lo + gr * (hi - lo); f2 = g(x2); } }
    return std::min(best, std::min(f1, f2));
}

// distance from point p to segment [A, B]
inline double point_segment(const double *p, const double *A, const double *B) {
    double d[3] = {B[0] - A[0], B[1] - A[1], B[2] - A[2]}, w[3] = {p[0] - A[0], p[1] - A[1], p[2] - A[2]};
    const double dd = d[0] * d[0] + d[1] * d[1] + d[2] * d[2]; double t = dd > 0 ? (w[0] * d[0] + w[1] * d[1] + w[2] * d[2]) / dd : 0; t = std::max(0.0, std::min(1.0, t));
    double r = 0; for (int k = 0; k < 3; ++k) { const double e = w[k] - t * d[k]; r += e * e; } return std::sqrt(r);
}

// support point of the ellipsoid (semi-axes a, frame E in some frame) in direction m, and its support value h(m) - c.m = |A R^T m|
inline double ell_support(const double *a, const Frame &E, const double *m, double *x) {
    double ml[3]; for (int k = 0; k < 3; ++k) { ml[k] = 0; for (int j = 0; j < 3; ++j) ml[k] += E.R[3 * j + k] * m[j]; }
    double w[3], nn = 0; for (int k = 0; k < 3; ++k) { w[k] = a[k] * ml[k]; nn += w[k] * w[k]; } nn = std::sqrt(nn);
    if (x) { double l[3]; for (int k = 0; k < 3; ++k) l[k] = a[k] * w[k] / nn; for (int k = 0; k < 3; ++k) { x[k] = E.c[k]; for (int j = 0; j < 3; ++j) x[k] += E.R[3 * k + j] * l[j]; } }
    return nn;
}

// a point given in frame P, expressed in the ellipsoid frame (P and E both in the same parent frame)
inline void into(const Frame &E, const double *p, double *y) { to_local(E, p, y); }

// gap from the ellipsoid (semi-axes a, frame E in the box frame) to the box of half-sizes s (axis-aligned at the origin)
inline double gap_ellipsoid_box(const double *a, const Frame &E, const double *s) {
    double best = 1e300, pg[3][2];   // pg[k][side]: signed distance of the ellipsoid from the plane of face (k, side)
    for (int k = 0; k < 3; ++k) for (int si = 0; si < 2; ++si) {
        const int sg = si ? 1 : -1; double m[3] = {0, 0, 0}; m[k] = -sg; double x[3];
        const double hs = ell_support(a, E, m, x);
        pg[k][si] = sg * E.c[k] - s[k] - hs;
        if (sg * E.c[k] <= s[k]) continue;   // faces whose plane has the ellipsoid's centre outside: the plane gap when the foot is inside
        bool in = true; for (int j = 0; j < 3; ++j) if (j != k && std::fabs(x[j]) > s[j]) in = false;
        if (in) best = std::min(best, pg[k][si]);
    }
    const double amax = std::max(a[0], std::max(a[1], a[2]));
    for (int k = 0; k < 3; ++k) for (int s1 = -1; s1 <= 1; s1 += 2) for (int s2 = -1; s2 <= 1; s2 += 2) {   // the 12 edges; an edge's gap is
        const int i = (k + 1) % 3, j = (k + 2) % 3; double P[3], Q[3];                                       // at least the ellipsoid's distance
        P[k] = -s[k]; Q[k] = s[k]; P[i] = Q[i] = s1 * s[i]; P[j] = Q[j] = s2 * s[j];                         // to either face plane through it,
        const double lb = std::max(std::max(pg[i][s1 > 0], pg[j][s2 > 0]), point_segment(E.c, P, Q) - amax); // and its centre distance less the
        if (lb >= best) continue;                                                                            // largest semi-axis
        double A[3], B[3]; into(E, P, A); into(E, Q, B);
        best = std::min(best, seg_min_ellipsoid(a, A, B));
    }
    return best;
}

// gap from the ellipsoid (frame E in the cylinder frame) to the cylinder (radius R, half-height h, axis local z)
inline double gap_ellipsoid_cylinder(const double *a, const Frame &E, double R, double h) {
    double best = 1e300;
    for (int sg = -1; sg <= 1; sg += 2) {
        if (sg * E.c[2] > h) {   // flat end facing the ellipsoid
            const double m[3] = {0, 0, static_cast<double>(-sg)}; double x[3];
            const double hs = ell_support(a, E, m, x);
            if (x[0] * x[0] + x[1] * x[1] <= R * R) best = std::min(best, sg * E.c[2] - h - hs);
        }
        const double rho = std::sqrt(E.c[0] * E.c[0] + E.c[1] * E.c[1]), dz = E.c[2] - sg * h, amax = std::max(a[0], std::max(a[1], a[2]));
        if (std::sqrt((rho - R) * (rho - R) + dz * dz) - amax < best)   // the rim's lower bound: centre to circle less the largest semi-axis
            if (GEOM_FAST) { double tw = NAN;   // pruned by the bound |y| - a_max (distance to the centre less the largest semi-axis)
                best = std::min(best, circle_min_lb(R, sg * h, [&](const double *p) { double y[3]; into(E, p, y); return sdf_ellipsoid(a, y, nullptr, &tw); },
                    [&](const double *p) { double y[3]; into(E, p, y); return std::sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]) - amax; })); }
            else best = std::min(best, circle_min(R, sg * h, [&](const double *p) { double y[3]; into(E, p, y); return sdf_ellipsoid(a, y); }));
    }
    {   // lateral surface: distance from the ellipsoid to the axis line, less R, when the closest axis point is within the height
        const double P[3] = {0, 0, -1}, Q[3] = {0, 0, 1}; double A[3], B[3]; into(E, P, A); into(E, Q, B);
        double ts; const double dl = seg_min_ellipsoid(a, A, B, &ts); const double zt = -1 + 2 * ts;
        if (zt >= -h && zt <= h && dl > 0.5 * R) best = std::min(best, dl - R);
    }
    return best;
}

// signed distance from a point to the solid cylinder (radius R, half-height h, axis local z)
inline double sdf_cylinder(const double *q, double R, double h) {
    const double dr = std::sqrt(q[0] * q[0] + q[1] * q[1]) - R, dz = std::fabs(q[2]) - h;
    return std::min(std::max(dr, dz), 0.0) + std::sqrt(std::max(dr, 0.0) * std::max(dr, 0.0) + std::max(dz, 0.0) * std::max(dz, 0.0));
}

// gap from the capsule (segment A-B in the cylinder frame, radius r) to the cylinder (radius R, half-height h, axis local z)
inline double gap_capsule_cylinder(const double *A, const double *B, double r, double R, double h) {
    double best = std::min(sdf_cylinder(A, R, h), sdf_cylinder(B, R, h));
    for (int sg = -1; sg <= 1; sg += 2) best = std::min(best, circle_min(R, sg * h, [&](const double *p) { return point_segment(p, A, B); }));
    {   // lateral surface: segment against the axis line
        const double d[3] = {B[0] - A[0], B[1] - A[1], B[2] - A[2]};
        const double a2 = d[0] * d[0] + d[1] * d[1];   // the axis direction is z: minimise the horizontal distance over the segment
        double t = a2 > 0 ? -(A[0] * d[0] + A[1] * d[1]) / a2 : 0; t = std::max(0.0, std::min(1.0, t));
        const double px = A[0] + t * d[0], py = A[1] + t * d[1], pz = A[2] + t * d[2], dl = std::sqrt(px * px + py * py);
        if (pz >= -h && pz <= h && dl > 0.5 * R) best = std::min(best, dl - R);
    }
    return best - r;
}

}  // namespace ms3
