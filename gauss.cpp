// Gauss circle problem: N(R) = #{(x,y) in Z^2 : x^2+y^2 <= R}
// fast(): convex-hull walk on the Stern-Brocot tree, ~O(R^{1/3} log R)
// brute(): O(sqrt R) two-pointer oracle
//
// build: g++ -O2 -std=c++17 -o gauss gauss.cpp
// usage: ./gauss fast R | brute R | bench | test
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <chrono>
#include <random>
#include <vector>
#include <map>
#include <utility>

typedef long long ll;
typedef __int128 i128;
typedef unsigned long long ull;

static ll isqrt_ll(ll n) {
    if (n <= 0) return 0;
    ll r = (ll)sqrtl((long double)n);
    while ((i128)r * r > n) --r;
    while ((i128)(r + 1) * (r + 1) <= n) ++r;
    return r;
}

static void print128(i128 v) {
    if (v == 0) { printf("0"); return; }
    char buf[64]; int n = 0;
    bool neg = v < 0; if (neg) v = -v;
    while (v > 0) { buf[n++] = '0' + (int)(v % 10); v /= 10; }
    if (neg) putchar('-');
    while (n) putchar(buf[--n]);
}

struct Stats { ull inside_calls = 0, edges = 0, pushes = 0, descent = 0; };

static i128 brute(ll R) {
    if (R < 0) return 0;
    ll r = isqrt_ll(R), y = r;
    i128 tot = 0;
    for (ll x = 0; x <= r; ++x) {
        while ((i128)x * x + (i128)y * y > R) --y;
        i128 w = 2 * (i128)y + 1;
        tot += (x == 0) ? w : 2 * w;
    }
    return tot;
}

struct Walker {
    ll R, m;
    Stats *st;
    inline bool inside(ll x, ll y) {
        st->inside_calls++;
        if (x > m) return false;
        return (i128)x * x + (i128)y * y <= R;
    }
    // largest k>=0 with (x+k dx, y-k dy) inside
    ll maxsteps(ll x, ll y, ll dx, ll dy) {
        if (!inside(x + dx, y - dy)) return 0;
        ll lo = 1, hi = 2;
        while (inside(x + hi * dx, y - hi * dy)) { lo = hi; hi *= 2; }
        while (hi - lo > 1) {
            ll mid = lo + (hi - lo) / 2;
            if (inside(x + mid * dx, y - mid * dy)) lo = mid; else hi = mid;
        }
        return lo;
    }
    // returns F = sum_{x=1}^m floor(sqrt(R-x^2))
    i128 run() {
        ll r = isqrt_ll(R);
        m = isqrt_ll(R / 2);
        if (m == 0) return 0;
        ll x = 0, y = r;
        i128 F = 0;
        std::vector<std::pair<ll, ll>> stack;
        stack.push_back({1, 1});
        stack.push_back({1, 0});
        while (x < m) {
            ll dx = stack.back().first, dy = stack.back().second;
            ll k = maxsteps(x, y, dx, dy);
            if (k) {
                i128 C = ((i128)dx * dy + dx + dy - 1) / 2;
                F += (i128)k * dx * y - (i128)dx * dy * ((i128)k * (k - 1) / 2) - (i128)k * C;
                x += k * dx; y -= k * dy;
                st->edges++;
            }
            std::pair<ll, ll> L = stack.back(); stack.pop_back();
            if (x >= m) break;
            while (!inside(x + stack.back().first, y - stack.back().second)) {
                L = stack.back(); stack.pop_back();
            }
            // Stern-Brocot descent
            ll cachedDx = -1, cachedDy = -1, cachedK = 0;
            while (true) {
                std::pair<ll, ll> Rt = stack.back();
                ll Mx = L.first + Rt.first, My = L.second + Rt.second;
                if (x + Mx > m) break;
                st->descent++;
                if (inside(x + Mx, y - My)) {
                    stack.push_back({Mx, My});
                    st->pushes++;
                } else {
                    if (cachedDx != Rt.first || cachedDy != Rt.second) {
                        cachedK = maxsteps(x, y, Rt.first, Rt.second);
                        cachedDx = Rt.first; cachedDy = Rt.second;
                    }
                    if (Mx + Rt.first > (cachedK + 1) * Rt.first - 1) break;
                    L = {Mx, My};
                }
            }
        }
        return F;
    }
};

static i128 fast(ll R, Stats *st) {
    if (R < 0) return 0;
    if (R == 0) return 1;
    Walker w; w.R = R; w.st = st;
    i128 F = w.run();
    ll r = isqrt_ll(R);
    i128 T = 2 * F - (i128)w.m * w.m;
    return 1 + 4 * (i128)r + 4 * T;
}

static double now() {
    using namespace std::chrono;
    return duration<double>(steady_clock::now().time_since_epoch()).count();
}

static bool eq(i128 a, i128 b) { return a == b; }

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: gauss fast R | brute R | bench | test\n"); return 1; }
    std::string cmd = argv[1];
    if (cmd == "fast" && argc >= 3) {
        Stats s; ll R = strtoll(argv[2], 0, 10);
        print128(fast(R, &s)); printf("\n"); return 0;
    }
    if (cmd == "brute" && argc >= 3) {
        ll R = strtoll(argv[2], 0, 10);
        print128(brute(R)); printf("\n"); return 0;
    }
    if (cmd == "test") {
        long bad = 0, n = 0;
        for (ll R = 0; R <= 2000000; ++R) {
            Stats s; ++n;
            if (!eq(fast(R, &s), brute(R))) { if (bad++ < 5) printf("MISMATCH R=%lld\n", R); }
        }
        printf("exhaustive R in [0, 2e6]: %ld tested, %ld mismatches\n", n, bad);
        // structured edge cases: perfect squares, R = a^2+b^2, neighbours
        bad = 0; n = 0;
        std::mt19937_64 rng(12345);
        for (int t = 0; t < 4000; ++t) {
            ll a = (ll)(rng() % 3000000) + 1, b = (ll)(rng() % 3000000);
            ll cands[6] = { a * a, a * a - 1, a * a + 1, a * a + b * b, a * a + b * b - 1, 2 * a * a };
            for (ll R : cands) { Stats s; ++n; if (!eq(fast(R, &s), brute(R))) { if (bad++ < 5) printf("MISMATCH R=%lld\n", R); } }
        }
        printf("structured (squares, sums of two squares, 2a^2, +-1) up to ~1.8e13: %ld tested, %ld mismatches\n", n, bad);
        // random large
        bad = 0; n = 0;
        for (int t = 0; t < 300; ++t) {
            ll R = (ll)(rng() % 900000000000000ULL) + 100000000000ULL; // up to ~9e14
            Stats s; ++n;
            if (!eq(fast(R, &s), brute(R))) { if (bad++ < 5) printf("MISMATCH R=%lld\n", R); }
        }
        printf("random R in [1e11, 9e14]: %ld tested, %ld mismatches\n", n, bad);
        return 0;
    }
    if (cmd == "bench") {
        printf("%-8s %-22s %10s %10s %10s %10s %10s\n", "R", "N(R)", "t_fast(s)", "t_brute(s)", "inside", "edges", "R^(1/3)");
        for (int e = 6; e <= 18; e += 2) {
            ll R = 1; for (int i = 0; i < e; ++i) R *= 10; R += 12345;
            Stats s; double t0 = now(); i128 a = fast(R, &s); double t1 = now();
            double tb = -1; bool ok = true;
            if (e <= 18) { double b0 = now(); i128 b = brute(R); tb = now() - b0; ok = eq(a, b); }
            printf("1e%-6d ", e); print128(a);
            printf("%*s %10.5f %10.3f %10llu %10llu %10.0f %s\n", 1, "", t1 - t0, tb,
                   s.inside_calls, s.edges, cbrt((double)R), ok ? "OK" : "MISMATCH");
        }
        // R near the int64 ceiling of the sum: 4e18
        {
            ll R = 4000000000000000000LL; Stats s; double t0 = now(); i128 a = fast(R, &s); double t1 = now();
            printf("R=4e18: N="); print128(a); printf("  t_fast=%.5f s  inside=%llu\n", t1 - t0, s.inside_calls);
        }
        return 0;
    }
    return 1;
}
