"""
Gauss circle problem: N(R) = #{(x,y) in Z^2 : x^2 + y^2 <= R}.

Three independent methods:
  brute(R)      O(sqrt R)  two-pointer sweep            (oracle A)
  leibniz(R)    O(R)       N = 1 + 4*sum floor(R/(4k+1)) - floor(R/(4k+3))  (oracle B)
  fast(R)       ~O(R^(1/3) log R) convex hull walk on the Stern-Brocot tree

Reduction used by fast():
  N = 1 + 4r + 4T,  r = isqrt(R),  T = #{x,y >= 1 : x^2+y^2 <= R}
  m = max{m : 2m^2 <= R},  T = 2*F - m^2,  F = sum_{x=1}^{m} floor(sqrt(R-x^2))
"""
from math import isqrt


def brute(R):
    if R < 0:
        return 0
    r = isqrt(R)
    y = r
    tot = 0
    for x in range(0, r + 1):
        while x * x + y * y > R:
            y -= 1
        tot += (2 * y + 1)          # y >= 0 always here
    # tot counts all (x>=0, y in [-y,y]); mirror x>0
    # recompute cleanly: sum over x in [-r,r]
    y = r
    tot = 0
    for x in range(0, r + 1):
        while x * x + y * y > R:
            y -= 1
        w = 2 * y + 1
        tot += w if x == 0 else 2 * w
    return tot


def leibniz(R):
    """N = 1 + 4 * sum_{k>=0} (floor(R/(4k+1)) - floor(R/(4k+3)))."""
    s = 0
    k = 0
    while 4 * k + 1 <= R:
        s += R // (4 * k + 1)
        if 4 * k + 3 <= R:
            s -= R // (4 * k + 3)
        k += 1
    return 1 + 4 * s


def F_bruteforce(R):
    m = isqrt(R // 2)
    return sum(isqrt(R - x * x) for x in range(1, m + 1))


class Stats:
    __slots__ = ("inside_calls", "edges", "pushes", "descent_steps")

    def __init__(self):
        self.inside_calls = 0
        self.edges = 0
        self.pushes = 0
        self.descent_steps = 0


def F_fast(R, stats=None):
    """sum_{x=1}^{m} floor(sqrt(R - x^2)) via hull walk. Returns (F, m)."""
    if stats is None:
        stats = Stats()
    r = isqrt(R)
    m = isqrt(R // 2)
    if m == 0:
        return 0, 0

    def inside(x, y):
        stats.inside_calls += 1
        return x <= m and x * x + y * y <= R

    def maxsteps(x, y, dx, dy):
        """largest k >= 0 with (x + k dx, y - k dy) inside (interval property)."""
        if not inside(x + dx, y - dy):
            return 0
        lo, hi = 1, 2
        while inside(x + hi * dx, y - hi * dy):
            lo, hi = hi, hi * 2
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if inside(x + mid * dx, y - mid * dy):
                lo = mid
            else:
                hi = mid
        return lo

    x, y = 0, r
    F = 0
    stack = [(1, 1), (1, 0)]  # top = flattest candidate
    while x < m:
        # 1. walk as far as possible along the flattest candidate
        dx, dy = stack[-1]
        k = maxsteps(x, y, dx, dy)
        if k:
            C = (dx * dy + dx + dy - 1) // 2
            F += k * dx * y - dx * dy * (k * (k - 1) // 2) - k * C
            x += k * dx
            y -= k * dy
            stats.edges += 1
        L = stack.pop()
        if x >= m:
            break
        # 2. absorb candidates that are now exhausted
        while not inside(x + stack[-1][0], y - stack[-1][1]):
            L = stack.pop()
        # 3. Stern-Brocot descent between L (outside) and Rt = stack top (inside)
        kcache = {}
        while True:
            Rt = stack[-1]
            Mx, My = L[0] + Rt[0], L[1] + Rt[1]
            if x + Mx > m:
                break
            stats.descent_steps += 1
            if inside(x + Mx, y - My):
                stack.append((Mx, My))
                stats.pushes += 1
            else:
                # reach bound: any candidate between L and Rt has dx <= (k+1)*dxRt - 1
                if Rt not in kcache:
                    kcache[Rt] = maxsteps(x, y, Rt[0], Rt[1])
                kk = kcache[Rt]
                if Mx + Rt[0] > (kk + 1) * Rt[0] - 1:
                    break
                L = (Mx, My)
    return F, m


def fast(R, stats=None):
    if R < 0:
        return 0
    if R == 0:
        return 1
    r = isqrt(R)
    F, m = F_fast(R, stats)
    T = 2 * F - m * m
    return 1 + 4 * r + 4 * T


if __name__ == "__main__":
    import sys
    R = int(sys.argv[1])
    print(fast(R))
