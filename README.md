# gauss-circle

Exact lattice-point counting in a disk in roughly O(R^(1/3)) time, using a Stern–Brocot convex-hull walk.

Computes the Gauss circle count

    N(R) = #{ (x, y) in Z^2 : x^2 + y^2 <= R }

exactly, for R up to about 4·10^18, with integer arithmetic only. At R = 10^18 it takes about 27 ms, versus 3.7 s for the standard O(√R) sweep.

## How it works

1. Reduce the problem to one octant: `N = 1 + 4r + 4(2F - m^2)`, where `F = sum_{x=1}^{m} floor(sqrt(R - x^2))`, `r = isqrt(R)`, and `m` is the largest integer with `2m^2 <= R`.
2. On that octant the boundary slope lies in [-1, 0]. Walk the upper convex hull of the lattice points under the arc, finding each hull edge by descending the Stern–Brocot tree.
3. Sum each edge in closed form (a floor-sum identity for coprime step vectors).

Full derivation, lemmas, and experiments are in [`whitepaper.pdf`](whitepaper.pdf).

## Files

| File | Purpose |
|---|---|
| `gauss.cpp` | C++17 implementation, `__int128` arithmetic, test suite and benchmark |
| `gauss.py` | Python reference implementation, plus two independent oracles |
| `whitepaper.tex` / `whitepaper.pdf` | Write-up: derivation, proofs (sketched), results, limitations |
| `bench.txt` | Benchmark output used in the paper |

## Usage

```bash
g++ -O2 -std=c++17 -o gauss gauss.cpp

./gauss fast  1000000000000000000    # fast algorithm
./gauss brute 1000000000000000000    # O(sqrt R) oracle
./gauss bench                        # timing table, R = 1e6 ... 1e18
./gauss test                         # correctness suite (~3.5 min)
```

Python:

```bash
python3 gauss.py 1000000000000
```

Build the paper:

```bash
pdflatex whitepaper.tex && pdflatex whitepaper.tex
```

## Results

| R | fast | O(√R) sweep | predicate calls | hull edges |
|---|---|---|---|---|
| 10^12 | 0.30 ms | 5 ms | 64,675 | 4,313 |
| 10^14 | 1.8 ms | 38 ms | 301,730 | 20,117 |
| 10^16 | 6.2 ms | 409 ms | 1,410,596 | 93,259 |
| 10^18 | 27 ms | 3.69 s | 6,553,646 | 432,570 |

Work fits R^0.338 (log–log fit), consistent with R^(1/3). Single core, `-O2`, single runs.

## Verification

Zero mismatches against the O(√R) oracle on:

- every R in [0, 2·10^6]
- 24,000 structured inputs (a^2, a^2±1, a^2+b^2, a^2+b^2−1, 2a^2)
- 300 random R in [10^11, 9·10^14]
- R = 10^k + 12345 for k = 6, 8, …, 18, and R = 4·10^18

## Limitations

- The O(R^(1/3) log R) bound on total descent work is measured, not proved.
- Correctness proofs in the paper are sketches, not machine-checked.
- Testing is exhaustive only up to 2·10^6; beyond that it is sampled.
- Input must fit in a signed 64-bit integer (R < 2^63).

## Prior art

This is not a new algorithm. Counting lattice points under a convex curve with a Stern–Brocot walk is known for the hyperbola (divisor summatory function); this project adapts it to the circle and documents it with tests. The literature search was not exhaustive.

## Author

Joseph Anointed — seth.cse@gmail.com
