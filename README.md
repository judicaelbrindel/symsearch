# SymSearch — Symbolic Search for Physical Constants

**Author:** Judicael Brindel, Independent Researcher, Cotonou, Benin  
**ORCID:** 0009-0007-4590-9874  
**License:** MIT  
**Computational assistance:** Claude (Anthropic)

---

## What is SymSearch?

SymSearch is an open-source tool that searches for the **simplest symbolic expression** matching any given number, using a formal grammar over fundamental physical and mathematical constants.

You enter a number and a tolerance. SymSearch searches millions of symbolic expressions and returns the simplest candidates that match.

---

## What SymSearch can do

- Find simple symbolic expressions that approximate a given number
- Suggest connections between constants (e.g. a new measurement expressed as a ratio of known masses)
- Rank results by symbolic complexity (simpler = ranked higher)
- Give a statistical rarity score
- Save results to JSON for further analysis

## What SymSearch cannot do

- **Guarantee the exact theoretical formula** — it finds candidates, not proofs
- **Cover the full expression space** — depth=3 and N=1,000,000 is a bounded sample
- **Replace physical reasoning** — results are suggestions, not conclusions

SymSearch is a tool for **exploration and hypothesis generation**, not for verification.

---

## Example

```
Enter target number: 137.036
Enter tolerance (sigma): 0.001
Choose [1] Quick

RESULTS
Valid expressions   : 91,888
Matches within sigma: 355
S_min               : 1.30

Rank  sigma     S        Value         Expression
1     0.0009   1.30   137.03599908   alpha_inv
2     0.0009   3.10   137.03599908   inv(alpha)
3     0.5437   3.60   137.03654370   (me_mp add alpha_inv)
4     0.5455   3.60   137.03545447   (alpha_inv sub me_mp)
```

---

## Grammar — 68 fundamental constants

**10 Mathematical constants:**
`pi, e, phi, sqrt2, sqrt3, sqrt5, sqrt7, ln2, zeta3, gamma_EM`

**58 Physical constants (PDG 2024 / CODATA 2018 / Planck 2018):**

| Category | Constants |
|---|---|
| Fine structure | `alpha, alpha_inv` |
| Mass ratios / electron | `mp_me, mn_me, mmu_me, mtau_me, mt_me, mb_me, mc_me, ms_me, mu_me, md_me, mW_me, mZ_me, mH_me, mpi0_me, mpipm_me` |
| Mass ratios / proton | `me_mp, mmu_mp, mtau_mp, mb_mp, mc_mp, mt_mp, mW_mp, mZ_mp, mH_mp, mn_mp` |
| Boson ratios | `mH_mW, mH_mZ, mZ_mW, mt_mW, mt_mZ, mt_mH` |
| Coupling constants | `alphas, sin2_thetaW, cos2_thetaW, sin_thetaW, cos_thetaW` |
| CKM matrix | `Vud, Vus, Vub, Vcd, Vcs, Vcb, Vtd, Vts, Vtb` |
| PMNS matrix | `sin2_theta12, sin2_theta23, sin2_theta13, sin_theta12, sin_theta23, sin_theta13` |
| Fermi constant | `GF_mp2` |
| Cosmological | `omega_b, omega_dm, omega_lambda, ns` |

**Operations:** `+, -, x, /, ^, sqrt, log, exp, sq, cb, neg, inv`

**Coefficients:** integers in [-50, 50], rationals a/b with b in [1, 20]

**Maximum depth:** 3 (configurable)

---

## Complexity Metric S

```
S(e) = 1.0 * N_op + 0.8 * depth + 0.5 * N_const + 0.3 * sum(log2(1 + |coeff|))
```

Lower S = simpler expression. Weights fixed before any computation.

---

## Installation

```bash
git clone https://github.com/judicaelbrindel/symsearch
cd symsearch
pip install -r requirements.txt
python search.py
```

---

## Google Colab

Upload `search.py` to Colab and run:

```python
!pip install numpy scipy matplotlib
exec(open('search.py').read())
results = symsearch(1836.15267343, sigma=0.001, n=100_000, verbose=True)
```

---

## Search sizes

| Option | N | Time (approx) |
|---|---|---|
| Quick | 100,000 | ~20 seconds |
| Standard | 1,000,000 | ~3 minutes |
| Deep | 5,000,000 | ~15 minutes |

---

## Scientific Background

> Brindel, J. (2026). *Algorithmic Compressibility of Dimensionless Physical Constants Under Bounded Symbolic Grammars: A Negative Empirical Study.* Zenodo.

At depth d=3, no statistically significant difference in compressibility was found between physical constants and random numbers (Mann-Whitney p=0.2754). Deeper searches may reveal structure not visible at depth 3.

---

## Contributing

- More physical constants
- Deeper grammar (d=4, d=5)
- Alternative complexity metrics
- Web interface
- Performance optimization (C++, Julia)

---

## Citation

```bibtex
@software{brindel2026symsearch,
  author    = {Brindel, Judicael},
  title     = {SymSearch: Symbolic Search for Physical Constants},
  year      = {2026},
  publisher = {GitHub},
  url       = doi = {10.5281/zenodo.18792523},
url = {https://github.com/judicaelbrindel/symsearch}
```

---

## License

MIT License. Free to use, modify, and distribute.
