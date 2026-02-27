---
title: 'SymSearch: A Python Tool for Symbolic Search of Physical Constants'
tags:
  - Python
  - physics
  - symbolic computation
  - physical constants
  - hypothesis generation
authors:
  - name: Judicael Brindel
    orcid: 0009-0007-4590-9874
    affiliation: 1
affiliations:
  - name: Independent Researcher, Cotonou, Benin
    index: 1
date: 27 February 2026
bibliography: paper.bib
---

# Summary

SymSearch is an open-source Python tool that searches for symbolic mathematical expressions matching a given numerical value, using a formal grammar over 68 fundamental physical and mathematical constants. It supports two complementary search modes: forward search (given a number, find candidate expressions) and inverse search (given a formula, identify which known physical constant it most closely approximates). The tool is designed to assist researchers in hypothesis generation and in identifying unexpected numerical connections between physical constants. All results include a symbolic complexity score and an empirical rarity measure. SymSearch is fully reproducible, archived on Zenodo [@brindel2026symsearch], and accessible via Google Colab without local installation.

# Statement of Need

The search for simple mathematical relationships among dimensionless physical constants has a long history in theoretical physics [@wyler1969; @chekanov2025]. However, most empirical searches suffer from methodological weaknesses: they are post-hoc, lack a formalized search space, and provide no statistical baseline for evaluating whether a proposed formula is genuinely rare or merely a numerical coincidence.

Existing tools such as the Inverse Symbolic Calculator (ISC) and Wolfram Alpha offer partial functionality but are closed-source, do not include physical constants as grammar atoms, and provide no statistical rarity measure. SymSearch addresses these gaps by providing a formally defined, reproducible, and statistically grounded symbolic search framework built specifically around the vocabulary of fundamental physics.

The tool is useful in two scenarios. First, when an experimentalist measures a new dimensionless quantity and wants to identify candidate symbolic representations. Second, when a theorist derives an expression and wants to verify whether it approximates a known physical constant. Both use cases are supported natively by SymSearch's two search modes.

# Description

SymSearch generates random free syntax trees over a grammar G. Terminal nodes are drawn from 68 fundamental constants (Table 1) and rational coefficients $a/b$ with $a \in [-50, 50]$ and $b \in [1, 20]$. Internal nodes are unary or binary operators: $\{+, -, \times, /, \hat{\;}, \sqrt{\cdot}, \log, \exp, \text{sq}, \text{cb}, \text{neg}, \text{inv}\}$. Tree depth is bounded by $d = 3$ by default.

For each expression $e$, the symbolic complexity is:

$$S(e) = 1.0 \cdot N_{op} + 0.8 \cdot \text{depth} + 0.5 \cdot N_{const} + 0.3 \sum \log_2(1 + |c|)$$

where $N_{op}$ is the number of operations, depth is the maximum tree depth, $N_{const}$ is the number of distinct constants used, and the sum runs over integer coefficients $c$. Lower $S$ indicates simpler expressions.

| Category | Constants |
|---|---|
| Mathematical (10) | pi, e, phi, sqrt2, sqrt3, sqrt5, sqrt7, ln2, zeta3, gamma_EM |
| Mass ratios / electron (15) | mp_me, mn_me, mmu_me, mtau_me, mt_me, mb_me, mc_me, ms_me, mu_me, md_me, mW_me, mZ_me, mH_me, mpi0_me, mpipm_me |
| Mass ratios / proton (10) | me_mp, mmu_mp, mtau_mp, mb_mp, mc_mp, mt_mp, mW_mp, mZ_mp, mH_mp, mn_mp |
| Boson ratios (6) | mH_mW, mH_mZ, mZ_mW, mt_mW, mt_mZ, mt_mH |
| Coupling constants (5) | alpha, alpha_inv, alphas, sin2_thetaW, cos2_thetaW |
| CKM matrix (9) | Vud, Vus, Vub, Vcd, Vcs, Vcb, Vtd, Vts, Vtb |
| PMNS matrix (6) | sin2_theta12, sin2_theta23, sin2_theta13, sin_theta12, sin_theta23, sin_theta13 |
| Cosmological (4) | omega_b, omega_dm, omega_lambda, ns |
| Other (3) | GF_mp2, sin_thetaW, cos_thetaW |

Table 1. Constants available as grammar atoms in SymSearch (PDG 2024, CODATA 2018, Planck 2018).

In **forward search** mode, the user provides a target value and a tolerance $\sigma$. SymSearch generates $N$ random expressions, evaluates each, and returns the top matches ranked by $S$. An empirical rarity probability $p = n_{matches} / n_{valid}$ is also reported. In **inverse search** mode, the user provides a formula string (e.g. `6*pi**5`), which is evaluated using the constants database, and the result is compared against all 68 known constants by relative difference.

# Usage Example

SymSearch can be used interactively or imported as a Python module:

```python
# Install and load
!wget -O search.py https://raw.githubusercontent.com/judicaelbrindel/symsearch/main/search.py
exec(open('search.py').read())

# Forward search: find formulas for the fine-structure constant
results = symsearch(0.0072973525693, sigma=1e-7, n=1_000_000,
                    exclude=['alpha', 'alpha_inv'])

# Inverse search: identify which constant 6*pi**5 approximates
results = inverse_search('6*pi**5', rel_tol=0.01)
# Output: mp_me = 1836.1527  (rel. diff = 1.88e-05)  YES ***
```

The inverse search correctly identifies that $6\pi^5 \approx 1836.12$ is within 0.002% of the proton-to-electron mass ratio $m_p/m_e = 1836.153$, confirming the formula proposed in @brindel2026mpme.

# Limitations

SymSearch performs a stochastic bounded search, not an exhaustive one. At depth $d = 3$ and $N = 1{,}000{,}000$, a small fraction of the full expression space is explored. A statistical study using SymSearch's batch protocol found no significant difference in compressibility between physical constants and random numbers at depth $d = 3$ (Mann-Whitney $p = 0.2754$) [@brindel2026compressibility], establishing a baseline for future deeper searches.

# Acknowledgements

The author thanks the open-source community for tools enabling independent research. Computational assistance for implementation and formal verification was provided by Claude (Anthropic).

# References
