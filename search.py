"""
=======================================================================
SymSearch — Symbolic Search for Physical Constants
=======================================================================

Author  : Judicael Brindel, Independent Researcher, Cotonou, Benin
ORCID   : 0009-0007-4590-9874
GitHub  : https://github.com/judicaelbrindel/symsearch
License : MIT

Description:
    Given any number, SymSearch searches for the simplest symbolic
    expression matching it, using a bounded grammar over fundamental
    physical and mathematical constants.

Usage:
    python search.py

    Or import as module:
        from search import symsearch
        results = symsearch(1836.15267343, sigma=0.001)

Citation:
    Brindel, J. (2026). SymSearch: Symbolic Search for Physical Constants.
    Zenodo. https://doi.org/10.5281/zenodo.18792523

Computational assistance: Claude (Anthropic)
=======================================================================
"""

import math
import random
import json
import time

# ===================================================================
# FROZEN PARAMETERS
# ===================================================================

SEED    = 42
N_DEFAULT = 1_000_000
DEPTH   = 3
W1, W2, W3, W4 = 1.0, 0.8, 0.5, 0.3
S_INF   = 1000.0

# ===================================================================
# FUNDAMENTAL CONSTANTS DATABASE
# ===================================================================
# All values from PDG 2024 and CODATA 2018
# These are used as ATOMS in the symbolic grammar

phi = (1 + math.sqrt(5)) / 2

MATH_CONSTANTS = {
    # Mathematical
    'pi':       math.pi,
    'e':        math.e,
    'phi':      phi,
    'sqrt2':    math.sqrt(2),
    'sqrt3':    math.sqrt(3),
    'sqrt5':    math.sqrt(5),
    'sqrt7':    math.sqrt(7),
    'ln2':      math.log(2),
    'zeta3':    1.2020569031595943,
    'gamma_EM': 0.5772156649015329,
}

PHYSICAL_CONSTANTS = {
    # Fine structure constant
    'alpha':        1/137.035999084,
    'alpha_inv':    137.035999084,

    # Mass ratios relative to electron mass (PDG 2024)
    'mp_me':        1836.15267343,
    'mn_me':        1838.68366173,
    'mmu_me':       206.7682830,
    'mtau_me':      3477.23,
    'mt_me':        337711.0,
    'mb_me':        8185.9,
    'mc_me':        2394.2,
    'ms_me':        182.5,
    'mu_me':        7.2958,
    'md_me':        9.5019,
    'mW_me':        157298.0,
    'mZ_me':        178450.0,
    'mH_me':        245280.0,
    'mpi0_me':      264.1426,
    'mpipm_me':     273.1320,

    # Mass ratios relative to proton mass
    'me_mp':        1/1836.15267343,
    'mmu_mp':       206.7682830/1836.15267343,
    'mtau_mp':      3477.23/1836.15267343,
    'mb_mp':        8185.9/1836.15267343,
    'mc_mp':        2394.2/1836.15267343,
    'mt_mp':        337711.0/1836.15267343,
    'mW_mp':        157298.0/1836.15267343,
    'mZ_mp':        178450.0/1836.15267343,
    'mH_mp':        245280.0/1836.15267343,
    'mn_mp':        1838.68366173/1836.15267343,

    # Mass ratios between bosons
    'mH_mW':        245280.0/157298.0,
    'mH_mZ':        245280.0/178450.0,
    'mZ_mW':        178450.0/157298.0,
    'mt_mW':        337711.0/157298.0,
    'mt_mZ':        337711.0/178450.0,
    'mt_mH':        337711.0/245280.0,

    # Coupling constants
    'alphas':       0.1180,
    'sin2_thetaW':  0.23122,
    'cos2_thetaW':  1 - 0.23122,
    'sin_thetaW':   math.sqrt(0.23122),
    'cos_thetaW':   math.sqrt(1 - 0.23122),

    # CKM matrix elements (quark mixing)
    'Vud':          0.97373,
    'Vus':          0.2243,
    'Vub':          0.00382,
    'Vcd':          0.221,
    'Vcs':          0.975,
    'Vcb':          0.0408,
    'Vtd':          0.0086,
    'Vts':          0.0415,
    'Vtb':          1.014,

    # PMNS matrix elements (neutrino mixing)
    'sin2_theta12': 0.307,
    'sin2_theta23': 0.546,
    'sin2_theta13': 0.02224,
    'sin_theta12':  math.sqrt(0.307),
    'sin_theta23':  math.sqrt(0.546),
    'sin_theta13':  math.sqrt(0.02224),

    # Fermi constant (dimensionless proxy)
    'GF_mp2':       1.1663788e-5 * (938.272046)**2,

    # Cosmological dimensionless parameters (Planck 2018)
    'omega_b':      0.0493,
    'omega_dm':     0.2607,
    'omega_lambda': 0.6847,
    'ns':           0.9665,
}


# Combined constants available as grammar atoms
ALL_CONSTANTS = {**MATH_CONSTANTS, **PHYSICAL_CONSTANTS}

# ===================================================================
# SYNTAX TREE
# ===================================================================

class Node:
    def __init__(self, kind, val=None, children=None, op=None):
        self.kind     = kind
        self.val      = val
        self.children = children or []
        self.op       = op

    def evaluate(self):
        if self.kind == 'int':   return float(self.val)
        if self.kind == 'rat':   return self.val[0] / self.val[1]
        if self.kind == 'const': return ALL_CONSTANTS[self.val]
        if self.kind == 'unary':
            x = self.children[0].evaluate()
            if self.op == 'neg':  return -x
            if self.op == 'inv':
                if abs(x) < 1e-15: raise ValueError
                return 1.0 / x
            if self.op == 'sqrt':
                if x < 0: raise ValueError
                return math.sqrt(x)
            if self.op == 'log':
                if x <= 0: raise ValueError
                return math.log(x)
            if self.op == 'exp':
                if x > 20: raise ValueError
                r = math.exp(x)
                if not math.isfinite(r): raise ValueError
                return r
            if self.op == 'sq':  return x * x
            if self.op == 'cb':  return x * x * x
        if self.kind == 'binary':
            x = self.children[0].evaluate()
            y = self.children[1].evaluate()
            if self.op == 'add': return x + y
            if self.op == 'sub': return x - y
            if self.op == 'mul': return x * y
            if self.op == 'div':
                if abs(y) < 1e-15: raise ValueError
                return x / y
            if self.op == 'pow':
                if abs(x) > 1e6: raise ValueError
                if x < 0 and not float(y).is_integer(): raise ValueError
                r = x ** float(y)
                if not math.isfinite(r): raise ValueError
                return r
        raise ValueError

    def _stats(self):
        if self.kind == 'int':
            return (0, 1, 0, math.log2(1 + abs(self.val)))
        if self.kind == 'rat':
            return (0, 1, 0, math.log2(1 + abs(self.val[0])))
        if self.kind == 'const':
            return (0, 1, 1, 0.0)
        if self.kind == 'unary':
            n, d, nc, cs = self.children[0]._stats()
            return (n + 1, d + 1, nc, cs)
        if self.kind == 'binary':
            n1, d1, nc1, cs1 = self.children[0]._stats()
            n2, d2, nc2, cs2 = self.children[1]._stats()
            return (n1 + n2 + 1, max(d1, d2) + 1, nc1 + nc2, cs1 + cs2)
        return (0, 0, 0, 0.0)

    def complexity(self):
        n_op, depth, n_const, coeff_sum = self._stats()
        return W1 * n_op + W2 * depth + W3 * n_const + W4 * coeff_sum

    def to_str(self):
        if self.kind == 'int':   return str(self.val)
        if self.kind == 'rat':   return f"{self.val[0]}/{self.val[1]}"
        if self.kind == 'const': return self.val
        if self.kind == 'unary':
            return f"{self.op}({self.children[0].to_str()})"
        if self.kind == 'binary':
            return f"({self.children[0].to_str()} {self.op} {self.children[1].to_str()})"


ACTIVE_CONSTANTS = dict(ALL_CONSTANTS)  # can be modified by exclude

def random_tree(max_depth=3, cur=0):
    p_leaf = 0.25 + 0.35 * cur
    const_keys = list(ACTIVE_CONSTANTS.keys())
    if cur >= max_depth or random.random() < p_leaf:
        c = random.random()
        if c < 0.40:
            return Node('const', val=random.choice(const_keys))
        elif c < 0.70:
            n = random.choice([x for x in range(-50, 51) if x != 0])
            return Node('int', val=n)
        else:
            a = random.choice([x for x in range(-50, 51) if x != 0])
            b = random.randint(1, 20)
            return Node('rat', val=(a, b))
    else:
        kind = random.choice(['unary', 'binary', 'binary'])
        if kind == 'unary':
            op = random.choice(['neg', 'inv', 'sqrt', 'log', 'exp', 'sq', 'cb'])
            return Node('unary', op=op, children=[random_tree(max_depth, cur + 1)])
        else:
            op = random.choice(['add', 'sub', 'mul', 'mul', 'mul', 'div', 'pow'])
            return Node('binary', op=op,
                        children=[random_tree(max_depth, cur + 1),
                                  random_tree(max_depth, cur + 1)])

# ===================================================================
# CORE SEARCH FUNCTION
# ===================================================================

def symsearch(target, sigma, n=N_DEFAULT, depth=DEPTH, seed=SEED, verbose=True, exclude=None):
    """
    Search for the simplest symbolic expression matching target within sigma.

    Parameters
    ----------
    target : float
    sigma : float
    n : int
    depth : int
    seed : int
    verbose : bool
    exclude : list of str
        Constants to exclude from the grammar (e.g. ['alpha', 'alpha_inv'])
    """
    global ACTIVE_CONSTANTS
    if exclude:
        ACTIVE_CONSTANTS = {k: v for k, v in ALL_CONSTANTS.items() if k not in exclude}
        if verbose:
            print(f"\n  Excluded: {exclude}")
            print(f"  Active constants: {len(ACTIVE_CONSTANTS)}")
    else:
        ACTIVE_CONSTANTS = dict(ALL_CONSTANTS)

    random.seed(seed)
    n_valid   = 0
    matches   = []
    S_min     = float('inf')

    if verbose:
        print(f"\n  Target  : {target}")
        print(f"  Sigma   : {sigma}")
        print(f"  N       : {n:,}")
        print(f"  Depth   : {depth}")
        print(f"  Grammar : {len(ALL_CONSTANTS)} constants + integers [-50,50]")
        print(f"  Searching", end='', flush=True)

    t0 = time.time()
    report_every = n // 10

    for i in range(n):
        if verbose and i % report_every == 0:
            print('.', end='', flush=True)
        try:
            tree = random_tree(max_depth=depth)
            val  = tree.evaluate()
            if not math.isfinite(val):
                continue
            S = tree.complexity()
            n_valid += 1

            dev = abs(val - target) / sigma
            if dev <= 1.0:
                if S < S_min:
                    S_min = S
                matches.append((round(dev, 6), round(S, 4), round(val, 10), tree.to_str()))
        except:
            continue

    elapsed = time.time() - t0
    S_final = S_min if math.isfinite(S_min) else S_INF

    # Deduplicate by expression string, keep best sigma per expression
    seen = {}
    for dev, S, val, expr in matches:
        if expr not in seen or S < seen[expr][1]:
            seen[expr] = (dev, S, val, expr)
    top10 = sorted(seen.values(), key=lambda x: x[1])[:10]

    if verbose:
        print(f" done ({elapsed:.1f}s)\n")
        print("=" * 60)
        print(f"  RESULTS")
        print("=" * 60)
        print(f"  Valid expressions   : {n_valid:,}")
        print(f"  Matches within sigma: {len(matches)}")
        print(f"  Empirical p         : {len(matches)/n_valid:.2e}" if n_valid > 0 else "  Empirical p : n/a")
        print(f"  S_min               : {S_final:.2f}" if math.isfinite(S_min) else "  S_min : no match found")
        print()

        if top10:
            print(f"  {'Rank':<5} {'sigma':>8}  {'S':>6}  {'Value':>14}  Expression")
            print(f"  {'-'*5} {'-'*8}  {'-'*6}  {'-'*14}  {'-'*30}")
            for rank, (dev, S, val, expr) in enumerate(top10, 1):
                print(f"  {rank:<5} {dev:>8.4f}  {S:>6.2f}  {val:>14.8f}  {expr[:50]}")
        else:
            print("  No match found in this search space.")
            print("  Try increasing N or depth.")
        print("=" * 60)

    return {
        'target':       target,
        'sigma':        sigma,
        'n_valid':      n_valid,
        'n_matches':    len(matches),
        'p_empirical':  len(matches) / n_valid if n_valid > 0 else 0,
        'S_min':        S_final,
        'has_match':    math.isfinite(S_min),
        'top_results':  [
            {'sigma': m[0], 'S': m[1], 'value': m[2], 'expression': m[3]}
            for m in top10
        ],
        'elapsed_s':    round(elapsed, 1),
        'parameters':   {'N': n, 'depth': depth, 'seed': seed,
                         'n_constants': len(ALL_CONSTANTS)}
    }

# ===================================================================
# COMMAND LINE INTERFACE
# ===================================================================

def main():
    print()
    print("=" * 60)
    print("  SymSearch — Symbolic Search for Physical Constants")
    print("  Author: Judicael Brindel (Cotonou, Benin)")
    print("  Computational assistance: Claude (Anthropic)")
    print("=" * 60)
    print()
    print(f"  Grammar contains {len(ALL_CONSTANTS)} fundamental constants:")
    print(f"  Mathematical : {len(MATH_CONSTANTS)} constants")
    print(f"  Physical     : {len(PHYSICAL_CONSTANTS)} constants")
    print()

    # Get target
    while True:
        try:
            target_str = input("  Enter target number: ").strip()
            target = float(target_str)
            break
        except ValueError:
            print("  Invalid number. Please try again.")

    # Get sigma
    while True:
        try:
            sigma_str = input("  Enter tolerance (sigma): ").strip()
            sigma = float(sigma_str)
            if sigma <= 0:
                raise ValueError
            break
        except ValueError:
            print("  Invalid tolerance. Must be positive.")

    # Get N
    print()
    print("  Search size:")
    print("  [1] Quick    — N=100,000   (~20 seconds)")
    print("  [2] Standard — N=1,000,000 (~3 minutes)")
    print("  [3] Deep     — N=5,000,000 (~15 minutes)")
    print("  [4] Custom")

    choice = input("  Choose [1/2/3/4]: ").strip()
    if choice == '1':
        n = 100_000
    elif choice == '2':
        n = 1_000_000
    elif choice == '3':
        n = 5_000_000
    elif choice == '4':
        n = int(input("  Enter N: ").strip().replace(',', '').replace(' ', ''))
    else:
        n = 100_000

    # Get exclude
    print()
    print("  Exclude constants? (useful to find formulas not using a constant directly)")
    excl_str = input("  Enter names to exclude, comma-separated (or press Enter to skip): ").strip()
    exclude = [x.strip() for x in excl_str.split(',') if x.strip()] if excl_str else None

    # Run search
    results = symsearch(target, sigma, n=n, verbose=True, exclude=exclude)

    # Save results
    save = input("\n  Save results to JSON? [y/n]: ").strip().lower()
    if save == 'y':
        filename = f"symsearch_{target}_{n}.json"
        with open(filename, 'w') as f:
            json.dump(results, f, indent=2)
        print(f"  Saved to {filename}")


if __name__ == "__main__":
    main()
