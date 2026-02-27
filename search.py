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
    Zenodo. https://doi.org/10.5281/zenodo.18805643

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
# COINCIDENCE CATALOGUE
# ===================================================================

CATALOGUE_FILE = "symsearch_catalogue.json"
CATALOGUE_P_THRESHOLD = 1e-4  # Auto-log if p < this value

def load_catalogue():
    try:
        with open(CATALOGUE_FILE) as f:
            return json.load(f)
    except:
        return {"metadata": {"tool": "SymSearch Coincidence Catalogue",
                             "author": "Judicael Brindel",
                             "github": "https://github.com/judicaelbrindel/symsearch"},
                "entries": []}

def save_catalogue(cat):
    with open(CATALOGUE_FILE, 'w') as f:
        json.dump(cat, f, indent=2)

def auto_log_catalogue(target, sigma, n_valid, n_matches, p_empirical, top_results, label=None):
    """Automatically log result to catalogue if p < threshold."""
    if p_empirical == 0 or p_empirical >= CATALOGUE_P_THRESHOLD:
        return False
    if not top_results:
        return False
    cat = load_catalogue()
    best = top_results[0]
    entry = {
        "date": time.strftime("%Y-%m-%d"),
        "label": label or f"target_{target}",
        "target_value": target,
        "sigma": sigma,
        "formula": best[3] if len(best) > 3 else str(best),
        "formula_value": best[2] if len(best) > 2 else None,
        "S": best[1] if len(best) > 1 else None,
        "p_empirical": round(p_empirical, 8),
        "n_valid": n_valid,
        "status": "candidate",
        "author": "Judicael Brindel",
    }
    existing = [e for e in cat["entries"]
                if e.get("target_value") == target and e.get("formula") == entry["formula"]]
    if existing:
        return False
    cat["entries"].append(entry)
    save_catalogue(cat)
    return True

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
    report_every = max(1, n // 10)

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

    p_emp = len(matches) / n_valid if n_valid > 0 else 0

    # Auto-log to catalogue if result is rare
    top_raw = [(m[0], m[1], m[2], m[3]) for m in top10]
    logged = auto_log_catalogue(target, sigma, n_valid, len(matches), p_emp, top_raw)
    if logged and verbose:
        print(f"  *** Rare result auto-logged to {CATALOGUE_FILE} ***")

    return {
        'target':       target,
        'sigma':        sigma,
        'n_valid':      n_valid,
        'n_matches':    len(matches),
        'p_empirical':  p_emp,
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

# ===================================================================
# INVERSE SEARCH
# ===================================================================

def safe_eval_formula(formula_str):
    """
    Safely evaluate a formula string using the constants database.
    Supported: +, -, *, /, **, sqrt(), log(), exp(), pi, e, phi, etc.
    Example: "6*pi**5" or "mp_me / (4*pi)"
    """
    # Build safe namespace
    safe_ns = dict(ALL_CONSTANTS)
    safe_ns.update({
        'sqrt': math.sqrt,
        'log':  math.log,
        'exp':  math.exp,
        'sin':  math.sin,
        'cos':  math.cos,
        'tan':  math.tan,
        'abs':  abs,
    })
    try:
        result = eval(formula_str, {"__builtins__": {}}, safe_ns)
        return float(result)
    except Exception as e:
        return None


def inverse_search(formula_str, top_n=10, rel_tol=0.01):
    """
    Given a formula string, find which known constants it is closest to.

    Parameters
    ----------
    formula_str : str
        A mathematical expression, e.g. "6*pi**5" or "mp_me/alpha"
    top_n : int
        Number of closest constants to return
    rel_tol : float
        Relative tolerance for flagging close matches (default 1%)

    Returns
    -------
    dict with computed value and ranked list of closest constants
    """
    value = safe_eval_formula(formula_str)
    if value is None:
        print(f"  Error: could not evaluate '{formula_str}'")
        print(f"  Check syntax. Use **, *, /, +, - and constant names.")
        return None

    print()
    print("=" * 60)
    print(f"  INVERSE SEARCH")
    print("=" * 60)
    print(f"  Formula : {formula_str}")
    print(f"  Value   : {value:.10g}")
    print()

    # Compare to all known constants
    distances = []
    for name, cval in ALL_CONSTANTS.items():
        if cval == 0:
            continue
        rel_diff = abs(value - cval) / abs(cval)
        distances.append((rel_diff, name, cval))

    distances.sort(key=lambda x: x[0])
    top = distances[:top_n]

    print(f"  {'Rank':<5} {'Constant':20s}  {'Known value':>14}  {'Rel. diff':>10}  {'Match?'}")
    print(f"  {'-'*5} {'-'*20}  {'-'*14}  {'-'*10}  {'-'*6}")
    for rank, (rel_diff, name, cval) in enumerate(top, 1):
        match = "YES ***" if rel_diff < rel_tol else ""
        print(f"  {rank:<5} {name:20s}  {cval:>14.8g}  {rel_diff:>10.2e}  {match}")

    print()
    close = [(n, v, d) for d, n, v in top if d < rel_tol]
    if close:
        print(f"  Close matches (rel. diff < {rel_tol*100:.0f}%):")
        for name, cval, rel_diff in close:
            print(f"  -> {name} = {cval:.10g}  (diff = {rel_diff:.2e})")
    else:
        print(f"  No known constant within {rel_tol*100:.0f}% of this value.")
        print(f"  Closest: {top[0][1]} = {top[0][2]:.6g}  (diff = {top[0][0]:.2e})")
    print("=" * 60)

    return {
        'formula':   formula_str,
        'value':     value,
        'top_matches': [
            {'constant': n, 'known_value': v, 'rel_diff': d}
            for d, n, v in top
        ]
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
    print(f"  Grammar: {len(ALL_CONSTANTS)} fundamental constants")
    print(f"  Mathematical: {len(MATH_CONSTANTS)}  |  Physical: {len(PHYSICAL_CONSTANTS)}")
    print()
    print("  Mode:")
    print("  [1] Forward search  — enter a number, get formulas")
    print("  [2] Inverse search  — enter a formula, get matching constants")
    print()
    mode = input("  Choose [1/2]: ").strip()

    # ── INVERSE MODE ─────────────────────────────────────────────
    if mode == '2':
        print()
        print("  Enter a formula using standard Python syntax.")
        print("  Available: +, -, *, /, ** (power), sqrt(), log(), exp()")
        print(f"  Constants: pi, e, phi, sqrt2, mp_me, alpha, alphas, ...")
        print(f"  Example: 6*pi**5")
        print(f"  Example: mp_me / (4*pi**2)")
        print(f"  Example: sqrt(alpha_inv) + phi")
        print()
        formula_str = input("  Enter formula: ").strip()

        tol_str = input("  Relative tolerance for match (default 0.01 = 1%): ").strip()
        try:
            rel_tol = float(tol_str) if tol_str else 0.01
        except:
            rel_tol = 0.01

        results = inverse_search(formula_str, rel_tol=rel_tol)

        if results:
            save = input("\n  Save results to JSON? [y/n]: ").strip().lower()
            if save == 'y':
                filename = f"symsearch_inverse_{formula_str[:20].replace('/', 'div').replace('*','x')}.json"
                with open(filename, 'w') as f:
                    json.dump(results, f, indent=2)
                print(f"  Saved to {filename}")
        return

    # ── FORWARD MODE ─────────────────────────────────────────────
    while True:
        try:
            target_str = input("  Enter target number: ").strip()
            target = float(target_str)
            break
        except ValueError:
            print("  Invalid number. Please try again.")

    while True:
        try:
            sigma_str = input("  Enter tolerance (sigma): ").strip()
            sigma = float(sigma_str)
            if sigma <= 0:
                raise ValueError
            break
        except ValueError:
            print("  Invalid tolerance. Must be positive.")

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
        n = int(input("  Enter custom N (e.g. 1000000): ").strip().replace(',', '').replace(' ', ''))
    else:
        n = 100_000

    print()
    depth_str = input("  Search depth (default 3, max 5): ").strip()
    try:
        depth = int(depth_str)
        if depth < 1 or depth > 5:
            depth = 3
    except:
        depth = 3

    print()
    print("  Exclude constants? (useful to find formulas not using a constant directly)")
    excl_str = input("  Enter names to exclude, comma-separated (or press Enter to skip): ").strip()
    exclude = [x.strip() for x in excl_str.split(',') if x.strip()] if excl_str else None

    results = symsearch(target, sigma, n=n, depth=depth, verbose=True, exclude=exclude)

    save = input("\n  Save results to JSON? [y/n]: ").strip().lower()
    if save == 'y':
        filename = f"symsearch_{target}_{n}.json"
        with open(filename, 'w') as f:
            json.dump(results, f, indent=2)
        print(f"  Saved to {filename}")


if __name__ == "__main__":
    main()
