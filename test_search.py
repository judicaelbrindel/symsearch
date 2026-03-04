"""
test_search.py — Automated tests for SymSearch
===============================================

Author  : Judicael Brindel, Independent Researcher, Cotonou, Benin
ORCID   : 0009-0007-4590-9874
License : MIT

Run with:
    pip install pytest
    pytest test_search.py -v

These tests verify the core functionality of the symsearch() function
and the inverse_search() function.
"""

import math
import pytest
from search import symsearch, inverse_search, safe_eval_formula


# ===================================================================
# FIXTURES — small N for fast tests
# ===================================================================

N_TEST = 100_000   # fast — enough for basic coverage
SEED   = 42


# ===================================================================
# TEST 1 — alpha_inv finds itself
# ===================================================================

def test_alpha_inv_finds_itself():
    """
    alpha_inv = 137.035999084 must be found at depth=3.
    It is a terminal atom in the grammar — S should be very low.
    """
    result = symsearch(
        target=137.035999084,
        sigma=0.001,
        n=N_TEST,
        depth=3,
        seed=SEED,
        verbose=False
    )
    assert result['has_match'], "alpha_inv should be found at depth=3"
    assert result['S_min'] <= 2.0, (
        f"alpha_inv should have S_min <= 2.0, got {result['S_min']}"
    )


# ===================================================================
# TEST 2 — mp_me finds itself
# ===================================================================

def test_mp_me_finds_itself():
    """
    mp_me = 1836.15267343 must be found — it is a terminal atom.
    """
    result = symsearch(
        target=1836.15267343,
        sigma=0.001,
        n=N_TEST,
        depth=3,
        seed=SEED,
        verbose=False
    )
    assert result['has_match'], "mp_me should be found at depth=3"
    assert result['S_min'] <= 2.0, (
        f"mp_me should have S_min <= 2.0, got {result['S_min']}"
    )


# ===================================================================
# TEST 3 — pi finds itself
# ===================================================================

def test_pi_finds_itself():
    """
    pi = 3.14159... must be found — it is a terminal atom.
    """
    result = symsearch(
        target=math.pi,
        sigma=1e-6,
        n=N_TEST,
        depth=3,
        seed=SEED,
        verbose=False
    )
    assert result['has_match'], "pi should be found at depth=3"
    assert result['S_min'] <= 2.0, (
        f"pi should have S_min <= 2.0, got {result['S_min']}"
    )


# ===================================================================
# TEST 4 — impossible target returns no match
# ===================================================================

def test_impossible_target_no_match():
    """
    A number with extremely tight sigma should produce no match.
    Uses a value far from any constant combination.
    """
    result = symsearch(
        target=0.123456789123456,
        sigma=1e-12,
        n=N_TEST,
        depth=3,
        seed=SEED,
        verbose=False
    )
    assert not result['has_match'], (
        "An impossible target with tight sigma should return no match"
    )


# ===================================================================
# TEST 5 — result structure is correct
# ===================================================================

def test_result_structure():
    """
    The result dictionary must contain all required keys.
    """
    result = symsearch(
        target=math.pi,
        sigma=1e-4,
        n=N_TEST,
        depth=3,
        seed=SEED,
        verbose=False
    )
    required_keys = [
        'target', 'sigma', 'n_valid', 'n_matches',
        'p_empirical', 'S_min', 'has_match',
        'top_results', 'elapsed_s', 'parameters'
    ]
    for key in required_keys:
        assert key in result, f"Missing key in result: {key}"


# ===================================================================
# TEST 6 — top_results are sorted by S (simplest first)
# ===================================================================

def test_results_sorted_by_complexity():
    """
    top_results must be sorted by S ascending (simplest first).
    """
    result = symsearch(
        target=137.035999084,
        sigma=0.01,
        n=N_TEST,
        depth=3,
        seed=SEED,
        verbose=False
    )
    if len(result['top_results']) > 1:
        S_values = [r['S'] for r in result['top_results']]
        assert S_values == sorted(S_values), (
            "top_results must be sorted by S ascending"
        )


# ===================================================================
# TEST 7 — exclude parameter works
# ===================================================================

def test_exclude_parameter():
    """
    When alpha_inv is excluded, the result for 137.036 should
    either not be found or have higher S than without exclusion.
    """
    result_with = symsearch(
        target=137.035999084,
        sigma=0.001,
        n=N_TEST,
        depth=3,
        seed=SEED,
        verbose=False
    )
    result_without = symsearch(
        target=137.035999084,
        sigma=0.001,
        n=N_TEST,
        depth=3,
        seed=SEED,
        verbose=False,
        exclude=['alpha_inv', 'alpha']
    )
    # Without the constant itself, S_min should be higher or no match
    if result_without['has_match']:
        assert result_without['S_min'] >= result_with['S_min'], (
            "Excluding alpha_inv should increase S_min for alpha_inv target"
        )


# ===================================================================
# TEST 8 — n_valid is positive
# ===================================================================

def test_n_valid_positive():
    """
    n_valid must be positive — the grammar must generate valid expressions.
    """
    result = symsearch(
        target=1.0,
        sigma=0.1,
        n=N_TEST,
        depth=3,
        seed=SEED,
        verbose=False
    )
    assert result['n_valid'] > 0, "n_valid must be positive"


# ===================================================================
# TEST 9 — inverse search works for a known formula
# ===================================================================

def test_inverse_search_pi_OLD():
    """
    The formula 'pi' should identify math.pi as the closest constant.
    """
    results = inverse_search('pi', top_n=5, rel_tol=0.01)
    assert len(results) > 0, "inverse_search('pi') should return results"
    # The closest match should be pi itself or a pi-related constant
    values = [r[1] for r in results]
    assert any(abs(v - math.pi) / math.pi < 0.001 for v in values), (
        "inverse_search('pi') should find a value close to pi"
    )


# ===================================================================
# TEST 10 — safe_eval_formula handles errors gracefully
# ===================================================================

def test_safe_eval_formula_valid():
    """safe_eval_formula should evaluate 'pi' correctly."""
    val = safe_eval_formula('pi')
    assert val is not None
    assert abs(val - math.pi) < 1e-10


def test_safe_eval_formula_invalid():
    """safe_eval_formula should return None for invalid input."""
    val = safe_eval_formula('not_a_formula_xyz_999')
    assert val is None


# ===================================================================
# TEST 11 — reproducibility with same seed
# ===================================================================

def test_reproducibility():
    """
    Two runs with the same seed must return identical results.
    """
    result_1 = symsearch(
        target=math.pi,
        sigma=1e-4,
        n=N_TEST,
        depth=3,
        seed=SEED,
        verbose=False
    )
    result_2 = symsearch(
        target=math.pi,
        sigma=1e-4,
        n=N_TEST,
        depth=3,
        seed=SEED,
        verbose=False
    )
    assert result_1['n_matches'] == result_2['n_matches'], (
        "Same seed must give same number of matches"
    )
    assert result_1['S_min'] == result_2['S_min'], (
        "Same seed must give same S_min"
    )


# ===================================================================
# TEST 12 — p_empirical is between 0 and 1
# ===================================================================

def test_p_empirical_range():
    """p_empirical must be in [0, 1]."""
    result = symsearch(
        target=math.pi,
        sigma=1e-4,
        n=N_TEST,
        depth=3,
        seed=SEED,
        verbose=False
    )
    assert 0.0 <= result['p_empirical'] <= 1.0, (
        f"p_empirical must be in [0,1], got {result['p_empirical']}"
    )


def test_inverse_search_pi():
    """
    The formula 'pi' should identify math.pi as the closest constant.
    inverse_search returns a dict with 'value' and 'top_matches'.
    """
    result = inverse_search('pi', top_n=5, rel_tol=0.01)
    assert result is not None, "inverse_search('pi') should return a result"
    assert 'value' in result, "Result must have 'value' key"
    assert abs(result['value'] - math.pi) < 1e-6, (
        "inverse_search('pi') should evaluate to pi"
    )
    assert len(result['top_matches']) > 0, "top_matches must not be empty"
    best = result['top_matches'][0]
    assert best['rel_diff'] < 1e-6, (
        f"Best match should be pi itself, got rel_diff={best['rel_diff']}"
    )
