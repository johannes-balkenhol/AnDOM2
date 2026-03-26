"""Unit tests for search/ensemble.py"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
sys.path.insert(0, str(Path(__file__).parent.parent))
import pandas as pd

def test_domain_bar_empty():
    from search.ensemble import domain_bar_html
    html = domain_bar_html(None, None, 100)
    assert "<div" in html

def test_evalue_to_score():
    from search.ensemble import _evalue_to_score
    assert _evalue_to_score(1e-30) == 1.0
    assert _evalue_to_score(1.0) == 0.0
    assert 0 < _evalue_to_score(1e-10) < 1.0

if __name__ == "__main__":
    test_domain_bar_empty()
    print("test_domain_bar_empty passed")
    test_evalue_to_score()
    print("test_evalue_to_score passed")
    print("All tests passed.")
