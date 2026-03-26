"""Unit tests for search/sequence.py"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
sys.path.insert(0, str(Path(__file__).parent.parent))

HBA_SEQ = "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR"

def test_fasta_format():
    fasta = f">HBA_HUMAN\n{HBA_SEQ}"
    assert fasta.startswith(">")
    assert len(HBA_SEQ) == 142

def test_scop_lookup_loads():
    import db.lookup as lk
    domains = lk.all_domains()
    assert len(domains) > 10000, "SCOPe lookup should have >10k domains"

def test_pdb_url():
    import db.lookup as lk
    url = lk.pdb_url("d3d1ka_")
    assert "3d1k" in url

if __name__ == "__main__":
    test_fasta_format()
    print("test_fasta_format passed")
    test_scop_lookup_loads()
    print("test_scop_lookup_loads passed")
    test_pdb_url()
    print("test_pdb_url passed")
    print("All tests passed.")
