import streamlit as st
import subprocess
import pandas as pd
import os
import requests
from pathlib import Path

MMSEQS = "mmseqs"
SCOPE_DB = "scopeSeqDB"
SCOPE_FA = "scopeseq_40.fa"
CATH_DB  = "cathDB"
TMP = "tmp"
ESMFOLD_API = "https://api.esmatlas.com/foldSequence/v1/pdb/"

st.set_page_config(page_title="AnDOM 2.0", layout="wide")

EXAMPLES = {
    "Ex 1": {
        "label": "Ex 1 — Haemoglobin (globin, all-alpha)",
        "seq": "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR",
        "desc": "Human haemoglobin alpha chain (142 aa) — classic all-alpha globin fold, well-characterised SCOP domain."
    },
    "Ex 2": {
        "label": "Ex 2 — HIV-1 protease (aspartyl, all-beta)",
        "seq": "PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF",
        "desc": "HIV-1 protease (99 aa) — all-beta aspartyl protease, major drug target, short enough for both search layers."
    },
    "Ex 3": {
        "label": "Ex 3 — TIM barrel (alpha/beta classic)",
        "seq": "MAPSRKFFVGGNWKMNGRKQSLGELIGTLNAAKVPADTEVVCAPPTAYIDFARQKLSQELGASNEILLGAQNVDLNLPKDKFVVLIIVYNKPKDILFKDNEALENGKPFKQNLKAKDLALAGVTPDKMKDLKAKGISGAFVPNIVNLHSQAPADCLMSKLVAGEFEGNIYMGLKPNPEELAAAKSSKLSELIQAAYATGNQVAFKPLTDAAQKAAQESSGKKSATIFAGQATVEDGDTVYL",
        "desc": "Triosephosphate isomerase (248 aa) — canonical TIM barrel, most common enzyme fold in nature. Sequence-only search (>400 aa for ESMFold)."
    },
    "Ex 4": {
        "label": "Ex 4 — Src SH2+SH3 (multi-domain)",
        "seq": "MGSNKSKPKDASQRRRSLEPAENVHGAGGGAFPASQTPSKPASADGHRGPSAAFAPAAAEKVLFGGFNSSDTVTSPQRAGPLAGGVTTFVALYDYESRTETDLSFKKGERLQIVNNTEGDWWLAHSLSTGQTGYIPSNYVAPSDSIQAEEWYFGKITRRESERLLLNAENPRGTFLVRESETSGY",
        "desc": "Src kinase SH2+SH3 domains (180 aa) — multi-domain signalling protein, tests domain boundary detection across two distinct folds."
    },
    "Ex 5": {
        "label": "Ex 5 — p53 DNA-binding domain (beta-sandwich)",
        "seq": "SVVRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGQMNRRPILTIITLEDSSGKLLGRNSFEVRVCACPGRDRRTEEENLRKKGEVVAPQHL",
        "desc": "p53 tumour suppressor DNA-binding domain (130 aa) — beta-sandwich with functional loops, central to cancer biology."
    },
}

SCOP_COLORS = {
    "a": "#e74c3c", "b": "#3498db", "c": "#2ecc71",
    "d": "#f39c12", "e": "#9b59b6", "f": "#1abc9c", "g": "#e67e22",
}
SCOP_CLASSES = {
    "a": "All alpha", "b": "All beta", "c": "Alpha/beta",
    "d": "Alpha+beta", "e": "Multi-domain", "f": "Membrane", "g": "Small proteins",
}

METHOD_SVG = """<svg width="100%" viewBox="0 0 680 740" xmlns="http://www.w3.org/2000/svg">
<defs>
<marker id="arr2" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
<path d="M2 1L8 5L2 9" fill="none" stroke="#888" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round"/>
</marker>
</defs>
<rect x="40" y="20" width="270" height="36" rx="8" fill="#e8e6df" stroke="#888" stroke-width="0.5"/>
<text font-size="14" font-weight="500" fill="#2c2c2a" x="175" y="38" text-anchor="middle" dominant-baseline="central">AnDOM 2002</text>
<rect x="370" y="20" width="270" height="36" rx="8" fill="#EEEDFE" stroke="#534AB7" stroke-width="0.5"/>
<text font-size="14" font-weight="500" fill="#3C3489" x="505" y="38" text-anchor="middle" dominant-baseline="central">AnDOM 2.0 (2026)</text>
<line x1="340" y1="20" x2="340" y2="720" stroke="#ccc" stroke-width="0.5" stroke-dasharray="4 3"/>
<text font-size="12" fill="#888" x="175" y="82" text-anchor="middle">Input</text>
<text font-size="12" fill="#888" x="505" y="82" text-anchor="middle">Input</text>
<rect x="60" y="90" width="230" height="44" rx="6" fill="#f1efe8" stroke="#888" stroke-width="0.5"/>
<text font-size="14" font-weight="500" fill="#2c2c2a" x="175" y="108" text-anchor="middle" dominant-baseline="central">Protein sequence</text>
<text font-size="12" fill="#5f5e5a" x="175" y="123" text-anchor="middle" dominant-baseline="central">FASTA or raw</text>
<rect x="390" y="90" width="230" height="44" rx="6" fill="#EEEDFE" stroke="#534AB7" stroke-width="0.5"/>
<text font-size="14" font-weight="500" fill="#3C3489" x="505" y="108" text-anchor="middle" dominant-baseline="central">Protein sequence</text>
<text font-size="12" fill="#534AB7" x="505" y="123" text-anchor="middle" dominant-baseline="central">FASTA or raw — any length</text>
<line x1="175" y1="134" x2="175" y2="165" stroke="#888" stroke-width="1.5" marker-end="url(#arr2)"/>
<line x1="505" y1="134" x2="505" y2="165" stroke="#534AB7" stroke-width="1.5" marker-end="url(#arr2)"/>
<text font-size="12" fill="#888" x="175" y="160" text-anchor="middle">Sequence search</text>
<text font-size="12" fill="#534AB7" x="505" y="160" text-anchor="middle">Sequence search</text>
<rect x="60" y="168" width="230" height="56" rx="6" fill="#f1efe8" stroke="#888" stroke-width="0.5"/>
<text font-size="14" font-weight="500" fill="#2c2c2a" x="175" y="188" text-anchor="middle" dominant-baseline="central">PSI-BLAST + IMPALA</text>
<text font-size="12" fill="#5f5e5a" x="175" y="205" text-anchor="middle" dominant-baseline="central">SCOP 1.50 · ~7k domains</text>
<text font-size="12" fill="#5f5e5a" x="175" y="217" text-anchor="middle" dominant-baseline="central">ASTRAL 40% · 10 iterations</text>
<rect x="390" y="168" width="230" height="56" rx="6" fill="#EEEDFE" stroke="#534AB7" stroke-width="0.5"/>
<text font-size="14" font-weight="500" fill="#3C3489" x="505" y="188" text-anchor="middle" dominant-baseline="central">MMseqs2 PSI-search</text>
<text font-size="12" fill="#534AB7" x="505" y="205" text-anchor="middle" dominant-baseline="central">SCOPe 2.08 · 15,177 domains</text>
<text font-size="12" fill="#534AB7" x="505" y="217" text-anchor="middle" dominant-baseline="central">ASTRAL 40% · configurable</text>
<line x1="175" y1="224" x2="175" y2="255" stroke="#888" stroke-width="1.5" marker-end="url(#arr2)"/>
<line x1="505" y1="224" x2="505" y2="255" stroke="#534AB7" stroke-width="1.5" marker-end="url(#arr2)"/>
<text font-size="12" fill="#888" x="175" y="251" text-anchor="middle">Structure search</text>
<text font-size="12" fill="#534AB7" x="505" y="251" text-anchor="middle">Structure search (NEW)</text>
<rect x="60" y="258" width="230" height="44" rx="6" fill="#f1efe8" stroke="#888" stroke-width="0.5"/>
<text font-size="14" font-weight="500" fill="#2c2c2a" x="175" y="276" text-anchor="middle" dominant-baseline="central">Not available</text>
<text font-size="12" fill="#5f5e5a" x="175" y="291" text-anchor="middle" dominant-baseline="central">Sequence similarity only</text>
<rect x="390" y="258" width="230" height="44" rx="6" fill="#E1F5EE" stroke="#0F6E56" stroke-width="0.5"/>
<text font-size="14" font-weight="500" fill="#085041" x="505" y="276" text-anchor="middle" dominant-baseline="central">ESMFold + Foldseek</text>
<text font-size="12" fill="#0F6E56" x="505" y="291" text-anchor="middle" dominant-baseline="central">CATH50 · lDDT confidence</text>
<line x1="175" y1="302" x2="175" y2="333" stroke="#888" stroke-width="1.5" marker-end="url(#arr2)"/>
<line x1="505" y1="302" x2="505" y2="333" stroke="#534AB7" stroke-width="1.5" marker-end="url(#arr2)"/>
<text font-size="12" fill="#888" x="175" y="329" text-anchor="middle">Output</text>
<text font-size="12" fill="#534AB7" x="505" y="329" text-anchor="middle">Ensemble output</text>
<rect x="60" y="336" width="230" height="56" rx="6" fill="#f1efe8" stroke="#888" stroke-width="0.5"/>
<text font-size="14" font-weight="500" fill="#2c2c2a" x="175" y="356" text-anchor="middle" dominant-baseline="central">Domain map</text>
<text font-size="12" fill="#5f5e5a" x="175" y="373" text-anchor="middle" dominant-baseline="central">SCOP class · e-value</text>
<text font-size="12" fill="#5f5e5a" x="175" y="385" text-anchor="middle" dominant-baseline="central">Web server (offline ~2010)</text>
<rect x="390" y="336" width="230" height="56" rx="6" fill="#EEEDFE" stroke="#534AB7" stroke-width="0.5"/>
<text font-size="14" font-weight="500" fill="#3C3489" x="505" y="356" text-anchor="middle" dominant-baseline="central">Ensemble domain map</text>
<text font-size="12" fill="#534AB7" x="505" y="373" text-anchor="middle" dominant-baseline="central">Seq + struct · metadata · links</text>
<text font-size="12" fill="#534AB7" x="505" y="385" text-anchor="middle" dominant-baseline="central">Local Streamlit · TSV download</text>
<rect x="40" y="428" width="600" height="280" rx="10" fill="none" stroke="#ccc" stroke-width="0.5"/>
<rect x="40" y="428" width="600" height="32" rx="8" fill="#f1efe8" stroke="#888" stroke-width="0.5"/>
<text font-size="14" font-weight="500" fill="#2c2c2a" x="160" y="444" text-anchor="middle" dominant-baseline="central">Component</text>
<text font-size="14" font-weight="500" fill="#2c2c2a" x="370" y="444" text-anchor="middle" dominant-baseline="central">AnDOM 2002</text>
<text font-size="14" font-weight="500" fill="#2c2c2a" x="550" y="444" text-anchor="middle" dominant-baseline="central">AnDOM 2.0</text>
<line x1="240" y1="428" x2="240" y2="708" stroke="#ccc" stroke-width="0.5"/>
<line x1="420" y1="428" x2="420" y2="708" stroke="#ccc" stroke-width="0.5"/>
<line x1="40" y1="460" x2="640" y2="460" stroke="#ccc" stroke-width="0.5"/>
<text font-size="12" fill="#444" x="60" y="476" dominant-baseline="central">Domain DB</text>
<text font-size="12" fill="#888" x="330" y="476" text-anchor="middle" dominant-baseline="central">SCOP 1.50</text>
<text font-size="12" fill="#3C3489" x="530" y="476" text-anchor="middle" dominant-baseline="central">SCOPe 2.08</text>
<line x1="40" y1="492" x2="640" y2="492" stroke="#ccc" stroke-width="0.5"/>
<text font-size="12" fill="#444" x="60" y="508" dominant-baseline="central">Domain count</text>
<text font-size="12" fill="#888" x="330" y="508" text-anchor="middle" dominant-baseline="central">~7,000</text>
<text font-size="12" fill="#3C3489" x="530" y="508" text-anchor="middle" dominant-baseline="central">15,177 (+117%)</text>
<line x1="40" y1="524" x2="640" y2="524" stroke="#ccc" stroke-width="0.5"/>
<text font-size="12" fill="#444" x="60" y="540" dominant-baseline="central">Search engine</text>
<text font-size="12" fill="#888" x="330" y="540" text-anchor="middle" dominant-baseline="central">PSI-BLAST + IMPALA</text>
<text font-size="12" fill="#3C3489" x="530" y="540" text-anchor="middle" dominant-baseline="central">MMseqs2 (100x faster)</text>
<line x1="40" y1="556" x2="640" y2="556" stroke="#ccc" stroke-width="0.5"/>
<text font-size="12" fill="#444" x="60" y="572" dominant-baseline="central">Structure search</text>
<text font-size="12" fill="#888" x="330" y="572" text-anchor="middle" dominant-baseline="central">None</text>
<text font-size="12" fill="#0F6E56" x="530" y="572" text-anchor="middle" dominant-baseline="central">ESMFold + Foldseek/CATH50</text>
<line x1="40" y1="588" x2="640" y2="588" stroke="#ccc" stroke-width="0.5"/>
<text font-size="12" fill="#444" x="60" y="604" dominant-baseline="central">Hit metadata</text>
<text font-size="12" fill="#888" x="330" y="604" text-anchor="middle" dominant-baseline="central">Domain ID only</text>
<text font-size="12" fill="#3C3489" x="530" y="604" text-anchor="middle" dominant-baseline="central">sccs, fold, organism, PDB link</text>
<line x1="40" y1="620" x2="640" y2="620" stroke="#ccc" stroke-width="0.5"/>
<text font-size="12" fill="#444" x="60" y="636" dominant-baseline="central">Dark proteome</text>
<text font-size="12" fill="#888" x="330" y="636" text-anchor="middle" dominant-baseline="central">Fails (sequence only)</text>
<text font-size="12" fill="#0F6E56" x="530" y="636" text-anchor="middle" dominant-baseline="central">Structural layer recovers hits</text>
<line x1="40" y1="652" x2="640" y2="652" stroke="#ccc" stroke-width="0.5"/>
<text font-size="12" fill="#444" x="60" y="668" dominant-baseline="central">Availability</text>
<text font-size="12" fill="#888" x="330" y="668" text-anchor="middle" dominant-baseline="central">Offline since ~2010</text>
<text font-size="12" fill="#3C3489" x="530" y="668" text-anchor="middle" dominant-baseline="central">Local deploy · open source</text>
<line x1="40" y1="684" x2="640" y2="684" stroke="#ccc" stroke-width="0.5"/>
<text font-size="12" fill="#444" x="60" y="700" dominant-baseline="central">Benchmark</text>
<text font-size="12" fill="#888" x="330" y="700" text-anchor="middle" dominant-baseline="central">Swiss-Prot 38.0 (77%)</text>
<text font-size="12" fill="#3C3489" x="530" y="700" text-anchor="middle" dominant-baseline="central">UniProt current (TBD)</text>
</svg>"""

@st.cache_data
def build_scop_lookup(fa_file):
    lookup = {}
    with open(fa_file) as f:
        for line in f:
            if line.startswith(">"):
                parts = line.strip().lstrip(">").split(None, 3)
                if len(parts) >= 2:
                    sid  = parts[0]
                    sccs = parts[1]
                    desc = parts[2] if len(parts) > 2 else ""
                    org  = line[line.index("{")+1:line.index("}")] if "{" in line and "}" in line else ""
                    lookup[sid] = {"sccs": sccs, "cls": sccs[0] if sccs else "?", "desc": desc, "org": org}
    return lookup

lookup = build_scop_lookup(SCOPE_FA)

def scop_url(sid):
    try:
        return f"https://www.rcsb.org/structure/{sid[1:5].lower()}"
    except:
        return "https://scop.berkeley.edu"

def cath_url(target):
    try:
        return f"https://www.rcsb.org/structure/{target[:4].lower()}"
    except:
        return "https://www.cathdb.info"

with st.sidebar:
    st.header("Search Parameters")
    evalue     = st.select_slider("E-value cutoff",
                    options=[1e-30,1e-20,1e-10,1e-5,1e-3,0.01,0.1,1.0], value=1e-3)
    iterations = st.slider("PSI-MMseqs iterations", 1, 5, 3)
    max_hits   = st.slider("Max hits shown", 5, 100, 30)
    use_struct = st.toggle("Enable structural search (ESMFold + Foldseek)", value=True)
    st.divider()
    st.markdown("**SCOP class colours**")
    for k, v in SCOP_CLASSES.items():
        st.markdown(
            f'<span style="background:{SCOP_COLORS[k]};padding:2px 8px;border-radius:3px;color:white;font-size:12px">{k}</span> {v}',
            unsafe_allow_html=True)
    st.divider()
    st.caption(f"SCOPe 2.08 — {len(lookup):,} domains")
    st.caption("CATH50 structural DB")

page = st.tabs(["Search", "Methods"])

# ── SEARCH TAB ────────────────────────────────────────────────────────────────
with page[0]:
    st.title("AnDOM 2.0 — Structural Domain Finder")
    st.caption("Sequence + Structure ensemble search | SCOPe 2.08 + CATH50 | MMseqs2 + ESMFold + Foldseek | Dandekar Lab Wuerzburg")

    st.markdown("**Try an example:**")
    ex_cols = st.columns(5)
    for i, (key, ex) in enumerate(EXAMPLES.items()):
        if ex_cols[i].button(ex["label"], use_container_width=True, key=f"btn_{key}"):
            st.session_state["injected_seq"] = ex["seq"]
            st.session_state["injected_desc"] = ex["desc"]
            st.rerun()

    if "injected_seq" in st.session_state:
        st.info(st.session_state.get("injected_desc", ""))

    seq_input = st.text_area(
        "Paste protein sequence (FASTA or raw):",
        height=130,
        value=st.session_state.get("injected_seq", "")
    )

    def run_sequence_search(fasta, evalue, iterations):
        import glob as _glob
        for _f in _glob.glob("queryDB*") + _glob.glob("resultDB*") + ["seq_results.tsv"]:
            try: os.remove(_f)
            except: pass
        Path("query_input.fasta").write_text(fasta)
        for cmd in [
            f"{MMSEQS} createdb query_input.fasta queryDB -v 0",
            f"{MMSEQS} search queryDB {SCOPE_DB} resultDB {TMP} --num-iterations {iterations} -e {evalue} -v 0 --threads 8",
            f"{MMSEQS} convertalis queryDB {SCOPE_DB} resultDB seq_results.tsv --format-output query,target,evalue,bits,qstart,qend,tstart,tend,pident -v 0",
        ]:
            r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            if r.returncode != 0:
                return None, r.stderr
        return "seq_results.tsv", None

    def run_structural_search(clean_seq):
        try:
            r = requests.post(ESMFOLD_API,
                headers={"Content-Type": "application/x-www-form-urlencoded"},
                data=clean_seq, timeout=90)
            if not r.text.startswith("HEADER"):
                return None, "ESMFold API error: " + r.text[:200]
            Path("query_struct.pdb").write_text(r.text)
        except Exception as e:
            return None, str(e)
        r2 = subprocess.run(
            f"foldseek easy-search query_struct.pdb {CATH_DB} struct_results.tsv {TMP} "
            f"--format-output query,target,evalue,bits,lddt,qstart,qend -e 0.001 -v 0",
            shell=True, capture_output=True, text=True)
        if r2.returncode != 0:
            return None, r2.stderr
        return "struct_results.tsv", None

    def domain_bar(df_seq, df_str, seq_len):
        bar = '<div style="position:relative;height:64px;background:#f0f2f6;border-radius:8px;width:100%;margin-bottom:6px">'
        if df_seq is not None:
            for _, row in df_seq.iterrows():
                cls = lookup.get(row["target"], {}).get("cls", "?")
                color = SCOP_COLORS.get(cls, "#888")
                left = (row["qstart"] / seq_len) * 100
                width = max(((row["qend"] - row["qstart"]) / seq_len) * 100, 2)
                tip = f"{row['target']} | {lookup.get(row['target'],{}).get('sccs','?')} | e={float(row['evalue']):.1e}"
                bar += f'<div title="{tip}" style="position:absolute;top:4px;left:{left:.1f}%;width:{width:.1f}%;height:24px;background:{color};opacity:0.9;border-radius:4px;border:1px solid rgba(255,255,255,0.6)"></div>'
        if df_str is not None:
            for _, row in df_str.iterrows():
                left = (row["qstart"] / seq_len) * 100
                width = max(((row["qend"] - row["qstart"]) / seq_len) * 100, 2)
                tip = f"{row['target']} | lddt={float(row['lddt']):.2f} | e={float(row['evalue']):.1e} | CATH"
                bar += f'<div title="{tip}" style="position:absolute;top:34px;left:{left:.1f}%;width:{width:.1f}%;height:24px;background:#2c3e50;opacity:0.75;border-radius:4px;border:1px solid rgba(255,255,255,0.4)"></div>'
        bar += "</div>"
        return bar

    if st.button("Run AnDOM 2.0 Search", type="primary"):
        if not seq_input.strip():
            st.warning("Please paste a sequence or click an example above.")
        else:
            lines     = seq_input.strip().splitlines()
            fasta     = seq_input if lines[0].startswith(">") else ">query\n" + seq_input
            clean_seq = "".join([l.strip() for l in lines if not l.startswith(">")])
            n = len(clean_seq)
            use_struct_run = use_struct and n <= 400
            if use_struct and n > 400:
                st.warning(f"Sequence is {n} aa — ESMFold API limit is 400 aa. Structural search disabled.")

            col_seq, col_str = st.columns(2)
            df_seq = df_str = sf = sf2 = None

            with col_seq:
                with st.spinner("Sequence search — SCOPe 2.08..."):
                    sf, err = run_sequence_search(fasta, evalue, iterations)
                if err:
                    st.error(f"Sequence search failed: {err}")
                elif not sf or os.path.getsize(sf) == 0:
                    st.warning("No sequence hits.")
                else:
                    df_seq = pd.read_csv(sf, sep="\t",
                        names=["query","target","evalue","bits","qstart","qend","tstart","tend","pident"])
                    df_seq = df_seq.sort_values("evalue").head(max_hits)
                    df_seq["sccs"]        = df_seq["target"].map(lambda x: lookup.get(x,{}).get("sccs","?"))
                    df_seq["cls"]         = df_seq["target"].map(lambda x: lookup.get(x,{}).get("cls","?"))
                    df_seq["class_name"]  = df_seq["cls"].map(lambda x: SCOP_CLASSES.get(x,"?"))
                    df_seq["description"] = df_seq["target"].map(lambda x: lookup.get(x,{}).get("desc",""))
                    df_seq["organism"]    = df_seq["target"].map(lambda x: lookup.get(x,{}).get("org",""))
                    df_seq["PDB link"]    = df_seq["target"].map(scop_url)
                    st.success(f"Sequence: {len(df_seq)} SCOPe hits")

            if use_struct_run:
                with col_str:
                    with st.spinner("Structure — ESMFold + Foldseek CATH50..."):
                        sf2, err2 = run_structural_search(clean_seq)
                    if err2:
                        st.error(f"Structural search failed: {err2}")
                    elif not sf2 or os.path.getsize(sf2) == 0:
                        st.warning("No structural hits.")
                    else:
                        df_str = pd.read_csv(sf2, sep="\t",
                            names=["query","target","evalue","bits","lddt","qstart","qend"])
                        df_str = df_str.sort_values("evalue").head(max_hits)
                        df_str["PDB link"] = df_str["target"].map(cath_url)
                        st.success(f"Structure: {len(df_str)} CATH hits")

            st.subheader("Domain Architecture Map")
            bar_len = max(
                df_seq["qend"].max() if df_seq is not None and len(df_seq)>0 else 1,
                df_str["qend"].max() if df_str is not None and len(df_str)>0 else 1,
                n,
            )
            st.markdown(
                '<p style="font-size:11px;color:#888">Top row: SCOPe sequence hits (coloured by SCOP class) &nbsp;|&nbsp; '
                f'Bottom row: CATH structural hits (dark) &nbsp;|&nbsp; Hover for details &nbsp;|&nbsp; 1 to {bar_len} aa</p>'
                + domain_bar(df_seq, df_str, bar_len), unsafe_allow_html=True)

            tab1, tab2 = st.tabs(["SCOPe sequence hits", "CATH structural hits"])
            with tab1:
                if df_seq is not None:
                    disp = df_seq[["target","sccs","class_name","description","organism","evalue","bits","qstart","qend","pident","PDB link"]].copy()
                    disp["evalue"] = disp["evalue"].apply(lambda x: f"{float(x):.2e}")
                    disp["pident"] = disp["pident"].apply(lambda x: f"{float(x):.1f}%")
                    st.dataframe(disp, use_container_width=True,
                        column_config={"PDB link": st.column_config.LinkColumn("PDB link")})
                    if sf:
                        st.download_button("Download SCOPe TSV", open(sf).read(), "AnDOM_scope.tsv")
                else:
                    st.info("No sequence hits.")
            with tab2:
                if df_str is not None:
                    disp2 = df_str[["target","evalue","bits","lddt","qstart","qend","PDB link"]].copy()
                    disp2["evalue"] = disp2["evalue"].apply(lambda x: f"{float(x):.2e}")
                    disp2["lddt"]   = disp2["lddt"].apply(lambda x: f"{float(x):.3f}")
                    st.dataframe(disp2, use_container_width=True,
                        column_config={"PDB link": st.column_config.LinkColumn("PDB link")})
                    if sf2:
                        st.download_button("Download CATH TSV", open(sf2).read(), "AnDOM_cath.tsv")
                else:
                    st.info("No structural hits or structural search disabled.")

    st.divider()
    st.caption("AnDOM 2.0 | SCOPe 2.08 ASTRAL 40% (MMseqs2 PSI) + CATH50 (ESMFold + Foldseek) | Dandekar Lab Wuerzburg")

# ── METHODS TAB ───────────────────────────────────────────────────────────────
with page[1]:
    st.title("Methods — AnDOM 2.0")
    st.markdown("""
AnDOM 2.0 is an updated and extended version of the original AnDOM server
(Schmidt, Bork & Dandekar, *J Chem Inf Comput Sci* 2002).
It combines sequence-based and structure-based domain assignment into a single ensemble pipeline.
    """)

    st.subheader("Pipeline comparison: AnDOM 2002 vs AnDOM 2.0")
    st.markdown(METHOD_SVG, unsafe_allow_html=True)

    st.subheader("What is new in AnDOM 2.0")
    col1, col2 = st.columns(2)
    with col1:
        st.markdown("""
**Database updates**
- SCOP 1.50 (2001) replaced by SCOPe 2.08 (2023)
- Domain count: ~7,000 → 15,177 (+117%)
- ASTRAL compendium now maintained within SCOPe at Berkeley
- 40% sequence identity cutoff preserved for comparability

**Search engine**
- PSI-BLAST + IMPALA replaced by MMseqs2
- ~100x faster at equivalent sensitivity
- Configurable iteration count and e-value threshold
- Result metadata enriched: sccs, fold description, organism, PDB link
        """)
    with col2:
        st.markdown("""
**New: structural search layer**
- ESMFold v1 (Meta AI) predicts 3D structure from sequence
- Foldseek searches predicted structure against CATH50
- Returns lDDT confidence score per hit
- Recovers domains where sequence similarity has decayed below detection
- Critical for dark proteome proteins

**Ensemble output**
- Two-row domain architecture map (sequence + structure)
- Hover tooltips with full annotation
- Clickable PDB links per hit
- TSV download for both layers
- Local deployment via Streamlit (original server offline ~2010)
        """)

    st.subheader("References")
    st.markdown("""
- Schmidt S, Bork P, Dandekar T (2002) A versatile structural domain analysis server using profile weight matrices. *J Chem Inf Comput Sci* 42:405–407
- van Kempen M et al. (2023) Fast and accurate protein structure search with Foldseek. *Nat Biotechnol* doi:10.1038/s41587-023-01773-0
- Lin Z et al. (2023) Evolutionary-scale prediction of atomic-level protein structure with a language model. *Science* 379:1123–1130
- Fox NK, Brenner SE, Chandonia JM (2014) SCOPe: Structural Classification of Proteins extended. *Nucleic Acids Res* 42:D304–309
- Sillitoe I et al. (2021) CATH: increased structural coverage of functional space. *Nucleic Acids Res* 49:D266–273
    """)
