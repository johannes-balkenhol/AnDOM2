"""
AnDOM 2.0 — Streamlit web application.

This file is intentionally thin: UI layout and user interaction only.
All domain search logic lives in search/, db/, and batch/ modules.

To add a new search method:
    1. Create search/your_method.py with a run() function
    2. Call it here alongside sequence.run() and structure.run()
    3. Add its results to ensemble.domain_bar_html()
"""
import sys
import os
import streamlit as st
import pandas as pd
from pathlib import Path

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent / 'src'))
sys.path.insert(0, str(Path(__file__).parent))

from config import SCOP_COLORS, SCOP_CLASSES, EXAMPLES, ESMFOLD_MAXLEN
import db.lookup as lookup
from search import sequence as seq_search
from search import structure as str_search
from search.ensemble import domain_bar_html, add_scores

st.set_page_config(page_title="AnDOM 2.0", layout="wide")

# ── sidebar ───────────────────────────────────────────────────────────────────
with st.sidebar:
    st.header("Search Parameters")
    evalue     = st.select_slider("E-value cutoff",
                    options=[1e-30,1e-20,1e-10,1e-5,1e-3,0.01,0.1,1.0], value=1e-3)
    iterations = st.slider("PSI-MMseqs2 iterations", 1, 5, 3)
    max_hits   = st.slider("Max hits shown", 5, 100, 30)
    use_struct = st.toggle("Structural search (ESMFold + Foldseek)", value=True)
    st.divider()
    st.markdown("**SCOP class colours**")
    for k, v in SCOP_CLASSES.items():
        st.markdown(
            f'<span style="background:{SCOP_COLORS[k]};padding:2px 8px;'
            f'border-radius:3px;color:white;font-size:12px">{k}</span> {v}',
            unsafe_allow_html=True,
        )
    st.divider()
    lk = lookup.all_domains()
    st.caption(f"SCOPe 2.08 — {len(lk):,} domains")
    st.caption("CATH50 structural DB")

page = st.tabs(["Search", "Methods"])

# ── search tab ────────────────────────────────────────────────────────────────
with page[0]:
    st.title("AnDOM 2.0 — Structural Domain Finder")
    st.caption(
        "Sequence + Structure ensemble | SCOPe 2.08 + CATH50 | "
        "MMseqs2 + ESMFold + Foldseek | Dandekar Lab Wuerzburg"
    )

    # example buttons
    st.markdown("**Try an example:**")
    ex_cols = st.columns(5)
    for i, (key, ex) in enumerate(EXAMPLES.items()):
        if ex_cols[i].button(ex["label"], use_container_width=True, key=f"btn_{key}"):
            st.session_state["injected_seq"]  = ex["seq"]
            st.session_state["injected_desc"] = ex["desc"]
            st.rerun()

    if "injected_seq" in st.session_state:
        st.info(st.session_state.get("injected_desc", ""))

    seq_input = st.text_area(
        "Paste protein sequence (FASTA or raw):",
        height=130,
        value=st.session_state.get("injected_seq", ""),
    )

    if st.button("Run AnDOM 2.0 Search", type="primary"):
        if not seq_input.strip():
            st.warning("Please paste a sequence or click an example.")
        else:
            lines     = seq_input.strip().splitlines()
            fasta     = seq_input if lines[0].startswith(">") else f">query\n{seq_input}"
            clean_seq = "".join(l.strip() for l in lines if not l.startswith(">"))
            n         = len(clean_seq)

            use_struct_run = use_struct and n <= ESMFOLD_MAXLEN
            if use_struct and n > ESMFOLD_MAXLEN:
                st.warning(
                    f"Sequence is {n} aa — ESMFold public API limit is "
                    f"{ESMFOLD_MAXLEN} aa. Structural search disabled."
                )

            col_seq, col_str = st.columns(2)
            df_seq = df_str = None

            with col_seq:
                with st.spinner("Sequence search — SCOPe 2.08..."):
                    df_seq, err = seq_search.run(
                        fasta, evalue=evalue, iterations=iterations
                    )
                if err:
                    st.error(f"Sequence search failed: {err}")
                elif df_seq is None or len(df_seq) == 0:
                    st.warning("No sequence hits.")
                else:
                    df_seq = df_seq.head(max_hits)
                    st.success(f"Sequence: {len(df_seq)} SCOPe hits")

            if use_struct_run:
                with col_str:
                    with st.spinner("Structure — ESMFold + Foldseek CATH50..."):
                        df_str, err2 = str_search.run(clean_seq)
                    if err2:
                        st.error(f"Structural search failed: {err2}")
                    elif df_str is None or len(df_str) == 0:
                        st.warning("No structural hits.")
                    else:
                        df_str = df_str.head(max_hits)
                        st.success(f"Structure: {len(df_str)} CATH hits")

            # add confidence scores
            df_seq, df_str = add_scores(df_seq, df_str)

            # domain bar
            st.subheader("Domain Architecture Map")
            bar_len = max(
                df_seq["qend"].max() if df_seq is not None and len(df_seq) > 0 else 1,
                df_str["qend"].max() if df_str is not None and len(df_str) > 0 else 1,
                n,
            )
            st.markdown(
                '<p style="font-size:11px;color:#888">'
                "Top: SCOPe sequence hits (coloured by SCOP class) &nbsp;|&nbsp; "
                f"Bottom: CATH structural hits (dark) &nbsp;|&nbsp; "
                f"Hover for details &nbsp;|&nbsp; 1 to {bar_len} aa</p>"
                + domain_bar_html(df_seq, df_str, bar_len),
                unsafe_allow_html=True,
            )

            # results tables
            tab1, tab2 = st.tabs(["SCOPe sequence hits", "CATH structural hits"])

            with tab1:
                if df_seq is not None and len(df_seq) > 0:
                    disp = df_seq[[
                        "target","sccs","class_name","description",
                        "organism","evalue","bits","qstart","qend","pident","PDB link"
                    ]].copy()
                    disp["evalue"] = disp["evalue"].apply(lambda x: f"{float(x):.2e}")
                    disp["pident"] = disp["pident"].apply(lambda x: f"{float(x):.1f}%")
                    st.dataframe(
                        disp, use_container_width=True,
                        column_config={"PDB link": st.column_config.LinkColumn("PDB link")},
                    )
                    if os.path.exists("seq_results.tsv"):
                        st.download_button(
                            "Download SCOPe TSV",
                            open("seq_results.tsv").read(),
                            "AnDOM_scope.tsv",
                        )
                else:
                    st.info("No sequence hits.")

            with tab2:
                if df_str is not None and len(df_str) > 0:
                    disp2 = df_str[[
                        "target","evalue","bits","lddt","qstart","qend","PDB link"
                    ]].copy()
                    disp2["evalue"] = disp2["evalue"].apply(lambda x: f"{float(x):.2e}")
                    disp2["lddt"]   = disp2["lddt"].apply(lambda x: f"{float(x):.3f}")
                    st.dataframe(
                        disp2, use_container_width=True,
                        column_config={"PDB link": st.column_config.LinkColumn("PDB link")},
                    )
                    if os.path.exists("struct_results.tsv"):
                        st.download_button(
                            "Download CATH TSV",
                            open("struct_results.tsv").read(),
                            "AnDOM_cath.tsv",
                        )
                else:
                    st.info("No structural hits or structural search disabled.")

    st.divider()
    st.caption(
        "AnDOM 2.0 | SCOPe 2.08 ASTRAL 40% (MMseqs2 PSI) + "
        "CATH50 (ESMFold + Foldseek) | Dandekar Lab Wuerzburg"
    )

# ── methods tab ───────────────────────────────────────────────────────────────
with page[1]:
    st.title("Methods — AnDOM 2.0")
    st.markdown(
        "AnDOM 2.0 updates the original AnDOM server "
        "(Schmidt, Bork & Dandekar, *J Chem Inf Comput Sci* 2002) "
        "with modern tools and a new structural search layer."
    )

    st.subheader("Pipeline comparison: AnDOM 2002 vs AnDOM 2.0")
    svg_path = Path(__file__).parent / "docs" / "AnDOM_old_vs_new_method.svg"
    if svg_path.exists():
        st.markdown(svg_path.read_text(), unsafe_allow_html=True)

    st.subheader("What is new in AnDOM 2.0")
    c1, c2 = st.columns(2)
    with c1:
        st.markdown("""
**Database updates**
- SCOP 1.50 → SCOPe 2.08 (2023)
- ~7,000 → 15,177 domains (+117%)
- ASTRAL compendium maintained within SCOPe

**Search engine**
- PSI-BLAST + IMPALA → MMseqs2
- ~100× faster at equivalent sensitivity
- Configurable iterations and e-value threshold
- Enriched metadata: sccs, fold, organism, PDB link
        """)
    with c2:
        st.markdown("""
**New: structural search layer**
- ESMFold v1 (Meta AI) predicts 3D structure
- Foldseek searches against CATH50
- lDDT confidence per hit
- Recovers domains invisible to sequence methods

**Ensemble output**
- Two-row domain architecture map
- Hover tooltips with full annotation
- Clickable PDB links per hit
- TSV download for both layers
- Local Streamlit (original server offline ~2010)
        """)

    st.subheader("References")
    st.markdown("""
- Schmidt S, Bork P, Dandekar T (2002) *J Chem Inf Comput Sci* 42:405–407
- van Kempen M et al. (2023) *Nat Biotechnol* doi:10.1038/s41587-023-01773-0
- Lin Z et al. (2023) *Science* 379:1123–1130
- Fox NK, Brenner SE, Chandonia JM (2014) *Nucleic Acids Res* 42:D304–309
- Sillitoe I et al. (2021) *Nucleic Acids Res* 49:D266–273
    """)
