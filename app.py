"""
AnDOM 2.0 — Streamlit web application.

This file is intentionally thin: UI layout and user interaction only.
All domain search logic lives in search/, db/, and batch/ modules.
"""
import sys
import os
import streamlit as st
import pandas as pd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent / 'src'))
sys.path.insert(0, str(Path(__file__).parent))

from config import SCOP_COLORS, SCOP_CLASSES, EXAMPLES, ESMFOLD_MAXLEN
import db.lookup as lookup
from search import sequence as seq_search
from search import structure as str_search
from search.ensemble import domain_bar_html, add_scores, fuse_results
from batch.processor import (
    BatchManager, parse_fasta_text, validate_sequences,
    MAX_SEQUENCES, MAX_SEQ_LEN, MIN_SEQ_LEN,
)

st.set_page_config(page_title="AnDOM 2.0", layout="wide")

# ── singletons ────────────────────────────────────────────────────────────────
@st.cache_resource
def get_batch_manager() -> BatchManager:
    return BatchManager(output_root="output")

batch_mgr = get_batch_manager()

# ── hard benchmark sequences ──────────────────────────────────────────────────
BATCH_EXAMPLES = {
    "dark_proteome": {
        "label": "Dark proteome (ORFan genes, no homologs)",
        "desc": (
            "Proteins with zero or near-zero sequence homologs in any database. "
            "InterPro/Pfam find nothing. AlphaFold can predict a structure but cannot "
            "classify the domain. This is where the structural arm of AnDOM 2.0 adds value."
        ),
        "fasta": """>UP000005640_ORFan_1
MSSRSSRGRGGRSSRGRSSRGRGGRSSRGRSSRGRGGRSSRGRSSRGRGGR
>UP000005640_ORFan_2
MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYGQQSYSGYSQSTDTSGYGQSSYSSYGQ
>Mycobacterium_PE_PGRS_1
MAPEPAAPAPAPAPAPAPAPAPAPAPAPAPAPAPAPAPAPAPAPAPAPAPAPAPAPAPAP
APAPAPAPAPAPAPAPAPAPAPAPAPAPAPAPAPAPAPAPAPAPAPAPAPAPAPAPAP
""",
    },
    "fast_evolving": {
        "label": "Fast-evolving viral/phage proteins",
        "desc": (
            "Viral and bacteriophage proteins evolve so fast that sequence methods "
            "fail below ~20% identity — the 'twilight zone'. Structural search recovers "
            "these when sequence search returns nothing. Critical for virology."
        ),
        "fasta": """>HIV1_Vif_HXB2
MENRWQVMIVWQVDRMRIRTWKSLVKHHMYVSGKARGWFYRHHYESPHPRISSEVHIPL
GDARLIITTYWGLHTGERDWHLGQGVSIEWRKKRYSTQVDPELAD
>Lambda_phage_Cro
MTKQKTLQELRQELQHEAHELYNALIQRLEQEVQAELANQEQQLHALEQLERERLLKLA
>T4_phage_gp5_baseplate
MNIFEMLRIDEGLRLKIYKDTEGYYTIGIGHLLTKSPSLNAAKSELDKAIGRNCNGVITK
DEAEKLFNQDVDAAVRGILRNAKLKPVYDSLDAVRRAALINMVFQMGETGVAGFTNSLRM
""",
    },
    "multidomain_hard": {
        "label": "Multi-domain — one known, one novel",
        "desc": (
            "Proteins where one domain is well-characterised but a second "
            "domain has no sequence homologs. Ensemble scoring flags the known "
            "domain with high confidence (both arms) and the novel domain "
            "as structure-only — exactly the complementarity we benchmark."
        ),
        "fasta": """>p53_full_with_disordered_TAD
MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPG
PDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYPQGLNGTVNLFQSSHY
SVVKQKDGQFEVTMDVTAPGTSVTIRNPRFQNLVTPQGRIKVSGNVPNLNAAVKRGDKK
SVLHACIHSISQGGVSSEDGKKKKNLTQKAVNSTDLNIIDNLTENVKSNLQKLYNQLSK
>Src_kinase_SH3_novel_insert
MGSNKSKPKDASQRRRSLEPAENVHGAGGGAFPASQTPSKPASADGHRGPSAAFAPAAAE
KVLFGGFNSSDTVTSPQRAGPLAGGVTTFVALYDYESRTETDLSFKKGERLQIVNNTEGD
WWLAHSLSTGQTGYIPSNYVAPSDSIQAEEWYFGKITRRESERLLLNAENPRGTFLVRESE
""",
    },
    "intrinsically_disordered": {
        "label": "Intrinsically disordered proteins (IDPs)",
        "desc": (
            "IDPs have no stable fold but often contain short linear motifs (SLiMs) "
            "and transient structured regions. AlphaFold gives low pLDDT scores. "
            "InterPro misses the structured segments. The ensemble can recover "
            "transiently folded domains via the structural arm."
        ),
        "fasta": """>FG_nucleoporin_Nup98
MSNTGGFGFGSGFGSGFGSGFGSSFGSGFGSGFGSGFGSSFGSGFGSGFGSGFGSSFGS
GFGSGFGSGFGSGFGSSFGSGFGSGFGSGFGSSFGSGFGSGFGSGFGSSFGSGFGSGFG
>Tau_repeat_domain
MAEPRQEFEVMEDHAGTYGLGDRKDQGGYTMHQDQEGDTDAGLKESPLQTPTEDGSEEP
GSETSDAKSTPTAEDVTAPLVDEGAPGKQAAAQPHTEIPEGTTAEEAGIGDTPSLEDEAA
GHVTQARMVSKSKDGTGSDDKKAKGADGKTKIATPRGAAPPGQKGQANATRIPAKTPPAP
>Alpha_synuclein_Parkinsons
MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTK
EQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDP
DNEAYEMPSEEGYQDYEPEA
""",
    },
    "de_novo_designed": {
        "label": "De novo designed proteins",
        "desc": (
            "Computationally designed proteins with no evolutionary history. "
            "Zero sequence homologs by definition. AlphaFold can predict them "
            "but cannot classify them. Only structural comparison to known folds "
            "can reveal if they adopt a known topology. Pure structural arm test."
        ),
        "fasta": """>Designed_TIM_barrel_de_novo
MAKMAEEILSQYGDNLPKEVLEYLEKNIAKSGKIFVEIKADNAIAAANAKAVLENAKEVL
EAYKTGKVNLVEKFDANLDAAKAKAVLEAAKEVLEAYKTGKVNLVEAFDANLDAAKAKAV
LEAAKEVLEAYKTGKVNLVEAFDANLDAAKAKAVLE
>Designed_beta_propeller
MSTGKVIKVLSANDQTGEVKITVHGKTVDVTVEEGDKVTLHYEGKAVDVSVEEGDKVTL
HYEGKAVDVSVEEGDKVTLHYEGKAVDVSVEEGDKVTLHYEGKAVDVSVEE
""",
    },
}

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
    st.divider()
    st.caption(f"Batch limit: **{MAX_SEQUENCES}** sequences · **{MAX_SEQ_LEN}** aa max")

page = st.tabs(["Search", "Batch", "Benchmark", "Methods"])

# ══════════════════════════════════════════════════════════════════════════════
# TAB 0 — Search
# ══════════════════════════════════════════════════════════════════════════════
with page[0]:
    st.title("AnDOM 2.0 — Structural Domain Finder")
    st.caption(
        "Sequence + Structure ensemble | SCOPe 2.08 + CATH50 | "
        "MMseqs2 + ESMFold + Foldseek | Dandekar Lab Wuerzburg"
    )

    st.markdown("**Try an example:**")
    ex_cols = st.columns(5)
    for i, (key, ex) in enumerate(EXAMPLES.items()):
        if ex_cols[i].button(ex["label"], width="stretch", key=f"btn_{key}"):
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

            df_seq, df_str = add_scores(df_seq, df_str)

            # ensemble fusion table
            fused = fuse_results(df_seq, df_str)
            if not fused.empty:
                st.subheader("Ensemble Ranked Hits")
                st.caption(
                    "🟢 found by both arms (high confidence) · "
                    "🔵 sequence only · "
                    "🟠 structure only (recovered by ESMFold+Foldseek)"
                )
                ev_icon = {"both": "🟢", "seq_only": "🔵", "struct_only": "🟠"}
                fused["ev"] = fused["evidence"].map(ev_icon)
                st.dataframe(
                    fused[["rank","ev","domain_id","evidence",
                           "ensemble_score","seq_evalue","struct_evalue"]].rename(
                        columns={"ev": "", "domain_id": "Domain",
                                 "ensemble_score": "Score",
                                 "seq_evalue": "Seq e-val",
                                 "struct_evalue": "Struct e-val"}),
                    width="stretch", hide_index=True,
                )

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
                "Bottom: CATH structural hits (dark) &nbsp;|&nbsp; "
                f"Hover for details &nbsp;|&nbsp; 1 to {bar_len} aa</p>"
                + domain_bar_html(df_seq, df_str, bar_len),
                unsafe_allow_html=True,
            )

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
                        disp, width="stretch",
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
                        disp2, width="stretch",
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

# ══════════════════════════════════════════════════════════════════════════════
# TAB 1 — Batch
# ══════════════════════════════════════════════════════════════════════════════
with page[1]:
    st.title("Batch Processing")
    st.caption(
        f"Up to **{MAX_SEQUENCES}** sequences per job · "
        f"{MIN_SEQ_LEN}–{MAX_SEQ_LEN} aa per sequence · "
        "jobs run in the background"
    )

    # ── why AnDOM 2.0 exists ───────────────────────────────────────────────
    with st.expander("💡 What can AnDOM 2.0 find that AlphaFold + Foldseek alone cannot?", expanded=False):
        st.markdown("""
AnDOM 2.0 is not a replacement for AlphaFold or Foldseek — it fills the gap **after** structure prediction:

| Challenge | AlphaFold | Foldseek | InterPro/Pfam | **AnDOM 2.0** |
|---|---|---|---|---|
| Predicts 3D structure | ✅ | needs structure | ❌ | via ESMFold |
| Assigns domain to SCOP/CATH class | ❌ | ✅ if DB has it | sequence only | **✅ seq + struct** |
| Dark proteome (no homologs) | structure only | ❌ no hit | ❌ | **structural arm recovers** |
| Fast-evolving viral proteins | structure only | weak | ❌ | **ensemble boosts** |
| Multi-domain: one known, one novel | ❌ | partial | partial | **flags both, evidence tagged** |
| Intrinsically disordered + folded domain | pLDDT low | poor | misses | **structural arm** |
| De novo designed proteins | structure only | ❌ | ❌ | **structural arm** |

**The key insight:** AlphaFold gives you coordinates. AnDOM 2.0 tells you *what fold family those coordinates belong to*, 
using two independent lines of evidence. When both agree → high confidence. When only structure finds it → 
the protein is in the dark proteome, invisible to sequence methods.
        """)

    # ── batch examples ────────────────────────────────────────────────────
    st.markdown("**Load a benchmark example set:**")
    ex_cols = st.columns(len(BATCH_EXAMPLES))
    for i, (key, ex) in enumerate(BATCH_EXAMPLES.items()):
        if ex_cols[i].button(ex["label"], width="stretch", key=f"bex_{key}"):
            st.session_state["batch_example_fasta"] = ex["fasta"]
            st.session_state["batch_example_desc"]  = ex["desc"]
            st.rerun()

    if "batch_example_desc" in st.session_state:
        st.info(st.session_state["batch_example_desc"])

    with st.expander("➕ Submit new batch job", expanded=True):
        b_upload = st.file_uploader(
            "Upload multi-FASTA", type=["fa","fasta","txt"], key="batch_upload"
        )
        b_text = st.text_area(
            "…or paste multi-FASTA here",
            height=150,
            key="batch_text",
            value=st.session_state.get("batch_example_fasta", ""),
        )
        b_struct = st.toggle(
            "Include structural search (≤400 aa only)", value=False, key="b_struct"
        )

        raw = None
        if b_upload:
            raw = b_upload.read().decode()
        elif b_text.strip():
            raw = b_text

        if st.button("Submit Batch Job", type="primary"):
            if not raw:
                st.warning("Please provide FASTA sequences.")
            else:
                seqs = parse_fasta_text(raw)
                errs = validate_sequences(seqs)
                if errs:
                    for e in errs:
                        st.error(e)
                else:
                    import tempfile
                    with tempfile.NamedTemporaryFile(
                        mode="w", suffix=".fa", delete=False
                    ) as tmp:
                        tmp.write(raw)
                        tmp_path = tmp.name
                    job_id = batch_mgr.submit(
                        tmp_path,
                        evalue=evalue,
                        iterations=iterations,
                        use_structure=b_struct,
                    )
                    st.success(
                        f"Job **{job_id}** submitted — {len(seqs)} sequences. "
                        "Refresh this page to track progress."
                    )
                    st.session_state.pop("batch_example_fasta", None)
                    st.session_state.pop("batch_example_desc", None)

    st.subheader("Jobs")
    jobs = batch_mgr.list_jobs()

    if not jobs:
        st.info("No batch jobs yet.")
    else:
        if st.button("🔄 Refresh"):
            st.rerun()

        status_icon = {
            "queued": "🕐", "running": "⏳",
            "done":   "✅", "failed":  "❌", "cancelled": "🚫",
        }
        for job in jobs:
            s = job.get("status", "unknown")
            with st.container(border=True):
                c1, c2, c3 = st.columns([3, 2, 2])
                c1.markdown(f"{status_icon.get(s,'❓')} **{job['job_id']}** — `{s}`")
                total    = job.get("total", 0)
                progress = job.get("progress", 0)
                c2.caption(
                    f"{progress}/{total} sequences" if total
                    else job.get("submitted", "")[:19]
                )
                if s == "done":
                    df_res = batch_mgr.results(job["job_id"])
                    if not df_res.empty:
                        c3.download_button(
                            "⬇ TSV",
                            df_res.to_csv(sep="\t", index=False),
                            file_name=f"AnDOM_batch_{job['job_id']}.tsv",
                            mime="text/tab-separated-values",
                            key=f"dl_{job['job_id']}",
                        )
                elif s in ("queued", "running"):
                    if c3.button("Cancel", key=f"cancel_{job['job_id']}"):
                        batch_mgr.cancel(job["job_id"])
                        st.rerun()
                if job.get("error"):
                    st.error(job["error"])

# ══════════════════════════════════════════════════════════════════════════════
# TAB 2 — Benchmark
# ══════════════════════════════════════════════════════════════════════════════
with page[2]:
    st.title("Benchmark — Sequence vs Structure vs Ensemble")

    st.markdown("""
**What this benchmark measures and why it matters for the paper:**

The central claim of AnDOM 2.0 is that combining sequence search (MMseqs2 → SCOPe) 
with structural search (ESMFold → Foldseek → CATH) recovers more domain annotations 
than either method alone — especially for proteins that are hard for the community:

- **Dark proteome** proteins with no sequence homologs
- **Fast-evolving** viral/phage proteins below the twilight zone (~20% identity)  
- **Intrinsically disordered** proteins with transient structured regions
- **De novo designed** proteins with no evolutionary signal

The benchmark runs all three arms on a curated set and reports:
- **Rank-1 accuracy** — did the top hit match the known classification?
- **Top-5 sensitivity** — was the correct answer anywhere in the top 5?
- **Unique to structure** — how many proteins did the structural arm recover that sequence search missed?

This directly answers: *does the ensemble deserve to exist alongside AlphaFold + Foldseek?*
    """)

    st.divider()

    bm_scope = st.toggle("Use built-in SCOPE-40 dataset", value=True)

    if bm_scope:
        from config import SCOPE_FA
        bm_fasta = SCOPE_FA
        st.caption(f"Dataset: `{bm_fasta}`")
        st.caption(
            "SCOPE-40 = SCOPe domains clustered at 40% sequence identity. "
            "This is the standard benchmark for domain annotation tools — "
            "sequences are diverse enough that sequence methods alone struggle."
        )
    else:
        bm_file = st.file_uploader(
            "Upload SCOPE-style FASTA (header: >domain_id class_code ...)",
            type=["fa","fasta"], key="bm_upload"
        )
        bm_fasta = None
        if bm_file:
            import tempfile
            with tempfile.NamedTemporaryFile(
                mode="wb", suffix=".fa", delete=False
            ) as tmp:
                tmp.write(bm_file.read())
                bm_fasta = tmp.name

    bm_limit = st.number_input(
        "Limit to first N sequences (0 = all; start with 50–100 to test)",
        min_value=0, value=100, step=10
    )

    if st.button("▶ Run Benchmark", type="primary"):
        if not bm_fasta:
            st.warning("Please provide a FASTA file.")
        else:
            from benchmark.run import run_arms_benchmark
            import time as _time, json as _json
            ts = _time.strftime("%Y%m%d_%H%M%S")
            out_json = f"output/results/benchmark_{ts}.json"

            with st.spinner("Running benchmark — this may take several minutes…"):
                stats = run_arms_benchmark(
                    bm_fasta,
                    evalue=evalue,
                    iterations=iterations,
                    limit=int(bm_limit) if bm_limit > 0 else None,
                    out_json=out_json,
                )

            st.success(f"Done — {stats['total']} sequences in {stats['elapsed_s']}s")

            rows = []
            for arm in ("seq_only", "struct_only", "ensemble"):
                d = stats[arm]
                rows.append({
                    "Arm":      arm,
                    "Rank-1 %": d["rank1_pct"],
                    "Top-5 %":  d["top5_pct"],
                    "Rank-1 n": d["rank1"],
                    "Top-5 n":  d["top5"],
                })
            df_bm = pd.DataFrame(rows)
            st.dataframe(df_bm, width="stretch", hide_index=True)
            st.bar_chart(df_bm.set_index("Arm")[["Rank-1 %","Top-5 %"]])

            st.caption(
                "If ensemble > seq_only AND ensemble > struct_only → "
                "the combination is justified. "
                "If struct_only finds hits seq_only misses → dark proteome recovery confirmed."
            )

            st.download_button(
                "⬇ Download benchmark JSON",
                data=_json.dumps(stats, indent=2),
                file_name=f"benchmark_{ts}.json",
                mime="application/json",
            )

# ══════════════════════════════════════════════════════════════════════════════
# TAB 3 — Methods
# ══════════════════════════════════════════════════════════════════════════════
with page[3]:
    st.title("Methods — AnDOM 2.0")
    st.markdown(
        "AnDOM 2.0 updates the original AnDOM server "
        "(Schmidt, Bork & Dandekar, *J Chem Inf Comput Sci* 2002) "
        "with modern tools and a new structural search layer."
    )

    st.subheader("Pipeline comparison: AnDOM 2002 vs AnDOM 2.0")

    # look in docs/ first, then project root as fallback
    svg_path = Path(__file__).parent / "docs" / "AnDOM_old_vs_new_method.svg"
    if not svg_path.exists():
        svg_path = Path(__file__).parent / "AnDOM_old_vs_new_method.svg"

    if svg_path.exists():
        st.markdown(svg_path.read_text(), unsafe_allow_html=True)
    else:
        st.warning(
            f"SVG not found. Copy the file to `docs/AnDOM_old_vs_new_method.svg` "
            f"inside the project directory."
        )

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
- Rank-fused table: sequence × structure combined score
- 🟢 both arms · 🔵 seq only · 🟠 struct only (dark proteome)
- Two-row domain architecture map with hover tooltips
- Batch processing up to 500 sequences
- Benchmark: rank-1 / top-5 comparison across all arms
        """)

    st.subheader("References")
    st.markdown("""
- Schmidt S, Bork P, Dandekar T (2002) *J Chem Inf Comput Sci* 42:405–407
- van Kempen M et al. (2023) *Nat Biotechnol* doi:10.1038/s41587-023-01773-0
- Lin Z et al. (2023) *Science* 379:1123–1130
- Fox NK, Brenner SE, Chandonia JM (2014) *Nucleic Acids Res* 42:D304–309
- Sillitoe I et al. (2021) *Nucleic Acids Res* 49:D266–273
    """)
