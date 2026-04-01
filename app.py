"""
AnDOM 2.0 — Streamlit web application.
"""
import sys
import os
import json
import time
import threading
import streamlit as st
import uuid as _uuid
import pandas as pd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent / 'src'))
sys.path.insert(0, str(Path(__file__).parent))

from config import SCOP_COLORS, SCOP_CLASSES, EXAMPLES, ESMFOLD_MAXLEN, SCOPE_FA
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
    return BatchManager(output_root=os.environ.get("ANDOM_OUTPUT_DIR", "output"))

batch_mgr = get_batch_manager()

# ── batch benchmark examples ──────────────────────────────────────────────────
BATCH_EXAMPLES = {
    "dark_proteome": {
        "label": "Dark proteome",
        "fasta": (
            ">SARS-CoV-2_nsp7\n"
            "MSKMVLSGGFGAGVQFMNLAKSTGPVVAAAVNGLMNLFGQKSPKTVNLNKDYIQRDTGEALNKILQDYINGVAPKALENRQLAGLRQLQEQR\n"
            ">Phage_T4_Soc\n"
            "MSDNKQIKAIVESVKDKLTSIKTSNDKLNQEADLIAKNGKNAISGVLENGKAEITKLQEELAKKAGVSTLSADDLAKKKNTDLVSIDQAKLKAAK\n"
            ">Methanogen_hypothetical\n"
            "MKILIVDDHPVVREGILEYLLSAEGYEVVCAEDGQEALDIYEDHPDLVLMDLMMPGMDGFELCRQIRQLDPRIPVLMLTAKDDEYDKVLGLEIGADDYVTKPFSTREELLARIRAHL\n"
        ),
        "desc": "Dark proteome proteins — SARS-CoV-2 nsp7, phage capsid protein, archaeal hypothetical protein. Sequence search finds nothing; structural arm recovers CATH domain.",
        "use_structure": True,
    },
    "viral_phage": {
        "label": "Viral / phage proteins",
        "fasta": (
            ">SARS-CoV-2_nsp7\n"
            "MSKMVLSGGFGAGVQFMNLAKSTGPVVAAAVNGLMNLFGQKSPKTVNLNKDYIQRDTGEALNKILQDYINGVAPKALENRQLAGLRQLQEQR\n"
            ">SARS-CoV-2_nsp8\n"
            "MIAGGHYVFKEIVMKDPEKFNEALKMLPIDGETVIAEQIAGLKNTLKYLRKLEKDLALKLNHITNDMSSEMAKQYKEYVNKVLPQLENFEDLTKLK\n"
        ),
        "desc": "SARS-CoV-2 RNA polymerase cofactors nsp7+nsp8. Low sequence identity to known proteins but conserved folds — structural arm recovers CATH class.",
        "use_structure": True,
    },
    "multidomain": {
        "label": "Multi-domain",
        "fasta": (
            ">Src_SH2_SH3\n"
            "MGSNKSKPKDASQRRRSLEPAENVHGAGGGAFPASQTPSKPASADGHRGPSAAFAPAAAEKVLFGGFNSS\n"
            "DTVTSPQRAGPLAGGVTTFVALYDYESRTETDLSFKKGERLQIVNNTEGDWWLAHSLSTGQTGYIPSNYW\n"
            ">p53_DBD\n"
            "SVVRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSC\n"
            "MGQMNRRPILTIITLEDSSGKLLGRNSFEVRVCACPGRDRRTEEENLRKKGEVVAPQHL\n"
        ),
        "desc": "Multi-domain proteins — ensemble shows both domains, evidence tagged per PDB match.",
    },
    "disordered": {
        "label": "Intrinsically disordered",
        "fasta": (
            ">p53_TAD\n"
            "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP\n"
            ">MBP_linker\n"
            "GKPIPNPLLGLDSTRTGHHHHHH\n"
        ),
        "desc": "IDPs with transient structured regions — low pLDDT from AlphaFold, structural arm finds fold family.",
    },
    "de_novo": {
        "label": "De novo designed",
        "fasta": (
            ">Designed_miniprotein\n"
            "MAKMAEEILSQYGDNLPKEVLEYLEKNIAKSGKIFVEIKADNAIAAANAKAVLENAKEVLEAYKTGKVNLV\n"
            ">Rosetta_designed_helix\n"
            "MAEAAAKEAAAKEAAAKEAAAKEAAAKEAAAKEAAAK\n"
        ),
        "desc": "Computationally designed proteins — zero evolutionary signal, only structural arm can classify.",
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
    st.caption(f"Batch limit: **{MAX_SEQUENCES}** seq · **{MAX_SEQ_LEN}** aa max")


def render_batch_cards(df_res, job_id):
    import pandas as pd
    from db.lookup import get_cath_code, cath_code_to_scop_class
    if df_res.empty:
        st.warning("No domain hits found — try enabling structural search for dark proteome proteins.")
        return
    queries = df_res["query_id"].unique() if "query_id" in df_res.columns else df_res["query"].unique()
    for qid in queries:
        q_df = df_res[df_res["query_id"]==qid] if "query_id" in df_res.columns else df_res[df_res["query"]==qid]
        seq_hits = q_df[q_df["source"]=="SCOPe_sequence"] if "source" in q_df.columns else q_df
        str_hits = q_df[q_df["source"]=="CATH_structure"] if "source" in q_df.columns else pd.DataFrame()
        has_seq = len(seq_hits) > 0
        has_str = len(str_hits) > 0
        if has_seq and has_str: icon,label = "🟢","Both arms — high confidence"
        elif has_seq: icon,label = "🔵","Sequence only"
        elif has_str: icon,label = "🟠","Structure only (dark proteome recovery)"
        else: icon,label = "⚪","No hits found"
        with st.container(border=True):
            seq_len = int(seq_hits.iloc[0]["qend"]) if has_seq and "qend" in seq_hits.columns else "?"
            st.markdown(f"### {icon} {qid}")
            st.caption(f"{label} · {seq_len} aa analysed")
            col1, col2 = st.columns(2)
            with col1:
                st.markdown("**🧬 Sequence arm (SCOPe)**")
                if has_seq:
                    top = seq_hits.iloc[0]
                    sccs = top.get("sccs","—"); cls_name = top.get("class_name","—")
                    pdb_link = top.get("PDB link",""); pdb_id = pdb_link.split("/")[-1] if pdb_link else "—"
                    try: ev = f"{float(top.get('evalue',0)):.1e}"
                    except: ev = "—"
                    st.markdown(f"**Fold:** `{sccs}` — {cls_name}")
                    if pdb_link: st.markdown(f"**Top hit:** [{pdb_id}]({pdb_link})")
                    st.caption(f"e-value: {ev} · {len(seq_hits)} SCOPe hits")
                    if len(seq_hits) > 1:
                        alts = seq_hits.iloc[1:5]
                        parts = []
                        for _, r in alts.iterrows():
                            pl = str(r.get("PDB link",""))
                            pid = pl.split("/")[-1] if pl else ""
                            if pid: parts.append(f"[{pid}]({pl})")
                        if parts: st.caption("Also: " + " · ".join(parts))
                else:
                    st.info("No sequence homologs in SCOPe ASTRAL 40%")
            with col2:
                st.markdown("**🏗️ Structure arm (CATH50)**")
                if has_str:
                    top_s = str_hits.iloc[0]
                    target = top_s.get("target","—")
                    cath_code = get_cath_code(str(target))
                    scop_cls = cath_code_to_scop_class(cath_code) if cath_code else "?"
                    pdb_link_s = top_s.get("PDB link",""); pdb_id_s = pdb_link_s.split("/")[-1] if pdb_link_s else "—"
                    try: lddt_s = f"{float(top_s.get('lddt',0)):.2f}"
                    except: lddt_s = "—"
                    try: ev_s = f"{float(top_s.get('evalue',0)):.1e}"
                    except: ev_s = "—"
                    st.markdown(f"**CATH:** `{cath_code}` → class `{scop_cls}`")
                    if pdb_link_s: st.markdown(f"**Top hit:** [{pdb_id_s}]({pdb_link_s})")
                    st.caption(f"lDDT: {lddt_s} · e-value: {ev_s} · {len(str_hits)} CATH hits")
                else:
                    st.info("Structural search not run or no CATH hits")
            if has_seq:
                pdb_id_enr = seq_hits.iloc[0].get("PDB link","").split("/")[-1]
                if pdb_id_enr and len(pdb_id_enr)==4:
                    with st.expander("🔬 Functional annotation"):
                        try:
                            import requests as _rq
                            r = _rq.get(f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id_enr}/1",timeout=5)
                            uids = r.json().get("rcsb_polymer_entity_container_identifiers",{}).get("uniprot_ids",[]) if r.status_code==200 else []
                            if uids:
                                uid=uids[0]
                                from search.ensemble import enrich_with_pfam,enrich_with_uniprot
                                uni=enrich_with_uniprot(uid); pfam=enrich_with_pfam(uid)
                                st.markdown(f"**UniProt:** [{uid}](https://www.uniprot.org/uniprot/{uid}) · **Gene:** {uni.get('gene','—')} · **Organism:** {uni.get('organism','—')}")
                                if uni.get("function"): st.markdown(f"**Function:** {uni['function'][:300]}")
                                if pfam: st.markdown("**Pfam:** "+" · ".join(f"`{p['pfam_id']}` {p['name']}" for p in pfam[:3]))
                            else: st.caption("No UniProt mapping found.")
                        except Exception as _e: st.caption(f"Enrichment unavailable: {_e}")

page = st.tabs(["Search", "Batch", "Methods"])

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
                    f"Sequence is {n} aa — ESMFold limit is {ESMFOLD_MAXLEN} aa. "
                    "Structural search disabled."
                )

            col_seq, col_str = st.columns(2)
            df_seq = df_str = None

            with col_seq:
                with st.spinner("Sequence search — SCOPe 2.08..."):
                    import uuid as _uuid; _tid = _uuid.uuid4().hex[:8]
                    df_seq, err = seq_search.run(fasta, evalue=evalue, iterations=iterations, tmp_dir=f"/output/tmp/{_tid}")
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

            # ── ensemble fusion (matched by PDB code) ──────────────────────
            fused = fuse_results(df_seq, df_str)
            if not fused.empty:
                st.subheader("Ensemble Ranked Hits")
                st.caption(
                    "🟢 found by **both arms** (high confidence) · "
                    "🔵 **sequence only** (SCOPe hit, no structure match) · "
                    "🟠 **structure only** (dark proteome — recovered by ESMFold+Foldseek)"
                )
                ev_icon = {"both": "🟢", "seq_only": "🔵", "struct_only": "🟠"}
                fused["ev"] = fused["evidence"].map(ev_icon)

                # count summary
                n_both   = (fused["evidence"] == "both").sum()
                n_seq    = (fused["evidence"] == "seq_only").sum()
                n_struct = (fused["evidence"] == "struct_only").sum()
                c1, c2, c3 = st.columns(3)
                c1.metric("🟢 Both arms", n_both)
                c2.metric("🔵 Seq only",  n_seq)
                c3.metric("🟠 Struct only (dark proteome)", n_struct)

                from db.lookup import get_cath_code
                fused["cath_code"] = fused["cath_domain"].apply(get_cath_code)
                st.dataframe(
                    fused[["rank","ev","scope_domain","sccs","cath_domain","cath_code",
                           "evidence","ensemble_score","seq_evalue","struct_evalue","lddt"]].rename(
                        columns={"ev": "",
                                 "scope_domain": "SCOPe domain", "sccs": "SCOP class",
                                 "cath_domain": "CATH domain", "cath_code": "CATH code",
                                 "ensemble_score": "Score",
                                 "seq_evalue": "Seq e-val",
                                 "struct_evalue": "Struct e-val"}),
                    use_container_width=True, hide_index=True,
                )
                # functional enrichment for top both hit
                both_hits = fused[fused["evidence"]=="both"]
                if not both_hits.empty:
                    top = both_hits.iloc[0]
                    pdb = str(top["scope_domain"])[1:5].lower()
                    with st.expander("🔬 Functional annotation (top hit)"):
                        try:
                            import requests as _req
                            r = _req.get(f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb}/1", timeout=5)
                            uids = r.json().get("rcsb_polymer_entity_container_identifiers",{}).get("uniprot_ids",[]) if r.status_code==200 else []
                            if uids:
                                uid = uids[0]
                                from search.ensemble import enrich_with_pfam, enrich_with_uniprot
                                uni = enrich_with_uniprot(uid)
                                pfam = enrich_with_pfam(uid)
                                st.markdown(f"**UniProt:** [{uid}](https://www.uniprot.org/uniprot/{uid}) | **Gene:** {uni.get('gene','—')} | **Organism:** {uni.get('organism','—')}")
                                if uni.get("function"):
                                    st.markdown(f"**Function:** {uni['function'][:300]}")
                                if pfam:
                                    st.markdown("**Pfam:** " + ", ".join(f"{p['pfam_id']} {p['name']}" for p in pfam))
                            else:
                                st.info("No UniProt mapping found.")
                        except Exception as _e:
                            st.warning(f"Enrichment unavailable: {_e}")
            st.subheader("Domain Architecture Map")
            bar_len = max(
                df_seq["qend"].max() if df_seq is not None and len(df_seq) > 0 else 1,
                df_str["qend"].max() if df_str is not None and len(df_str) > 0 else 1,
                n,
            )
            st.markdown(
                '<p style="font-size:11px;color:#888">'
                "Top: SCOPe sequence hits (coloured by SCOP class) · "
                "Bottom: CATH structural hits (dark) · "
                f"Hover for details · 1 to {bar_len} aa</p>"
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
                    st.dataframe(disp, use_container_width=True,
                        column_config={"PDB link": st.column_config.LinkColumn("PDB link")})
                    if os.path.exists("seq_results.tsv"):
                        st.download_button("Download SCOPe TSV",
                            open("seq_results.tsv").read(), "AnDOM_scope.tsv")
                else:
                    st.info("No sequence hits.")

            with tab2:
                if df_str is not None and len(df_str) > 0:
                    disp2 = df_str[[
                        "target","evalue","bits","lddt","qstart","qend","PDB link"
                    ]].copy()
                    disp2["evalue"] = disp2["evalue"].apply(lambda x: f"{float(x):.2e}")
                    disp2["lddt"]   = disp2["lddt"].apply(lambda x: f"{float(x):.3f}")
                    st.dataframe(disp2, use_container_width=True,
                        column_config={"PDB link": st.column_config.LinkColumn("PDB link")})
                    if os.path.exists("struct_results.tsv"):
                        st.download_button("Download CATH TSV",
                            open("struct_results.tsv").read(), "AnDOM_cath.tsv")
                else:
                    st.info("No structural hits or structural search disabled.")

    st.divider()
    st.caption("AnDOM 2.0 | SCOPe 2.08 ASTRAL 40% (MMseqs2 PSI) + CATH50 (ESMFold + Foldseek) | Dandekar Lab Wuerzburg")

# ══════════════════════════════════════════════════════════════════════════════
# TAB 1 — Batch
# ══════════════════════════════════════════════════════════════════════════════
with page[1]:
    st.title("Batch Processing")
    st.caption(
        f"Up to **{MAX_SEQUENCES}** sequences per job · "
        f"{MIN_SEQ_LEN}–{MAX_SEQ_LEN} aa · jobs run in background"
    )

    with st.expander("💡 Why AnDOM 2.0 alongside AlphaFold + Foldseek?", expanded=False):
        st.markdown("""
| Challenge | AlphaFold | Foldseek | InterPro | **AnDOM 2.0** |
|---|---|---|---|---|
| Assigns domain to SCOP/CATH class | ❌ | ✅ if DB hit | seq only | **seq + struct** |
| Dark proteome (no homologs) | structure only | ❌ | ❌ | **structural arm** |
| Fast-evolving viral proteins | structure only | weak | ❌ | **ensemble** |
| Multi-domain, one novel | ❌ | partial | partial | **both flagged** |
| Disordered + folded domain | low pLDDT | poor | misses | **structural arm** |

🟠 `struct_only` hits = dark proteome recovery.
🟢 `both` hits = high-confidence annotation.
        """)

    # ── example buttons ───────────────────────────────────────────────────
    st.markdown("**Load a benchmark example set:**")
    ex_cols = st.columns(len(BATCH_EXAMPLES))
    for i, (key, ex) in enumerate(BATCH_EXAMPLES.items()):
        if ex_cols[i].button(ex["label"], use_container_width=True, key=f"bex_{key}"):
            st.session_state["batch_text"] = ex["fasta"]
            st.session_state["_batch_desc"]  = ex["desc"]
            if ex.get("use_structure", False):
                st.session_state["b_struct"] = True

    if "_batch_desc" in st.session_state:
        st.info(st.session_state["_batch_desc"])

    # ── submit form ───────────────────────────────────────────────────────
    with st.expander("➕ Submit new batch job", expanded=True):
        b_upload = st.file_uploader(
            "Upload multi-FASTA", type=["fa","fasta","txt"], key="batch_upload"
        )

        # pre-fill from example button if set
        b_text  = st.text_area(
            "…or paste multi-FASTA here", height=180,
            key="batch_text"
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
                        tmp_path, evalue=evalue,
                        iterations=iterations, use_structure=b_struct,
                    )
                    st.success(f"Job **{job_id}** submitted — {len(seqs)} sequences.")
                    st.session_state.pop("_batch_fasta", None)
                    st.session_state.pop("_batch_desc", None)

    # ── job list ──────────────────────────────────────────────────────────
    st.subheader("Jobs")
    jobs = batch_mgr.list_jobs()
    if not jobs:
        st.info("No batch jobs yet.")
    else:
        if st.button("🔄 Refresh"):
            st.rerun()
        status_icon = {"queued":"🕐","running":"⏳","done":"✅","failed":"❌","cancelled":"🚫"}
        for job in jobs:
            s = job.get("status","unknown")
            with st.container(border=True):
                c1, c2, c3 = st.columns([3,2,2])
                c1.markdown(f"{status_icon.get(s,'❓')} **{job['job_id']}** — `{s}`")
                total = job.get("total",0); prog = job.get("progress",0)
                c2.caption(f"{prog}/{total} sequences" if total else job.get("submitted","")[:19])
                if s == "done":
                    df_res = batch_mgr.results(job["job_id"])
                    if not df_res.empty:
                        c3.download_button("⬇ TSV",
                            df_res.to_csv(sep="\t",index=False),
                            file_name=f"AnDOM_batch_{job['job_id']}.tsv",
                            mime="text/tab-separated-values",
                            key=f"dl_{job['job_id']}")
                        render_batch_cards(df_res, job["job_id"])
                    if s in ("queued","running"):
                        if c3.button("Cancel", key="cancel_"+job["job_id"]):
                            batch_mgr.cancel(job["job_id"]); st.rerun()
                if job.get("error"):
                    st.error(job["error"])

# ══════════════════════════════════════════════════════════════════════════════
# TAB 2 — Benchmark (background thread, non-blocking) — HIDDEN FROM PUBLIC UI
# ══════════════════════════════════════════════════════════════════════════════
if False:
    st.title("Benchmark — Sequence vs Structure vs Ensemble")
    st.markdown(
        "Compares all three arms on **SCOPE-40** (domains clustered at 40% sequence identity — "
        "diverse enough that sequence methods alone struggle). Reports rank-1 accuracy and "
        "top-5 sensitivity. This is the evidence for the paper claim: "
        "*the ensemble recovers domains that neither arm finds alone.*"
    )
    st.divider()

    bm_scope = st.toggle("Use built-in SCOPE-40 dataset", value=True)
    if bm_scope:
        bm_fasta = SCOPE_FA
        st.caption(f"Dataset: `{bm_fasta}`")
    else:
        bm_file = st.file_uploader("Upload SCOPE-style FASTA", type=["fa","fasta"], key="bm_upload")
        bm_fasta = None
        if bm_file:
            import tempfile
            with tempfile.NamedTemporaryFile(mode="wb", suffix=".fa", delete=False) as tmp:
                tmp.write(bm_file.read()); bm_fasta = tmp.name

    bm_limit = st.number_input("Limit to first N sequences (0 = all)", min_value=0, value=50, step=10)

    # results stored in session state so they survive reruns
    if "bm_result" not in st.session_state:
        st.session_state["bm_result"] = None
    if "bm_running" not in st.session_state:
        st.session_state["bm_running"] = False
    if "bm_error" not in st.session_state:
        st.session_state["bm_error"] = None

    def _run_benchmark_thread(fasta_path, ev, iters, limit):
        """Runs in background thread, stores result in session state."""
        try:
            from benchmark.run import run_arms_benchmark
            ts = time.strftime("%Y%m%d_%H%M%S")
            stats = run_arms_benchmark(
                fasta_path, evalue=ev, iterations=iters,
                limit=int(limit) if limit > 0 else None,
                out_json=f"output/results/benchmark_{ts}.json",
            )
            st.session_state["bm_result"]  = stats
            st.session_state["bm_error"]   = None
        except Exception as exc:
            st.session_state["bm_error"]   = str(exc)
        finally:
            st.session_state["bm_running"] = False

    if st.button("▶ Run Benchmark", type="primary", disabled=st.session_state["bm_running"]):
        if not bm_fasta:
            st.warning("Please provide a FASTA file.")
        else:
            st.session_state["bm_running"] = True
            st.session_state["bm_result"]  = None
            st.session_state["bm_error"]   = None
            t = threading.Thread(
                target=_run_benchmark_thread,
                args=(bm_fasta, evalue, iterations, bm_limit),
                daemon=True,
            )
            t.start()

    if st.session_state["bm_running"]:
        st.info("⏳ Benchmark running in background — come back to this tab in a few minutes and refresh.")

    if st.session_state["bm_error"]:
        st.error(f"Benchmark failed: {st.session_state['bm_error']}")

    stats = st.session_state["bm_result"]
    if stats:
        st.success(f"Done — {stats['total']} sequences in {stats.get('elapsed_s','?')}s")
        rows = []
        for arm in ("seq_only","struct_only","ensemble"):
            d = stats[arm]
            rows.append({"Arm":arm,"Rank-1 %":d["rank1_pct"],"Top-5 %":d["top5_pct"],
                          "Rank-1 n":d["rank1"],"Top-5 n":d["top5"]})
        df_bm = pd.DataFrame(rows)
        st.dataframe(df_bm, use_container_width=True, hide_index=True)
        st.bar_chart(df_bm.set_index("Arm")[["Rank-1 %","Top-5 %"]])
        st.caption(
            "ensemble > seq_only AND ensemble > struct_only → combination justified. "
            "struct_only > 0 → dark proteome recovery confirmed."
        )
        st.download_button("⬇ Download JSON", data=json.dumps(stats,indent=2),
            file_name="benchmark.json", mime="application/json")

# ══════════════════════════════════════════════════════════════════════════════
# TAB 3 — Methods
# ══════════════════════════════════════════════════════════════════════════════
with page[2]:
    st.title("Methods — AnDOM 2.0")
    st.markdown(
        "AnDOM 2.0 updates the original AnDOM server "
        "(Schmidt, Bork & Dandekar, *J Chem Inf Comput Sci* 2002) "
        "with modern tools and a new structural search layer."
    )

    st.subheader("Pipeline comparison: AnDOM 2002 vs AnDOM 2.0")

    # try both locations
    for candidate in [
        Path(__file__).parent / "docs" / "AnDOM_old_vs_new_method.svg",
        Path(__file__).parent / "AnDOM_old_vs_new_method.svg",
        Path("/app/docs/AnDOM_old_vs_new_method.svg"),
    ]:
        if candidate.exists():
            svg_text = candidate.read_text()
            st.image(str(candidate))
            break
    else:
        st.warning("SVG diagram not found. Expected at `docs/AnDOM_old_vs_new_method.svg`.")

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

**Ensemble fusion (new)**
- Matched by 4-char PDB code across both arms
- 🟢 both arms = high-confidence annotation
- 🟠 struct_only = dark proteome recovery
- Rank-fused score: 40% seq + 60% struct
- Batch up to 500 sequences, background jobs
- Benchmark: rank-1 / top-5 per arm on SCOPE-40
        """)

    st.subheader("References")
    st.markdown("""
- Schmidt S, Bork P, Dandekar T (2002) *J Chem Inf Comput Sci* 42:405–407
- van Kempen M et al. (2023) *Nat Biotechnol* doi:10.1038/s41587-023-01773-0
- Lin Z et al. (2023) *Science* 379:1123–1130
- Fox NK, Brenner SE, Chandonia JM (2014) *Nucleic Acids Res* 42:D304–309
- Sillitoe I et al. (2021) *Nucleic Acids Res* 49:D266–273
    """)
