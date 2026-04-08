"""
AnDOM 2.0 — Streamlit web application.
Three-arm parallel display: Sequence · Structure · Profile
Each arm shows: SCOPe top 7 · CATH top 7 · PDB top 7 with functional lookup
Ensemble voting fuses all three arms.
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
from search.ensemble import (
    domain_bar_html, add_scores, fuse_results_three,
    fetch_pdb_function, enrich_with_pfam, enrich_with_uniprot,
)
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

# ── batch examples ────────────────────────────────────────────────────────────
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
        "desc": "Dark proteome — sequence search finds nothing; structural arm recovers CATH domain.",
        "use_structure": True,
    },
    "viral_phage": {
        "label": "Viral / phage",
        "fasta": (
            ">SARS-CoV-2_nsp7\n"
            "MSKMVLSGGFGAGVQFMNLAKSTGPVVAAAVNGLMNLFGQKSPKTVNLNKDYIQRDTGEALNKILQDYINGVAPKALENRQLAGLRQLQEQR\n"
            ">SARS-CoV-2_nsp8\n"
            "MIAGGHYVFKEIVMKDPEKFNEALKMLPIDGETVIAEQIAGLKNTLKYLRKLEKDLALKLNHITNDMSSEMAKQYKEYVNKVLPQLENFEDLTKLK\n"
        ),
        "desc": "SARS-CoV-2 RNA polymerase cofactors nsp7+nsp8.",
        "use_structure": True,
    },
    "multidomain": {
        "label": "Multi-domain",
        "fasta": (
            ">Src_SH2_SH3\n"
            "MGSNKSKPKDASQRRRSLEPAENVHGAGGGAFPASQTPSKPASADGHRGPSAAFAPAAAEKVLFGGFNSS"
            "DTVTSPQRAGPLAGGVTTFVALYDYESRTETDLSFKKGERLQIVNNTEGDWWLAHSLSTGQTGYIPSNYW\n"
            ">p53_DBD\n"
            "SVVRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSC"
            "MGQMNRRPILTIITLEDSSGKLLGRNSFEVRVCACPGRDRRTEEENLRKKGEVVAPQHL\n"
        ),
        "desc": "Multi-domain proteins — ensemble shows both domains.",
    },
}

# ── helper: clean CATH domain display ─────────────────────────────────────────
def clean_cath_domain(d: str) -> str:
    if str(d).startswith("af_"):
        parts = str(d).split("_")
        return f"AlphaFold:{parts[1]}" if len(parts) > 1 else str(d)
    return str(d)

# ── helper: format evalue ──────────────────────────────────────────────────────
def fmt_e(v) -> str:
    try:
        return f"{float(v):.2e}"
    except Exception:
        return "—"

# ── helper: pdb id from domain ─────────────────────────────────────────────────
def pdb_from_scope(domain: str) -> str:
    try:
        return domain[1:5].lower()
    except Exception:
        return ""

# ── per-hit functional expander ────────────────────────────────────────────────
def render_func_expander(pdb_id: str, label: str = "🔬 Functional info"):
    if not pdb_id or len(pdb_id) != 4:
        return
    with st.expander(label):
        with st.spinner("Fetching from RCSB + UniProt…"):
            info = fetch_pdb_function(pdb_id)
        if info["uniprot_id"]:
            uid = info["uniprot_id"]
            st.markdown(
                f"**UniProt:** [{uid}](https://www.uniprot.org/uniprot/{uid}) &nbsp;|&nbsp; "
                f"**Gene:** {info.get('gene','—')} &nbsp;|&nbsp; "
                f"**Organism:** {info.get('organism','—')}"
            )
            if info.get("function"):
                st.markdown(f"**Function:** {info['function'][:350]}")
            if info.get("pfam"):
                pfam_str = " · ".join(
                    f"`{p['pfam_id']}` {p['name']}" for p in info["pfam"][:4]
                )
                st.markdown(f"**Pfam:** {pfam_str}")
        else:
            st.caption("No UniProt mapping found for this PDB entry.")

# ── arm result panel ─────────────────────────────────────────────────────────
def render_arm_panel(
    title: str,
    color: str,
    df: pd.DataFrame | None,
    arm: str,   # "seq" | "str" | "hh"
    top_n: int = 7,
):
    """Render one arm's results: SCOPe class · CATH code · top PDB hits."""
    lk = lookup.all_domains()

    st.markdown(
        f'<div style="font-weight:600;font-size:15px;color:{color};'
        f'border-left:3px solid {color};padding-left:8px;margin-bottom:8px">'
        f'{title}</div>',
        unsafe_allow_html=True,
    )

    if df is None or len(df) == 0:
        st.info("No hits found.")
        return

    hits = df.head(top_n)

    # ── SCOPe class ──
    st.markdown("**SCOPe classes**")
    if arm == "seq":
        sccs_rows = []
        for _, r in hits.iterrows():
            sccs = lk.get(str(r["target"]), {}).get("sccs", "—")
            cls  = lk.get(str(r["target"]), {}).get("cls", "?")
            color_dot = SCOP_COLORS.get(cls, "#888")
            sccs_rows.append({
                "●": f'<span style="color:{color_dot}">●</span>',
                "SCOPe": sccs,
                "class": lk.get(str(r["target"]), {}).get("class_name", cls),
                "e-val": fmt_e(r["evalue"]),
            })
        df_sccs = pd.DataFrame(sccs_rows)
        st.write(df_sccs[["SCOPe", "class", "e-val"]].to_html(index=False, escape=False), unsafe_allow_html=True)
    elif arm == "str":
        from db.lookup import get_cath_code, cath_code_to_scop_class
        for _, r in hits.iterrows():
            cc   = get_cath_code(str(r["target"]))
            cls  = cath_code_to_scop_class(cc) if cc else "?"
            lbl  = SCOP_CLASSES.get(cls, cls)
            col  = SCOP_COLORS.get(cls, "#888")
            st.markdown(
                f'<span style="color:{col}">●</span> '
                f'`{cc}` → **{cls}** {lbl} &nbsp; lDDT={float(r.get("lddt",0)):.2f}',
                unsafe_allow_html=True,
            )
    elif arm == "hh":
        for _, r in hits.iterrows():
            sccs = str(r.get("sccs", "—"))
            cls  = sccs[0] if sccs not in ("—","?","") else "?"
            col  = SCOP_COLORS.get(cls, "#888")
            lbl  = SCOP_CLASSES.get(cls, cls)
            st.markdown(
                f'<span style="color:{col}">●</span> '
                f'`{sccs}` → **{cls}** {lbl} &nbsp; prob={float(r.get("prob",0)):.0f}%',
                unsafe_allow_html=True,
            )

    st.divider()

    # ── CATH codes ──
    st.markdown("**CATH codes**")
    from db.lookup import get_cath_code
    cath_seen: list[str] = []
    for _, r in hits.iterrows():
        if arm == "hh":
            pdb = str(r.get("pdb", ""))
            cc  = get_cath_code(pdb)
        else:
            cc = get_cath_code(str(r["target"]))
        if cc and cc not in cath_seen:
            cath_seen.append(cc)
            cath_url = f"https://www.cathdb.info/version/v4_3_0/superfamily/{cc}"
            st.markdown(f"[`{cc}`]({cath_url})")
        if len(cath_seen) >= top_n:
            break

    st.divider()

    # ── top PDB hits ──
    st.markdown("**Top PDB hits**")
    for i, (_, r) in enumerate(hits.iterrows()):
        if arm == "seq":
            pdb_id  = pdb_from_scope(str(r["target"]))
            pdb_url = str(r.get("PDB link", f"https://www.rcsb.org/structure/{pdb_id}"))
            ev_str  = fmt_e(r["evalue"])
            extra   = f"e={ev_str} · {float(r.get('pident',0)):.0f}%id"
            desc    = lk.get(str(r["target"]), {}).get("desc", "")[:50]
        elif arm == "str":
            pdb_id  = str(r["target"])[:4].lower()
            pdb_url = str(r.get("PDB link", f"https://www.rcsb.org/structure/{pdb_id}"))
            ev_str  = fmt_e(r["evalue"])
            extra   = f"e={ev_str} · lDDT={float(r.get('lddt',0)):.2f}"
            desc    = clean_cath_domain(str(r["target"]))
        else:
            pdb_id  = str(r.get("pdb", ""))
            pdb_url = str(r.get("PDB link", f"https://www.rcsb.org/structure/{pdb_id}"))
            ev_str  = fmt_e(r["evalue"])
            extra   = f"e={ev_str} · prob={float(r.get('prob',0)):.0f}%"
            desc    = str(r.get("hit_name", ""))

        with st.container():
            c1, c2 = st.columns([3, 1])
            with c1:
                st.markdown(f"**{i+1}.** [{pdb_id.upper()}]({pdb_url}) — {desc}")
                st.caption(extra)
            with c2:
                if st.button("func", key=f"func_{arm}_{i}_{pdb_id}", help="Functional info"):
                    st.session_state[f"func_open_{arm}_{i}"] = True

            if st.session_state.get(f"func_open_{arm}_{i}"):
                render_func_expander(pdb_id, label=f"Function: {pdb_id.upper()}")


# ── sidebar ───────────────────────────────────────────────────────────────────
with st.sidebar:
    st.header("Search Parameters")
    evalue     = st.select_slider("E-value cutoff",
                    options=[1e-30,1e-20,1e-10,1e-5,1e-3,0.01,0.1,1.0], value=1e-3)
    iterations = st.slider("PSI-MMseqs2 iterations", 1, 5, 3)
    max_hits   = st.slider("Max hits shown", 5, 100, 30)
    use_struct = st.toggle("Structural search (ESMFold + Foldseek)", value=True)
    use_hhblits = st.toggle(
        "🔬 Deep search (HHblits twilight zone)", value=False, key="use_hhblits",
        help="Profile-profile search via HHblits+UniClust30. Finds remote homologs at 10-15% identity. Adds ~2 min.",
    )
    if use_hhblits:
        st.caption("HHblits arm: PDB70 database loading...")
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


# ── batch card renderer ───────────────────────────────────────────────────────
def render_batch_cards(df_res, job_id):
    from db.lookup import get_cath_code, cath_code_to_scop_class
    if df_res.empty:
        st.warning("No domain hits found.")
        return
    queries = [q for q in (
        df_res["query_id"].unique() if "query_id" in df_res.columns
        else df_res["query"].unique()
    ) if str(q) != "nan"]
    for qid in queries:
        q_df = (df_res[df_res["query_id"] == qid] if "query_id" in df_res.columns
                else df_res[df_res["query"] == qid])
        seq_hits = q_df[q_df["source"] == "SCOPe_sequence"] if "source" in q_df.columns else q_df
        str_hits = q_df[q_df["source"] == "CATH_structure"] if "source" in q_df.columns else pd.DataFrame()
        has_seq = len(seq_hits) > 0
        has_str = len(str_hits) > 0
        if has_seq and has_str:   icon, label = "🟢", "Both arms — high confidence"
        elif has_seq:             icon, label = "🔵", "Sequence only"
        elif has_str:             icon, label = "🟠", "Structure only (dark proteome recovery)"
        else:                     icon, label = "⚪", "No hits found"

        with st.container(border=True):
            seq_len = (int(seq_hits.iloc[0]["qend"]) if has_seq and "qend" in seq_hits.columns
                       else (int(str_hits.iloc[0]["qend"]) if has_str and "qend" in str_hits.columns else "?"))
            st.markdown(f"### {icon} {qid}")
            st.caption(f"{label} · {seq_len} aa analysed")
            col1, col2 = st.columns(2)
            with col1:
                st.markdown("**🧬 Sequence arm (SCOPe)**")
                if has_seq:
                    top = seq_hits.iloc[0]
                    sccs = top.get("sccs", "—")
                    cls_name = top.get("class_name", "—")
                    pdb_link = top.get("PDB link", "")
                    pdb_id = pdb_link.split("/")[-1] if pdb_link else "—"
                    st.markdown(f"**Fold:** `{sccs}` — {cls_name}")
                    if pdb_link:
                        st.markdown(f"**Top hit:** [{pdb_id}]({pdb_link})")
                    st.caption(f"e-value: {fmt_e(top.get('evalue',0))} · {len(seq_hits)} SCOPe hits")
                else:
                    st.info("No sequence homologs in SCOPe ASTRAL 95%")
            with col2:
                st.markdown("**🏗️ Structure arm (CATH50)**")
                if has_str:
                    top_s = str_hits.iloc[0]
                    target = top_s.get("target", "—")
                    cath_code = get_cath_code(str(target))
                    scop_cls = cath_code_to_scop_class(cath_code) if cath_code else "?"
                    pdb_link_s = top_s.get("PDB link", "")
                    pdb_id_s = pdb_link_s.split("/")[-1] if pdb_link_s else "—"
                    st.markdown(f"**CATH:** `{cath_code}` → class `{scop_cls}`")
                    if pdb_link_s:
                        st.markdown(f"**Top hit:** [{pdb_id_s}]({pdb_link_s})")
                    st.caption(
                        f"lDDT: {float(top_s.get('lddt', 0)):.2f} · "
                        f"e-value: {fmt_e(top_s.get('evalue', 0))} · {len(str_hits)} CATH hits"
                    )
                else:
                    st.info("Structural search not run or no CATH hits")


page = st.tabs(["Search", "Batch", "Methods"])

# ══════════════════════════════════════════════════════════════════════════════
# TAB 0 — Search
# ══════════════════════════════════════════════════════════════════════════════
with page[0]:
    st.title("AnDOM 2.0 — Structural Domain Finder")
    st.caption(
        "Sequence + Structure + Profile ensemble | SCOPe 2.08 + CATH50 + PDB70 | "
        "MMseqs2 + ESMFold + Foldseek + HHblits | Dandekar Lab Wuerzburg"
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

            df_seq = df_str = df_hh = None
            _tid = _uuid.uuid4().hex[:8]

            # ── run the three arms ────────────────────────────────────────────
            col_s, col_t, col_h = st.columns(3)

            with col_s:
                with st.spinner("🧬 Sequence search…"):
                    df_seq, err = seq_search.run(
                        fasta, evalue=evalue, iterations=iterations,
                        tmp_dir=f"/output/tmp/{_tid}"
                    )
                if err:
                    st.error(f"Sequence search failed: {err}")
                elif df_seq is None or len(df_seq) == 0:
                    st.warning("No sequence hits.")
                else:
                    df_seq = df_seq.head(max_hits)
                    st.success(f"Sequence: {len(df_seq)} SCOPe hits")

            with col_t:
                if use_struct_run:
                    with st.spinner("🏗 Structure search…"):
                        df_str, err2 = str_search.run(clean_seq)
                    if err2:
                        st.error(f"Structural search failed: {err2}")
                    elif df_str is None or len(df_str) == 0:
                        st.warning("No structural hits.")
                    else:
                        df_str = df_str.head(max_hits)
                        st.success(f"Structure: {len(df_str)} CATH hits")
                else:
                    st.info("Structural search disabled.")

            with col_h:
                if use_hhblits:
                    with st.spinner("🔬 Deep search — HHblits + PDB70…"):
                        from search.sequence_hhblits import run_hhblits
                        _hid = _uuid.uuid4().hex[:8]
                        df_hh, err_hh = run_hhblits(
                            fasta, iterations=3, threads=8,
                            tmp_dir=f"/output/tmp/hh_{_hid}"
                        )
                    if err_hh:
                        st.warning(f"HHblits: {err_hh}")
                    elif df_hh is not None and len(df_hh) > 0:
                        df_hh["source"] = "HHblits_PDB70"
                        st.success(f"HHblits: {len(df_hh)} PDB70 hits")
                    else:
                        st.info("HHblits: no hits found")
                else:
                    st.info("Deep search disabled.")

            df_seq, df_str = add_scores(df_seq, df_str)

            # ── THREE-ARM COLUMN DISPLAY ──────────────────────────────────────
            st.subheader("Results by arm")
            st.caption(
                "Each arm independently predicts SCOPe class · CATH code · top PDB hits. "
                "Click **func** on any PDB hit for functional annotation."
            )

            arm_col1, arm_col2, arm_col3 = st.columns(3)
            with arm_col1:
                render_arm_panel(
                    "🧬 Sequence arm — MMseqs2 / SCOPe",
                    "#3B82F6", df_seq, "seq", top_n=7,
                )
            with arm_col2:
                render_arm_panel(
                    "🏗 Structure arm — ESMFold / CATH50",
                    "#0F6E56", df_str, "str", top_n=7,
                )
            with arm_col3:
                render_arm_panel(
                    "🔬 Profile arm — HHblits / PDB70",
                    "#7F77DD", df_hh, "hh", top_n=7,
                )

            # ── DOMAIN ARCHITECTURE MAP (3 rows) ─────────────────────────────
            st.subheader("Domain Architecture Map")
            bar_len = max(
                df_seq["qend"].max() if df_seq is not None and len(df_seq) > 0 else 1,
                df_str["qend"].max() if df_str is not None and len(df_str) > 0 else 1,
                df_hh["qend"].max()  if df_hh  is not None and len(df_hh)  > 0 else 1,
                n,
            )
            st.markdown(
                domain_bar_html(df_seq, df_str, bar_len, df_hh),
                unsafe_allow_html=True,
            )

            # ── ENSEMBLE VOTING TABLE ─────────────────────────────────────────
            fused = fuse_results_three(df_seq, df_str, df_hh)
            if not fused.empty:
                st.subheader("Ensemble Ranked Hits")

                ev_icon  = {
                    "all_three":   "🟢",
                    "two_arms":    "🟡",
                    "seq_only":    "🔵",
                    "struct_only": "🟠",
                    "hhblits_only":"🟣",
                }
                ev_label = {
                    "all_three":   "All three arms — highest confidence",
                    "two_arms":    "Two arms agree — high confidence",
                    "seq_only":    "Sequence arm only",
                    "struct_only": "Structure arm only (dark proteome)",
                    "hhblits_only":"Profile arm only (twilight zone)",
                }
                st.caption(
                    " · ".join(f"{v} {ev_label[k]}" for k, v in ev_icon.items())
                )

                n_all    = (fused["evidence"] == "all_three").sum()
                n_two    = (fused["evidence"] == "two_arms").sum()
                n_seq    = (fused["evidence"] == "seq_only").sum()
                n_struct = (fused["evidence"] == "struct_only").sum()
                n_hh     = (fused["evidence"] == "hhblits_only").sum()
                c1, c2, c3, c4, c5 = st.columns(5)
                c1.metric("🟢 All three", n_all)
                c2.metric("🟡 Two arms",  n_two)
                c3.metric("🔵 Seq only",  n_seq)
                c4.metric("🟠 Struct only", n_struct)
                c5.metric("🟣 Profile only", n_hh)

                fused["ev"] = fused["evidence"].map(ev_icon)
                fused["cath_display"] = fused["cath_domain"].apply(clean_cath_domain)
                for col in ["seq_evalue", "struct_evalue", "hh_evalue"]:
                    if col in fused.columns:
                        fused[col] = fused[col].apply(
                            lambda x: fmt_e(x) if x is not None and str(x) not in ("None","nan") else "—"
                        )
                fused["hh_prob"] = fused["hh_prob"].apply(
                    lambda x: f"{float(x):.0f}%" if x is not None and str(x) not in ("None","nan") else "—"
                )

                st.dataframe(
                    fused[[
                        "rank","ev","scope_domain","sccs","cath_display","cath_code",
                        "hh_hit","votes","ensemble_score",
                        "seq_evalue","struct_evalue","hh_evalue","hh_prob","lddt"
                    ]].rename(columns={
                        "ev":             "",
                        "scope_domain":   "SCOPe domain",
                        "sccs":           "SCOP class",
                        "cath_display":   "CATH domain",
                        "cath_code":      "CATH code",
                        "hh_hit":         "HHblits hit",
                        "votes":          "Votes",
                        "ensemble_score": "Score",
                        "seq_evalue":     "Seq e-val",
                        "struct_evalue":  "Struct e-val",
                        "hh_evalue":      "HH e-val",
                        "hh_prob":        "HH prob",
                        "lddt":           "lDDT",
                    }),
                    use_container_width=True,
                    hide_index=True,
                )

                # functional annotation for top hit
                top_row = fused.iloc[0]
                top_pdb = pdb_from_scope(str(top_row.get("scope_domain", "")))
                if not top_pdb:
                    top_pdb = str(top_row.get("hh_hit", ""))[:4].lower()
                if top_pdb:
                    render_func_expander(top_pdb, "🔬 Functional annotation (top ensemble hit)")



    st.divider()
    st.caption(
        "AnDOM 2.0 | SCOPe 2.08 ASTRAL 95% (MMseqs2 PSI) + CATH50 (ESMFold + Foldseek) + PDB70 (HHblits) | "
        "Dandekar Lab Wuerzburg"
    )

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
| Assigns domain to SCOP/CATH class | ❌ | ✅ if DB hit | seq only | **seq + struct + profile** |
| Dark proteome (no homologs) | structure only | ❌ | ❌ | **structural arm** |
| Twilight zone (10-15% identity) | ❌ | ❌ | ❌ | **HHblits profile arm** |
| Fast-evolving viral proteins | structure only | weak | ❌ | **ensemble** |
| Multi-domain, one novel | ❌ | partial | partial | **both flagged** |

🟣 `hhblits_only` = twilight zone recovery.
🟠 `struct_only` = dark proteome recovery.
🟢 `all_three` = highest-confidence annotation.
        """)

    st.markdown("**Load a benchmark example set:**")
    ex_cols = st.columns(len(BATCH_EXAMPLES))
    for i, (key, ex) in enumerate(BATCH_EXAMPLES.items()):
        if ex_cols[i].button(ex["label"], use_container_width=True, key=f"bex_{key}"):
            st.session_state["batch_text"]  = ex["fasta"]
            st.session_state["_batch_desc"] = ex["desc"]
            if ex.get("use_structure", False):
                st.session_state["b_struct"] = True

    if "_batch_desc" in st.session_state:
        st.info(st.session_state["_batch_desc"])

    with st.expander("➕ Submit new batch job", expanded=True):
        b_upload = st.file_uploader(
            "Upload multi-FASTA", type=["fa","fasta","txt"], key="batch_upload"
        )
        b_text = st.text_area("…or paste multi-FASTA here", height=180, key="batch_text")
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
                    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as tmp:
                        tmp.write(raw)
                        tmp_path = tmp.name
                    job_id = batch_mgr.submit(
                        tmp_path, evalue=evalue,
                        iterations=iterations, use_structure=b_struct,
                    )
                    st.success(f"Job **{job_id}** submitted — {len(seqs)} sequences.")
                    st.session_state.pop("_batch_desc", None)

    st.subheader("Jobs")
    jobs = batch_mgr.list_jobs()
    if not jobs:
        st.info("No batch jobs yet.")
    else:
        if st.button("🔄 Refresh"):
            st.rerun()
        status_icon = {"queued":"🕐","running":"⏳","done":"✅","failed":"❌","cancelled":"🚫"}
        for job in jobs:
            s = job.get("status", "unknown")
            with st.container(border=True):
                c1, c2, c3 = st.columns([3,2,2])
                c1.markdown(f"{status_icon.get(s,'❓')} **{job['job_id']}** — `{s}`")
                total = job.get("total",0); prog = job.get("progress",0)
                c2.caption(f"{prog}/{total} sequences" if total else job.get("submitted","")[:19])
                if s == "done":
                    df_res = batch_mgr.results(job["job_id"])
                    if not df_res.empty:
                        c3.download_button("⬇ TSV",
                            df_res.to_csv(sep="\t", index=False),
                            file_name=f"AnDOM_batch_{job['job_id']}.tsv",
                            mime="text/tab-separated-values",
                            key=f"dl_{job['job_id']}")
                        render_batch_cards(df_res, job["job_id"])
                if s in ("queued","running"):
                    if c3.button("Cancel", key="cancel_" + job["job_id"]):
                        batch_mgr.cancel(job["job_id"]); st.rerun()
                if job.get("error"):
                    st.error(job["error"])

# ══════════════════════════════════════════════════════════════════════════════
# TAB 2 — Methods
# ══════════════════════════════════════════════════════════════════════════════
with page[2]:
    st.title("Methods — AnDOM 2.0")
    st.markdown(
        "AnDOM 2.0 updates the original AnDOM server "
        "(Schmidt, Bork & Dandekar, *J Chem Inf Comput Sci* 2002) "
        "with modern tools, a structural search layer, and a profile-profile twilight zone arm."
    )

    st.subheader("Pipeline comparison: AnDOM 2002 vs AnDOM 2.0")
    for candidate in [
        Path(__file__).parent / "docs" / "AnDOM_old_vs_new_method.svg",
        Path(__file__).parent / "AnDOM_old_vs_new_method.svg",
        Path("/app/docs/AnDOM_old_vs_new_method.svg"),
    ]:
        if candidate.exists():
            st.image(str(candidate))
            break
    else:
        st.warning("SVG diagram not found. Expected at `docs/AnDOM_old_vs_new_method.svg`.")

    st.subheader("Three-arm ensemble")
    c1, c2, c3 = st.columns(3)
    with c1:
        st.markdown("""
**🧬 Sequence arm**
- MMseqs2 PSI-search vs SCOPe ASTRAL 95%
- 35,580 representative domains
- Returns sccs + PDB links + organism
- 3 configurable iterations
        """)
    with c2:
        st.markdown("""
**🏗 Structure arm**
- ESMFold v1 predicts 3D structure
- Foldseek searches vs CATH50 v4.3.0
- Returns CATH code + lDDT confidence
- Recovers dark proteome domains
        """)
    with c3:
        st.markdown("""
**🔬 Profile arm (new)**
- HHblits builds MSA vs UniClust30
- HHsearch profiles vs PDB70 (92k HMMs)
- Finds homologs at 10–15% identity
- Twilight zone recovery
        """)

    st.subheader("Ensemble voting")
    st.markdown("""
Each arm independently predicts a SCOPe class letter (a/b/c/d/e/f).
Hits from different arms are matched by query region overlap (≥50% of shorter region).

| Votes | Confidence | Evidence tag |
|---|---|---|
| 3 arms agree | 🟢 Highest | `all_three` |
| 2 arms agree | 🟡 High | `two_arms` |
| Seq only | 🔵 Medium | `seq_only` |
| Struct only | 🟠 Medium | `struct_only` — dark proteome |
| Profile only | 🟣 Low | `hhblits_only` — twilight zone |

Score = 0.35 × seq_norm + 0.40 × struct_norm + 0.25 × HH_norm (normalised −log₁₀ e-value).
    """)

    st.subheader("References")
    st.markdown("""
- Schmidt S, Bork P, Dandekar T (2002) *J Chem Inf Comput Sci* 42:405–407
- van Kempen M et al. (2023) *Nat Biotechnol* doi:10.1038/s41587-023-01773-0
- Lin Z et al. (2023) *Science* 379:1123–1130
- Fox NK, Brenner SE, Chandonia JM (2014) *Nucleic Acids Res* 42:D304–309
- Sillitoe I et al. (2021) *Nucleic Acids Res* 49:D266–273
- Steinegger M & Söding J (2017) *Nat Methods* 14:1101–1102
    """)
