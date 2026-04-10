"""
AnDOM 2.0 — Streamlit web application.
Default view: compact domain summary (Prof. Dandekar request).
Expert mode: full three-arm parallel columns.
"""
import sys
import os
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
    fetch_pdb_function,
)
from batch.processor import (
    BatchManager, parse_fasta_text, validate_sequences,
    MAX_SEQUENCES, MAX_SEQ_LEN, MIN_SEQ_LEN,
)

st.set_page_config(page_title="AnDOM 2.0", layout="wide")

@st.cache_resource
def get_batch_manager():
    return BatchManager(output_root=os.environ.get("ANDOM_OUTPUT_DIR", "output"))

batch_mgr = get_batch_manager()

BATCH_EXAMPLES = {
    "dark_proteome": {
        "label": "Dark proteome",
        "fasta": (
            ">SARS-CoV-2_nsp7\nMSKMVLSGGFGAGVQFMNLAKSTGPVVAAAVNGLMNLFGQKSPKTVNLNKDYIQRDTGEALNKILQDYINGVAPKALENRQLAGLRQLQEQR\n"
            ">Phage_T4_Soc\nMSDNKQIKAIVESVKDKLTSIKTSNDKLNQEADLIAKNGKNAISGVLENGKAEITKLQEELAKKAGVSTLSADDLAKKKNTDLVSIDQAKLKAAK\n"
            ">Methanogen_hypothetical\nMKILIVDDHPVVREGILEYLLSAEGYEVVCAEDGQEALDIYEDHPDLVLMDLMMPGMDGFELCRQIRQLDPRIPVLMLTAKDDEYDKVLGLEIGADDYVTKPFSTREELLARIRAHL\n"
        ),
        "desc": "Dark proteome — sequence search finds nothing; structural arm recovers CATH domain.",
        "use_structure": True,
    },
    "viral_phage": {
        "label": "Viral / phage",
        "fasta": (
            ">SARS-CoV-2_nsp7\nMSKMVLSGGFGAGVQFMNLAKSTGPVVAAAVNGLMNLFGQKSPKTVNLNKDYIQRDTGEALNKILQDYINGVAPKALENRQLAGLRQLQEQR\n"
            ">SARS-CoV-2_nsp8\nMIAGGHYVFKEIVMKDPEKFNEALKMLPIDGETVIAEQIAGLKNTLKYLRKLEKDLALKLNHITNDMSSEMAKQYKEYVNKVLPQLENFEDLTKLK\n"
        ),
        "desc": "SARS-CoV-2 RNA polymerase cofactors nsp7+nsp8.",
        "use_structure": True,
    },
    "multidomain": {
        "label": "Multi-domain",
        "fasta": (
            ">Src_SH2_SH3\nMGSNKSKPKDASQRRRSLEPAENVHGAGGGAFPASQTPSKPASADGHRGPSAAFAPAAAEKVLFGGFNSSDTVTSPQRAGPLAGGVTTFVALYDYESRTETDLSFKKGERLQIVNNTEGDWWLAHSLSTGQTGYIPSNYW\n"
            ">p53_DBD\nSVVRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGQMNRRPILTIITLEDSSGKLLGRNSFEVRVCACPGRDRRTEEENLRKKGEVVAPQHL\n"
        ),
        "desc": "Multi-domain proteins — ensemble shows both domains.",
    },
}

# ── helpers ───────────────────────────────────────────────────────────────────
def fmt_e(v) -> str:
    try:
        return f"{float(v):.2e}"
    except Exception:
        return "—"

def pdb_from_scope(domain: str) -> str:
    try:
        return domain[1:5].lower()
    except Exception:
        return ""

def clean_cath_domain(d: str) -> str:
    if str(d).startswith("af_"):
        parts = str(d).split("_")
        return f"AlphaFold:{parts[1]}" if len(parts) > 1 else str(d)
    return str(d)

def lddt_color(lddt: float) -> str:
    if lddt >= 0.9:  return "#1D9E75"
    if lddt >= 0.7:  return "#639922"
    if lddt >= 0.5:  return "#BA7517"
    return "#A32D2D"


# ── helper: single map row ───────────────────────────────────────────────────
def _map_row(label, color, segs, seq_len, sublabel=""):
    content = ""
    for seg in segs:
        lt  = seg.get("label","")
        tip = seg.get("tip","").replace('"',"'").replace("\n"," ")
        lbl = (f'<span style="position:absolute;left:50%;top:50%;transform:translate(-50%,-50%);'
               f'font-size:10px;color:white;font-weight:500;white-space:nowrap;'
               f'overflow:hidden;max-width:90%;text-overflow:ellipsis">{lt}</span>'
               if lt and seg["w"] > 10 else "")
        tid = f"tip_{id(seg) % 99999}"
        content += (
            f'<div style="position:absolute;top:0;left:{seg["l"]:.1f}%;width:{seg["w"]:.1f}%;'
            f'height:100%;background:{seg["c"]};opacity:{seg.get("o",0.85)};'
            f'border-radius:3px;border:1px solid rgba(255,255,255,0.5);overflow:hidden;cursor:pointer"'
            f' onmouseover="var t=document.getElementById(\'tipbar\');if(t){{t.innerHTML=\'{tip}\';t.style.display=\'block\';}}"'
            f' onmouseout="var t=document.getElementById(\'tipbar\');if(t)t.style.display=\'none\';">'
            f'{lbl}</div>'
        )
    bg = "#f0f2f6" if content else "#f8f9fa"
    sub = (f'<br><span style="font-size:9px;font-weight:400;color:#bbb">{sublabel}</span>'
           if sublabel else "")
    return (
        f'<div style="display:flex;align-items:center;margin-bottom:4px">'
        f'<span style="width:62px;font-size:11px;color:{color};font-weight:500;flex-shrink:0;line-height:1.2">{label}{sub}</span>'
        f'<div style="position:relative;height:28px;background:{bg};border-radius:6px;'
        f'flex:1;border:1px solid #e0e0e0">{content}</div></div>'
    )

def _map_tooltip_bar():
    return (
        '<div id="tipbar" style="display:none;background:#333;color:#fff;font-size:11px;'
        'padding:3px 10px;border-radius:4px;margin-bottom:5px;'
        'white-space:nowrap;overflow:hidden;text-overflow:ellipsis;max-width:100%"></div>'
    )

CATH_NAMES = {"1":"mainly alpha","2":"mainly beta","3":"alpha/beta","4":"few secondary"}

def _seg(df, seq_len, arm, mode, row_start=0, row_end=None):
    segs = []
    if df is None or len(df) == 0:
        return segs
    df_use = df.iloc[row_start:row_end] if row_end is not None else df.iloc[row_start:]
    lk = lookup.all_domains()
    from db.lookup import get_cath_code, cath_code_to_scop_class
    CATH_C = {"1":"#E85D24","2":"#1D6FAE","3":"#1D9E75","4":"#BA7517"}
    ARM_C  = {"seq":"#3B82F6","str":"#0F6E56","hh":"#7F77DD"}
    for _, row in df_use.iterrows():
        qs = int(row.get("qstart",0)); qe = int(row.get("qend",0))
        l  = (qs / seq_len) * 100
        w  = max(((qe - qs) / seq_len) * 100, 2)
        rng = f"{qs}-{qe}"
        if mode == "scope":
            if arm == "seq":
                cls  = lk.get(str(row["target"]),{}).get("cls","?")
                sccs = lk.get(str(row["target"]),{}).get("sccs","?")
                c = SCOP_COLORS.get(cls,"#888"); o = 0.9
                lbl = f"{sccs} {rng}"
                tip = f"{row['target']} | {sccs} | e={fmt_e(row['evalue'])} | {rng} aa"
            elif arm == "str":
                cc  = get_cath_code(str(row["target"]))
                cls = cath_code_to_scop_class(cc) if cc else "?"
                c = SCOP_COLORS.get(cls,"#888"); o = 0.85
                lbl = f"class {cls} {rng}"
                tip = f"class {cls} via CATH:{cc} | lDDT={float(row.get('lddt',0)):.2f} | {rng} aa"
            else:
                sccs = str(row.get("sccs","?"))
                cls  = sccs[0] if sccs not in ("—","?","") else "?"
                # Fallback for unknown: use CATH class name
                if cls == "?":
                    cc = str(row.get("cath_code",""))
                    if cc and cc not in ("—","?",""):
                        c1 = cc.split(".")[0]
                        sccs = f"CATH-{CATH_NAMES.get(c1,c1)}"
                        cls  = {"1":"a","2":"b","3":"c","4":"g"}.get(c1,"?")
                prob = float(row.get("prob",50)) / 100
                c = SCOP_COLORS.get(cls,"#888"); o = 0.35 + 0.65*prob
                lbl = f"{sccs} {rng}"
                tip = f"{row.get('hit_name','')} | {sccs} | prob={float(row.get('prob',0)):.0f}% | {rng} aa"
        elif mode == "cath":
            if arm == "str":
                cc   = get_cath_code(str(row["target"]))
                lddt = float(row.get("lddt",0.7))
                c = CATH_C.get(cc.split(".")[0] if cc else "?","#888"); o = 0.5+0.5*lddt
            elif arm == "hh":
                cc   = str(row.get("cath_code","?"))
                prob = float(row.get("prob",50)) / 100
                c = CATH_C.get(cc.split(".")[0] if cc and cc not in ("—","?") else "?","#888"); o = 0.35+0.65*prob
            else:
                cc = get_cath_code(str(row["target"]))
                c  = CATH_C.get(cc.split(".")[0] if cc else "?","#888"); o = 0.9
            lbl = f"{cc} {rng}" if cc and cc not in ("—","?") else rng
            tip = f"CATH:{cc} | {rng} aa"
        else:  # pdb
            if arm == "seq":
                pdb = pdb_from_scope(str(row["target"])); o = 0.85
                lbl = f"{pdb.upper()} {rng}"
                tip = f"{pdb.upper()} | e={fmt_e(row['evalue'])} | {float(row.get('pident',0)):.0f}%id | {rng} aa"
            elif arm == "str":
                pdb  = str(row["target"])[:4].lower()
                lddt = float(row.get("lddt",0.7)); o = 0.5+0.5*lddt
                lbl  = f"{pdb.upper()} {rng}"
                tip  = f"{pdb.upper()} | lDDT={lddt:.2f} | {rng} aa"
            else:
                pdb  = str(row.get("pdb",""))
                prob = float(row.get("prob",50)) / 100; o = 0.35+0.65*prob
                lbl  = f"{str(row.get('hit_name',pdb))[:8]} {rng}"
                tip  = f"{row.get('hit_name',pdb)} | prob={float(row.get('prob',0)):.0f}% | {rng} aa"
            c = ARM_C[arm]
        segs.append({"l":l,"w":w,"c":c,"o":o,"label":lbl,"tip":tip})
    return segs


def render_three_domain_maps(df_seq, df_str, df_hh, seq_len):
    """Three maps (SCOPe/CATH/PDB). Profile arm split into 3 sub-rows of 7 to avoid overlap."""
    n_hh = len(df_hh) if df_hh is not None else 0
    for mode, title, legend in [
        ("scope","SCOPe classification per arm",
         " ".join(f'<span style="color:{SCOP_COLORS.get(k,"#888")}">&#9632;</span> {k}' for k in SCOP_CLASSES)),
        ("cath","CATH classification per arm",
         '<span style="color:#E85D24">&#9632;</span> mainly alpha &nbsp;'
         '<span style="color:#1D6FAE">&#9632;</span> mainly beta &nbsp;'
         '<span style="color:#1D9E75">&#9632;</span> alpha/beta &nbsp;'
         '<span style="color:#BA7517">&#9632;</span> few secondary &nbsp;&middot; struct opacity = lDDT'),
        ("pdb","PDB hits per arm",
         '<span style="color:#3B82F6">&#9632;</span> seq &nbsp;'
         '<span style="color:#0F6E56">&#9632;</span> struct &nbsp;'
         '<span style="color:#7F77DD">&#9632;</span> profile (top 21, 3 rows of 7) &nbsp;&middot; hover for details'),
    ]:
        st.markdown(
            f'<div style="font-size:13px;font-weight:500;margin:10px 0 2px">{title}</div>'
            f'<div style="font-size:10px;color:#999;margin-bottom:4px">{legend}</div>',
            unsafe_allow_html=True
        )
        html  = _map_tooltip_bar()
        html += _map_row("Seq",    "#3B82F6", _seg(df_seq, seq_len, "seq", mode), seq_len)
        html += _map_row("Struct", "#0F6E56", _seg(df_str, seq_len, "str", mode), seq_len)
        # Profile: 3 sub-rows of 3 hits each (9 total) — avoids label overlap
        html += _map_row("Prof", "#7F77DD",
                         _seg(df_hh, seq_len, "hh", mode, 0, 3), seq_len, "1-3")
        html += _map_row("Prof", "#7F77DD",
                         _seg(df_hh, seq_len, "hh", mode, 3, 6) if n_hh > 3 else [],
                         seq_len, "4-6")
        html += _map_row("Prof", "#7F77DD",
                         _seg(df_hh, seq_len, "hh", mode, 6, 9) if n_hh > 6 else [],
                         seq_len, "7-9")
        st.markdown(html, unsafe_allow_html=True)


def render_compact_summary(df_seq, df_str, df_hh, fused, seq_len):
    lk = lookup.all_domains()
    from db.lookup import get_cath_code, cath_code_to_scop_class
    render_three_domain_maps(df_seq, df_str, df_hh, seq_len)
    # ── Top-1 per arm and category ────────────────────────────────────────────
    st.markdown('<div style="font-size:13px;font-weight:500;margin-bottom:10px">Top-1 per arm and category</div>', unsafe_allow_html=True)

    from db.lookup import get_cath_code, cath_code_to_scop_class
    CATH_NAMES_FULL = {
        "1": "Mainly alpha",
        "2": "Mainly beta",
        "3": "Alpha/beta (mixed barrel)",
        "4": "Few secondary structures",
    }

    def _hh_sccs_resolved(r):
        """Return (display_label, class_letter) for HHblits hit."""
        sccs   = str(r.get("sccs","?"))
        source = str(r.get("sccs_source",""))
        cc     = str(r.get("cath_code",""))
        # Direct SCOPe hit from HHsearch vs ASTRAL
        if sccs not in ("—","?","") and not sccs.startswith("~"):
            return sccs, sccs[0]
        # CATH-inferred (~c)
        if sccs.startswith("~"):
            cls = sccs[1:2]
            c1  = cc.split(".")[0] if cc and cc not in ("—","?") else "?"
            return f"~{cls} (CATH-inferred)", cls
        # Try CATH code fallback
        if cc and cc not in ("—","?",""):
            c1  = cc.split(".")[0]
            cls = {"1":"a","2":"b","3":"c","4":"g"}.get(c1,"?")
            return f"~{cls} (CATH-inferred)", cls
        return "? (novel/unannotated)", "?"

    # Column headers
    hdr, sc, tc, pc = st.columns([1.1, 2.3, 2.3, 2.3])
    hdr.markdown("")
    sc.markdown('<b style="color:#3B82F6;font-size:12px">Sequence arm</b><br><span style="font-size:10px;color:#999">MMseqs2 vs SCOPe ASTRAL</span>', unsafe_allow_html=True)
    tc.markdown('<b style="color:#0F6E56;font-size:12px">Structure arm</b><br><span style="font-size:10px;color:#999">ESMFold + Foldseek vs CATH50</span>', unsafe_allow_html=True)
    pc.markdown('<b style="color:#7F77DD;font-size:12px">Profile arm</b><br><span style="font-size:10px;color:#999">HHblits + HHsearch vs PDB70</span>', unsafe_allow_html=True)

    st.markdown('<div style="height:4px"></div>', unsafe_allow_html=True)

    # SCOPe row
    lc1, sc1, tc1, pc1 = st.columns([1.1, 2.3, 2.3, 2.3])
    lc1.markdown('<div style="font-size:11px;color:#888;padding-top:6px">SCOPe class</div>', unsafe_allow_html=True)
    with sc1:
        if df_seq is not None and len(df_seq) > 0:
            r = df_seq.iloc[0]
            sccs = lk.get(str(r["target"]),{}).get("sccs","—")
            cls  = lk.get(str(r["target"]),{}).get("cls","?")
            col  = SCOP_COLORS.get(cls,"#888")
            rng  = f"{int(r.get('qstart',0))}–{int(r.get('qend',0))} aa"
            scope_url = f"https://scop.berkeley.edu/sunid={sccs}" if sccs != "—" else "#"
            st.markdown(
                f'<span style="color:{col};font-size:16px">●</span> '
                f'<a href="{scope_url}" style="font-weight:600;font-size:13px">{sccs}</a><br>'
                f'<span style="font-size:11px;color:#666">{SCOP_CLASSES.get(cls,cls)}</span><br>'
                f'<span style="font-size:11px;color:#aaa">{rng} · e={fmt_e(r["evalue"])}</span>',
                unsafe_allow_html=True
            )
        else:
            st.caption("no sequence hits")
    with tc1:
        if df_str is not None and len(df_str) > 0:
            r   = df_str.iloc[0]
            cc  = get_cath_code(str(r["target"]))
            cls = cath_code_to_scop_class(cc) if cc else "?"
            col = SCOP_COLORS.get(cls,"#888")
            rng = f"{int(r.get('qstart',0))}–{int(r.get('qend',0))} aa"
            lc  = lddt_color(float(r.get("lddt",0)))
            st.markdown(
                f'<span style="color:{col};font-size:16px">●</span> '
                f'<span style="font-weight:600;font-size:13px">class {cls}</span><br>'
                f'<span style="font-size:11px;color:#666">{SCOP_CLASSES.get(cls,cls)} (via CATH crosswalk)</span><br>'
                f'<span style="font-size:11px;color:#aaa">{rng} · </span>'
                f'<span style="font-size:11px;color:{lc}">lDDT={float(r.get("lddt",0)):.2f}</span>',
                unsafe_allow_html=True
            )
        else:
            st.caption("no structural hits")
    with pc1:
        if df_hh is not None and len(df_hh) > 0:
            r            = df_hh.iloc[0]
            sccs_r, cls  = _hh_sccs_resolved(r)
            col          = SCOP_COLORS.get(cls,"#888")
            rng          = f"{int(r.get('qstart',0))}–{int(r.get('qend',0))} aa"
            source       = str(r.get("sccs_source",""))
            if source == "cath_inferred":
                note = '<br><span style="font-size:10px;color:#BA7517">~ class inferred from CATH (not in SCOPe ASTRAL)</span>'
            elif source == "unknown":
                note = '<br><span style="font-size:10px;color:#A32D2D">not in SCOPe or CATH — novel fold?</span>'
            else:
                note = ""
            st.markdown(
                f'<span style="color:{col};font-size:16px">●</span> '
                f'<span style="font-weight:600;font-size:13px">{sccs_r}</span><br>'
                f'<span style="font-size:11px;color:#666">{SCOP_CLASSES.get(cls,cls)}</span>'
                f'{note}<br>'
                f'<span style="font-size:11px;color:#aaa">{rng} · prob={float(r.get("prob",0)):.0f}%</span>',
                unsafe_allow_html=True
            )
        else:
            st.caption("deep search off or no hits")

    st.markdown('<div style="height:8px"></div>', unsafe_allow_html=True)

    # CATH row
    lc2, sc2, tc2, pc2 = st.columns([1.1, 2.3, 2.3, 2.3])
    lc2.markdown('<div style="font-size:11px;color:#888;padding-top:6px">CATH code</div>', unsafe_allow_html=True)
    with sc2:
        if df_seq is not None and len(df_seq) > 0:
            r   = df_seq.iloc[0]
            cc  = get_cath_code(str(r["target"]))
            rng = f"{int(r.get('qstart',0))}–{int(r.get('qend',0))} aa"
            url = f"https://www.cathdb.info/version/v4_3_0/superfamily/{cc}" if cc else "#"
            c1  = cc.split(".")[0] if cc else "?"
            st.markdown(
                f'<a href="{url}" style="font-weight:600;font-size:13px">{cc or "—"}</a><br>'
                f'<span style="font-size:11px;color:#666">{CATH_NAMES_FULL.get(c1,"—")}</span><br>'
                f'<span style="font-size:11px;color:#aaa">{rng}</span>',
                unsafe_allow_html=True
            )
        else:
            st.caption("—")
    with tc2:
        if df_str is not None and len(df_str) > 0:
            r   = df_str.iloc[0]
            cc  = get_cath_code(str(r["target"]))
            rng = f"{int(r.get('qstart',0))}–{int(r.get('qend',0))} aa"
            url = f"https://www.cathdb.info/version/v4_3_0/superfamily/{cc}" if cc else "#"
            c1  = cc.split(".")[0] if cc else "?"
            st.markdown(
                f'<a href="{url}" style="font-weight:600;font-size:13px">{cc or "—"}</a><br>'
                f'<span style="font-size:11px;color:#666">{CATH_NAMES_FULL.get(c1,"—")}</span><br>'
                f'<span style="font-size:11px;color:#aaa">{rng}</span>',
                unsafe_allow_html=True
            )
        else:
            st.caption("—")
    with pc2:
        if df_hh is not None and len(df_hh) > 0:
            r   = df_hh.iloc[0]
            cc  = str(r.get("cath_code","—"))
            rng = f"{int(r.get('qstart',0))}–{int(r.get('qend',0))} aa"
            url = f"https://www.cathdb.info/version/v4_3_0/superfamily/{cc}" if cc and cc not in ("—","?") else "#"
            c1  = cc.split(".")[0] if cc and cc not in ("—","?") else "?"
            st.markdown(
                f'<a href="{url}" style="font-weight:600;font-size:13px">{cc}</a><br>'
                f'<span style="font-size:11px;color:#666">{CATH_NAMES_FULL.get(c1,"—")}</span><br>'
                f'<span style="font-size:11px;color:#aaa">{rng}</span>',
                unsafe_allow_html=True
            )
        else:
            st.caption("—")

    st.markdown('<div style="height:8px"></div>', unsafe_allow_html=True)

    # PDB row
    lc3, sc3, tc3, pc3 = st.columns([1.1, 2.3, 2.3, 2.3])
    lc3.markdown('<div style="font-size:11px;color:#888;padding-top:6px">PDB top hit</div>', unsafe_allow_html=True)
    with sc3:
        if df_seq is not None and len(df_seq) > 0:
            r   = df_seq.iloc[0]
            pdb = pdb_from_scope(str(r["target"]))
            purl= str(r.get("PDB link", f"https://www.rcsb.org/structure/{pdb}"))
            rng = f"{int(r.get('qstart',0))}–{int(r.get('qend',0))} aa"
            desc= lk.get(str(r["target"]),{}).get("desc","")[:50]
            st.markdown(
                f'<a href="{purl}" style="font-weight:600;font-size:14px">{pdb.upper()}</a><br>'
                f'<span style="font-size:11px;color:#666">{desc}</span><br>'
                f'<span style="font-size:11px;color:#aaa">{rng} · {float(r.get("pident",0)):.0f}%id · e={fmt_e(r["evalue"])}</span>',
                unsafe_allow_html=True
            )
        else:
            st.caption("—")
    with tc3:
        if df_str is not None and len(df_str) > 0:
            r    = df_str.iloc[0]
            pdb  = str(r["target"])[:4].lower()
            purl = str(r.get("PDB link", f"https://www.rcsb.org/structure/{pdb}"))
            rng  = f"{int(r.get('qstart',0))}–{int(r.get('qend',0))} aa"
            lddt = float(r.get("lddt",0))
            lc   = lddt_color(lddt)
            desc = str(r.get("target",""))
            st.markdown(
                f'<a href="{purl}" style="font-weight:600;font-size:14px">{pdb.upper()}</a><br>'
                f'<span style="font-size:11px;color:#666">{desc}</span><br>'
                f'<span style="font-size:11px;color:#aaa">{rng} · e={fmt_e(r["evalue"])} · </span>'
                f'<span style="font-size:11px;color:{lc}">lDDT={lddt:.2f}</span>',
                unsafe_allow_html=True
            )
        else:
            st.caption("—")
    with pc3:
        if df_hh is not None and len(df_hh) > 0:
            r    = df_hh.iloc[0]
            pdb  = str(r.get("pdb",""))
            purl = str(r.get("PDB link", f"https://www.rcsb.org/structure/{pdb}"))
            rng  = f"{int(r.get('qstart',0))}–{int(r.get('qend',0))} aa"
            name = str(r.get("hit_name", pdb.upper()))
            st.markdown(
                f'<a href="{purl}" style="font-weight:600;font-size:14px">{name}</a><br>'
                f'<span style="font-size:11px;color:#666">PDB70 hit</span><br>'
                f'<span style="font-size:11px;color:#aaa">{rng} · prob={float(r.get("prob",0)):.0f}% · e={fmt_e(r["evalue"])}</span>',
                unsafe_allow_html=True
            )
        else:
            st.caption("—")

    st.markdown('<div style="height:4px"></div>', unsafe_allow_html=True)
    st.caption("Each arm searches a different database — PDB IDs may differ but represent the same fold. Ensemble matches by query region overlap.")


    # Ensemble verdict — prominent card with links
    if not fused.empty:
        st.markdown("---")
        top  = fused.iloc[0]
        ev   = top.get("evidence",""); votes = int(top.get("votes",0))
        conf_color = {"all_three":"#1D9E75","two_arms":"#BA7517"}.get(ev,"#888780")
        conf_label = {
            "all_three":    "All three arms agree — highest confidence",
            "two_arms":     "Two arms agree — high confidence",
            "seq_only":     "Sequence arm only",
            "struct_only":  "Structure arm only (dark proteome)",
            "hhblits_only": "Profile arm only (twilight zone)",
        }.get(ev, ev)

        # Clean consensus values — strip ~ prefix for display
        agreed_sccs_raw = top.get("agreed_sccs", top.get("sccs","—"))
        agreed_sccs = agreed_sccs_raw.lstrip("~") if agreed_sccs_raw else "—"
        agreed_cath = top.get("agreed_cath", top.get("cath_code","—"))
        score       = float(top.get("ensemble_score",0))
        arms_cls    = top.get("arms_classes","")

        # Links
        sccs_url = f"https://scop.berkeley.edu/search?sortField=10&val={agreed_sccs}" if agreed_sccs not in ("—","?") else "#"
        cath_url = f"https://www.cathdb.info/version/v4_3_0/superfamily/{agreed_cath}" if agreed_cath not in ("—","?","") else "#"

        # Top PDB hit for link
        top_pdb_scope = pdb_from_scope(str(top.get("scope_domain","")))
        top_pdb_hh    = str(top.get("hh_hit",""))[:4].lower()
        top_pdb_str   = str(top.get("cath_domain",""))[:4].lower()
        top_pdb       = top_pdb_scope or top_pdb_str or top_pdb_hh
        pdb_url_top   = f"https://www.rcsb.org/structure/{top_pdb.upper()}" if top_pdb else "#"

        st.markdown(
            f'<div style="background:var(--color-background-secondary);border-left:5px solid {conf_color};'
            f'border-radius:0 12px 12px 0;padding:16px 20px;margin-top:8px">'
            f'<div style="font-weight:700;font-size:16px;color:{conf_color}">'
            f'Ensemble: {votes}/3 votes &nbsp;·&nbsp; {conf_label}</div>'
            f'<div style="display:flex;gap:32px;margin-top:10px;flex-wrap:wrap">'
            f'<div><div style="font-size:10px;color:#aaa;text-transform:uppercase;letter-spacing:0.5px">SCOPe class</div>'
            f'<a href="{sccs_url}" style="font-size:16px;font-weight:600;text-decoration:none">{agreed_sccs}</a><br>'
            f'<span style="font-size:11px;color:#888">{SCOP_CLASSES.get(agreed_sccs[0] if agreed_sccs and agreed_sccs not in ("—","?") else "?","")}</span></div>'
            f'<div><div style="font-size:10px;color:#aaa;text-transform:uppercase;letter-spacing:0.5px">CATH superfamily</div>'
            f'<a href="{cath_url}" style="font-size:16px;font-weight:600;text-decoration:none">{agreed_cath}</a></div>'
            f'<div><div style="font-size:10px;color:#aaa;text-transform:uppercase;letter-spacing:0.5px">Top PDB hit</div>'
            f'<a href="{pdb_url_top}" style="font-size:16px;font-weight:600;text-decoration:none">{top_pdb.upper() if top_pdb else "—"}</a></div>'
            f'<div><div style="font-size:10px;color:#aaa;text-transform:uppercase;letter-spacing:0.5px">Score</div>'
            f'<span style="font-size:16px;font-weight:600">{score:.3f}</span></div>'
            f'</div>'
            f'<div style="font-size:11px;color:#aaa;margin-top:8px">'
            f'Arm classes: {arms_cls} &nbsp;·&nbsp; '
            f'Score = weighted e-value norm + class agreement bonus when ≥2 arms agree'
            f'</div></div>',
            unsafe_allow_html=True
        )

    # Functional annotation — fixed AlphaFold URL (needs UniProt ID not PDB ID)
    top_pdb = ""
    if not fused.empty:
        top_pdb = pdb_from_scope(str(fused.iloc[0].get("scope_domain","")))
        if not top_pdb: top_pdb = str(fused.iloc[0].get("hh_hit",""))[:4].lower()
        if not top_pdb: top_pdb = str(fused.iloc[0].get("cath_domain",""))[:4].lower()
    if top_pdb and len(top_pdb) == 4:
        st.markdown(
            f'<div style="font-size:11px;color:#888;margin:10px 0 2px">'
            f'Functional annotation for top ensemble hit ({top_pdb.upper()})</div>',
            unsafe_allow_html=True
        )
        with st.expander(f"Load UniProt · Pfam · AlphaFold structure — {top_pdb.upper()}", expanded=False):
            with st.spinner(f"Fetching {top_pdb.upper()}…"):
                info = fetch_pdb_function(top_pdb)
            if info.get("uniprot_id"):
                uid = info["uniprot_id"]
                st.markdown(
                    f'**UniProt:** [{uid}](https://www.uniprot.org/uniprot/{uid}) &nbsp;·&nbsp; '
                    f'**Gene:** {info.get("gene","—")} &nbsp;·&nbsp; '
                    f'**Organism:** {info.get("organism","—")}'
                )
                if info.get("function"):
                    st.markdown(f'_{info["function"][:400]}_')
                if info.get("pfam"):
                    st.markdown("**Pfam:** " + " · ".join(f'`{p["pfam_id"]}` {p["name"]}' for p in info["pfam"][:4]))
                # AlphaFold uses UniProt accession — this is the correct URL
                st.markdown(
                    f'**AlphaFold structure** ([{uid}](https://alphafold.ebi.ac.uk/entry/{uid})):'
                )
                st.markdown(
                    f'<iframe src="https://alphafold.ebi.ac.uk/entry/{uid}" '
                    f'width="100%" height="520" style="border:none;border-radius:8px;margin-top:4px" '
                    f'title="AlphaFold structure for {uid}"></iframe>',
                    unsafe_allow_html=True
                )
                st.markdown(
                    f'[View {top_pdb.upper()} on RCSB PDB](https://www.rcsb.org/structure/{top_pdb.upper()})'
                )
            else:
                st.caption(f"No UniProt mapping found for {top_pdb.upper()}.")
                st.markdown(f"[View {top_pdb.upper()} on RCSB PDB](https://www.rcsb.org/structure/{top_pdb.upper()})")



# ── expert arm panel ──────────────────────────────────────────────────────────
def render_arm_panel(title, color, df, arm, top_n=7):
    lk = lookup.all_domains()
    st.markdown(f'<div style="font-weight:600;font-size:14px;color:{color};border-left:3px solid {color};padding-left:8px;margin-bottom:8px">{title}</div>', unsafe_allow_html=True)
    if df is None or len(df) == 0:
        st.info("No hits found.")
        return
    hits = df.head(top_n)

    st.markdown("**SCOPe classes**")
    if arm == "seq":
        for _, r in hits.iterrows():
            sccs = lk.get(str(r["target"]),{}).get("sccs","—")
            cls  = lk.get(str(r["target"]),{}).get("cls","?")
            col  = SCOP_COLORS.get(cls,"#888")
            st.markdown(f'<span style="color:{col}">●</span> `{sccs}` {SCOP_CLASSES.get(cls,cls)} · e={fmt_e(r["evalue"])}', unsafe_allow_html=True)
    elif arm == "str":
        from db.lookup import get_cath_code, cath_code_to_scop_class
        for _, r in hits.iterrows():
            cc   = get_cath_code(str(r["target"]))
            cls  = cath_code_to_scop_class(cc) if cc else "?"
            col  = SCOP_COLORS.get(cls,"#888")
            lddt = float(r.get("lddt",0))
            lc   = lddt_color(lddt)
            st.markdown(f'<span style="color:{col}">●</span> `{cc}` → {SCOP_CLASSES.get(cls,cls)} · <span style="color:{lc}">lDDT={lddt:.2f}</span>', unsafe_allow_html=True)
    else:
        for _, r in hits.iterrows():
            sccs = str(r.get("sccs","—"))
            cls  = sccs[0] if sccs not in ("—","?","") else "?"
            col  = SCOP_COLORS.get(cls,"#888")
            st.markdown(f'<span style="color:{col}">●</span> `{sccs}` {SCOP_CLASSES.get(cls,cls)} · prob={float(r.get("prob",0)):.0f}%', unsafe_allow_html=True)

    st.divider()
    st.markdown("**CATH codes**")
    from db.lookup import get_cath_code
    seen = []
    for _, r in hits.iterrows():
        pdb = str(r.get("pdb","")) if arm == "hh" else str(r.get("target",""))
        cc  = get_cath_code(pdb)
        if cc and cc not in seen:
            seen.append(cc)
            st.markdown(f"[`{cc}`](https://www.cathdb.info/version/v4_3_0/superfamily/{cc})")
        if len(seen) >= top_n: break

    st.divider()
    st.markdown("**Top PDB hits**")
    for i, (_, r) in enumerate(hits.iterrows()):
        if arm == "seq":
            pdb_id = pdb_from_scope(str(r["target"]))
            purl   = str(r.get("PDB link", f"https://www.rcsb.org/structure/{pdb_id}"))
            extra  = f"e={fmt_e(r['evalue'])} · {float(r.get('pident',0)):.0f}%id · {int(r.get('qstart',0))}–{int(r.get('qend',0))} aa"
            desc   = lk.get(str(r["target"]),{}).get("desc","")[:40]
        elif arm == "str":
            pdb_id = str(r["target"])[:4].lower()
            purl   = str(r.get("PDB link", f"https://www.rcsb.org/structure/{pdb_id}"))
            lddt   = float(r.get("lddt",0))
            extra  = f"e={fmt_e(r['evalue'])} · lDDT={lddt:.2f} · {int(r.get('qstart',0))}–{int(r.get('qend',0))} aa"
            desc   = clean_cath_domain(str(r["target"]))
        else:
            pdb_id = str(r.get("pdb",""))
            purl   = str(r.get("PDB link", f"https://www.rcsb.org/structure/{pdb_id}"))
            extra  = f"e={fmt_e(r['evalue'])} · prob={float(r.get('prob',0)):.0f}% · {int(r.get('qstart',0))}–{int(r.get('qend',0))} aa"
            desc   = str(r.get("hit_name",""))
        c1, c2 = st.columns([4,1])
        with c1:
            st.markdown(f"**{i+1}.** [{pdb_id.upper()}]({purl}) {desc}")
            st.caption(extra)
        with c2:
            pass  # expander below handles func
        with st.expander(f"Function: {pdb_id.upper()}", expanded=False):
            with st.spinner(f"Fetching {pdb_id.upper()}…"):
                info = fetch_pdb_function(pdb_id)
            if info.get("uniprot_id"):
                uid = info["uniprot_id"]
                st.markdown(f"**[{uid}](https://www.uniprot.org/uniprot/{uid})** · {info.get('gene','—')} · {info.get('organism','—')}")
                if info.get("function"): st.caption(info["function"][:250])
                if info.get("pfam"): st.caption("Pfam: " + " · ".join(f'{p["pfam_id"]} {p["name"]}' for p in info["pfam"][:3]))
            else:
                st.caption("No UniProt mapping.")


# ── sidebar ───────────────────────────────────────────────────────────────────
with st.sidebar:
    st.header("Search Parameters")
    evalue      = st.select_slider("E-value cutoff", options=[1e-30,1e-20,1e-10,1e-5,1e-3,0.01,0.1,1.0], value=1e-3)
    iterations  = st.slider("PSI-MMseqs2 iterations", 1, 5, 3)
    max_hits    = st.slider("Max hits shown", 5, 100, 30)
    use_struct  = st.toggle("Structural search (ESMFold + Foldseek)", value=True)
    use_hhblits = st.toggle("🔬 Deep search (HHblits twilight zone)", value=False, key="use_hhblits",
                    help="Profile-profile search via HHblits+UniClust30. 10-15% identity. ~2 min extra.")
    if use_hhblits: st.caption("HHblits arm: PDB70 database loading...")
    st.divider()
    st.markdown("**SCOP class colours**")
    for k, v in SCOP_CLASSES.items():
        st.markdown(f'<span style="background:{SCOP_COLORS[k]};padding:2px 8px;border-radius:3px;color:white;font-size:12px">{k}</span> {v}', unsafe_allow_html=True)
    st.divider()
    st.caption(f"SCOPe 2.08 — {len(lookup.all_domains()):,} domains")
    st.caption("CATH50 structural DB")
    st.divider()
    st.caption(f"Batch limit: **{MAX_SEQUENCES}** seq · **{MAX_SEQ_LEN}** aa max")


# ── batch cards ───────────────────────────────────────────────────────────────
def render_batch_cards(df_res, job_id):
    from db.lookup import get_cath_code, cath_code_to_scop_class
    if df_res.empty:
        st.warning("No domain hits found.")
        return
    col_id = "query_id" if "query_id" in df_res.columns else "query"
    for qid in [q for q in df_res[col_id].unique() if str(q) != "nan"]:
        q_df     = df_res[df_res[col_id] == qid]
        seq_hits = q_df[q_df["source"]=="SCOPe_sequence"] if "source" in q_df.columns else q_df
        str_hits = q_df[q_df["source"]=="CATH_structure"]  if "source" in q_df.columns else pd.DataFrame()
        has_seq  = len(seq_hits) > 0
        has_str  = len(str_hits) > 0
        if has_seq and has_str:   icon, label = "🟢", "Both arms"
        elif has_seq:             icon, label = "🔵", "Sequence only"
        elif has_str:             icon, label = "🟠", "Structure only (dark proteome)"
        else:                     icon, label = "⚪", "No hits"
        with st.container(border=True):
            seq_len = int(seq_hits.iloc[0]["qend"]) if has_seq and "qend" in seq_hits.columns else (int(str_hits.iloc[0]["qend"]) if has_str and "qend" in str_hits.columns else "?")
            st.markdown(f"### {icon} {qid}")
            st.caption(f"{label} · {seq_len} aa")
            col1, col2 = st.columns(2)
            with col1:
                st.markdown("**Sequence arm (SCOPe)**")
                if has_seq:
                    top = seq_hits.iloc[0]
                    pdb_link = top.get("PDB link",""); pdb_id = pdb_link.split("/")[-1] if pdb_link else "—"
                    st.markdown(f"`{top.get('sccs','—')}` — {top.get('class_name','—')}")
                    if pdb_link: st.caption(f"[{pdb_id}]({pdb_link}) · e={fmt_e(top.get('evalue',0))}")
                else: st.info("No sequence homologs")
            with col2:
                st.markdown("**Structure arm (CATH50)**")
                if has_str:
                    top_s = str_hits.iloc[0]
                    cc    = get_cath_code(str(top_s.get("target","—")))
                    cls   = cath_code_to_scop_class(cc) if cc else "?"
                    pdb_link_s = top_s.get("PDB link",""); pdb_id_s = pdb_link_s.split("/")[-1] if pdb_link_s else "—"
                    st.markdown(f"`{cc}` → class {cls}")
                    if pdb_link_s: st.caption(f"[{pdb_id_s}]({pdb_link_s}) · lDDT={float(top_s.get('lddt',0)):.2f}")
                else: st.info("Structural search not run or no hits")


page = st.tabs(["Search", "Batch", "Methods"])

# ══════════════════════════════════════════════════════════════════════════════
# TAB 0 — Search
# ══════════════════════════════════════════════════════════════════════════════
with page[0]:
    st.title("AnDOM 2.0 — Structural Domain Finder")
    st.caption("Sequence + Structure + Profile ensemble | SCOPe 2.08 + CATH50 + PDB70 | MMseqs2 + ESMFold + Foldseek + HHblits | Dandekar Lab Wuerzburg")

    st.markdown("**Try an example:**")
    ex_cols = st.columns(5)
    for i, (key, ex) in enumerate(EXAMPLES.items()):
        if ex_cols[i].button(ex["label"], use_container_width=True, key=f"btn_{key}"):
            st.session_state["injected_seq"]  = ex["seq"]
            st.session_state["injected_desc"] = ex["desc"]
            st.rerun()

    if "injected_seq" in st.session_state:
        st.info(st.session_state.get("injected_desc",""))

    seq_input = st.text_area("Paste protein sequence (FASTA or raw):", height=120, value=st.session_state.get("injected_seq",""))

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
                st.warning(f"Sequence {n} aa — ESMFold limit {ESMFOLD_MAXLEN} aa. Structural search disabled.")

            df_seq = df_str = df_hh = None
            _tid = _uuid.uuid4().hex[:8]

            col_s, col_t, col_h = st.columns(3)
            with col_s:
                with st.spinner("Sequence search…"):
                    df_seq, err = seq_search.run(fasta, evalue=evalue, iterations=iterations, tmp_dir=f"/output/tmp/{_tid}")
                if err: st.error(f"Sequence search failed: {err}")
                elif df_seq is None or len(df_seq)==0: st.warning("No sequence hits.")
                else:
                    df_seq = df_seq.head(max_hits)
                    st.success(f"Sequence: {len(df_seq)} SCOPe hits")

            with col_t:
                if use_struct_run:
                    with st.spinner("Structure search…"):
                        df_str, err2 = str_search.run(clean_seq)
                    if err2: st.error(f"Structure search failed: {err2}")
                    elif df_str is None or len(df_str)==0: st.warning("No structural hits.")
                    else:
                        df_str = df_str.head(max_hits)
                        st.success(f"Structure: {len(df_str)} CATH hits")
                else:
                    st.info("Structural search disabled.")

            with col_h:
                if use_hhblits:
                    with st.spinner("Deep search — HHblits + PDB70…"):
                        from search.sequence_hhblits import run_hhblits
                        _hid = _uuid.uuid4().hex[:8]
                        df_hh, err_hh = run_hhblits(fasta, iterations=3, threads=8, tmp_dir=f"/output/tmp/hh_{_hid}")
                    if err_hh: st.warning(f"HHblits: {err_hh}")
                    elif df_hh is not None and len(df_hh)>0:
                        df_hh["source"] = "HHblits_PDB70"
                        st.success(f"HHblits: {len(df_hh)} PDB70 hits")
                    else: st.info("HHblits: no hits found")
                else:
                    st.info("Deep search disabled.")

            df_seq, df_str = add_scores(df_seq, df_str)
            fused = fuse_results_three(df_seq, df_str, df_hh)

            bar_len = max(
                df_seq["qend"].max() if df_seq is not None and len(df_seq)>0 else 1,
                df_str["qend"].max() if df_str is not None and len(df_str)>0 else 1,
                df_hh["qend"].max()  if df_hh  is not None and len(df_hh)>0  else 1,
                n,
            )

            # ── DEFAULT: COMPACT SUMMARY ──────────────────────────────────────
            st.subheader("Domain annotation")
            render_compact_summary(df_seq, df_str, df_hh, fused, bar_len)

            # ── EXPERT MODE ───────────────────────────────────────────────────
            with st.expander("Expert mode — detailed three-arm results", expanded=False):
                st.caption("Each arm independently: SCOPe class · CATH code · top 7 PDB hits · functional lookup")
                arm_col1, arm_col2, arm_col3 = st.columns(3)
                with arm_col1: render_arm_panel("Sequence arm — MMseqs2 / SCOPe", "#3B82F6", df_seq, "seq")
                with arm_col2: render_arm_panel("Structure arm — ESMFold / CATH50", "#0F6E56", df_str, "str")
                with arm_col3: render_arm_panel("Profile arm — HHblits / PDB70", "#7F77DD", df_hh, "hh")

            # ── ENSEMBLE TABLE ────────────────────────────────────────────────
            if not fused.empty:
                with st.expander("Ensemble ranked hits table", expanded=False):
                    ev_icon = {"all_three":"🟢","two_arms":"🟡","seq_only":"🔵","struct_only":"🟠","hhblits_only":"🟣"}
                    fused["ev"] = fused["evidence"].map(ev_icon)
                    fused["cath_d"] = fused["cath_domain"].apply(clean_cath_domain)
                    for col in ["seq_evalue","struct_evalue","hh_evalue"]:
                        if col in fused.columns:
                            fused[col] = fused[col].apply(lambda x: fmt_e(x) if x is not None and str(x) not in ("None","nan") else "—")
                    fused["hh_prob"] = fused["hh_prob"].apply(lambda x: f"{float(x):.0f}%" if x is not None and str(x) not in ("None","nan") else "—")
                    n_all=(fused["evidence"]=="all_three").sum(); n_two=(fused["evidence"]=="two_arms").sum()
                    n_seq=(fused["evidence"]=="seq_only").sum();  n_st=(fused["evidence"]=="struct_only").sum()
                    n_hh=(fused["evidence"]=="hhblits_only").sum()
                    cc1,cc2,cc3,cc4,cc5 = st.columns(5)
                    cc1.metric("🟢 All three",n_all); cc2.metric("🟡 Two arms",n_two)
                    cc3.metric("🔵 Seq only",n_seq);  cc4.metric("🟠 Struct only",n_st); cc5.metric("🟣 Profile only",n_hh)
                    st.dataframe(
                        fused[["rank","ev","scope_domain","sccs","cath_d","cath_code","hh_hit","votes","ensemble_score","qstart","qend","seq_evalue","struct_evalue","hh_evalue","hh_prob","lddt"]
                        ].rename(columns={"ev":"","scope_domain":"SCOPe domain","sccs":"SCOP class","cath_d":"CATH domain","cath_code":"CATH code","hh_hit":"HHblits hit","votes":"Votes","ensemble_score":"Score","qstart":"Start aa","qend":"End aa","seq_evalue":"Seq e-val","struct_evalue":"Struct e-val","hh_evalue":"HH e-val","hh_prob":"HH prob","lddt":"lDDT"}),
                        use_container_width=True, hide_index=True,
                    )
                    st.download_button("Download TSV", fused.to_csv(sep="\t",index=False), file_name="AnDOM_results.tsv", mime="text/tab-separated-values")

    st.divider()
    st.caption("AnDOM 2.0 | SCOPe 2.08 ASTRAL 95% (MMseqs2 PSI) + CATH50 (ESMFold + Foldseek) + PDB70 (HHblits) | Dandekar Lab Wuerzburg")

# ══════════════════════════════════════════════════════════════════════════════
# TAB 1 — Batch
# ══════════════════════════════════════════════════════════════════════════════
with page[1]:
    st.title("Batch Processing")
    st.caption(f"Up to **{MAX_SEQUENCES}** sequences · {MIN_SEQ_LEN}–{MAX_SEQ_LEN} aa · background jobs")

    with st.expander("Why AnDOM 2.0 alongside AlphaFold + Foldseek?", expanded=False):
        st.markdown("""
| Challenge | AlphaFold | Foldseek | InterPro | **AnDOM 2.0** |
|---|---|---|---|---|
| Assigns domain to SCOP/CATH | ❌ | ✅ if DB hit | seq only | **seq + struct + profile** |
| Dark proteome | structure only | ❌ | ❌ | **structural arm** |
| Twilight zone 10–15% id | ❌ | ❌ | ❌ | **HHblits arm** |
| Fast-evolving viral | structure only | weak | ❌ | **ensemble** |
        """)

    st.markdown("**Load example set:**")
    ex_cols = st.columns(len(BATCH_EXAMPLES))
    for i, (key, ex) in enumerate(BATCH_EXAMPLES.items()):
        if ex_cols[i].button(ex["label"], use_container_width=True, key=f"bex_{key}"):
            st.session_state["batch_text"]  = ex["fasta"]
            st.session_state["_batch_desc"] = ex["desc"]
            if ex.get("use_structure", False): st.session_state["b_struct"] = True
    if "_batch_desc" in st.session_state: st.info(st.session_state["_batch_desc"])

    with st.expander("Submit new batch job", expanded=True):
        b_upload = st.file_uploader("Upload multi-FASTA", type=["fa","fasta","txt"], key="batch_upload")
        b_text   = st.text_area("…or paste multi-FASTA here", height=180, key="batch_text")
        b_struct = st.toggle("Include structural search (≤400 aa only)", value=False, key="b_struct")
        raw = b_upload.read().decode() if b_upload else (b_text if b_text.strip() else None)
        if st.button("Submit Batch Job", type="primary"):
            if not raw: st.warning("Please provide FASTA sequences.")
            else:
                seqs = parse_fasta_text(raw); errs = validate_sequences(seqs)
                if errs:
                    for e in errs: st.error(e)
                else:
                    import tempfile
                    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as tmp:
                        tmp.write(raw); tmp_path = tmp.name
                    job_id = batch_mgr.submit(tmp_path, evalue=evalue, iterations=iterations, use_structure=b_struct)
                    st.success(f"Job **{job_id}** submitted — {len(seqs)} sequences.")
                    st.session_state.pop("_batch_desc", None)

    st.subheader("Jobs")
    jobs = batch_mgr.list_jobs()
    if not jobs: st.info("No batch jobs yet.")
    else:
        if st.button("Refresh"): st.rerun()
        status_icon = {"queued":"🕐","running":"⏳","done":"✅","failed":"❌","cancelled":"🚫"}
        for job in jobs:
            s = job.get("status","unknown")
            with st.container(border=True):
                c1,c2,c3 = st.columns([3,2,2])
                c1.markdown(f"{status_icon.get(s,'❓')} **{job['job_id']}** — `{s}`")
                total=job.get("total",0); prog=job.get("progress",0)
                c2.caption(f"{prog}/{total} sequences" if total else job.get("submitted","")[:19])
                if s == "done":
                    df_res = batch_mgr.results(job["job_id"])
                    if not df_res.empty:
                        c3.download_button("⬇ TSV", df_res.to_csv(sep="\t",index=False),
                            file_name=f"AnDOM_batch_{job['job_id']}.tsv", mime="text/tab-separated-values", key=f"dl_{job['job_id']}")
                        render_batch_cards(df_res, job["job_id"])
                if s in ("queued","running"):
                    if c3.button("Cancel", key="cancel_"+job["job_id"]): batch_mgr.cancel(job["job_id"]); st.rerun()
                if job.get("error"): st.error(job["error"])

# ══════════════════════════════════════════════════════════════════════════════
# TAB 2 — Methods
# ══════════════════════════════════════════════════════════════════════════════
with page[2]:
    st.title("Methods — AnDOM 2.0")
    st.markdown("AnDOM 2.0 updates the original AnDOM server (Schmidt, Bork & Dandekar, *J Chem Inf Comput Sci* 2002) with modern tools, a structural search layer, and a profile-profile twilight zone arm.")
    st.subheader("Pipeline comparison: AnDOM 2002 vs AnDOM 2.0")
    for candidate in [Path(__file__).parent/"docs"/"AnDOM_old_vs_new_method.svg", Path(__file__).parent/"AnDOM_old_vs_new_method.svg", Path("/app/docs/AnDOM_old_vs_new_method.svg")]:
        if candidate.exists(): st.image(str(candidate)); break
    else: st.warning("SVG diagram not found.")

    st.subheader("Three-arm ensemble")
    c1,c2,c3 = st.columns(3)
    with c1: st.markdown("**Sequence arm**\n- MMseqs2 PSI vs SCOPe ASTRAL 95%\n- 35,580 domains · sccs + PDB + organism")
    with c2: st.markdown("**Structure arm**\n- ESMFold v1 + Foldseek vs CATH50\n- CATH code + lDDT · dark proteome recovery")
    with c3: st.markdown("**Profile arm**\n- HHblits vs UniClust30 + HHsearch PDB70\n- 10–15% identity · twilight zone recovery")

    st.subheader("Ensemble voting")
    st.markdown("""
| Votes | Confidence | Tag |
|---|---|---|
| 3 arms | 🟢 Highest | `all_three` |
| 2 arms | 🟡 High | `two_arms` |
| Seq only | 🔵 Medium | `seq_only` |
| Struct only | 🟠 Medium | `struct_only` |
| Profile only | 🟣 Low | `hhblits_only` |

Score = 0.35 × seq + 0.40 × struct + 0.25 × profile (normalised −log₁₀ e-value)
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
