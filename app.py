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
import hashlib, json, time
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
            ">Src_SH2_SH3\nMGSNKSKPKDASQRRRSLEPAENVHGAGGGAFPASQTPSKPASADGHRGPSAAFAPAAAEKVLFGGFNSSDTVTSPQRAGPLAGGVTTFVALYDYESRTETDLSFKKGERLQIVNNTEGDWWLAHSLSTGQTGYIPSNYWAGNEKINGQEPIPKAKIKALRQLRISDDAHERVEFDQENPAQVAREAFVSELDMHLNYGQSSANFYFKQLFKKSGESNEEVAENYLEHLGDSWRRFQPSFDPTKNPGYEIVNQIKTCTDDFRGMECGKEPYASEMTRIIREALQKLQKQRGGIEYMQKQGKDHSMSHSV\n"
            ">p53_DBD\nSVVRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGQMNRRPILTIITLEDSSGKLLGRNSFEVRVCACPGRDRRTEEENLRKKGEVVAPQHL\n"
        ),
        "desc": "Multi-domain proteins — Src SH3 (all-alpha, a.118.1) + SH2 (all-beta, b.34.2.1) + p53 DNA-binding domain (beta-sandwich).",
    },
}

# ── helpers ───────────────────────────────────────────────────────────────────
def fmt_e(v) -> str:
    try:
        return f"{float(v):.2e}"
    except Exception:
        return "—"

CACHE_DIR = Path(os.environ.get('ANDOM_OUTPUT_DIR', 'output')) / 'cache'
CACHE_TTL = 48 * 3600

def _cache_key(seq): return hashlib.md5(seq.strip().upper().encode()).hexdigest()

def _load_cache(seq):
    try:
        p = CACHE_DIR / (_cache_key(seq) + '.json')
        if not p.exists(): return None
        data = json.loads(p.read_text())
        if time.time() - data.get('ts', 0) > CACHE_TTL: return None
        return {k: pd.DataFrame(data[k]) if data.get(k) else None for k in ('df_seq','df_str','df_hh')}
    except: return None

def _save_cache(seq, df_seq, df_str, df_hh):
    try:
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        p = CACHE_DIR / (_cache_key(seq) + '.json')
        data = {'ts': time.time()}
        for k, df in (('df_seq',df_seq),('df_str',df_str),('df_hh',df_hh)):
            data[k] = df.to_dict(orient='records') if df is not None and len(df)>0 else []
        p.write_text(json.dumps(data))
    except Exception as e:
        import logging; logging.warning(f"Cache save failed: {e}")

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


def render_three_domain_maps(df_seq, df_str, df_hh, seq_len, domains=None):
    import streamlit.components.v1 as components
    import math
    from db.lookup import get_cath_code, cath_code_to_scop_class
    from search.ensemble import get_domain_clusters
    lk = lookup.all_domains()

    if domains is None:
        domains = get_domain_clusters(df_seq, df_str, df_hh, min_overlap=0.6)

    # ── segment builder ───────────────────────────────────────────────────
    def _seg(qs, qe, color, opacity, label, tip, url=''):
        l = (qs / seq_len) * 100
        w = max(((qe - qs) / seq_len) * 100, 1.0)
        ts = tip.replace('"', '&quot;').replace("'", '&#39;')
        inner = ''
        if w > 4:
            inner = ('<span style="position:absolute;left:4px;right:4px;top:50%;'
                     'transform:translateY(-50%);font-size:10px;font-weight:500;'
                     'color:rgba(255,255,255,0.95);white-space:nowrap;overflow:hidden;'
                     'text-overflow:ellipsis">' + label + '</span>')
        sty = ('position:absolute;top:1px;left:%.2f%%;width:%.2f%%;height:20px;'
               'background:%s;opacity:%.2f;border-radius:3px;overflow:hidden;'
               'cursor:default') % (l, w, color, opacity)
        return ('<div onmouseover="T(\'' + ts + '\')" onmouseout="T()"'
                ' style="' + sty + '">' + inner + '</div>')

    # ── single labeled bar row ────────────────────────────────────────────
    def _bar(arm_label, arm_color, annotation, segs_html, has_hits):
        bg = '#f0f3f8' if has_hits else '#f7f8fa'
        no_h = ('<span style="font-size:10px;color:#bbb;padding:0 8px;'
                'line-height:22px">no hits</span>') if not has_hits else ''
        return (
            '<div style="display:flex;align-items:center;margin-bottom:5px">'
            '<div style="width:110px;flex-shrink:0;padding-right:8px;text-align:right">'
            '<span style="font-size:10px;font-weight:600;color:' + arm_color + '">' + arm_label + '</span>'
            '<br><span style="font-size:9px;color:#aaa">' + annotation + '</span>'
            '</div>'
            '<div style="position:relative;height:22px;background:' + bg + ';'
            'border-radius:4px;flex:1;border:0.5px solid #e0e3ea">'
            + segs_html + no_h +
            '</div></div>'
        )

    # ── build segments for each arm ───────────────────────────────────────
    def _seq_bars(df):
        if df is None or len(df) == 0:
            return [('SCOPe', '#3B82F6', 'class · sccs', '', False),
                    ('CATH',  '#3B82F6', 'via lookup',   '', False),
                    ('PDB',   '#3B82F6', 'top hits',     '', False)]
        rows = []
        sc = ca = pb = ''
        for _hi, (_, r) in enumerate(df.head(1).iterrows()):
            qs = int(r.get('qstart',0)); qe = int(r.get('qend',0))
            tgt  = str(r['target'])
            cls  = lk.get(tgt,{}).get('cls','?')
            sccs = lk.get(tgt,{}).get('sccs','?')
            desc = lk.get(tgt,{}).get('desc','')[:30]
            ev   = float(r['evalue'])
            cc   = get_cath_code(tgt) or '?'
            pdb  = tgt[1:5].lower() if len(tgt)>=5 else ''
            col  = SCOP_COLORS.get(cls,'#888')
            op   = max(0.3, min(0.95, -math.log10(max(ev,1e-300))/30))
            _scope_url = 'https://scop.berkeley.edu/search/?key=' + sccs if sccs not in ('?','--','') else ''
            sc += _seg(qs,qe,col,op,sccs if _hi==0 else '',
                       'SCOPe: %s | e=%.1e | %d-%d aa | %s — click to open SCOPe'%(sccs,ev,qs,qe,desc), _scope_url)
            _cath_url = 'https://www.cathdb.info/version/v4_3_0/superfamily/' + cc if cc not in ('?','--','') else ''
            ca += _seg(qs,qe,col,op*0.8,cc if _hi==0 else '',
                       'CATH: %s | %s | %d-%d aa — click to open CATH'%(cc,tgt,qs,qe), _cath_url)
            _pdb_url = 'https://www.rcsb.org/structure/' + pdb.upper() if pdb else ''
            pb += _seg(qs,qe,'#3B82F6',op,pdb.upper() if _hi==0 else '',
                       'PDB: %s | e=%.1e | %d-%d aa — click to open RCSB'%(pdb.upper(),ev,qs,qe), _pdb_url)
        return [('SCOPe','#3B82F6','class · sccs',sc,True),
                ('CATH', '#3B82F6','via lookup',  ca,True),
                ('PDB',  '#3B82F6','top hits',    pb,True)]

    def _str_bars(df):
        if df is None or len(df) == 0:
            return [('SCOPe','#0F6E56','crosswalk','',False),
                    ('CATH', '#0F6E56','direct',   '',False),
                    ('PDB',  '#0F6E56','top hits', '',False)]
        sc = ca = pb = ''
        for _hi, (_, r) in enumerate(df.head(1).iterrows()):
            qs = int(r.get('qstart',0)); qe = int(r.get('qend',0))
            tgt  = str(r['target'])
            cc   = get_cath_code(tgt) or '?'
            cls  = cath_code_to_scop_class(cc) if cc and cc!='?' else '?'
            lddt = float(r.get('lddt',0.7))
            pdb  = tgt[:4].lower()
            col  = SCOP_COLORS.get(cls,'#888')
            op   = 0.4+0.6*lddt
            sc += _seg(qs,qe,col,op,'cls:'+cls if _hi==0 else '',
                       'CATH→SCOPe: class %s | lDDT=%.2f | %d-%d aa'%(cls,lddt,qs,qe))
            _cath_url2 = 'https://www.cathdb.info/version/v4_3_0/superfamily/' + cc if cc not in ('?','--','') else ''
            ca += _seg(qs,qe,'#0F6E56',op,cc if _hi==0 else '',
                       'CATH: %s | lDDT=%.2f | %d-%d aa — click to open CATH'%(cc,lddt,qs,qe), _cath_url2)
            pb += _seg(qs,qe,'#0F6E56',op*0.85,pdb.upper() if _hi==0 else '',
                       'PDB: %s | lDDT=%.2f | %d-%d aa'%(pdb.upper(),lddt,qs,qe))
        return [('SCOPe','#0F6E56','crosswalk',sc,True),
                ('CATH', '#0F6E56','direct',   ca,True),
                ('PDB',  '#0F6E56','top hits', pb,True)]

    def _hh_bars(df):
        if df is None or len(df) == 0:
            return [('SCOPe','#7F77DD','direct sccs','',False),
                    ('CATH', '#7F77DD','via PDB hit','',False),
                    ('PDB',  '#7F77DD','top hits',   '',False)]
        sc = ca = pb = ''
        for _hi, (_, r) in enumerate(df.head(1).iterrows()):
            qs   = int(r.get('qstart',0)); qe = int(r.get('qend',0))
            sccs = str(r.get('sccs','?'))
            cls  = sccs[0] if sccs not in ('--','?','') and not sccs.startswith('~') else '?'
            prob = float(r.get('prob',50))/100
            pdb  = str(r.get('pdb',''))[:4].lower()
            name = str(r.get('hit_name',pdb))[:12]
            cc   = str(r.get('cath_code','?')) or '?'
            col  = SCOP_COLORS.get(cls,'#888')
            op   = 0.3+0.7*prob
            sc += _seg(qs,qe,col,op,sccs if _hi==0 else '',
                       'HH: %s | SCOPe: %s | prob=%.0f%% | %d-%d aa'%(name,sccs,prob*100,qs,qe))
            ca += _seg(qs,qe,'#7F77DD',op*0.85,cc if _hi==0 else '',
                       'CATH: %s | %s | %d-%d aa'%(cc,name,qs,qe))
            pb += _seg(qs,qe,'#7F77DD',op,pdb.upper() if _hi==0 else '',
                       'PDB: %s | prob=%.0f%% | %d-%d aa'%(pdb.upper(),prob*100,qs,qe))
        return [('SCOPe','#7F77DD','direct sccs',sc,True),
                ('CATH', '#7F77DD','via PDB hit', ca,True),
                ('PDB',  '#7F77DD','top hits',    pb,True)]

    # ── build full HTML ───────────────────────────────────────────────────
    ticks = ''.join(
        '<span style="position:absolute;left:%.0f%%;font-size:9px;color:#bbb;'
        'transform:translateX(-50%%)">%d</span>' % (p, int(p*seq_len/100))
        for p in [0,25,50,75,100]
    )
    legend = ' &nbsp;'.join(
        '<span style="color:'+SCOP_COLORS.get(k,'#888')+'">&#9632;</span>'
        '<span style="font-size:9px;color:#666">'+k+'</span>'
        for k in SCOP_CLASSES
    )

    body = ''
    if domains:
        for dom in domains:
            qs = dom['qstart']; qe = dom['qend']; idx = dom['domain_idx']
            n  = len(domains)
            # Domain section header
            if n > 1:
                body += ('<div style="font-size:11px;font-weight:600;color:#333;'
                         'margin:14px 0 6px 0;padding:5px 0 5px 110px;'
                         'border-bottom:1.5px solid #c8d0e0">'
                         'Domain %d &nbsp;·&nbsp; aa %d–%d</div>' % (idx, qs, qe))
            # Arm bars
            for lbl,col,ann,segs,ok in _seq_bars(dom['df_seq_sub']):
                body += _bar('Seq · '+lbl, col, ann, segs, ok)
            body += '<div style="height:4px"></div>'
            for lbl,col,ann,segs,ok in _str_bars(dom['df_str_sub']):
                body += _bar('Str · '+lbl, col, ann, segs, ok)
            body += '<div style="height:4px"></div>'
            for lbl,col,ann,segs,ok in _hh_bars(dom['df_hh_sub']):
                body += _bar('HH · '+lbl, col, ann, segs, ok)
    else:
        # No domain clustering — show global bars
        for lbl,col,ann,segs,ok in _seq_bars(df_seq):
            body += _bar('Seq · '+lbl, col, ann, segs, ok)
        body += '<div style="height:4px"></div>'
        for lbl,col,ann,segs,ok in _str_bars(df_str):
            body += _bar('Str · '+lbl, col, ann, segs, ok)
        body += '<div style="height:4px"></div>'
        for lbl,col,ann,segs,ok in _hh_bars(df_hh):
            body += _bar('HH · '+lbl, col, ann, segs, ok)

    n_dom = len(domains) if domains else 1
    # 9 bars × 27px + 3 gaps × 4px + domain headers × 40px + ruler/tip/legend = overhead
    height = 80 + n_dom * (9 * 27 + 3 * 4 + 40)

    html = (
        '<!DOCTYPE html><html><head><meta charset="utf-8">'
        '<style>body{margin:0;padding:4px 0;font-family:-apple-system,sans-serif;'
        'background:transparent;font-size:12px}</style></head><body>'
        '<div id="tip" style="height:18px;font-size:11px;font-family:monospace;'
        'color:#555;white-space:nowrap;overflow:hidden;text-overflow:ellipsis;'
        'margin-bottom:6px;padding-left:110px"></div>'
        '<div style="position:relative;height:14px;margin-bottom:8px;padding-left:110px">'
        '<div style="position:absolute;top:7px;left:110px;right:0;height:0.5px;background:#d8dce4"></div>'
        + ticks + '</div>'
        + body
        + '<div style="font-size:9px;color:#bbb;margin-top:8px;padding-left:110px">'
        + legend + ' &nbsp;·&nbsp; opacity = confidence</div>'
        '<script>function T(s){document.getElementById("tip").textContent=s||"";}</script>'
        '</body></html>'
    )

    components.html(html, height=height, scrolling=True)


def render_compact_summary(df_seq, df_str, df_hh, fused, seq_len, domains=None):
    from db.lookup import get_cath_code, cath_code_to_scop_class
    from search.ensemble import get_domain_clusters
    lk = lookup.all_domains()

    render_three_domain_maps(df_seq, df_str, df_hh, seq_len, domains)

    # ── Per-domain ensemble cards ─────────────────────────────────────────────
    if domains is None:
        domains = get_domain_clusters(df_seq, df_str, df_hh, min_overlap=0.6)

    # Fallback: if clustering gives 0 domains (all arms empty), skip
    if not domains:
        st.caption("No domain hits found.")
        return

    n_domains = len(domains)
    multi = n_domains > 1

    SCOP_FULL = {
        "a":"All alpha","b":"All beta","c":"Alpha/beta","d":"Alpha+beta",
        "e":"Multi-domain","f":"Membrane","g":"Small proteins"
    }
    CATH_FULL = {"1":"Mainly alpha","2":"Mainly beta","3":"Alpha/beta","4":"Few secondary"}

    # Check for uncovered N/C terminal regions
    if domains:
        min_start = min(d["qstart"] for d in domains)
        max_end   = max(d["qend"]   for d in domains)
        if min_start > 30:
            st.caption(f"⚠️ aa 1–{min_start-1} has no hits — possible undetected domain or disordered region. Try raising e-value cutoff to 0.01.")
        if seq_len - max_end > 30:
            st.caption(f"⚠️ aa {max_end+1}–{seq_len} has no hits — possible undetected domain or disordered region.")
    if multi:
        st.markdown(
            f'<div style="font-size:13px;font-weight:500;margin:14px 0 8px">' +
            f'{n_domains} domain regions detected</div>',
            unsafe_allow_html=True
        )

    for dom in domains:
        idx       = dom["domain_idx"]
        qs        = dom["qstart"]
        qe        = dom["qend"]
        sub_seq   = dom["df_seq_sub"]
        sub_str   = dom["df_str_sub"]
        sub_hh    = dom["df_hh_sub"]
        dom_fused = dom["fused"]

        top         = dom_fused.iloc[0] if dom_fused is not None and not dom_fused.empty else None
        votes       = int(top.get("votes",0)) if top is not None else 0
        ev          = top.get("evidence","") if top is not None else ""
        agreed_sccs = ((top.get("agreed_sccs","—") or "—").lstrip("~")) if top is not None else "—"
        agreed_cath = top.get("agreed_cath","—") if top is not None else "—"
        # Override: prefer seq arm sccs/cath (more specific than HH arm)
        if sub_seq is not None and len(sub_seq) > 0 and 'target' in sub_seq.columns:
            _sq0 = sub_seq.iloc[0]
            _sq_sccs = lk.get(str(_sq0['target']), {}).get('sccs', '')
            if _sq_sccs and _sq_sccs not in ('--', '?', ''):
                agreed_sccs = _sq_sccs.lstrip('~')
            _sq_cc = get_cath_code(str(_sq0['target']))
            if _sq_cc and _sq_cc not in ('--', '?', ''):
                agreed_cath = _sq_cc
        score       = float(top.get("ensemble_score",0)) if top is not None else 0.0
        arms_cls    = top.get("arms_classes","") if top is not None else ""
        top_pdb     = ""
        if top is not None:
            top_pdb = (pdb_from_scope(str(top.get("scope_domain",""))) or
                       str(top.get("hh_hit",""))[:4].lower() or
                       str(top.get("cath_domain",""))[:4].lower())

        conf_color = {"all_three":"#1D9E75","two_arms":"#BA7517"}.get(ev,"#888")
        conf_label = {
            "all_three":"All 3 arms agree","two_arms":"2 arms agree",
            "seq_only":"Sequence arm only","struct_only":"Structure arm only",
            "hhblits_only":"Profile arm only (twilight zone)",
        }.get(ev, ev or "—")
        vote_bg  = {"all_three":"#D1FAE5","two_arms":"#FEF3C7"}.get(ev,"#F1EFE8")
        vote_col = {"all_three":"#065F46","two_arms":"#92400E"}.get(ev,"#5F5E5A")

        sccs_cls  = agreed_sccs[0] if agreed_sccs and agreed_sccs not in ("—","?") else "?"
        sccs_name = SCOP_FULL.get(sccs_cls,"")
        cath_c    = agreed_cath.split(".")[0] if agreed_cath and agreed_cath not in ("—","?") else "?"
        cath_name = CATH_FULL.get(cath_c,"")
        sccs_url  = f"https://scop.berkeley.edu/search/?key={agreed_sccs}" if agreed_sccs not in ("—","?") else "#"
        cath_url  = f"https://www.cathdb.info/version/v4_3_0/superfamily/{agreed_cath}" if agreed_cath not in ("—","?","") else "#"

        # Collect top PDB per arm (deduplicated)
        pdb_hits: list[tuple] = []
        if sub_seq is not None and len(sub_seq) > 0:
            p = pdb_from_scope(str(sub_seq.iloc[0].get("target","")))
            if p: pdb_hits.append((p, "seq"))
        if sub_str is not None and len(sub_str) > 0:
            p = str(sub_str.iloc[0].get("target",""))[:4].lower()
            if p and p not in [x[0] for x in pdb_hits]: pdb_hits.append((p, "struct"))
        if sub_hh is not None and len(sub_hh) > 0:
            p = str(sub_hh.iloc[0].get("pdb",""))[:4].lower()
            if p and p not in [x[0] for x in pdb_hits]: pdb_hits.append((p, "profile"))

        pdb_btns = " ".join(
            f'<a href="https://www.rcsb.org/structure/{p.upper()}" target="_blank" ' +
            f'style="font-size:11px;padding:3px 10px;border-radius:6px;border:0.5px solid #ccc;' +
            f'text-decoration:none;color:#185FA5;margin-right:4px">{p.upper()} ↗</a>'
            for p, _ in pdb_hits[:3]
        )
        af_btn = ""
        if top_pdb:
            af_btn = (f'<a href="https://alphafold.ebi.ac.uk/search/text/{top_pdb.upper()}" target="_blank" ' +
                      f'style="font-size:11px;padding:3px 10px;border-radius:6px;border:0.5px solid #ccc;' +
                      f'text-decoration:none;color:#0F6E56;margin-right:4px">AlphaFold ↗</a>')

        dom_label = f"Domain {idx} &nbsp;·&nbsp; aa {qs}–{qe}" if multi else f"aa {qs}–{qe}"

        st.markdown(
            f'<div style="border:0.5px solid #e0e0e0;border-left:4px solid {conf_color};' +
            f'border-radius:0 12px 12px 0;padding:14px 18px;margin:8px 0">' +
            f'<div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:10px">' +
            f'<span style="font-size:13px;font-weight:600">{dom_label}</span>' +
            f'<span style="font-size:11px;padding:2px 10px;border-radius:10px;font-weight:500;' +
            f'background:{vote_bg};color:{vote_col}">{votes}/3 &nbsp;·&nbsp; {conf_label}</span>' +
            f'</div>' +
            f'<div style="display:flex;gap:28px;flex-wrap:wrap;margin-bottom:10px">' +
            f'<div><div style="font-size:10px;color:#aaa;letter-spacing:0.4px;text-transform:uppercase">SCOPe</div>' +
            f'<a href="{sccs_url}" style="font-size:15px;font-weight:600;text-decoration:none;color:#185FA5">{agreed_sccs}</a>' +
            f'<div style="font-size:11px;color:#888">{sccs_name}</div></div>' +
            f'<div><div style="font-size:10px;color:#aaa;letter-spacing:0.4px;text-transform:uppercase">CATH</div>' +
            f'<a href="{cath_url}" style="font-size:15px;font-weight:600;text-decoration:none;color:#185FA5">{agreed_cath}</a>' +
            f'<div style="font-size:11px;color:#888">{cath_name}</div></div>' +
            f'<div><div style="font-size:10px;color:#aaa;letter-spacing:0.4px;text-transform:uppercase">Score</div>' +
            f'<span style="font-size:15px;font-weight:600">{score:.3f}</span></div>' +
            f'</div>' +
            f'<div style="margin-bottom:6px">{pdb_btns}{af_btn}</div>' +
            f'<div style="font-size:10px;color:#bbb">Arm classes: {arms_cls}</div>' +
            f'</div>',
            unsafe_allow_html=True
        )

        if multi and idx < n_domains:
            st.markdown('<hr style="border:none;border-top:0.5px solid #eee;margin:12px 0"/>', unsafe_allow_html=True)

    if not multi:
        st.caption("Each arm searches a different database — PDB IDs may differ but represent the same fold.")



# ── expert arm panel ──────────────────────────────────────────────────────────
def render_arm_panel(title, color, df, arm, top_n=7):
    lk = lookup.all_domains()
    st.markdown(f'<div style="font-weight:600;font-size:14px;color:{color};border-left:3px solid {color};padding-left:8px;margin-bottom:8px">{title}</div>', unsafe_allow_html=True)
    if df is None or len(df) == 0:
        st.info("No hits found.")
        return
    hits = df.head(top_n)

    st.markdown("**SCOPe classes**" if arm != "str" else "**CATH → SCOPe class** *(structure arm gives class letter only)*")
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
            st.markdown(f'<span style="color:{col}">●</span> CATH:`{cc}` → class **{cls}** ({SCOP_CLASSES.get(cls,cls)}) · <span style="color:{lc}">lDDT={lddt:.2f}</span>', unsafe_allow_html=True)
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
            _cached = _load_cache(clean_seq)
            if _cached:
                df_seq = _cached['df_seq']; df_str = _cached['df_str']; df_hh = _cached['df_hh']
                st.info('⚡ Cached result (< 48h). Run again to refresh.')

            # Progress indicator
            if not _cached:
                _est = '~20s'
                if use_hhblits: _est = '~5-10 min (HHblits enabled)'
                st.info(f'🔍 Running search — estimated time: {_est}')

            col_s, col_t, col_h = st.columns(3)
            with col_s:
                if _cached and df_seq is not None:
                    st.success(f"Sequence: {len(df_seq)} SCOPe hits (cached)")
                else:
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
                    pass  # HHblits disabled — no message shown

            if not _cached: _save_cache(clean_seq, df_seq, df_str, df_hh)
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
            from search.ensemble import get_domain_clusters as _gdc_main
            _domains_main = _gdc_main(df_seq, df_str, df_hh, min_overlap=0.6)
            render_compact_summary(df_seq, df_str, df_hh, fused, bar_len, _domains_main)

            # ── EXPERT MODE ───────────────────────────────────────────────────
            with st.expander("Expert mode — detailed three-arm results", expanded=False):
                st.caption("Each arm independently: SCOPe class · CATH code · top 7 PDB hits · functional lookup")
                from search.ensemble import get_domain_clusters as _gdc
                _domains_exp = _gdc(df_seq, df_str, df_hh, min_overlap=0.6)
                if len(_domains_exp) > 1:
                    _exp_tab_labels = [f"Domain {d['domain_idx']} · aa {d['qstart']}–{d['qend']}" for d in _domains_exp]
                    _exp_tabs = st.tabs(_exp_tab_labels)
                    for _dom, _tab in zip(_domains_exp, _exp_tabs):
                        with _tab:
                            _c1, _c2, _c3 = st.columns(3)
                            with _c1: render_arm_panel("Sequence arm — MMseqs2 / SCOPe", "#3B82F6", _dom["df_seq_sub"], "seq")
                            with _c2: render_arm_panel("Structure arm — ESMFold / CATH50", "#0F6E56", _dom["df_str_sub"], "str")
                            with _c3: render_arm_panel("Profile arm — HHblits / PDB70", "#7F77DD", _dom["df_hh_sub"], "hh")
                else:
                    arm_col1, arm_col2, arm_col3 = st.columns(3)
                    with arm_col1: render_arm_panel("Sequence arm — MMseqs2 / SCOPe", "#3B82F6", df_seq, "seq")
                    with arm_col2: render_arm_panel("Structure arm — ESMFold / CATH50", "#0F6E56", df_str, "str")
                    with arm_col3: render_arm_panel("Profile arm — HHblits / PDB70", "#7F77DD", df_hh, "hh")

            # ── ENSEMBLE TABLE ────────────────────────────────────────────────
            if not fused.empty:
                with st.expander("Ensemble ranked hits table", expanded=False):
                    from search.ensemble import get_domain_clusters as _gdc_ens
                    _domains_ens = _gdc_ens(df_seq, df_str, df_hh, min_overlap=0.6)
                    ev_icon = {"all_three":"🟢","two_arms":"🟡","seq_only":"🔵","struct_only":"🟠","hhblits_only":"🟣"}
                    def _prep_df(df):
                        df = df.copy()
                        df["ev"] = df["evidence"].map(ev_icon)
                        df["cath_d"] = df["cath_domain"].apply(clean_cath_domain)
                        for col in ["seq_evalue","struct_evalue","hh_evalue"]:
                            if col in df.columns:
                                df[col] = df[col].apply(lambda x: fmt_e(x) if x is not None and str(x) not in ("None","nan") else "—")
                        df["hh_prob"] = df["hh_prob"].apply(lambda x: f"{float(x):.0f}%" if x is not None and str(x) not in ("None","nan") else "—")
                        # Cap hhblits_only flood
                        df = pd.concat([
                            df[df["evidence"] != "hhblits_only"],
                            df[df["evidence"] == "hhblits_only"].head(10)
                        ]).sort_values("ensemble_score", ascending=False).reset_index(drop=True)
                        if "rank" in df.columns:
                            df["rank"] = range(1, len(df)+1)
                        else:
                            df.insert(0, "rank", range(1, len(df)+1))
                        return df
                    _cols = ["rank","ev","scope_domain","sccs","cath_d","cath_code","hh_hit","votes","ensemble_score","qstart","qend","seq_evalue","struct_evalue","hh_evalue","hh_prob","lddt"]
                    _rename = {"ev":"","scope_domain":"SCOPe domain","sccs":"SCOP class","cath_d":"CATH domain","cath_code":"CATH code","hh_hit":"HHblits hit","votes":"Votes","ensemble_score":"Score","qstart":"Start aa","qend":"End aa","seq_evalue":"Seq e-val","struct_evalue":"Struct e-val","hh_evalue":"HH e-val","hh_prob":"HH prob","lddt":"lDDT"}
                    if len(_domains_ens) > 1:
                        _ens_tab_labels = [f"Domain {d['domain_idx']} · aa {d['qstart']}–{d['qend']}" for d in _domains_ens]
                        _ens_tabs = st.tabs(_ens_tab_labels)
                        for _dom_e, _tab_e in zip(_domains_ens, _ens_tabs):
                            with _tab_e:
                                _df_e = _dom_e["fused"]
                                if _df_e is not None and not _df_e.empty:
                                    _df_e = _prep_df(_df_e)
                                    n_all=(_df_e["evidence"]=="all_three").sum(); n_two=(_df_e["evidence"]=="two_arms").sum()
                                    n_hh=(_df_e["evidence"]=="hhblits_only").sum()
                                    _c1,_c2,_c3 = st.columns(3)
                                    _c1.metric("🟢 All three",n_all); _c2.metric("🟡 Two arms",n_two); _c3.metric("🟣 Profile only (shown: max 10)",n_hh)
                                    st.dataframe(_df_e[[c for c in _cols if c in _df_e.columns]].rename(columns=_rename), use_container_width=True, hide_index=True)
                                else:
                                    st.info("No hits for this domain.")
                    else:
                        _df_g = _prep_df(fused)
                        n_all=(_df_g["evidence"]=="all_three").sum(); n_two=(_df_g["evidence"]=="two_arms").sum()
                        n_seq=(_df_g["evidence"]=="seq_only").sum(); n_st=(_df_g["evidence"]=="struct_only").sum()
                        n_hh=(_df_g["evidence"]=="hhblits_only").sum()
                        cc1,cc2,cc3,cc4,cc5 = st.columns(5)
                        cc1.metric("🟢 All three",n_all); cc2.metric("🟡 Two arms",n_two)
                        cc3.metric("🔵 Seq only",n_seq); cc4.metric("🟠 Struct only",n_st); cc5.metric("🟣 Profile only (max 10)",n_hh)
                        st.dataframe(_df_g[[c for c in _cols if c in _df_g.columns]].rename(columns=_rename), use_container_width=True, hide_index=True)
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