#!/usr/bin/env python3
"""Inject explanatory Text cells into 06_response_defs.nb (TES & MKID variants).

Run AFTER split_TES.py / split_MKID.py. Idempotent: if marker already present,
this script does nothing.
"""

import re
import uuid
from pathlib import Path

BASE = Path(__file__).parent

# Sentinel string we embed in our first injected cell, so reruns can detect
# that annotations have already been added.
SENTINEL = "QSENS-DOC-MARKER-v1"


def _mma_escape(s):
    """Encode non-ASCII characters to Mathematica's \\:XXXX notebook escape.

    Preserves existing escape sequences (\\n, \\t, etc.) and ASCII chars.
    Mathematica's front end stores Unicode in notebook files as \\:XXXX,
    not raw UTF-8; raw UTF-8 displays as garbled text (mojibake).
    """
    out = []
    i = 0
    while i < len(s):
        c = s[i]
        # Preserve backslash escapes
        if c == "\\" and i + 1 < len(s):
            out.append(s[i:i+2])
            i += 2
            continue
        cp = ord(c)
        if cp < 0x80:
            out.append(c)
        elif cp <= 0xFFFF:
            out.append(f"\\:{cp:04x}")
        else:
            # Surrogate pair for chars beyond BMP
            cp -= 0x10000
            hi = 0xD800 + (cp >> 10)
            lo = 0xDC00 + (cp & 0x3FF)
            out.append(f"\\:{hi:04x}\\:{lo:04x}")
        i += 1
    return "".join(out)


def _text_cell(content, marker=False):
    """Build a Cell[TextData[...], "Text", ...] string from a list of segments.

    Each segment is either a plain string (rendered as normal text) or a
    2-tuple (style, text), e.g. ("Code", "BoxPDF[Ep, sig]").
    """
    parts = []
    for seg in content:
        if isinstance(seg, tuple):
            style, txt = seg
            parts.append(f'StyleBox["{_mma_escape(txt)}", "{style}"]')
        else:
            parts.append(f'"{_mma_escape(seg)}"')
    body = ", ".join(parts)
    extra = f' (* {SENTINEL} *)' if marker else ""
    return (
        f'Cell[TextData[{{{body}}}], "Text",\n'
        f' CellTags->"doc",\n'
        f' Background->RGBColor[0.95, 0.97, 1.0],\n'
        f' ExpressionUUID->"{uuid.uuid4()}"],{extra}\n\n'
    )


# ============================================================================
# Documentation content. Keys are the function name (or a unique tag for a
# combined Code cell that defines multiple functions). Values are lists of
# TextData segments. Use ("Code", "...") for monospace, ("StandardForm", ...)
# for emphasis, etc.
# ============================================================================

DOCS = {
    # ---- Subsection: defined functions ----
    "BoxPDF": [
        ("Text", "■ BoxPDF[Ep, sig][ER]   — 矩形 (box) 型エネルギー分解能関数\\n\\n"),
        "中心 ", ("Code", "Ep"), " 半幅 ", ("Code", "sig"),
        " の uniform PDF。観測エネルギー Ep に対し、真の反跳エネルギー ER が ",
        "|ER−Ep| ≤ sig の範囲なら 1/(2·sig)、外側は 0。\\n\\n",
        "Gaussian の代わりに高速計算したい時のための単純化版。",
        "実際の畳み込みでは IntkerRAl* 側で NormalDistribution[Ep, sig] を使っているので、",
        "この BoxPDF は IntkerRAl* のコメントアウトされた代替版として温存されている。",
    ],

    # ---- Subsubsection: Al ----
    "kerRAl_LR": [
        ("Text", "■ kerRAll[E][mχ][n][vmin], kerRAlr[E][mχ][n][vmin]   — 微分散乱率カーネル (Left/Right ブランチ)\\n\\n"),
        "DM-電子散乱の反跳エネルギー E あたりの ",
        ("Code", "dR/dE"), " 寄与。kinematic に運動量 q が 2 解 ",
        "(left: ", ("Code", "ql"), ", right: ", ("Code", "qr"),
        ") を持つので、左右別関数化。\\n\\n",
        "式:\\n",
        "    kerRAlα = (1/ρ_Al) × (1/8π²α) × (1/μ_χe²) × Jacobα × qα³ × FDM(qα)[n]² × Im[ε_Al(E, qα)]\\n",
        "    (α ∈ {l, r})\\n\\n",
        "引数:\\n",
        "  • ", ("Code", "E"), "    反跳エネルギー (eV 単位の数値)\\n",
        "  • ", ("Code", "mχ"), "   DM 質量 (MeV/GeV)\\n",
        "  • ", ("Code", "n"), "    媒介子モデル番号 (1=heavy, 2=interm, 3=light) — FDM(q)[n] の q 依存性を決める\\n",
        "  • ", ("Code", "vmin"), " DM 最小速度 (km/s スケール内部単位)\\n\\n",
        "依存:\\n",
        "  • ", ("Code", "ql, qr, Jacobl, Jacobr, FDM"), " (03_functions_response)\\n",
        "  • ", ("Code", "ImepsAlf"), " (04_material — Al の Im ε 補間)\\n",
        "  • ", ("Code", "rhoAl, alpha, me, μχt"), " (01_setup の定数)\\n\\n",
        "数値判定: ql/qr が NumericQ で実数なら式を返す、それ以外は 0。",
    ],

    "kerRAlCut": [
        ("Text", "■ kerRAllCut[ER][mχ][n][vmin], kerRAlrCut[ER][mχ][n][vmin]   — 閾値カット付きカーネル\\n\\n"),
        "kinematic 上限 ", ("Code", "vmin ≥ √(2·ER/mχ)·kps"),
        " より速い DM は ER を出せないため 0、それ以外は内側カーネルを呼ぶ。\\n\\n",
        "⚠ 注意: 元コードのバグで kerRTiNl/kerRTiNr (TiN 関数) を呼び出している。",
        "TES として正しく動かすには kerRAll/kerRAlr に修正が必要。",
        "詳細は README 参照。",
    ],

    "RangeAlGen": [
        ("Text", "■ RangeAlGen[vmin][md][Ep][Sig]   — 数値積分の積分域生成\\n\\n"),
        "観測エネルギー Ep ± Sig と kinematic 上限から、真の反跳エネルギーの積分域 ",
        ("Code", "{emin, emax}"), " を返す。\\n\\n",
        "式:\\n",
        "    emin = max(0.1 eV, Ep − Sig)\\n",
        "    emax = min(½·md·vmin², Ep + Sig)\\n\\n",
        "emin ≥ emax (積分域が空) の場合は ", ("Code", "{0, 0}"), " を返す。",
        "上位の KerRAll/KerRAlr がこれをチェックして 0 を返す目印として使う。",
    ],

    "IntkerRAl_LR": [
        ("Text", "■ IntkerRAll[E1, E2][mχ][n][Ep, sig][vmin], IntkerRAlr[...]   — 分解能畳み込み カーネル\\n\\n"),
        "kerRAlα を Gaussian 分解能関数 ", ("Code", "N(Ep, sig)"),
        " で ER ∈ [E1, E2] にわたって畳み込んだ「観測される微分率」。\\n\\n",
        "式:\\n",
        "    IntkerRAlα(E1, E2)[...][Ep, sig][vmin]\\n",
        "      ≈ Σ_{ER = E1..E2, Δ=0.001 eV} PDF[N(Ep,sig), ER] × kerRAlα[ER][...] × Δ\\n\\n",
        "Riemann 和で離散化 (Δ = 0.001 eV)。NIntegrate より高速だが粗い。\\n",
        "コードには BoxPDF を使う代替版がコメントアウトで残されている。",
    ],

    "KerRAl_LR": [
        ("Text", "■ KerRAll[md][n][Ep, sig][vmin], KerRAlr[md][n][Ep, sig][vmin]   — ユーザ向け観測カーネル\\n\\n"),
        "RangeAlGen で自動的に積分域 {r1, r2} を決め、IntkerRAlα を呼ぶラッパー。\\n",
        "積分域が空 ({0,0}) なら 0 を返す。\\n\\n",
        "外部呼び出し時は基本これを使う。CRvmTES* の中で各 vm 値について繰り返し呼ばれる。",
    ],

    # ---- Subsubsection: defined (Matrix data から移動) ----
    "CRvmTES_group": [
        ("Text", "■ CRvmTES[md, n][E1, E2], CRvmTESLeft[md, n][E1, E2], CRvmTESRight[md, n][E1, E2]   — vm をパラメータとする観測エネルギー積分のルックアップテーブル\\n\\n"),
        "重要: この関数は ", ("Code", "w"), " (= 観測エネルギー = 反跳エネルギー ER) 方向の積分しか行わない。",
        ("Code", "vm"), " (DM 最小速度) は While ループでサンプリングするだけで、",
        "vm 積分は行わない。vm 積分は後段 (10_data.nb) で halo speed distribution ",
        ("Code", "η(vm)"), " を掛けて行う。\\n\\n",
        "数学的にやっていること:\\n",
        "    Rlist[i] = ∫_{E1}^{E2} (KerRAll[md][n][w, TESsig·w][vm[i]·kps] + KerRAlr[...]) dw\\n",
        "    ≈ Σ_{w=E1, E1+0.01, ..., E2}  KerR(w, vm[i]) × 0.01 eV   (Riemann 和、Δw = 0.01 eV)\\n\\n",
        "vm グリッド: i = 1..200 (or 800)、各 vm[i] が決まった値の Rlist[i] を返す。\\n\\n",
        "3 つの変種:\\n",
        "  • ", ("Code", "CRvmTES"),
        "       両ブランチ和 (KerRAll + KerRAlr)。vm: 5 → 108 km/s、刻み 0.51758、200 点\\n",
        "  • ", ("Code", "CRvmTESLeft"),
        "   左ブランチのみ (KerRAll)。vm: 5 → 200 km/s、刻み 196/800 ≈ 0.245、800 点\\n",
        "  • ", ("Code", "CRvmTESRight"),
        "  右ブランチのみ (KerRAlr)。同じ vm グリッド (5 → 200 km/s, 800 点)\\n\\n",
        "引数:\\n",
        "  • ", ("Code", "md"), "  DM 質量\\n",
        "  • ", ("Code", "n"), "   媒介子モデル (1/2/3)\\n",
        "  • ", ("Code", "E1, E2"), "  観測エネルギー bin の下端・上端 (eV)\\n\\n",
        "戻り値: 長さ 200 (or 800) のリスト、各要素は vm[i] における w 積分値。\\n",
        "TESsig は分解能パラメータ (≈ 0.036)、KerR の Gauss 幅は TESsig·w (= 観測エネルギーに比例)。",
    ],

    # ---- (MKID side) ----
    "kerRTiN_LR": [
        ("Text", "■ kerRTiNl[E][mχ][n][vmin], kerRTiNr[E][mχ][n][vmin]   — TiN 用微分散乱率カーネル (Left/Right)\\n\\n"),
        "kerRAl* の TiN 版。誘電関数は Lindhard 解析モデル (", ("Code", "ImepsLTiN"), ") を使う。\\n",
        "式は kerRAl* と同形、ρ_Al → ρ_TiN、ImepsAlf → ImepsLTiN に置換。",
    ],

    "kerRTiNCut": [
        ("Text", "■ kerRTiNlCut, kerRTiNrCut   — TiN 用の閾値カット付きカーネル\\n\\n"),
        "kerRAl*Cut の TiN 版。",
    ],

    "RangeTiNGen": [
        ("Text", "■ RangeTiNGen[vmin][md][Ep][Sig]   — TiN 用積分域生成\\n\\n"),
        "RangeAlGen と同構造の TiN 版。",
    ],

    "IntkerRTiN_LR": [
        ("Text", "■ IntkerRTiNl[...], IntkerRTiNr[...]   — TiN 用分解能畳み込み\\n\\n"),
        "IntkerRAl* の TiN 版。",
    ],

    "KerRTiN_LR": [
        ("Text", "■ KerRTiNl[md][n][Ep, sig][vmin], KerRTiNr[...]   — TiN 用ユーザ向け観測カーネル\\n\\n"),
        "KerRAl* の TiN 版。",
    ],

    "CRvmTiN_group": [
        ("Text", "■ CRvmTiN[md, n][E1, E2], CRvmTiNLeft[md, n][E1, E2], CRvmTiNRight[md, n][E1, E2]   — TiN/MKID 応答ルックアップテーブル群\\n\\n"),
        "CRvmTES* の TiN 版。w (観測エネルギー) 方向の積分のみ行い、vm はサンプリングするだけ。",
        "vm 積分は後段で η(vm) と畳み込んで行う。MKIDsig (≈ 0.127) を分解能として使う。\\n\\n",
        "  • ", ("Code", "CRvmTiN"),
        "       両ブランチ和。vm: 6 → 108 km/s、刻み 0.51256、200 点\\n",
        "  • ", ("Code", "CRvmTiNLeft"),
        "   左ブランチのみ。vm: 6 → 200 km/s、刻み 195/800、800 点\\n",
        "  • ", ("Code", "CRvmTiNRight"),
        "  右ブランチのみ。同じ vm グリッド",
    ],
}


# ============================================================================
# Where to insert each doc cell. We anchor on RowBox patterns that appear at
# the START of the function's definition (LHS of :=). We insert the doc cell
# immediately BEFORE the enclosing Cell[BoxData[ that contains this anchor.
# ============================================================================

# (doc_key, anchor_regex). The regex matches a single line in the notebook
# (after splitlines), and we walk back to find the nearest Cell[BoxData[
# at column 0 to use as the insertion point.
ANCHORS = [
    ("BoxPDF",        re.compile(r'RowBox\[\{"BoxPDF", "\["')),
    ("kerRAl_LR",     re.compile(r'RowBox\[\{"kerRAll", "\["')),
    ("kerRAlCut",     re.compile(r'RowBox\[\{"kerRAllCut", "\["')),
    ("RangeAlGen",    re.compile(r'RowBox\[\{"RangeAlGen", "\["')),
    ("IntkerRAl_LR",  re.compile(r'RowBox\[\{"IntkerRAll", "\["')),
    ("KerRAl_LR",     re.compile(r'RowBox\[\{"KerRAll", "\["')),
    ("CRvmTES_group", re.compile(r'RowBox\[\{"CRvmTES", "\["')),
    # MKID side
    ("kerRTiN_LR",    re.compile(r'RowBox\[\{"kerRTiNl", "\["')),
    ("kerRTiNCut",    re.compile(r'RowBox\[\{"kerRTiNlCut", "\["')),
    ("RangeTiNGen",   re.compile(r'RowBox\[\{"RangeTiNGen", "\["')),
    ("IntkerRTiN_LR", re.compile(r'RowBox\[\{"IntkerRTiNl", "\["')),
    ("KerRTiN_LR",    re.compile(r'RowBox\[\{"KerRTiNl", "\["')),
    ("CRvmTiN_group", re.compile(r'RowBox\[\{"CRvmTiN", "\["')),
]

CELL_BOXDATA_RE = re.compile(r'^Cell\[BoxData\[\{?$')
# How many lines past the start of a Cell[BoxData[ block the LHS of := must appear.
# Filters out matches inside function bodies.
LHS_WINDOW = 8


def annotate(file_path: Path):
    text = file_path.read_text()
    if SENTINEL in text:
        print(f"  {file_path.name}: already annotated, skipping")
        return
    lines = text.splitlines(keepends=True)

    insertions = []  # list of (line_idx, doc_text)
    used_cells = set()  # Cell[BoxData[ line indices already targeted
    first = True
    for doc_key, anchor_re in ANCHORS:
        # Search lines for first hit of this anchor that's near a Cell[BoxData[
        for i, line in enumerate(lines):
            if not anchor_re.search(line):
                continue
            # Walk back to enclosing Cell[BoxData[
            cell_idx = None
            for j in range(i, -1, -1):
                if CELL_BOXDATA_RE.match(lines[j]):
                    cell_idx = j
                    break
            if cell_idx is None:
                continue
            # Skip if match is too deep (likely inside a function body, not LHS)
            if i - cell_idx > LHS_WINDOW:
                continue
            # Skip if this cell already has a doc
            if cell_idx in used_cells:
                continue
            insertions.append((cell_idx, _text_cell(DOCS[doc_key], marker=first)))
            used_cells.add(cell_idx)
            first = False
            break

    # Apply insertions: sort by idx descending; for ties, keep collection order
    # (we want first-collected to end up TOPMOST, so insert it LAST at that idx).
    insertions.sort(key=lambda x: x[0], reverse=True)
    for idx, doc in insertions:
        lines.insert(idx, doc)

    file_path.write_text("".join(lines))
    print(f"  {file_path.name}: inserted {len(insertions)} doc cells")


def main():
    for detector in ("TES", "MKID"):
        f = BASE / detector / "06_response_defs.nb"
        if not f.exists():
            print(f"  {f}: not found, skipping")
            continue
        print(f"== {detector} ==")
        annotate(f)


if __name__ == "__main__":
    main()
