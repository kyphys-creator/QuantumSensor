#!/usr/bin/env python3
"""Build a clean 06_response_defs.nb from scratch.

Replaces the prior split-extract + annotate + add_response_kernels pipeline
with a single source of truth. Overwrites the 06_response_defs.nb that
split_TES.py / split_MKID.py produced.

Six sections per notebook:

  1. Setup                     — Direc, Get prev, SetDirectory
  2. Helpers and resolution    — energySum, midpointSum, BoxPDF
  3. Differential rate kernels — kerR{All,Alr|TiNl,TiNr}, RangeGen,
                                 IntkerR*, KerR*  (clean rewrite)
  4. Energy-bin response       — CRTES*/CRTiN* (continuous in vm)
  5. Velocity integration       — CRintTES*/CRintTiN*
  6. Legacy table API           — CRvmTES*/CRvmTiN* (Table wrappers for
                                  backward compatibility with 08_response_matrix)

Dropped from the original:
  • kerR{All,Alr}Cut — buggy (called kerRTiNl/r inside Al subsubsection)
  • CRvmTiN* in the TES file (TiN code leaked into Al subsubsection)

Run order:
  split_TES.py / split_MKID.py        (generate other notebooks)
  build_06_response_defs.py           (this script — overwrites 06)
"""

import re
import uuid
from pathlib import Path

BASE = Path(__file__).parent


# ============================================================================
# String encoding helpers
# ============================================================================

def mma_escape(s: str) -> str:
    """Encode non-ASCII characters to Mathematica's \\:XXXX notebook escape.

    Existing backslash escapes (\\n, \\t, etc.) are left alone. We need this
    because Mathematica's front end expects Unicode characters in .nb files
    as \\:XXXX hex escapes; raw UTF-8 displays as garbled text.
    """
    out = []
    i = 0
    while i < len(s):
        c = s[i]
        if c == "\\" and i + 1 < len(s):
            out.append(s[i:i + 2])
            i += 2
            continue
        cp = ord(c)
        if cp < 0x80:
            out.append(c)
        elif cp <= 0xFFFF:
            out.append(f"\\:{cp:04x}")
        else:
            cp -= 0x10000
            hi = 0xD800 + (cp >> 10)
            lo = 0xDC00 + (cp & 0x3FF)
            out.append(f"\\:{hi:04x}\\:{lo:04x}")
        i += 1
    return "".join(out)


def cell_string(s: str) -> str:
    r"""Escape a multi-line Mathematica code/text snippet for use as a Cell content
    string. Doubles backslashes and escapes quotes, then applies mma_escape.
    """
    s = s.replace("\\", "\\\\").replace('"', '\\"')
    return mma_escape(s)


def fresh_uuid() -> str:
    return str(uuid.uuid4())


# ============================================================================
# Cell builders
# ============================================================================

def section_header(title: str) -> str:
    """A Section-style header cell."""
    return f'Cell["{mma_escape(title)}", "Section",ExpressionUUID->"{fresh_uuid()}"],\n\n'


def text_cell(body: str) -> str:
    """A Text cell with multi-line body content (newlines preserved)."""
    return (
        f'Cell["{cell_string(body)}", "Text",\n'
        f' Background->RGBColor[0.96, 0.96, 0.92],\n'
        f' CellTags->"doc",\n'
        f' ExpressionUUID->"{fresh_uuid()}"],\n\n'
    )


def code_cell(code: str, bg_rgb: str = "0.96, 0.96, 0.96") -> str:
    """An Input cell with Mathematica code, marked as InitializationCell."""
    return (
        f'Cell["{cell_string(code)}", "Input",\n'
        f' InitializationCell->True,\n'
        f' Background->RGBColor[{bg_rgb}],\n'
        f' ExpressionUUID->"{fresh_uuid()}"],\n\n'
    )


def section(title: str, doc: str, code: str, code_bg: str = "0.96, 0.96, 0.96") -> str:
    """Wrap a section title + doc cell + code cell in a CellGroupData."""
    return (
        'Cell[CellGroupData[{\n\n'
        + section_header(title)
        + text_cell(doc)
        + code_cell(code, code_bg)
        + '}, Open  ]],\n\n'
    )


# ============================================================================
# Notebook header / footer templates
# ============================================================================

NB_HEADER = '''(* Content-type: application/vnd.wolfram.mathematica *)

(*** Wolfram Notebook File ***)
(* http://www.wolfram.com/nb *)

(* CreatedBy='Mathematica 14.3' *)

(*CacheID: 234*)
(* Internal cache information:
NotebookFileLineBreakTest
NotebookFileLineBreakTest
NotebookDataPosition[         0,          0]
NotebookDataLength[         0,          0]
NotebookOptionsPosition[         0,          0]
NotebookOutlinePosition[         0,          0]
CellTagsIndexPosition[         0,          0]
WindowFrame->Normal*)

(* Beginning of Notebook Content *)
Notebook[{
Cell[CellGroupData[{

Cell["{TITLE}", "Title",
 ExpressionUUID->"{UUID_TITLE}"],

'''

NB_FOOTER = '''}, Open  ]]
},
WindowSize->{1920, 954},
WindowMargins->{{0, Automatic}, {Automatic, 0}},
FrontEndVersion->"14.3 for Mac OS X ARM (64-bit) (July 8, 2025)",
StyleDefinitions->"Default.nb",
ExpressionUUID->"{UUID_NB}"
]
(* End of Notebook Content *)
'''


def assemble_notebook(title: str, sections: list[str]) -> str:
    head = (
        NB_HEADER
        .replace("{TITLE}", mma_escape(title))
        .replace("{UUID_TITLE}", fresh_uuid())
    )
    foot = NB_FOOTER.replace("{UUID_NB}", fresh_uuid())
    return head + "".join(sections) + foot


# ============================================================================
# Per-detector configuration — names of material-specific symbols and grid
# parameters for the legacy CRvmTES*/CRvmTiN* wrappers.
# ============================================================================

CONFIG = {
    "TES": {
        "material":  "Al",
        "rho":       "rhoAl",
        "imeps":     "ImepsAlf",
        "sig":       "TESsig",
        "kerR_l":    "kerRAll",
        "kerR_r":    "kerRAlr",
        "intkerR_l": "IntkerRAll",
        "intkerR_r": "IntkerRAlr",
        "KerR_l":    "KerRAll",
        "KerR_r":    "KerRAlr",
        "RangeGen":  "RangeAlGen",
        "CR":        "CRTES",
        "CRleft":    "CRTESLeft",
        "CRright":   "CRTESRight",
        "CRint":     "CRintTES",
        "CRvm":      "CRvmTES",
        # Legacy CRvm grid:  CRvmTES uses vm 5 → 108 step 0.51758 (200 pts)
        # CRvmTESLeft/Right use vm 5 → 200 step 196/800 (~796 pts)
        "vm_total":  ("5", "108", "0.51758"),
        "vm_branch": ("5", "200", "196/800"),
        "prev_nb":   "05_parameters.nb",
    },
    "MKID": {
        "material":  "TiN",
        "rho":       "rhoTiN",
        "imeps":     "ImepsLTiN",
        "sig":       "MKIDsig",
        "kerR_l":    "kerRTiNl",
        "kerR_r":    "kerRTiNr",
        "intkerR_l": "IntkerRTiNl",
        "intkerR_r": "IntkerRTiNr",
        "KerR_l":    "KerRTiNl",
        "KerR_r":    "KerRTiNr",
        "RangeGen":  "RangeTiNGen",
        "CR":        "CRTiN",
        "CRleft":    "CRTiNLeft",
        "CRright":   "CRTiNRight",
        "CRint":     "CRintTiN",
        "CRvm":      "CRvmTiN",
        # MKID grids: CRvmTiN uses vm 6 → 108 step 0.51256 (200 pts)
        # CRvmTiNLeft/Right: vm 6 → 200 step 195/800
        "vm_total":  ("6", "108", "0.51256"),
        "vm_branch": ("6", "200", "195/800"),
        "prev_nb":   "05_parameters.nb",
    },
}


# ============================================================================
# Section content (doc + code)
# ============================================================================

def setup_section(det: str) -> str:
    """Section 1: Setup — directory, Get upstream, SetDirectory."""
    cfg = CONFIG[det]
    doc = (
        "■ Setup\n\n"
        f"Direc を ../../input/{det} に指す。上流ノートブック ({cfg['prev_nb']}) を Get で読み込み、"
        "依存する定数・関数を取り込む。SetDirectory[Direc] で cwd をデータディレクトリに合わせる。"
    )
    code = (
        f'nbDir = NotebookDirectory[];\n'
        f'Direc = FileNameJoin[{{ParentDirectory[ParentDirectory[nbDir]], "input", "{det}"}}];\n'
        f'Get[FileNameJoin[{{nbDir, "{cfg["prev_nb"]}"}}]];\n'
        f'SetDirectory[Direc];\n'
        f'Print["Direc = ", Direc];'
    )
    return section("Setup", doc, code)


def helpers_section() -> str:
    """Section 2: Helpers and resolution — energySum, midpointSum, BoxPDF."""
    doc = (
        "■ Helpers and resolution\n\n"
        "以降の節で繰り返し使う汎用数値積分ヘルパーと、エネルギー分解能の矩形 PDF。\n\n"
        "  • energySum[f, E1, E2, dE]    観測エネルギー w を E1..E2 で dE 刻みリーマン和 (× eV)\n"
        "  • midpointSum[f, vmin, vmax, n]   速度 vm を [vmin, vmax] の n 中点でサンプルし中点則積分\n"
        "  • BoxPDF[Ep, sig][ER]    中心 Ep 半幅 sig の矩形 PDF (高速計算用、本流では Gaussian)"
    )
    code = (
        '(* ---- Generic numerical integrators ---- *)\n'
        '\n'
        '(* energySum[f, E1, E2, dE]\n'
        '   観測エネルギー w 方向のリーマン和。\n'
        '     ・f は単引数の純関数 (例: kerRAll[md][n][#, sig #][vm kps] &)\n'
        '     ・戻り値 = (Σ f(w_i)) × dE × eV  (エネルギー次元を付与)\n'
        '     ・dE はオプション (デフォルト 0.01)\n'
        '*)\n'
        'energySum[f_, E1_, E2_, dE_:0.01] := dE eV Total[f /@ Range[E1, E2, dE]];\n'
        '\n'
        '(* midpointSum[f, vmin, vmax, n]\n'
        '   速度 vm を [vmin, vmax] 区間 n 個の中点で評価する中点則。\n'
        '     ・f は単引数の純関数 (例: CRTES[md, n][E1, E2][#] &)\n'
        '     ・戻り値 = (Σ f(vm_i)) × dv,  vm_i = vmin + (i - 0.5) dv\n'
        '*)\n'
        'midpointSum[f_, vmin_, vmax_, n_Integer] := With[\n'
        '  {dv = (vmax - vmin)/n},\n'
        '  dv Total[f /@ Table[vmin + (i - 0.5) dv, {i, n}]]\n'
        '];\n'
        '\n'
        '(* ---- Resolution function ---- *)\n'
        '\n'
        '(* BoxPDF[Ep, sig][ER]\n'
        '   観測エネルギー Ep ± sig の uniform PDF。\n'
        '     ・|ER - Ep| ≤ sig で 1/(2 sig)、外で 0\n'
        '     ・本流の IntkerR* は Gaussian (NormalDistribution) を使うが、\n'
        '       単純化テスト用にこの矩形 PDF も用意。\n'
        '*)\n'
        'BoxPDF[Ep_, sig_][ER_] := UnitStep[sig - Abs[ER - Ep]] / (2 sig);'
    )
    return section("Helpers and resolution", doc, code, code_bg="0.94, 0.94, 0.94")


def kernels_section(det: str) -> str:
    """Section 3: Differential rate kernels — clean rewrite."""
    cfg = CONFIG[det]
    mat = cfg["material"]
    rho, imeps, sig = cfg["rho"], cfg["imeps"], cfg["sig"]
    kL, kR = cfg["kerR_l"], cfg["kerR_r"]
    iL, iR = cfg["intkerR_l"], cfg["intkerR_r"]
    KL, KR = cfg["KerR_l"], cfg["KerR_r"]
    rg = cfg["RangeGen"]
    doc = (
        "■ Differential rate kernels\n\n"
        f"DM-電子散乱の反跳エネルギー E に対する dR/dE 寄与 ({mat}, left/right 2 ブランチ)。\n"
        "原コード (kerR*Cut) のバグ (内部参照ミス) は除去。階層は以下:\n\n"
        f"  Level 0: {kL}, {kR}                 raw 微分散乱率カーネル (関数 of E, mχ, n, vmin)\n"
        f"  Level 1: {rg}                       数値積分の積分域 [emin, emax] 生成\n"
        f"  Level 2: {iL}, {iR}                 Gaussian 分解能で畳み込んだ観測量 dR/dE'\n"
        f"  Level 3: {KL}, {KR}                 自動 range 付き user-facing kernel\n\n"
        "全体共通の式:\n"
        f"    kerR{mat}α(E, mχ, n, vmin) = (1/ρ_{mat}) (1/8π² α) (1/μ_χe²)\n"
        "                                  × Jacobα(E, mχ, vmin)\n"
        "                                  × qα(E, mχ, vmin)³\n"
        "                                  × FDM(qα(E, mχ, vmin))[n]²\n"
        f"                                  × Im[ε_{mat}(E, qα(E, mχ, vmin))]\n"
        "    (α ∈ {l, r}; qα が実数値でなければ 0)"
    )
    code = (
        '(* ---- Level 0: Raw differential rate kernels ---- *)\n'
        '\n'
        '(* Left branch (qα = ql) *)\n'
        f'{kL}[E_][m\\[Chi]_][n_][vmin_] := With[\n'
        '  {q = ql[E][m\\[Chi]][vmin]},\n'
        '  If[Internal`RealValuedNumericQ[q],\n'
        f'    (1/{rho}) (1/(8 \\[Pi]^2 alpha)) (1/(\\[Mu]\\[Chi]t[m\\[Chi], me])^2) *\n'
        '      Jacobl[E][m\\[Chi]][vmin] q^3 *\n'
        '      FDM[q][n]^2 *\n'
        f'      {imeps}[E, q],\n'
        '    0\n'
        '  ]\n'
        '];\n'
        '\n'
        '(* Right branch (qα = qr) *)\n'
        f'{kR}[E_][m\\[Chi]_][n_][vmin_] := With[\n'
        '  {q = qr[E][m\\[Chi]][vmin]},\n'
        '  If[Internal`RealValuedNumericQ[q],\n'
        f'    (1/{rho}) (1/(8 \\[Pi]^2 alpha)) (1/(\\[Mu]\\[Chi]t[m\\[Chi], me])^2) *\n'
        '      Jacobr[E][m\\[Chi]][vmin] q^3 *\n'
        '      FDM[q][n]^2 *\n'
        f'      {imeps}[E, q],\n'
        '    0\n'
        '  ]\n'
        '];\n'
        '\n'
        '(* ---- Level 1: Integration domain ---- *)\n'
        '\n'
        f'(* {rg}[vmin][md][Ep][Sig]\n'
        '   観測 Ep ± Sig と kinematic 上限 (1/2 md vmin²) から、真の反跳エネルギー\n'
        '   の積分域 {emin, emax} を返す。emin ≥ emax (空) のときは {0, 0}。\n'
        '*)\n'
        f'{rg}[vmin_][md_][Ep_][Sig_] := Module[\n'
        '  {emin, emax},\n'
        '  emax = Min[(md vmin^2)/2, Ep + Sig];\n'
        '  emin = Max[0.1 eV, Ep - Sig];\n'
        '  If[emin < emax, {emin, emax}, {0, 0}]\n'
        '];\n'
        '\n'
        '(* ---- Level 2: Resolution-convolved kernels ---- *)\n'
        '\n'
        f'(* {iL}, {iR}\n'
        '   kerR* を Gaussian N(Ep, sig) で ER∈[E1, E2] にわたって畳み込んだ観測量。\n'
        '   energySum を流用し、Riemann 和 (ΔE = 0.001 eV) で離散化。\n'
        '*)\n'
        f'{iL}[E1_, E2_][m\\[Chi]_][n_][Ep_, sig_][vmin_] :=\n'
        '  energySum[\n'
        f'    PDF[NormalDistribution[Ep, sig], #] {kL}[#][m\\[Chi]][n][vmin] &,\n'
        '    E1, E2, 0.001\n'
        '  ];\n'
        '\n'
        f'{iR}[E1_, E2_][m\\[Chi]_][n_][Ep_, sig_][vmin_] :=\n'
        '  energySum[\n'
        f'    PDF[NormalDistribution[Ep, sig], #] {kR}[#][m\\[Chi]][n][vmin] &,\n'
        '    E1, E2, 0.001\n'
        '  ];\n'
        '\n'
        '(* ---- Level 3: User-facing kernels (auto range) ---- *)\n'
        '\n'
        f'(* {KL}, {KR}\n'
        f'   {rg} で積分域を自動決定して IntkerR* を呼ぶラッパー。\n'
        '   range が {0, 0} (空) なら 0 を返す。\n'
        '*)\n'
        f'{KL}[md_][n_][Ep_, sig_][vmin_] := Module[\n'
        '  {r1, r2},\n'
        f'  {{r1, r2}} = {rg}[vmin][md][Ep][sig];\n'
        f'  If[r1 == 0 && r2 == 0, 0, {iL}[r1, r2][md][n][Ep, sig][vmin]]\n'
        '];\n'
        '\n'
        f'{KR}[md_][n_][Ep_, sig_][vmin_] := Module[\n'
        '  {r1, r2},\n'
        f'  {{r1, r2}} = {rg}[vmin][md][Ep][sig];\n'
        f'  If[r1 == 0 && r2 == 0, 0, {iR}[r1, r2][md][n][Ep, sig][vmin]]\n'
        '];'
    )
    return section("Differential rate kernels", doc, code, code_bg="0.98, 0.94, 0.94")


def cr_section(det: str) -> str:
    """Section 4: Energy-bin response (continuous in vm)."""
    cfg = CONFIG[det]
    KL, KR, sig = cfg["KerR_l"], cfg["KerR_r"], cfg["sig"]
    CR, CRl, CRr = cfg["CR"], cfg["CRleft"], cfg["CRright"]
    doc = (
        "■ Energy-bin response (continuous in vm)\n\n"
        f"任意の数値 vm に対し、観測エネルギー w を E1..E2 で w 方向にリーマン和を取って返す。\n"
        f"vm は ?NumericQ パターンで連続変数として扱える (Plot, NIntegrate, D 可)。\n\n"
        f"  {CRl}[md, n][E1, E2][vm]  = ∫_E1^E2 {KL}[md][n][w, {sig}·w][vm·kps] dw\n"
        f"  {CRr}[md, n][E1, E2][vm]  = ∫_E1^E2 {KR}[md][n][w, {sig}·w][vm·kps] dw\n"
        f"  {CR}[md, n][E1, E2][vm]   = Left + Right  (energySum を 1 回で済ます最適化)"
    )
    code = (
        f'{CRl}[md_, n_][E1_, E2_][vm_?NumericQ] :=\n'
        f'  energySum[{KL}[md][n][#, {sig} #][vm kps] &, E1, E2];\n'
        '\n'
        f'{CRr}[md_, n_][E1_, E2_][vm_?NumericQ] :=\n'
        f'  energySum[{KR}[md][n][#, {sig} #][vm kps] &, E1, E2];\n'
        '\n'
        '(* 両ブランチ和は Left + Right を別個呼ぶより、energySum を 1 回で済ます方が速い。 *)\n'
        f'{CR}[md_, n_][E1_, E2_][vm_?NumericQ] :=\n'
        f'  energySum[({KL}[md][n][#, {sig} #][vm kps] +\n'
        f'             {KR}[md][n][#, {sig} #][vm kps]) &, E1, E2];'
    )
    return section("Energy-bin response (continuous in vm)", doc, code, code_bg="0.94, 0.96, 0.98")


def crint_section(det: str) -> str:
    """Section 5: Velocity-integrated response."""
    cfg = CONFIG[det]
    CR, CRl, CRr = cfg["CR"], cfg["CRleft"], cfg["CRright"]
    CRi = cfg["CRint"]
    doc = (
        "■ Velocity-integrated response\n\n"
        f"{CR}* を vm 方向に [vmin, vmax] で ns 点の中点則で積分。midpointSum を 1 行呼ぶだけ。\n\n"
        f"  {CRi}[md, n][E1, E2][vmin, vmax, ns]            = ∫_vmin^vmax {CR}[…][vm] dvm\n"
        f"  {CRi}Weighted[md, n][E1, E2][vmin, vmax, ns, eta] = ∫_vmin^vmax {CR}[…][vm] · eta[vm] dvm\n\n"
        "Weighted 版は halo 速度分布 η(vm) との畳み込み (= イベント率の vm 積分項)。"
    )
    code = (
        f'{CRi}Left[md_, n_][E1_, E2_][vmin_?NumericQ, vmax_?NumericQ, ns_Integer] :=\n'
        f'  midpointSum[{CRl}[md, n][E1, E2][#] &, vmin, vmax, ns];\n'
        '\n'
        f'{CRi}Right[md_, n_][E1_, E2_][vmin_?NumericQ, vmax_?NumericQ, ns_Integer] :=\n'
        f'  midpointSum[{CRr}[md, n][E1, E2][#] &, vmin, vmax, ns];\n'
        '\n'
        f'{CRi}[md_, n_][E1_, E2_][vmin_?NumericQ, vmax_?NumericQ, ns_Integer] :=\n'
        f'  midpointSum[{CR}[md, n][E1, E2][#] &, vmin, vmax, ns];\n'
        '\n'
        '(* halo 速度分布 η(vm) を掛けた畳み込み版。 *)\n'
        f'{CRi}Weighted[md_, n_][E1_, E2_][vmin_?NumericQ, vmax_?NumericQ, ns_Integer, eta_] :=\n'
        f'  midpointSum[{CR}[md, n][E1, E2][#] eta[#] &, vmin, vmax, ns];'
    )
    return section("Velocity-integrated response", doc, code, code_bg="0.94, 0.98, 0.94")


def crvm_section(det: str) -> str:
    """Section 6: Legacy table API."""
    cfg = CONFIG[det]
    CR, CRl, CRr = cfg["CR"], cfg["CRleft"], cfg["CRright"]
    CRv = cfg["CRvm"]
    vmin_t, vmax_t, step_t = cfg["vm_total"]
    vmin_b, vmax_b, step_b = cfg["vm_branch"]
    doc = (
        "■ Legacy table API\n\n"
        "原コードの While ループ実装を Table で書き直した薄いラッパー。\n"
        f"08_response_matrix.nb との後方互換のため残置。{CR}* の連続関数版を vm グリッドでサンプルする。\n\n"
        "原実装との違い:\n"
        "  • 戻り値の Length が実際のサンプル点数と一致 (原コードは 800 確保し末尾を 0 埋めしていた)\n"
        f"  • {CRv}[md, n][E1, E2]      vm = {vmin_t} → {vmax_t} km/s, ステップ {step_t} (両ブランチ和)\n"
        f"  • {CRv}Left/Right[md, n][E1, E2]  vm = {vmin_b} → {vmax_b} km/s, ステップ {step_b}"
    )
    code = (
        f'{CRv}[md_, n_][E1_, E2_] :=\n'
        f'  Table[{CR}[md, n][E1, E2][vm], {{vm, {vmin_t}, {vmax_t}, {step_t}}}];\n'
        '\n'
        f'{CRv}Left[md_, n_][E1_, E2_] :=\n'
        f'  Table[{CRl}[md, n][E1, E2][vm], {{vm, {vmin_b}, {vmax_b}, {step_b}}}];\n'
        '\n'
        f'{CRv}Right[md_, n_][E1_, E2_] :=\n'
        f'  Table[{CRr}[md, n][E1, E2][vm], {{vm, {vmin_b}, {vmax_b}, {step_b}}}];'
    )
    return section("Legacy table API", doc, code, code_bg="0.97, 0.95, 0.97")


# ============================================================================
# Build and write
# ============================================================================

def build_for(det: str) -> str:
    """Assemble the full notebook string for one detector."""
    title = f"{det} - Response Function Definitions"
    sections = [
        setup_section(det),
        helpers_section(),
        kernels_section(det),
        cr_section(det),
        crint_section(det),
        crvm_section(det),
    ]
    return assemble_notebook(title, sections)


def main():
    for det in ("TES", "MKID"):
        out = BASE / det / "06_response_defs.nb"
        if not out.parent.exists():
            print(f"  {out.parent}: not found, skipping")
            continue
        nb = build_for(det)
        out.write_text(nb)
        # Basic sanity check
        op = nb.count("Cell[CellGroupData[{")
        cl = len(re.findall(r"^\}, Open\s*\]\],?\s*$", nb, re.MULTILINE))
        print(f"  {det}: wrote {out.name}  {len(nb)//1024}KB  CellGroupOpen={op} Close={cl}")


if __name__ == "__main__":
    main()
