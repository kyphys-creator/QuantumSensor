#!/usr/bin/env python3
"""Inject continuous-vm response kernels (CRTES*, CRTiN*) into 06_response_defs.nb.

These are continuous-in-vm versions of CRvmTES*: they integrate over the
observed energy w from E1 to E2 (same Riemann sum, Δw = 0.01 eV) but do
not sample over vm — vm is a NumericQ argument. The result is a function
of a continuous vm rather than a precomputed lookup table.

Run after split_*.py and annotate_06_response_defs.py.
Idempotent: detects a sentinel comment and skips on re-run.
"""

import re
import uuid
from pathlib import Path

BASE = Path(__file__).parent
SENTINEL = "QSENS-CONT-KERNELS-v1"


TES_CODE = """\
(* Continuous-vm versions of CRvmTES* *)
(* vm をサンプリングせず、numeric な vm 値での w 積分を返す。 *)
(* CRTES[md, n][E1, E2][vm] is a continuous function of vm (km/s value). *)
CRTES[md_, n_][E1_, E2_][vm_?NumericQ] := Total[Table[
  0.01 eV * (
    KerRAll[md][n][w, TESsig*w][vm*kps] +
    KerRAlr[md][n][w, TESsig*w][vm*kps]
  ),
  {w, E1, E2, 0.01}
]];

CRTESLeft[md_, n_][E1_, E2_][vm_?NumericQ] := Total[Table[
  0.01 eV * KerRAll[md][n][w, TESsig*w][vm*kps],
  {w, E1, E2, 0.01}
]];

CRTESRight[md_, n_][E1_, E2_][vm_?NumericQ] := Total[Table[
  0.01 eV * KerRAlr[md][n][w, TESsig*w][vm*kps],
  {w, E1, E2, 0.01}
]];"""


MKID_CODE = """\
(* Continuous-vm versions of CRvmTiN* *)
(* vm をサンプリングせず、numeric な vm 値での w 積分を返す。 *)
CRTiN[md_, n_][E1_, E2_][vm_?NumericQ] := Total[Table[
  0.01 eV * (
    KerRTiNl[md][n][w, MKIDsig*w][vm*kps] +
    KerRTiNr[md][n][w, MKIDsig*w][vm*kps]
  ),
  {w, E1, E2, 0.01}
]];

CRTiNLeft[md_, n_][E1_, E2_][vm_?NumericQ] := Total[Table[
  0.01 eV * KerRTiNl[md][n][w, MKIDsig*w][vm*kps],
  {w, E1, E2, 0.01}
]];

CRTiNRight[md_, n_][E1_, E2_][vm_?NumericQ] := Total[Table[
  0.01 eV * KerRTiNr[md][n][w, MKIDsig*w][vm*kps],
  {w, E1, E2, 0.01}
]];"""


def _mma_escape(s):
    """Encode non-ASCII to \\:XXXX notebook escapes; preserve \\-escapes."""
    out = []
    i = 0
    while i < len(s):
        c = s[i]
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
            cp -= 0x10000
            hi = 0xD800 + (cp >> 10)
            lo = 0xDC00 + (cp & 0x3FF)
            out.append(f"\\:{hi:04x}\\:{lo:04x}")
        i += 1
    return "".join(out)


def _escape_for_cell_string(s):
    """Escape a multi-line code string for Cell["...", "Input"] content."""
    s = s.replace("\\", "\\\\")
    s = s.replace('"', '\\"')
    s = _mma_escape(s)
    return s


def _make_code_cell(code, marker=True):
    content = _escape_for_cell_string(code)
    marker_text = f' (* {SENTINEL} *)' if marker else ''
    return (
        f'Cell["{content}", "Input",\n'
        f' InitializationCell->True,\n'
        f' Background->RGBColor[0.98, 0.95, 0.95],\n'
        f' ExpressionUUID->"{uuid.uuid4()}"],{marker_text}\n\n'
    )


def _make_doc_cell(detector):
    label = "CRTES" if detector == "TES" else "CRTiN"
    cr_tbl = "CRvmTES" if detector == "TES" else "CRvmTiN"
    kerL = "KerRAll" if detector == "TES" else "KerRTiNl"
    kerR = "KerRAlr" if detector == "TES" else "KerRTiNr"
    sig = "TESsig" if detector == "TES" else "MKIDsig"
    body = (
        f'\\:25a0 {label}, {label}Left, {label}Right   '
        f'\\:2014 vm \\:306e\\:9023\\:7d9a\\:95a2\\:6570\\:7248\\:5fdc\\:7b54\\:30ab\\:30fc\\:30cd\\:30eb\\n\\n'
        f'{cr_tbl}* \\:3068\\:540c\\:3058 w \\:7a4d\\:5206\\:3092\\:884c\\:3046\\:304c\\:3001'
        f'vm \\:3092 ?NumericQ \\:30d1\\:30bf\\:30fc\\:30f3\\:3067\\:53d7\\:3051\\:53d6\\:308a\\:3001\\:30b5\\:30f3\\:30d7\\:30ea\\:30f3\\:30b0\\:305b\\:305a\\:306b'
        f'\\:9023\\:7d9a\\:5024\\:3092\\:8fd4\\:3059\\:3002\\n\\n'
        f'\\:7528\\:9014: vm \\:7a4d\\:5206\\:3092 NIntegrate \\:3067\\:884c\\:3044\\:305f\\:3044\\:3068\\:304d\\:3001'
        f'\\:307e\\:305f\\:306f \\:03b7(vm) \\:3092 \\:7573\\:307f\\:8fbc\\:3080\\:3068\\:304d\\:306b\\:9023\\:7d9a\\:5024\\:304c\\:5fc5\\:8981\\:306a\\:30b1\\:30fc\\:30b9\\:3002\\n\\n'
        f'\\:5b9a\\:7fa9:\\n'
        f'  {label}[md, n][E1, E2][vm] = \\:222b_E1^E2 ({kerL}[md][n][w, {sig}\\:00b7w][vm\\:00b7kps] + {kerR}[\\:2026]) dw\\n'
        f'  {label}Left  = {kerL} \\:306e\\:307f\\n'
        f'  {label}Right = {kerR} \\:306e\\:307f\\n\\n'
        f'\\:5f15\\:6570\\:306e vm \\:306f km/s \\:5358\\:4f4d\\:306e\\:6570\\:5024 (\\:5185\\:90e8\\:3067 vm\\:00b7kps \\:306b\\:5909\\:63db)\\:3002'
    )
    return (
        f'Cell[TextData[{{"{body}"}}], "Text",\n'
        f' Background->RGBColor[0.95, 0.97, 1.],\n'
        f' CellTags->"doc",\n'
        f' ExpressionUUID->"{uuid.uuid4()}"],\n\n'
    )


OPEN_RE = re.compile(r'^Cell\[CellGroupData\[\{\s*$')
CLOSE_RE = re.compile(r'^\}, Open\s*\]\],?\s*$')


def _find_subsubsection_close(lines, sub_name):
    """Locate `Cell["<sub_name>", "Subsubsection"`, then find its CellGroupData
    closing }, Open ]] line by depth tracking. Returns the close line index."""
    header_re = re.compile(rf'^Cell\["{re.escape(sub_name)}", "Subsubsection",')
    header_idx = None
    for i, line in enumerate(lines):
        if header_re.match(line):
            header_idx = i
            break
    if header_idx is None:
        return None
    # Walk back to opening CellGroupData
    op_idx = None
    for j in range(header_idx - 1, -1, -1):
        s = lines[j].rstrip("\n")
        if not s.strip():
            continue
        if OPEN_RE.match(s):
            op_idx = j
            break
        return None
    if op_idx is None:
        return None
    # Depth track forward
    depth = 0
    for k in range(op_idx, len(lines)):
        s = lines[k].rstrip("\n")
        if OPEN_RE.match(s):
            depth += 1
        elif CLOSE_RE.match(s):
            depth -= 1
            if depth == 0:
                return k
    return None


def inject(file_path: Path, code: str, detector: str):
    text = file_path.read_text()
    if SENTINEL in text:
        print(f"  {file_path.name}: continuous kernels already present, skipping")
        return
    lines = text.splitlines(keepends=True)
    close_idx = _find_subsubsection_close(lines, "defined")
    if close_idx is None:
        print(f"  {file_path.name}: could not find 'defined' Subsubsection, skipping")
        return
    doc_cell = _make_doc_cell(detector)
    code_cell = _make_code_cell(code)
    # Insert doc + code just BEFORE the closing }, Open ]]
    lines.insert(close_idx, doc_cell + code_cell)
    file_path.write_text("".join(lines))
    print(f"  {file_path.name}: inserted continuous kernel cell at line {close_idx + 1}")


def main():
    for detector, code in [("TES", TES_CODE), ("MKID", MKID_CODE)]:
        f = BASE / detector / "06_response_defs.nb"
        if not f.exists():
            print(f"  {f}: not found, skipping")
            continue
        print(f"== {detector} ==")
        inject(f, code, detector)


if __name__ == "__main__":
    main()
