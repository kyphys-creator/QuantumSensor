#!/usr/bin/env python3
"""Inject velocity-integrated response functions into 06_response_defs.nb.

Adds CRintTES*, CRintTiN* — these wrap CRTES*/CRTiN* with a Riemann-sum
integration over vm from vmin to vmax, sampled at `nsamples` midpoints.

Run after add_continuous_kernels.py (depends on CRTES*/CRTiN*).
Idempotent via SENTINEL comment.
"""

import re
import uuid
from pathlib import Path

BASE = Path(__file__).parent
SENTINEL = "QSENS-VEL-INT-v1"


TES_CODE = """\
(* Velocity-integrated versions of CRTES* *)
(* Midpoint rule: nsamples の点を [vmin, vmax] の中点でサンプリング、dvm = (vmax-vmin)/nsamples を掛けて和を取る。 *)
CRintTES[md_, n_][E1_, E2_][vmin_?NumericQ, vmax_?NumericQ, nsamples_Integer] := Module[
  {dvm, vmList},
  dvm = (vmax - vmin) / nsamples;
  vmList = Table[vmin + (i - 0.5) dvm, {i, 1, nsamples}];
  Total[CRTES[md, n][E1, E2][#] & /@ vmList] * dvm
];

CRintTESLeft[md_, n_][E1_, E2_][vmin_?NumericQ, vmax_?NumericQ, nsamples_Integer] := Module[
  {dvm, vmList},
  dvm = (vmax - vmin) / nsamples;
  vmList = Table[vmin + (i - 0.5) dvm, {i, 1, nsamples}];
  Total[CRTESLeft[md, n][E1, E2][#] & /@ vmList] * dvm
];

CRintTESRight[md_, n_][E1_, E2_][vmin_?NumericQ, vmax_?NumericQ, nsamples_Integer] := Module[
  {dvm, vmList},
  dvm = (vmax - vmin) / nsamples;
  vmList = Table[vmin + (i - 0.5) dvm, {i, 1, nsamples}];
  Total[CRTESRight[md, n][E1, E2][#] & /@ vmList] * dvm
];

(* Halo-weighted version: 速度分布関数 etaFunc[vm] を掛けて積分 *)
CRintTESWeighted[md_, n_][E1_, E2_][vmin_?NumericQ, vmax_?NumericQ, nsamples_Integer, etaFunc_] := Module[
  {dvm, vmList},
  dvm = (vmax - vmin) / nsamples;
  vmList = Table[vmin + (i - 0.5) dvm, {i, 1, nsamples}];
  Total[(CRTES[md, n][E1, E2][#] * etaFunc[#]) & /@ vmList] * dvm
];"""


MKID_CODE = """\
(* Velocity-integrated versions of CRTiN* *)
(* Midpoint rule: nsamples の点を [vmin, vmax] の中点でサンプリング、dvm = (vmax-vmin)/nsamples を掛けて和を取る。 *)
CRintTiN[md_, n_][E1_, E2_][vmin_?NumericQ, vmax_?NumericQ, nsamples_Integer] := Module[
  {dvm, vmList},
  dvm = (vmax - vmin) / nsamples;
  vmList = Table[vmin + (i - 0.5) dvm, {i, 1, nsamples}];
  Total[CRTiN[md, n][E1, E2][#] & /@ vmList] * dvm
];

CRintTiNLeft[md_, n_][E1_, E2_][vmin_?NumericQ, vmax_?NumericQ, nsamples_Integer] := Module[
  {dvm, vmList},
  dvm = (vmax - vmin) / nsamples;
  vmList = Table[vmin + (i - 0.5) dvm, {i, 1, nsamples}];
  Total[CRTiNLeft[md, n][E1, E2][#] & /@ vmList] * dvm
];

CRintTiNRight[md_, n_][E1_, E2_][vmin_?NumericQ, vmax_?NumericQ, nsamples_Integer] := Module[
  {dvm, vmList},
  dvm = (vmax - vmin) / nsamples;
  vmList = Table[vmin + (i - 0.5) dvm, {i, 1, nsamples}];
  Total[CRTiNRight[md, n][E1, E2][#] & /@ vmList] * dvm
];

(* Halo-weighted version *)
CRintTiNWeighted[md_, n_][E1_, E2_][vmin_?NumericQ, vmax_?NumericQ, nsamples_Integer, etaFunc_] := Module[
  {dvm, vmList},
  dvm = (vmax - vmin) / nsamples;
  vmList = Table[vmin + (i - 0.5) dvm, {i, 1, nsamples}];
  Total[(CRTiN[md, n][E1, E2][#] * etaFunc[#]) & /@ vmList] * dvm
];"""


def _mma_escape(s):
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
    s = s.replace("\\", "\\\\")
    s = s.replace('"', '\\"')
    return _mma_escape(s)


def _make_code_cell(code, marker=True):
    content = _escape_for_cell_string(code)
    marker_text = f' (* {SENTINEL} *)' if marker else ''
    return (
        f'Cell["{content}", "Input",\n'
        f' InitializationCell->True,\n'
        f' Background->RGBColor[0.95, 0.98, 0.95],\n'
        f' ExpressionUUID->"{uuid.uuid4()}"],{marker_text}\n\n'
    )


def _make_doc_cell(detector):
    label = "CRintTES" if detector == "TES" else "CRintTiN"
    cont = "CRTES" if detector == "TES" else "CRTiN"
    body = (
        f'\\:25a0 {label}, {label}Left, {label}Right, {label}Weighted   '
        f'\\:2014 vm \\:65b9\\:5411\\:7a4d\\:5206 (midpoint rule)\\n\\n'
        f'\\:5b9a\\:7fa9:\\n'
        f'  {label}[md, n][E1, E2][vmin, vmax, nsamples]\\n'
        f'    = \\:222b_vmin^vmax {cont}[md, n][E1, E2][vm] dvm\\n'
        f'    \\:2248 dvm \\:00d7 \\:03a3_{{i=1..nsamples}} {cont}[md, n][E1, E2][vm_i]\\n'
        f'    \\:3053\\:3053\\:3067 vm_i = vmin + (i - 0.5) \\:00b7 dvm,  dvm = (vmax - vmin) / nsamples\\n\\n'
        f'\\:5f15\\:6570:\\n'
        f'  \\:2022 md           DM \\:8cea\\:91cf\\n'
        f'  \\:2022 n            \\:5a92\\:4ecb\\:5b50\\:30e2\\:30c7\\:30eb (1/2/3)\\n'
        f'  \\:2022 E1, E2       \\:89b3\\:6e2c\\:30a8\\:30cd\\:30eb\\:30ae\\:30fc bin (eV)\\n'
        f'  \\:2022 vmin, vmax   vm \\:7a4d\\:5206\\:7bc4\\:56f2 (km/s, \\:6570\\:5024)\\n'
        f'  \\:2022 nsamples     \\:30b5\\:30f3\\:30d7\\:30ea\\:30f3\\:30b0\\:70b9\\:6570 (\\:6574\\:6570)\\n\\n'
        f'\\:7a2e\\:985e:\\n'
        f'  \\:2022 {label}        \\:4e21\\:30d6\\:30e9\\:30f3\\:30c1\\:548c\\:306e vm \\:7a4d\\:5206\\n'
        f'  \\:2022 {label}Left    \\:5de6\\:30d6\\:30e9\\:30f3\\:30c1\\:306e\\:307f\\n'
        f'  \\:2022 {label}Right   \\:53f3\\:30d6\\:30e9\\:30f3\\:30c1\\:306e\\:307f\\n'
        f'  \\:2022 {label}Weighted   etaFunc[vm] \\:3092\\:639b\\:3051\\:3066\\:7a4d\\:5206 (\\:30cf\\:30ed\\:901f\\:5ea6\\:5206\\:5e03 \\:03b7(vm) \\:7528)\\n\\n'
        f'\\:4f8b:\\n'
        f'  result = {label}[100 MeV, 1][1.0, 1.5][5, 200, 100]\\n'
        f'  result = {label}Weighted[100 MeV, 1][1.0, 1.5][5, 200, 100, etaSHM]'
    )
    return (
        f'Cell[TextData[{{"{body}"}}], "Text",\n'
        f' Background->RGBColor[0.95, 0.98, 0.95],\n'
        f' CellTags->"doc",\n'
        f' ExpressionUUID->"{uuid.uuid4()}"],\n\n'
    )


OPEN_RE = re.compile(r'^Cell\[CellGroupData\[\{\s*$')
CLOSE_RE = re.compile(r'^\}, Open\s*\]\],?\s*$')


def _find_subsubsection_close(lines, sub_name):
    header_re = re.compile(rf'^Cell\["{re.escape(sub_name)}", "Subsubsection",')
    header_idx = None
    for i, line in enumerate(lines):
        if header_re.match(line):
            header_idx = i
            break
    if header_idx is None:
        return None
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
        print(f"  {file_path.name}: velocity integrals already present, skipping")
        return
    lines = text.splitlines(keepends=True)
    close_idx = _find_subsubsection_close(lines, "defined")
    if close_idx is None:
        print(f"  {file_path.name}: could not find 'defined' Subsubsection, skipping")
        return
    doc_cell = _make_doc_cell(detector)
    code_cell = _make_code_cell(code)
    lines.insert(close_idx, doc_cell + code_cell)
    file_path.write_text("".join(lines))
    print(f"  {file_path.name}: inserted velocity-integral cell at line {close_idx + 1}")


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
