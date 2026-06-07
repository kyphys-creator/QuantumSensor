"""Section-level Mathematica notebook splitter."""

import re
import uuid
from pathlib import Path

SECTIONS_TEMPLATE = [
    ("01_setup",              "{D} - Setup (Constants/Units)",      "Constants"),
    ("02_functions_math",     "{D} - Math/Utility Functions",        "Functions",            0),
    ("03_functions_response", "{D} - Domain/FDM/eta/Resolution",     "Functions",            1),
    ("04_material",           "{D} - Material ({M})",                "Material Input"),
    ("05_parameters",         "{D} - Calculation Parameters",        "Parameters"),
    # Response function is split into 3 notebooks: defs (incl. Matrix data > defined), plots, matrices
    ("06_response_defs",      "{D} - Response Function Defs",        "Response function",    0,
        {"only": "defined functions", "append_subsubs": ["defined"]}),
    ("07_dRdE_plot",          "{D} - dR/dE kernel plots",            "Response function",    0,
        {"only": "plot"}),
    ("08_response_matrix",    "{D} - Response (Matrix data)",        "Response function",    0,
        {"only": "Matrix data", "remove_subsubs": ["defined"]}),
    ("09_data",               "{D} - Data for minimization",         "Data for minimization"),
    ("10_minimization",       "{D} - Minimization",                  "Minimization"),
]

OPEN_RE    = re.compile(r'^Cell\[CellGroupData\[\{\s*$')
CLOSE_RE   = re.compile(r'^\}, Open\s*\]\],?\s*$')
SECTION_RE = re.compile(r'^Cell\["([^"]+)", "Section",')


def _find_section_lines(lines):
    return [(i, m.group(1))
            for i, line in enumerate(lines)
            if (m := SECTION_RE.match(line))]


def _find_open_before(lines, header_idx):
    for j in range(header_idx - 1, -1, -1):
        s = lines[j].rstrip("\n")
        if not s.strip():
            continue
        if OPEN_RE.match(s):
            return j
        return None
    return None


def _find_close_from(lines, open_idx):
    depth = 0
    for k in range(open_idx, len(lines)):
        s = lines[k].rstrip("\n")
        if OPEN_RE.match(s):
            depth += 1
        elif CLOSE_RE.match(s):
            depth -= 1
            if depth == 0:
                return k
    return None


def _extract_subblock(lines, name, level):
    """Find Cell["<name>", "<level>"... and return its enclosing CellGroupData lines."""
    pat = re.compile(rf'^Cell\["{re.escape(name)}", "{level}",')
    for i, line in enumerate(lines):
        if pat.match(line):
            op = _find_open_before(lines, i)
            cl = _find_close_from(lines, op) if op is not None else None
            if op is not None and cl is not None:
                return lines[op:cl+1]
    raise ValueError(f'{level} {name!r} not found')


def _list_subsections(lines, level="Subsection"):
    """Return names of cells matching `Cell["<name>", "<level>"`."""
    pat = re.compile(rf'^Cell\["([^"]+)", "{level}",')
    out = []
    for line in lines:
        m = pat.match(line)
        if m:
            out.append(m.group(1))
    return out


def _extract_section(lines, name, occurrence=0):
    matches = [i for i, n in _find_section_lines(lines) if n == name]
    if occurrence >= len(matches):
        raise ValueError(f"Section {name!r}#{occurrence} not found")
    op = _find_open_before(lines, matches[occurrence])
    cl = _find_close_from(lines, op)
    return op, cl


def _remove_subsections(lines, labels, levels=("Subsection", "Subsubsection")):
    """Remove every CellGroupData containing `Cell["<label>", "<level>"`
    for any label in `labels`. Operates on a list of lines, returns new list.

    Walks back/forward exactly like _find_open_before / _find_close_from to
    locate the enclosing group, then drops those lines. Adjusts trailing comma
    on the preceding sibling close if needed.
    """
    if not labels:
        return lines
    header_re = re.compile(
        r'^Cell\["(' + '|'.join(re.escape(l) for l in labels) + r')", "('
        + '|'.join(levels) + r')",'
    )
    # Find all matching headers
    drop_ranges = []
    for i, line in enumerate(lines):
        if header_re.match(line):
            op = _find_open_before(lines, i)
            cl = _find_close_from(lines, op) if op is not None else None
            if op is not None and cl is not None:
                drop_ranges.append((op, cl))
    # Sort descending so deletions don't invalidate earlier indices
    drop_ranges.sort(reverse=True)
    out = lines[:]
    for op, cl in drop_ranges:
        # Check if the close has a trailing comma — meaning a sibling follows.
        # If we drop this and it was the LAST sibling, the previous close needs
        # its trailing comma removed to keep syntax valid.
        was_not_last = lines[cl].rstrip("\n").endswith(",")
        del out[op:cl+1]
        # Also drop a single trailing blank line if present (cosmetic)
        if op < len(out) and out[op].strip() == "":
            del out[op]
        if not was_not_last:
            # The dropped block was the last sibling: the previous sibling's
            # close (now immediately before op) probably has a trailing comma
            # to drop. Find it.
            for k in range(op - 1, -1, -1):
                s = out[k].rstrip("\n")
                if not s.strip():
                    continue
                if CLOSE_RE.match(s) and s.endswith(","):
                    out[k] = s[:-1] + "\n"
                break
    return out


HEADER = '''(* Content-type: application/vnd.wolfram.mathematica *)

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

Cell[{TITLE!}, "Title",
 ExpressionUUID->"{UUID_TITLE}"],

Cell[CellGroupData[{

Cell["Setup", "Section",
 ExpressionUUID->"{UUID_SETUP_SEC}"],

Cell[TextData[{
"Loads upstream notebooks and points ",
StyleBox["Direc", "Code"],
" to ",
StyleBox["../../input/{DETECTOR}", "Code"],
"."
}], "Text",
 ExpressionUUID->"{UUID_NOTE}"],

Cell[
 "{SETUP_BODY}",
 "Input",
 InitializationCell->True,
 ExpressionUUID->"{UUID_DIREC}"],

}, Open  ]],

'''

FOOTER = '''
}, Open  ]]
},
WindowSize->{1920, 954},
WindowMargins->{{0, Automatic}, {Automatic, 0}},
FrontEndVersion->"14.3 for Mac OS X ARM (64-bit) (July 8, 2025)",
StyleDefinitions->"Default.nb",
ExpressionUUID->"{UUID_NB}"
]
(* End of Notebook Content *)
'''


def _setup_body(detector, prev_file):
    """Atomic Setup cell body: compute Direc, Get prev (if any), SetDirectory."""
    lines = [
        'nbDir = NotebookDirectory[];',
        f'Direc = FileNameJoin[{{ParentDirectory[ParentDirectory[nbDir]], \\"input\\", \\"{detector}\\"}}];',
    ]
    if prev_file:
        lines.append(f'Get[FileNameJoin[{{nbDir, \\"{prev_file}\\"}}]];')
    lines += [
        'SetDirectory[Direc];',
        'Print[\\"Direc = \\", Direc];',
    ]
    return '\\n'.join(lines)


def _make_notebook(title, detector, prev_file, section_block):
    head = (HEADER
            .replace('{TITLE!}', f'"{title}"')
            .replace('{DETECTOR}', detector)
            .replace('{SETUP_BODY}', _setup_body(detector, prev_file))
            .replace('{UUID_TITLE}',     str(uuid.uuid4()))
            .replace('{UUID_SETUP_SEC}', str(uuid.uuid4()))
            .replace('{UUID_NOTE}',      str(uuid.uuid4()))
            .replace('{UUID_DIREC}',     str(uuid.uuid4()))
            .replace('{UUID_SETDIR}',    str(uuid.uuid4())))
    foot = FOOTER.replace('{UUID_NB}', str(uuid.uuid4()))
    sb = section_block.rstrip()
    if sb.endswith(','):
        sb = sb[:-1]
    return head + sb + foot


def split_notebook(detector, src_name, material_label, remove_labels=None):
    """remove_labels: dict mapping output-stem (e.g. "04_material") to
    list of Subsection/Subsubsection labels to strip from that file."""
    remove_labels = remove_labels or {}
    base = Path(__file__).parent
    src = base / detector / src_name
    out_dir = base / detector

    # Source might already be moved to backup; try backup too
    if not src.exists():
        bk = base / "backup" / src_name.replace(".nb", "_original.nb")
        if bk.exists():
            src = bk
        else:
            raise FileNotFoundError(f"Could not find {src_name}")

    text = src.read_text()
    lines = text.splitlines(keepends=True)

    sections = []
    for entry in SECTIONS_TEMPLATE:
        opts = {}
        if len(entry) >= 5:
            stem, title_tmpl, name, occ, opts = entry[:5]
        elif len(entry) == 4:
            stem, title_tmpl, name, occ = entry
        else:
            stem, title_tmpl, name = entry
            occ = 0
        title = title_tmpl.format(D=detector, M=material_label)
        sections.append((stem, title, name, occ, opts))

    plans = []
    for stem, title, name, occ, opts in sections:
        op, cl = _extract_section(lines, name, occ)
        plans.append((stem, title, op, cl, opts))
        print(f"{stem}: {name!r}#{occ}  lines {op+1}-{cl+1}  ({cl-op+1} lines)")

    for idx, (stem, title, op, cl, opts) in enumerate(plans):
        prev = plans[idx-1][0] + ".nb" if idx > 0 else None
        section_lines = lines[op:cl+1]
        original_section = section_lines[:]  # snapshot before filtering, for append lookups

        # Apply only/skip filtering on Subsection level
        if "only" in opts:
            all_subs = _list_subsections(section_lines, "Subsection")
            to_remove = [s for s in all_subs if s != opts["only"]]
            if to_remove:
                before = len(section_lines)
                section_lines = _remove_subsections(section_lines, to_remove, levels=("Subsection",))
                print(f"  {stem}: kept only {opts['only']!r}, removed {to_remove} ({before}→{len(section_lines)} lines)")
        if "skip" in opts:
            before = len(section_lines)
            section_lines = _remove_subsections(section_lines, [opts["skip"]], levels=("Subsection",))
            print(f"  {stem}: skipped {opts['skip']!r} ({before}→{len(section_lines)} lines)")

        # Remove specified Subsubsections (anywhere in current section_lines)
        if "remove_subsubs" in opts:
            before = len(section_lines)
            section_lines = _remove_subsections(section_lines, opts["remove_subsubs"], levels=("Subsubsection",))
            print(f"  {stem}: removed Subsubsections {opts['remove_subsubs']} ({before}→{len(section_lines)} lines)")

        # Append Subsubsections from elsewhere in the original Section
        if "append_subsubs" in opts:
            insert_idx = _find_close_from(section_lines, 0)
            blocks = []
            for name in opts["append_subsubs"]:
                blocks.extend(_extract_subblock(original_section, name, "Subsubsection"))
            section_lines = section_lines[:insert_idx] + blocks + section_lines[insert_idx:]
            print(f"  {stem}: appended Subsubsections {opts['append_subsubs']} (+{len(blocks)} lines)")

        if stem in remove_labels:
            before = len(section_lines)
            section_lines = _remove_subsections(section_lines, remove_labels[stem])
            print(f"  {stem}: stripped material subsections {remove_labels[stem]} ({before}→{len(section_lines)} lines)")
        section_block = ''.join(section_lines)
        nb = _make_notebook(title, detector, prev, section_block)
        out_path = out_dir / f"{stem}.nb"
        out_path.write_text(nb)
        print(f"  wrote {out_path.name}  ({len(nb)} bytes)")
