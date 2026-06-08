#!/usr/bin/env wolframscript

(* Convert a single .nb file to .wl by extracting Code/Input cell contents.

   String cells are emitted verbatim. BoxData cells are converted to InputForm
   text via MakeExpression + HoldForm.

   Usage:  wolframscript -file convert_nb_to_wl.wls <input.nb> <output.wl>

   Replaces NotebookDirectory[] with DirectoryName[$InputFileName] in the
   Setup cell. Updates Get[FileNameJoin[{nbDir, "XX.nb"}]] → ...wl path.
*)

args = Rest[$ScriptCommandLine];
If[Length[args] < 2,
  Print["usage: convert_nb_to_wl.wls <input.nb> <output.wl>"];
  Exit[1]
];

inputPath  = args[[1]];
outputPath = args[[2]];

stripHoldForm[s_String] := StringReplace[
  s,
  {StartOfString ~~ "HoldForm[" -> "", "]" ~~ EndOfString -> ""}
];

(* Convert one box statement to InputForm text *)
boxStmtToText[box_] := Module[{expr, text},
  expr = Quiet@MakeExpression[box, StandardForm];
  If[Head[expr] =!= HoldComplete, Return[""]];
  text = ToString[HoldForm @@ expr, InputForm, PageWidth -> 120];
  StringTrim @ stripHoldForm[text]
];

(* Convert a cell's content to text *)
cellToText[s_String] := s;
cellToText[BoxData[boxes_]] := Module[{stmts, results},
  stmts = If[Head[boxes] === List,
    DeleteCases[boxes, _String?(StringMatchQ[#, Whitespace | "\n"] &)],
    {boxes}
  ];
  results = Select[boxStmtToText /@ stmts, # =!= "" && # =!= "Null" &];
  StringRiffle[results, "\n"]
];
cellToText[_] := "";  (* fallback *)

(* Update Setup pattern: NotebookDirectory[] → DirectoryName[$InputFileName] *)
(* Also rewrite Get[FileNameJoin[{nbDir, "XX.nb"}]] → ...wl              *)
postProcess[code_String] := StringReplace[code,
  {
    (* MakeExpression artifacts: cells with embedded comments produce
       (EXPR; )*Null  or  Null*((s1; )*Null*(s2; ))  — clean these up *)
    RegularExpression["\\(([^()]+?;)\\s*\\)(\\*Null)+"] :> "$1",
    RegularExpression["Null\\*\\(([^()]+?;)\\s*\\)"]    :> "$1",
    RegularExpression["\\)\\*\\("]                      :> "); (",
    "*Null"                                              -> "",
    "Null*"                                              -> "",

    (* Setup pattern *)
    "NotebookDirectory[]"            -> "DirectoryName[$InputFileName]",
    "nbDir = DirectoryName[$InputFileName]" -> "fileDir = DirectoryName[$InputFileName]",
    RegularExpression["nbDir"]       -> "fileDir",

    (* Get-chain *)
    RegularExpression["\"(\\d\\d_[A-Za-z_]+)\\.nb\""] -> "\"$1.wl\""
  }
];

Print["Reading ", inputPath, "..."];
nb = Get[inputPath];
If[Head[nb] =!= Notebook,
  Print["Error: not a notebook expression"];
  Exit[1]
];

cells = Cases[nb, Cell[content_, "Code" | "Input", ___] :> content, Infinity];
Print["  ", Length[cells], " code/input cells"];

texts = cellToText /@ cells;
combined = StringRiffle[Select[texts, # =!= "" &], "\n\n"];
combined = postProcess[combined];

(* Add a header comment *)
header = StringTemplate[
  "(* ========================================================================== *)
(*  `name`\n*)
(*  Auto-converted from `srcname` by convert_nb_to_wl.wls\n*)
(* ========================================================================== *)\n\n"
][<|"name" -> FileBaseName[outputPath], "srcname" -> FileNameTake[inputPath]|>];

Export[outputPath, header <> combined, "Text"];
Print["Wrote ", outputPath, " (", FileByteCount[outputPath], " bytes)"];
