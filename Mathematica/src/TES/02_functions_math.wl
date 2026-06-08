(* ========================================================================== *)
(*  TES - Math / Utility Functions                                          *)
(* ========================================================================== *)

fileDir = DirectoryName[$InputFileName];
Direc = FileNameJoin[{ParentDirectory[ParentDirectory[fileDir]], "input", "TES"}];
Get[FileNameJoin[{fileDir, "01_setup.wl"}]];
SetDirectory[Direc];
Print["Direc = ", Direc];


(* ============================ Math ============================ *)

(* LIntpl1D[x, y] — Build a piecewise-linear interpolation from parallel
   vectors x and y, returning an InterpolatingFunction (order 1). *)
LIntpl1D[x1_, y1_] := Module[{l1, Leng, lv},
  Leng = Length[x1];
  l1 = ConstantArray[0, Leng];
  For[lv = 1, lv <= Leng, lv++,
    l1[[lv]] = {{x1[[lv]]}, y1[[lv]]};
  ];
  ListInterpolation[l1, InterpolationOrder -> 1]
];

(* IntplArray[m] — Build an array of linear interpolations.
   m[[1]] is the common x grid; m[[k+1]] (k = 1..len-1) are y vectors. *)
IntplArray[m_] := Module[{len, intpl, lv},
  len = Length[m];
  intpl = ConstantArray[0, len - 1];
  For[lv = 1, lv <= len - 1, lv++,
    intpl[[lv]] = LIntpl1D[m[[1]], m[[lv + 1]]];
  ];
  intpl
];


(* ============================ Utility ============================ *)

(* Inword[word][size][x, y][color] — text annotation inset for figures.
   Places `word` at fractional position (x, y) of the image, with given
   font size and color, using a Times font. *)
Inword[word_][a_][x_, y_][color_] := Inset[
  Style[word, a, color, FontFamily -> "Times"],
  ImageScaled[{x, y}],
  {Center, Center}
];
