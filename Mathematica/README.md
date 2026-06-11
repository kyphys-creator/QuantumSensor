# QuantumSensor / Mathematica — Response-calculation pipeline (01–12)

📖 **[English](#english) · [日本語](#日本語)**

---

<a id="english"></a>

## English

A set of Wolfram Language packages for computing the **detector response** in
dark-matter (DM) direct detection. It evaluates the DM scattering rate in the
electron-excitation channel using the material's dielectric response
`Im[-1/ε(ω, q)]`, convolves it with the energy resolution, and prepares it up to
the point of integrating over the halo velocity distribution.

Two detector systems:

| System | Material | Dielectric function | Input data |
|------|------|----------|-----------|
| **TES**  | Al   | Mermin model (numerical table) | `input/TES/*.dat` |
| **MKID** | TiN  | Lindhard model (analytic)      | none (analytic only) |

`src/TES/` and `src/MKID/` share the **same structure**; only the kernel names
in 04 (material) and 06 differ per material.

---

### File layout and dependencies

Each file is a **linear dependency chain** in numerical order: the file with the
higher number `Get`s the preceding stage at its top. Loading
`06_response_defs.wl` alone pulls in 01–05 transitively.

```
01_setup ─▶ 02_functions_math ─▶ 03_functions_response ─▶ 04_material ─▶ 05_parameters ─▶ 06_response_defs
 units/consts   interpolation utils    kinematics/FDM/η          dielectric/DOS         xsec/resolution/exposure   response kernel (core)
```

- **01–06 are `.wl`** (hand-written packages; the build script `src/build_wl.py` generates both TES/MKID).
- **07, 07b, 08, 09, 09b, 10, 12 are `.wl`** (kernel/response plots and their CSV exports, saving response functions, building the response matrix, η. Run from the CLI. Both TES/MKID).
- **10_data, 11_minimization are `.nb`** (legacy data/fit notebooks).

---

### What each stage prepares

#### 01_setup — units and physical constants
- Uses **natural units** with every quantity expressed in **GeV** (`GeV = 10^9`).
  `eV`, `keV`, `cm`, `sec`, `Kel` (Kelvin), `grams`, `kg`, etc. are defined as GeV-equivalent factors.
- **DM halo velocities (SHM)**: `v0` (=238 km/s), `ve` (=250 km/s), `vesc` (=544 km/s). Also a **pure dark-disk** model (`v0DD`=70, `veDD`=100, `vescDD`=694) and an **Earth-bound** thermal population (`\[Rho]b`, `v0EB`=√(2kT/mχ), `vescEB`=11.2 km/s).
- **Physics constants / DM density**: `alpha` (fine-structure constant), `\[Rho]DM` (local DM density 0.4 GeV/cm³), `me` (electron mass).
- Material densities `rhoAl`, `rhoTiN`.

> The unit conversions chosen here propagate consistently through every later
> formula. When reading "the order of magnitude" of a number, always come back to
> these unit definitions.

#### 02_functions_math — interpolation utilities
- `LIntpl1D[x, y]` — build a first-order (piecewise-linear) `InterpolatingFunction` from parallel vectors.
- `IntplArray[m]` — interpolate multiple y-columns against a common x grid (`m[[1]]`). Used in 04 to turn `.dat` into functions.
- `Inword` — text inset for figure annotations (used by the plotting stages).

#### 03_functions_response — kinematics, form factor, velocity integral
Defines the **kinematics** of DM–electron scattering and the analytic building
blocks needed for the rate.
- **Momentum transfer q**: two branches `ql` (left) / `qr` (right) for energy transfer `Ee`, with their Jacobians `Jacobl/Jacobr`.
- `vmin` — minimum DM velocity for a given `(Ee, q)`.
- **DM form factor `FDM[q][n] = (α·me / q)^n`** — encodes the mediator nature (`n=0` heavy mediator, `n=2` light mediator). `FmedH/FmedL` are alternative parameterizations.
- **Halo velocity integral η** (`\[Eta]th`, `\[Eta]td`) — `η(vmin)` from integrating a truncated Maxwellian above `vmin`. Normalization is `KKf`.

#### 04_material — material dielectric response (**the one stage that differs most between TES and MKID**)
Prepares the "material-side" input `Im[-1/ε(ω, q)]` of the scattering rate.

- **TES (Al, Mermin)**: reads ε₁, ε₂ from `input/TES/Al_mermin.dat`, interpolates in 2D, and
  builds `ImepsAlf[w, q]`, `eps1Alf`, `eps2Alf` (0 outside the domain `0≤w≤99.3`, `0≤q≤37289.5`).
  It also prepares the phonon density of states `Al_pDoS.dat → Dw[w]` and the form-factor table `Al_Fn.dat → Fnw10GeV1..5`.
- **MKID (TiN, Lindhard)**: no data needed. From the Fermi wavenumber `kFTiN`, plasma frequency `wpTiN`, and damping `GammaTiN`,
  it directly defines the **analytic Lindhard dielectric function `ImepsLTiN[w, q]`**.

> TES is a measured numerical table, MKID an analytic formula. The `input/MKID`
> directory does not exist, but since TiN is analytic this does not affect loading.

#### 05_parameters — calculation parameters
- **Reference cross sections** `\[Sigma]e` (electron 10⁻³⁰ cm²), `\[Sigma]N`.
- **Energy resolutions** `TESsig`, `MKIDsig` (relative resolution, FWHM converted to a Gaussian σ).
- **Exposures** (active mass × time): `Alexp = 8200 μg·month` (TES), `TiNexp = 1e7 × 0.42 ng·yr` (MKID, design per arXiv:2404.10785).

#### 06_response_defs — response kernel (the body of the pipeline)
Assembles the parts of 01–05 and defines the functions that return the response
**already integrated over velocity**. TES uses `…Al…` names and MKID `…TiN…`
names, **fully in parallel** (described below with TES names).

The rate kernel is built as a 4-level nesting:

1. **Raw differential rate kernel** `kerRAll / kerRAlr` (left/right branches):
   `(1/ρ)·(1/8π²α)·(1/μχe²)·Jacob·q³·FDM²·Im[-1/ε]` evaluated at `q = ql/qr`.
2. **Resolution convolution** `IntkerRAll/r` — convolves the true energy `E` with a Gaussian around the observed energy `Ep` (`NormalDistribution[Ep, sig]`) as a Riemann sum (`energySum`).
3. **Auto-range version** `KerRAll/r` — picks a kinematically/resolution-valid integration window `[emin, emax]` via `RangeAlGen`, then runs step 2.
4. **Energy-bin response (continuous vm)** `CRTESLeft/Right/CRTES` — accumulates over the observed energy band `[E1, E2]` with `energySum`.
5. **Velocity-integrated version** `CRintTES…` — midpoint-rule integral over `vmin..vmax` with `midpointSum`. `CRintTESWeighted` is η-weighted.
6. **Legacy Table API** `CRvmTES…` — wrappers tabulating vm, kept for compatibility with old code.

Helpers: `energySum` (energy-direction Riemann sum), `midpointSum` (velocity-direction midpoint rule), `BoxPDF` (rectangular resolution).

---

### Flow of the calculation (what physically happens)

```
   DM flux (ρDM/mχ, velocity distribution)
        │
        ▼   ┌─ kinematics (03: q, vmin, Jacob)
   scattering-rate kernel ┤
        │      └─ material response (04: Im[-1/ε], FDM²)
        ▼
   convolve with energy resolution (06 step2, 05: σ)
        │
        ▼
   accumulate over observed energy band (06 step4)
        │
        ▼
   integrate over halo velocity distribution (06 step5, 03: η)
        │
        ▼
   detector response  →  (07–11: matrix form, fitting, sensitivity)
```

---

### How to load

```wolfram
(* Get-ing 06 chain-loads 01–05 *)
Get["src/TES/06_response_defs.wl"];      (* or src/MKID/... *)

(* e.g. velocity-integrated response at mχ=1 GeV, heavy mediator (n=0) *)
CRintTES[1 GeV, 0][E1, E2][vmin, vmax, ns]
```

---

### Smoke tests

A bundled test verifies that 01–06 are all defined and return sane numbers.

```bash
cd QuantumSensor/Mathematica
./src/tests/run_tests.sh          # both TES and MKID
./src/tests/run_tests.sh TES      # one detector only
```

- Checks three layers: load success, presence of key symbol definitions, and a numerical smoke at representative arguments.
- Reports go to `output/<DET>/test/report.txt` (human) and `report.m` (machine-readable).
- Current status: **TES 20/20, MKID 18/18 pass**.

---

### Post-pipeline (07, 07b, 08, 09, 09b, 10, 12)

These **save and discretize** the response functions defined in 01–06, turning
them into a matrix form usable for data analysis.

The post-pipeline (07–12) is **fully parallel `.wl` across TES and MKID**. The
examples below use TES names (`CRTES`, `Al_…` tags, `output/TES/…`); for MKID
just substitute `CRTiN`, `TiN_…` tags, and `output/MKID/…` — the CLI arguments,
output structure, and physics are identical (e.g. `wolframscript -file
src/MKID/08_response_functions.wl bin5 M3 q0` → `output/MKID/response_functions/TiN_q0M3_R5.wdx`).
The dependency chain (07→06→…→01) is the same for both.

#### 07_dcurlyRdEprime_plot.wl / 07b_dcurlyRdEprime_export.wl — differential kernel dℛ/dE′

**07** plots the differential response kernel `dℛ/dE′(v_min)` (the summed
left+right kinematic branches `KerRAll+KerRAlr`) for the three DM masses and two
mediators, at the bold observed energy E′ (0.1 eV TES / 0.2 eV MKID) and the
1 eV reference, as PDFs to `output/<DET>/dcurlyRdEprime/`.
**07b** writes the same curves to CSV over the full v_min range, in kg⁻¹ eV⁻¹
(`output/<DET>/dcurlyRdEprime_csv/{Al,TiN}_{heavy,light}.csv`), so the Python
`final_figures.ipynb` can plot them.

```bash
wolframscript -file 07_dcurlyRdEprime_plot.wl       # PDFs
wolframscript -file 07b_dcurlyRdEprime_export.wl    # CSVs (full domain)
```

#### 08_response_functions.wl — saving response functions

Evaluates 06's `CRTES` on a v_min grid and saves the per-energy-bin response as
an `InterpolatingFunction` (piecewise-linear) to `.wdx`.

```bash
wolframscript -file 08_response_functions.wl bin5             # 5 bins, all masses, heavy mediator
wolframscript -file 08_response_functions.wl bin10 M3 q2      # (TES) 10 bins, 1 GeV, light mediator
```

- **Bins** (TES, threshold 0.1 eV): `bin5` = 0.1–1.1 eV in 0.2 eV steps × 5, `bin10` = 0.1 eV steps × 10.
- **Bins** (MKID, threshold 0.2 eV): `bin5` = `[0.2,0.3],[0.3,0.5],[0.5,0.7],[0.7,0.9],[0.9,1.1]` (5), `bin9` = 0.2–1.1 eV in 0.1 eV steps × 9 (MKID uses `bin9`, not `bin10`).
- **Mass**: `M1`(10 MeV) / `M2`(100 MeV) / `M3`(1 GeV) / `ALL`.
- **Mediator**: `q0`(heavy, n=0) / `q2`(light, n=2).
- **Output**: `output/TES/response_functions/<name>.wdx` (MKID: `output/MKID/...`).

#### 09_response_function_plot.wl — plotting response functions

Reads the `.wdx` saved by 08 and outputs a PDF overlaying each energy bin's
response function R_bin(v_min). Does not need the 01–06 pipeline loaded.

```bash
wolframscript -file 09_response_function_plot.wl          # all files
wolframscript -file 09_response_function_plot.wl M3        # only names containing "M3"
```

- **Output**: `output/TES/response_function_plots/<name>.pdf`

#### 09b_response_function_export.wl — response functions to CSV (full domain)

Samples each `.wdx`'s per-bin `R_bin(v_min)` on a log-spaced grid over its full
domain (1–800 km/s for q0, 1–2000 for q2) and writes one CSV per config
(`output/<DET>/response_functions_csv/<name>.csv`, the [kg⁻¹] quantity
`R_bin·kg` that 09 plots) — for the Python `final_figures.ipynb`, which needs
the full domain the windowed matrix CSVs cannot provide.

```bash
wolframscript -file 09b_response_function_export.wl        # all .wdx
wolframscript -file 09b_response_function_export.wl M3     # filter by substring
```

#### 10_response_matrix.wl — building the response matrix

Integrates the response functions saved by 08 over v_min intervals to build a
response matrix `M[bin_i, interval_j]`. The integration measure is in natural
units (`dv` multiplied by `kps`).

```bash
wolframscript -file 10_response_matrix.wl Al_q0M3_R5 1 800 1000 0.3
wolframscript -file 10_response_matrix.wl ALL 1 800 1000 0.01     # process all .wdx at once
```

- `<name|ALL>` — base name of a `.wdx`, or `ALL` for all of them.
- `<vminLo> <vminHi>` — integration range [km/s] (auto-clipped to the function's domain).
- `<N>` — number of equal v_min intervals (matrix columns).
- `[alpha]` — per-row window tail-cut tolerance: keep each row from its kinematic threshold up to where its cumulative integral reaches `1−alpha` (default 0.5), zeroing the long tail and trimming the matrix/`vmin.csv` to the populated column window. **Currently TES q0 = 0.3, everything else 0.01.**

**Output structure** (per mass → per-run subfolder):

```
output/TES/response_matrix/
  M1/
    Al_q0M1_R5_v1-800_N1000/
      matrix.csv   # response matrix (nBins × N)
      vmin.csv     # per-column v_min interval {lo, hi, mid} [km/s]
      bins.csv     # per-row energy bin {lo, hi} [eV]
  M2/
    ...
  M3/
    ...
```

#### 12_eta.wl — natural-units η on the matrix grid

Samples `etaSHM[dmMass][v_mid]` (`\[Eta]th` from 03, the SHM speed integral) on
the same v_min interval mid-points that label the matrix columns (read from each
matrix folder's `vmin.csv`), and writes `eta_<model>.csv` next to `matrix.csv` so
the Python loader picks it up with no alignment or unit conversion. Models:
**Halo** (SHM), **Disk** (pure dark disk), **Bound** (Earth-bound thermal — an
isotropic dense population nonzero only below `vescEB = 11.2 km/s`).

```bash
wolframscript -file 12_eta.wl ALL                 # every matrix folder (Halo, default)
wolframscript -file 12_eta.wl ALL Disk            # pure dark-disk model
wolframscript -file 12_eta.wl ALL Bound           # Earth-bound model
wolframscript -file 12_eta.wl Al_q0M1_R5 Halo     # one folder
```

---

<a id="日本語"></a>

## 日本語

ダークマター（DM）直接検出における**検出器応答**を計算するための Wolfram Language パッケージ群です。
電子励起チャネルでの DM 散乱レートを、物質の誘電応答 `Im[-1/ε(ω, q)]` を用いて評価し、
エネルギー分解能で畳み込み、ハロー速度分布で積分するところまでを準備します。

検出器は 2 系統：

| 系統 | 物質 | 誘電関数 | 入力データ |
|------|------|----------|-----------|
| **TES**  | Al   | Mermin モデル（数値テーブル） | `input/TES/*.dat` |
| **MKID** | TiN  | Lindhard モデル（解析式）     | なし（解析式のみ） |

`src/TES/` と `src/MKID/` は**同一構造**で、04（物質）と 06 の核（kernel）名だけが物質ごとに異なります。

---

### ファイル構成と依存関係

各ファイルは番号順に、末尾の番号のファイルが先頭で `Get` して前段を読み込む**直線的な依存チェーン**です。
`06_response_defs.wl` を 1 つロードすれば 01–05 が芋づる式に読み込まれます。

```
01_setup ─▶ 02_functions_math ─▶ 03_functions_response ─▶ 04_material ─▶ 05_parameters ─▶ 06_response_defs
 単位・定数      補間ユーティリティ      運動学・FDM・η            誘電関数・DOS          断面積・分解能・露光      応答核（本体）
```

- **01–06 は `.wl`**（手書きパッケージ、ビルドスクリプト `src/build_wl.py` で TES/MKID 両方を生成）
- **07, 07b, 08, 09, 09b, 10, 12 は `.wl`**（カーネル/応答のプロットと CSV 出力・応答関数の保存・応答行列の構築・η。CLI で実行。TES/MKID 両系統）
- **10_data, 11_minimization は `.nb`**（データ・フィット等のレガシー）

---

### 各段階が「何を準備しているか」

#### 01_setup — 単位系と物理定数
- **自然単位系**を採用し、すべての量を **GeV 建て**で表現します（`GeV = 10^9`）。
  `eV`, `keV`, `cm`, `sec`, `Kel`（ケルビン）, `grams`, `kg` などは GeV 換算の係数として定義。
- **DM ハロー速度（SHM）**：`v0`（=238 km/s）, `ve`（=250 km/s）, `vesc`（脱出速度=544 km/s）。さらに**純ダークディスク**（`v0DD`=70, `veDD`=100, `vescDD`=694）と**地球束縛**熱的成分（`\[Rho]b`, `v0EB`=√(2kT/mχ), `vescEB`=11.2 km/s）も定義。
- **物理定数・DM 密度**：`alpha`（微細構造定数）, `\[Rho]DM`（局所 DM 密度 0.4 GeV/cm³）, `me`（電子質量）。
- 物質密度 `rhoAl`, `rhoTiN`。

> ここで決めた単位換算が以降すべての数式に一貫して効きます。「数値の桁」を読むときは必ずこの単位定義に立ち返ってください。

#### 02_functions_math — 補間ユーティリティ
- `LIntpl1D[x, y]` — 並行ベクトルから 1 次（区分線形）の `InterpolatingFunction` を構築。
- `IntplArray[m]` — 共通 x グリッド（`m[[1]]`）に対する複数 y 列を一括補間。04 で `.dat` を関数化するのに使用。
- `Inword` — 図注釈用のテキスト inset（プロット系で使用）。

#### 03_functions_response — 運動学・形状因子・速度積分
DM–電子散乱の**運動学**と、レート計算に必要な解析的部品を定義します。
- **運動量移行 q**：エネルギー移行 `Ee` に対する 2 つの分岐 `ql`（左）/ `qr`（右）と、そのヤコビアン `Jacobl/Jacobr`。
- `vmin` — 与えられた `(Ee, q)` に対する DM 最小速度。
- **DM 形状因子 `FDM[q][n] = (α·me / q)^n`** — 媒介子の性質を表す（`n=0` 重い媒介子、`n=2` 軽い媒介子）。`FmedH/FmedL` は別パラメータ化。
- **ハロー速度積分 η**（`\[Eta]th`, `\[Eta]td`）— 切り詰めマクスウェル分布を `vmin` 以上で積分した `η(vmin)`。規格化は `KKf`。

#### 04_material — 物質の誘電応答（**TES と MKID で唯一大きく異なる段階**）
散乱レートの「物質側」の入力 `Im[-1/ε(ω, q)]` を準備します。

- **TES（Al, Mermin）**：`input/TES/Al_mermin.dat` から ε₁, ε₂ を読み込み 2 次元補間し、
  `ImepsAlf[w, q]`, `eps1Alf`, `eps2Alf` を構築（定義域 `0≤w≤99.3`, `0≤q≤37289.5` 外は 0）。
  さらにフォノン状態密度 `Al_pDoS.dat → Dw[w]`、フォームファクター表 `Al_Fn.dat → Fnw10GeV1..5` も準備。
- **MKID（TiN, Lindhard）**：データ不要。フェルミ波数 `kFTiN`、プラズマ振動数 `wpTiN`、減衰 `GammaTiN` から
  **解析的 Lindhard 誘電関数 `ImepsLTiN[w, q]`** を直接定義。

> TES は実測ベースの数値テーブル、MKID は解析式という違い。`input/MKID` ディレクトリは存在しませんが、
> TiN は解析式なのでロードに影響しません。

#### 05_parameters — 計算パラメータ
- **参照断面積** `\[Sigma]e`（電子 10⁻³⁰ cm²）, `\[Sigma]N`。
- **エネルギー分解能** `TESsig`, `MKIDsig`（FWHM をガウス σ に換算した相対分解能）。
- **露光量**（有効質量 × 時間）：`Alexp = 8200 μg·month`（TES）, `TiNexp = 1e7 × 0.42 ng·yr`（MKID, arXiv:2404.10785 準拠）。

#### 06_response_defs — 応答核（パイプラインの本体）
01–05 の部品を組み上げ、最終的に**速度積分まで済んだ応答**を返す関数群を定義します。
TES は `…Al…`、MKID は `…TiN…` という名前で**完全に並行**しています（以下 TES 名で記述）。

レート核は 4 段の入れ子で組み立てられます：

1. **生の微分レート核** `kerRAll / kerRAlr`（左右分岐）
   `(1/ρ)·(1/8π²α)·(1/μχe²)·Jacob·q³·FDM²·Im[-1/ε]` を `q = ql/qr` で評価。
2. **分解能畳み込み** `IntkerRAll/r` — 真のエネルギー `E` を観測エネルギー `Ep` 周りのガウス（`NormalDistribution[Ep, sig]`）で畳み込み（`energySum` によるリーマン和）。
3. **自動レンジ版** `KerRAll/r` — `RangeAlGen` で運動学的・分解能的に有効な積分区間 `[emin, emax]` を決めてから 2 を実行。
4. **エネルギービン応答（vm 連続）** `CRTESLeft/Right/CRTES` — 観測エネルギー帯 `[E1, E2]` を `energySum` で積算。
5. **速度積分版** `CRintTES…` — `midpointSum` で `vmin..vmax` を中点則積分。`CRintTESWeighted` は η 重み付き。
6. **レガシー Table API** `CRvmTES…` — 旧コードとの互換用に vm をテーブル化するラッパー。

補助ヘルパ：`energySum`（エネルギー方向リーマン和）, `midpointSum`（速度方向中点則）, `BoxPDF`（矩形分解能）。

---

### 計算の流れ（物理的に何が起きているか）

```
   DM 流束（ρDM/mχ, 速度分布）
        │
        ▼   ┌─ 運動学（03: q, vmin, Jacob）
   散乱レート核  ┤
        │      └─ 物質応答（04: Im[-1/ε], FDM²）
        ▼
   エネルギー分解能で畳み込み（06 step2, 05: σ）
        │
        ▼
   観測エネルギー帯で積算（06 step4）
        │
        ▼
   ハロー速度分布で積分（06 step5, 03: η）
        │
        ▼
   検出器応答  →（07–11 で行列化・フィット・感度導出へ）
```

---

### ロード方法

```wolfram
(* 06 を Get すれば 01–05 も連鎖ロードされる *)
Get["src/TES/06_response_defs.wl"];      (* または src/MKID/... *)

(* 例：質量 mχ=1 GeV, 重い媒介子 (n=0) での速度積分応答 *)
CRintTES[1 GeV, 0][E1, E2][vmin, vmax, ns]
```

---

### 動作確認（スモークテスト）

01–06 がすべて定義され、数値が健全に返ることを検証するテストを同梱しています。

```bash
cd QuantumSensor/Mathematica
./src/tests/run_tests.sh          # TES と MKID 両方
./src/tests/run_tests.sh TES      # 片方だけ
```

- ロード成否・主要シンボルの定義有無・代表引数での数値スモークの 3 層を検査。
- レポートは `output/<DET>/test/report.txt`（人間用）と `report.m`（機械可読）に出力。
- 現状：**TES 20/20・MKID 18/18 パス**。

---

### 後段パイプライン (07, 07b, 08, 09, 09b, 10, 12)

01–06 で定義した応答関数を**保存・離散化**して、データ解析に使える行列形式に変換します。

後段（07–12）は **TES と MKID で完全に並行**な `.wl` です。下記の例は TES 名（`CRTES`,
`Al_…` タグ, `output/TES/…`）で記述しますが、MKID では `CRTiN`・`TiN_…` タグ・`output/MKID/…`
に置き換わるだけで、CLI 引数・出力構造・物理は同一です（例: `wolframscript -file
src/MKID/08_response_functions.wl bin5 M3 q0` → `output/MKID/response_functions/TiN_q0M3_R5.wdx`）。
依存チェーン（07→06→…→01）も両系統で同じです。

#### 07_dcurlyRdEprime_plot.wl / 07b_dcurlyRdEprime_export.wl — 微分カーネル dℛ/dE′

**07** は微分応答カーネル `dℛ/dE′(v_min)`（左右分岐の和 `KerRAll+KerRAlr`）を、3 質量・2 媒介子について、太線の観測エネルギー E′（TES 0.1 eV / MKID 0.2 eV）と参照の 1 eV で PDF 出力します（`output/<DET>/dcurlyRdEprime/`）。
**07b** は同じ曲線を全 v_min 域で kg⁻¹ eV⁻¹ の CSV に出力し（`output/<DET>/dcurlyRdEprime_csv/{Al,TiN}_{heavy,light}.csv`）、Python の `final_figures.ipynb` で描けるようにします。

```bash
wolframscript -file 07_dcurlyRdEprime_plot.wl       # PDF
wolframscript -file 07b_dcurlyRdEprime_export.wl    # CSV（全定義域）
```

#### 08_response_functions.wl — 応答関数の保存

06 の `CRTES` を v_min グリッド上で評価し、各エネルギービンごとの応答を `InterpolatingFunction`（区分線形）として `.wdx` に保存します。

```bash
wolframscript -file 08_response_functions.wl bin5             # 5 bins, 全質量, heavy mediator
wolframscript -file 08_response_functions.wl bin10 M3 q2      # (TES) 10 bins, 1 GeV, light mediator
```

- **ビン**（TES, threshold 0.1 eV）: `bin5` = 0.1–1.1 eV を 0.2 eV 幅 × 5 本、`bin10` = 0.1 eV 幅 × 10 本
- **ビン**（MKID, threshold 0.2 eV）: `bin5` = `[0.2,0.3],[0.3,0.5],[0.5,0.7],[0.7,0.9],[0.9,1.1]`（5 本）、`bin9` = 0.2–1.1 eV を 0.1 eV 幅 × 9 本（MKID は `bin10` ではなく `bin9`）
- **質量**: `M1`(10 MeV) / `M2`(100 MeV) / `M3`(1 GeV) / `ALL`
- **媒介子**: `q0`(heavy, n=0) / `q2`(light, n=2)
- **出力**: `output/TES/response_functions/<name>.wdx`（MKID は `output/MKID/...`）

#### 09_response_function_plot.wl — 応答関数のプロット

08 で保存した `.wdx` を読み込み、各エネルギービンの応答関数 R_bin(v_min) を重ね描きした PDF を出力します。
01–06 パイプラインのロードは不要です。

```bash
wolframscript -file 09_response_function_plot.wl          # 全ファイル
wolframscript -file 09_response_function_plot.wl M3        # 名前に "M3" を含むもののみ
```

- **出力**: `output/TES/response_function_plots/<name>.pdf`

#### 09b_response_function_export.wl — 応答関数を CSV 出力（全定義域）

各 `.wdx` のビンごとの `R_bin(v_min)` を全定義域（q0 は 1–800 km/s、q2 は 1–2000）で対数グリッド上にサンプリングし、設定ごとに 1 つの CSV を出力します（`output/<DET>/response_functions_csv/<name>.csv`、09 が描く [kg⁻¹] 量 `R_bin·kg`）。窓掛け行列 CSV では得られない全定義域が必要な Python の `final_figures.ipynb` 用です。

```bash
wolframscript -file 09b_response_function_export.wl        # 全 .wdx
wolframscript -file 09b_response_function_export.wl M3     # 部分文字列でフィルタ
```

#### 10_response_matrix.wl — 応答行列の構築

08 で保存した応答関数を v_min 区間上で積分し、レスポンス行列 `M[bin_i, interval_j]` を構築します。
積分の測度は自然単位（`dv` に `kps` を乗算）で計算されます。

```bash
wolframscript -file 10_response_matrix.wl Al_q0M3_R5 1 800 1000 0.3
wolframscript -file 10_response_matrix.wl ALL 1 800 1000 0.01     # 全 .wdx を一括処理
```

- `<name|ALL>` — `.wdx` のベース名、または `ALL` で一括
- `<vminLo> <vminHi>` — 積分範囲 [km/s]（関数の定義域に自動クリップ）
- `<N>` — v_min の等分割数（行列の列数）
- `[alpha]` — 行ごとのウィンドウ裾カット許容度：各行をその運動学的閾値から累積積分が `1−alpha` に達する点まで残し（既定 0.5）、長い裾を 0 にして行列・`vmin.csv` を実際に値のある列範囲にトリム。**現状 TES q0 = 0.3、それ以外は 0.01。**

**出力構造**（質量別 → 実行ごとのサブフォルダ）:

```
output/TES/response_matrix/
  M1/
    Al_q0M1_R5_v1-800_N1000/
      matrix.csv   # レスポンス行列 (nBins × N)
      vmin.csv     # 各列の v_min 区間 {下限, 上限, 中央値} [km/s]
      bins.csv     # 各行のエネルギービン {下限, 上限} [eV]
  M2/
    ...
  M3/
    ...
```

#### 12_eta.wl — 行列グリッド上の自然単位 η

`etaSHM[dmMass][v_mid]`（03 の `\[Eta]th`、SHM 速度積分）を、行列の列を規定する v_min 区間の中央値
（各行列フォルダの `vmin.csv` から読み込み）と同じ点で評価し、`matrix.csv` の隣に `eta_<model>.csv` を
出力します。Python ローダがアライメント・単位変換なしで読めるようにするためです。モデルは **Halo**（SHM）, **Disk**（純ダークディスク）, **Bound**（地球束縛熱的 — 等方的で密な成分。`vescEB = 11.2 km/s` 未満でのみ非ゼロ）。

```bash
wolframscript -file 12_eta.wl ALL                 # 全行列フォルダ（Halo, 既定）
wolframscript -file 12_eta.wl ALL Disk            # 純ダークディスク
wolframscript -file 12_eta.wl ALL Bound           # 地球束縛
wolframscript -file 12_eta.wl Al_q0M1_R5 Halo     # 1 フォルダ
```
