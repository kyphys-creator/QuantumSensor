# QuantumSensor

**[English](#english) | [日本語](#日本語)**

Dark-matter direct-detection sensitivity studies for quantum sensors (TES and
MKID). The project has two stages: a **Mathematica** stage that computes the
detector *response* from material physics, and a **Python** stage that uses
that response to recover the dark-matter velocity distribution and derive
sensitivities.

---

## English

### What this is

Quantum sensors such as **TES** (Transition-Edge Sensors, here aluminium) and
**MKID** (Microwave Kinetic Inductance Detectors, here TiN) can detect tiny
energy depositions from dark-matter scattering. This repository:

1. **Computes the detector response** — the per-energy response kernel
   `dℛ/dω'` and the response matrix — from first principles, using the
   material's dielectric function `Im[-1/ε(ω, q)]`, the dark-matter kinematics,
   the halo velocity distribution, and the detector energy resolution.
2. **Optimises and analyses** — feeds that response into a Python pipeline that
   fits / optimises the signal against background scenarios and produces
   sensitivity results.

All physics quantities in the Mathematica stage are in **natural units (GeV)**;
results are converted to physical units (e.g. `kg⁻¹ eV⁻¹`) for plotting.

### Repository layout

| Folder | Stage | Contents |
|--------|-------|----------|
| [`Mathematica/`](Mathematica/) | Response (Wolfram) | The `01`–`12` pipeline that builds the response kernels, response functions and matrices for TES (Al) and MKID (TiN). See **[Mathematica/README.md](Mathematica/README.md)** for the detailed per-file description. |
| [`qsensor_analysis/`](qsensor_analysis/) | Analysis (Python) | The analysis package (`quantum_sensor`), material-parameterised for **both** detectors (`material="Al"`/`"TiN"`): `DarkMatterQuantumAnalysis` runs the self-consistent forward+inverse model on the response matrices. See **[qsensor_analysis/README.md](qsensor_analysis/README.md)**. |
| [`Obsolete/`](Obsolete/) | Archive | Earlier Python versions (`TES`, `TES_clean`, `TES_clean10`), kept for reference. |

### The Mathematica pipeline (summary)

The `01`–`12` files form a single dependency chain — loading one file `Get`s
all the earlier ones. `01`–`06` are clean hand-written `.wl` packages; the
plotting / CSV-export / response-function / matrix / η stages (`07`, `07b`, `08`,
`09`, `09b`, `10`, `12`) are `.wl` for **both** TES and MKID. Only the legacy
`10_data` / `11_minimization` remain as `.nb`.

```
01 setup        units & constants (natural units, GeV)
02 functions    interpolation utilities
03 functions    kinematics, DM form factor F_DM, halo speed integral η
04 material     dielectric function  (Al: Mermin table / TiN: analytic Lindhard)
05 parameters   cross sections, resolutions, exposures
06 response     response-kernel definitions (the core)
07 dℛ/dω' plots kernel plots → output/<DET>/dcurlyRdEprime/*.pdf
07b dℛ/dω' csv  full-domain kernel curves → output/<DET>/dcurlyRdEprime_csv/*.csv
08 response fn  binned CRTES/CRTiN saved as .wdx
09 plots        response-function R_bin(v_min) plots
09b response csv full-domain R_bin(v_min) → output/<DET>/response_functions_csv/*.csv
10 matrix       v_min-integrated response matrix (matrix.csv, vmin.csv, bins.csv;
                per-row window kept up to 1−alpha, TES q0=0.3 else 0.01)
12 eta          natural-units η(v_min) on the matrix grid: Halo / Disk / Bound
   (10_data, 11_minimization remain as legacy .nb)
```

`<DET>` = `TES` (Al) or `MKID` (TiN); the `07`–`12` stages exist in both
`src/TES/` and `src/MKID/`. See **[Mathematica/README.md](Mathematica/README.md)**
for what each stage prepares.

### Running it

**Mathematica stage** (needs `wolframscript`; on macOS it ships inside
`Wolfram.app`):

```bash
cd Mathematica
# generate the dℛ/dω' figures (07 chains in 06 → 01):
wolframscript -file src/TES/07_dcurlyRdEprime_plot.wl
# smoke-test that 01–06 are all defined (TES & MKID):
./src/tests/run_tests.sh
```

Tip: if `wolframscript` is not on your `PATH`, either call it by full path
(`/Applications/Wolfram.app/Contents/MacOS/wolframscript`) or symlink it once
into `/usr/local/bin`.

**Python stage** (the packaged version):

```bash
cd qsensor_analysis
pip install -e .            # installs the `quantum_sensor` package
jupyter notebook main_Refined.ipynb    # run the analysis (TES/Al, MKID/TiN)
jupyter notebook final_figures.ipynb   # publication figures: eta recovery, per-bin signal, response functions, dℛ/dE'
```

The Mathematica stage writes response matrices
(`Mathematica/output/TES/response_matrix/` for Al, `…/MKID/…` for TiN) that the
Python stage consumes; `material="Al"`/`"TiN"` selects which. `final_figures.ipynb`
then builds paper-ready, colour-blind- and projector-safe figures from the
stored results (under `qsensor_analysis/results/final/`).

---

## 日本語

### 概要

**TES**（超伝導転移端センサー、ここではアルミ）や **MKID**（マイクロ波力学
インダクタンス検出器、ここでは TiN）といった量子センサーは、ダークマター
散乱による微小なエネルギー付与を検出できます。本リポジトリは:

1. **検出器応答を計算する** — 物質の誘電関数 `Im[-1/ε(ω, q)]`、ダークマター
   の運動学、ハロー速度分布、検出器のエネルギー分解能から、エネルギー応答核
   `dℛ/dω'` と応答行列を第一原理的に計算します。
2. **最適化・解析する** — その応答を Python パイプラインに渡し、背景事象に対して
   信号をフィット／最適化し、ダークマターの速度分布回復・感度を導出します。

Mathematica 段の物理量はすべて**自然単位系（GeV建て）**で、プロット時に
物理単位（例: `kg⁻¹ eV⁻¹`）へ変換します。

### リポジトリ構成

| フォルダ | 段階 | 内容 |
|----------|------|------|
| [`Mathematica/`](Mathematica/) | 応答（Wolfram） | TES(Al)・MKID(TiN) の応答核・行列を構築する `01`–`12` パイプライン。各ファイルの詳細は **[Mathematica/README.md](Mathematica/README.md)** を参照。 |
| [`qsensor_analysis/`](qsensor_analysis/) | 解析（Python） | 解析パッケージ（`quantum_sensor`）。`material="Al"`/`"TiN"` で**両検出器**に対応し、`DarkMatterQuantumAnalysis` が応答行列に対して自己整合な順方向＋逆問題を実行。詳細は **[qsensor_analysis/README.md](qsensor_analysis/README.md)**。 |
| [`Obsolete/`](Obsolete/) | アーカイブ | 旧Python版（`TES`, `TES_clean`, `TES_clean10`）。参照用に保管。 |

### Mathematica パイプライン（要約）

`01`–`12` は一直線の依存チェーンで、1つを `Get` すると前段がすべて読み込まれ
ます。`01`–`06` は手書きのクリーンな `.wl` パッケージ。プロット／CSV出力／応答関数／行列／η
段（`07`, `07b`, `08`, `09`, `09b`, `10`, `12`）は **TES・MKID 両系統とも `.wl`** です。レガシーの
`10_data` / `11_minimization` のみ `.nb` のまま残しています。

```
01 setup        単位・定数（自然単位系, GeV）
02 functions    補間ユーティリティ
03 functions    運動学・DMフォームファクター F_DM・ハロー速度積分 η
04 material      誘電関数（Al: Merminテーブル / TiN: 解析的Lindhard）
05 parameters   断面積・分解能・露光量
06 response     応答核の定義（本体）
07 dℛ/dω' プロット 核のプロット → output/<DET>/dcurlyRdEprime/*.pdf
07b dℛ/dω' csv  全域の核曲線 → output/<DET>/dcurlyRdEprime_csv/*.csv
08 応答関数      ビン分けした CRTES/CRTiN を .wdx 保存
09 プロット      応答関数 R_bin(v_min) のプロット
09b 応答 csv    全域の R_bin(v_min) → output/<DET>/response_functions_csv/*.csv
10 行列         v_min 積分した応答行列（matrix.csv, vmin.csv, bins.csv;
                行ごとに 1−alpha まで残す窓, TES q0=0.3 他 0.01）
12 eta          行列グリッド上の自然単位 η(v_min): Halo / Disk / Bound
   （10_data, 11_minimization はレガシー .nb のまま）
```

`<DET>` = `TES`(Al) または `MKID`(TiN)。`07`–`12` 段は `src/TES/`・`src/MKID/`
の両方に存在します。各段が何を準備しているかは **[Mathematica/README.md](Mathematica/README.md)** を参照。

### 実行方法

**Mathematica 段**（`wolframscript` が必要。macOS では `Wolfram.app` に同梱）:

```bash
cd Mathematica
# dℛ/dω' の図を生成（07 が 06→01 を芋づる式に読み込む）:
wolframscript -file src/TES/07_dcurlyRdEprime_plot.wl
# 01–06 がすべて定義されているかスモークテスト（TES・MKID）:
./src/tests/run_tests.sh
```

補足: `wolframscript` が `PATH` に無い場合は、フルパス
（`/Applications/Wolfram.app/Contents/MacOS/wolframscript`）で呼ぶか、
`/usr/local/bin` に一度シンボリックリンクを張ってください。

**Python 段**（パッケージ版）:

```bash
cd qsensor_analysis
pip install -e .            # `quantum_sensor` パッケージをインストール
jupyter notebook main_Refined.ipynb    # 解析を実行（TES/Al・MKID/TiN）
jupyter notebook final_figures.ipynb   # 論文用図: eta 回復・bin ごと signal・応答関数・dℛ/dE'
```

Mathematica 段が書き出す応答行列（Al は `Mathematica/output/TES/response_matrix/`、
TiN は `…/MKID/…`）を Python 段が読み込みます。`material="Al"`/`"TiN"` で切り替えます。
`final_figures.ipynb` は保存済み結果から、色弱（CVD）・プロジェクター対応の論文用図を
生成します（`qsensor_analysis/results/final/` 配下）。
