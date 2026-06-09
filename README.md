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
| [`Mathematica/`](Mathematica/) | Response (Wolfram) | The `01`–`11` pipeline that builds the response kernels, response functions and matrices for TES (Al) and MKID (TiN). See **[Mathematica/README.md](Mathematica/README.md)** for the detailed per-file description. |
| [`TES_clean2/`](TES_clean2/) | Analysis (Python) | The current analysis package (`quantum_sensor`): `DarkMatterQuantumAnalysis` runs the self-consistent forward+inverse model on the new response matrices. See **[TES_clean2/README.md](TES_clean2/README.md)**. |
| [`Obsolete/`](Obsolete/) | Archive | Earlier Python versions (`TES`, `TES_clean`, `TES_clean10`), kept for reference. |

### The Mathematica pipeline (summary)

The `01`–`11` files form a single dependency chain — loading one file `Get`s
all the earlier ones. `01`–`06` are clean hand-written `.wl` packages; the
plotting / matrix / fit stages (`07`–`11`) are being migrated from `.nb`.

```
01 setup        units & constants (natural units, GeV)
02 functions    interpolation utilities
03 functions    kinematics, DM form factor F_DM, halo speed integral η
04 material     dielectric function  (Al: Mermin table / TiN: analytic Lindhard)
05 parameters   cross sections, resolutions, exposures
06 response     response-kernel definitions (the core)
07 dℛ/dω' plots kernel plots → output/TES/dcurlyRdEprime/*.pdf   (✅ migrated to .wl)
08–11           response matrix, plots, data, minimisation         (.nb)
```

See **[Mathematica/README.md](Mathematica/README.md)** for what each stage prepares.

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
cd TES_clean2
pip install -e .            # installs the `quantum_sensor` package
jupyter notebook main_Refined.ipynb
```

The Mathematica stage writes response matrices
(`Mathematica/output/TES/response_matrix/`) that the Python stage consumes.

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
| [`Mathematica/`](Mathematica/) | 応答（Wolfram） | TES(Al)・MKID(TiN) の応答核・行列を構築する `01`–`11` パイプライン。各ファイルの詳細は **[Mathematica/README.md](Mathematica/README.md)** を参照。 |
| [`TES_clean2/`](TES_clean2/) | 解析（Python） | 現行の解析パッケージ（`quantum_sensor`）。`DarkMatterQuantumAnalysis` が新応答行列に対して自己整合な順方向＋逆問題を実行。詳細は **[TES_clean2/README.md](TES_clean2/README.md)**。 |
| [`Obsolete/`](Obsolete/) | アーカイブ | 旧Python版（`TES`, `TES_clean`, `TES_clean10`）。参照用に保管。 |

### Mathematica パイプライン（要約）

`01`–`11` は一直線の依存チェーンで、1つを `Get` すると前段がすべて読み込まれ
ます。`01`–`06` は手書きのクリーンな `.wl` パッケージ、プロット／行列／フィット
段（`07`–`11`）は `.nb` から移行中です。

```
01 setup        単位・定数（自然単位系, GeV）
02 functions    補間ユーティリティ
03 functions    運動学・DMフォームファクター F_DM・ハロー速度積分 η
04 material      誘電関数（Al: Merminテーブル / TiN: 解析的Lindhard）
05 parameters   断面積・分解能・露光量
06 response     応答核の定義（本体）
07 dℛ/dω' プロット 核のプロット → output/TES/dcurlyRdEprime/*.pdf   （✅ .wl 移行済み）
08–11           応答行列・プロット・データ・最小化                  （.nb）
```

各段が何を準備しているかは **[Mathematica/README.md](Mathematica/README.md)** を参照。

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
cd TES_clean2
pip install -e .            # `quantum_sensor` パッケージをインストール
jupyter notebook main_Refined.ipynb
```

Mathematica 段が書き出す応答行列（`Mathematica/output/TES/response_matrix/`）を Python 段が読み込みます。
