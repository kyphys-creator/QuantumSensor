# quantum_sensor (TES_clean2)

ダークマター量子センサー（TES / Al）の **応答 (response) 解析パッケージ**。
Mathematica 側で作った応答行列 `M` と、ハローモデルの速度分布 `eta(v_min)` を入力に、

1. **順方向 (forward)**: 観測されるエネルギービンごとのイベント数を `観測 = 露光 × M @ eta (+ 背景)` で生成し、
2. **逆方向 (inverse)**: その観測カウントを再現する **単調・非負のフラックス `x(v_min)`** を復元する、

という **自己無撞着 (self-consistent)** な forward + inverse 解析を行う。観測カウントは行列と `eta` から自分で作るので、逆問題は必ずその `eta` を再現する解を持つ。

このドキュメントは「**入力は何か / どこから来るか / 途中で何がいつどこで計算されるか / 最終的に何が得られるか**」を余すことなく説明することを目的とする。

---

## 0. 用語と設定パラメータ

1回の解析は `RunConfig`（[`src/quantum_sensor/config.py`](src/quantum_sensor/config.py)）1個で完全に決まる。

| フィールド | 意味 | 取りうる値 |
|---|---|---|
| `material` | 検出器素材 | `"Al"`（TES） |
| `q` | 媒介子（FDM）タグ | `"0"` = 重い媒介子 (n=0) / `"2"` = 軽い媒介子 (n=2) |
| `mass` | DM 質量タグ | `"1"` = 10 MeV / `"2"` = 100 MeV / `"3"` = 1 GeV |
| `nbins` | 観測エネルギービン数 | `5`（0.2 eV 幅）/ `10`（0.1 eV 幅） |
| `eta` | ハロー速度分布モデル | `"Halo"`（SHM）/ `"Disk"` / `"Bound"` |
| `background` | 背景シナリオ名 | `"none"`,`"a"`,`"c"`,`"b"`,`"b2"`,`"flat"` |
| `run` | 行列フォルダを一意に選ぶための部分文字列（同一設定で複数 run があるとき） | 省略可 |

---

## 1. 入力 — 何で、どこから持ってくるか

入力は大きく **4種類**ある。それぞれ「実体」「出所」「読み込むコード」を示す。

### 1-A. 応答行列 `M` とその軸グリッド（**主入力**）

- **実体**: 1つの設定につき 3 つの CSV が入ったフォルダ。
  - `matrix.csv` — 純粋な数値行列。形は **(n_ebins, n_vmin)**。行 = 観測エネルギービン、列 = `v_min` 区間。各要素は「その区間にわたる応答の `v_min` 積分」（自然単位、露光は未適用）。
  - `vmin.csv` — 各列に対応する `{v_low, v_high, v_mid}`（km/s）。
  - `bins.csv` — 各行に対応するエネルギービン端 `{E_low, E_high}`（eV）。
- **出所**: Mathematica パイプラインの出力。
  ```
  ../Mathematica/output/TES/response_matrix/M<mass>/<material>_q<q>M<mass>_R<nbins>_v<lo>-<hi>_N<N>/
  例: Mathematica/output/TES/response_matrix/M1/Al_q0M1_R5_v1-800_N1000/
  ```
  - これは Mathematica の **stage 08 → stage 10** で生成される（[`Mathematica/src/TES/08_response_functions.wl`](../Mathematica/src/TES/08_response_functions.wl), [`10_response_matrix.wl`](../Mathematica/src/TES/10_response_matrix.wl)）。
    - **08**: 各エネルギービン `[E1,E2]` について応答 `CRTES[mass,n][E1,E2](v_min)` を `v_min` の関数（補間関数）として評価し `.wdx` に保存。
    - **10**: その `.wdx` を読み込み、`v_min` を `[lo,hi]` で N 等分し、各区間で応答を台形則で**厳密積分** → `M[i,j]`。`km/s` を自然単位に直すため `kps` を掛ける。先頭の全ゼロ列（最低エネルギービンが閾値を越える前の列）は削られ、`vmin.csv` も同じ数だけ頭が削られて列と整合する。
  - **重要**: 形・エネルギービン・`v_min` グリッドはすべてこれらのファイルから読む。**Python 側にハードコードは一切ない**。
- **読み込むコード**: [`data_loader.py`](src/quantum_sensor/data_loader.py) の `find_matrix_dir()` → `load_response_matrix()`。結果は `ResponseMatrix` データクラス（`matrix`, `vmin_low/high/mid`, `ebin_low/high`, `name`, `path`）。

### 1-B. 速度分布 `eta(v_min)`（ハローモデル）

- **実体**: 800 サンプルの 1 列 CSV（`v_min ∈ [1, 800] km/s` を等間隔にサンプリングした `eta` の値）。
- **出所**: パッケージ同梱のレガシーデータ。
  ```
  data/Eta_data/Eta<eta>M<mass>_Ko.csv
  例: data/Eta_data/EtaHaloM1_Ko.csv
  ```
- **読み込むコード**: [`data_loader.py`](src/quantum_sensor/data_loader.py) の `load_eta(eta, mass)`。

### 1-C. 物理定数

- **実体**: 自然単位（GeV 基準）の物理定数群（単位換算、`alpha`、DM 密度 `RHO_DM`、Al 密度、参照断面積、エネルギー分解能、そして **Al 露光 `AL_EXP = 8200 ng·month`**、質量タグ→DM 質量 `DM_MASS` など）。
- **出所**: [`constants.py`](src/quantum_sensor/constants.py)。**Python 側の唯一の真実の源**で、各値は Mathematica の [`01_setup.wl`](../Mathematica/src/TES/01_setup.wl) / [`05_parameters.wl`](../Mathematica/src/TES/05_parameters.wl) から一字一句写してあり、コメントに出所が書いてある。手調整値（マジックナンバー）は無い。

### 1-D. 背景シナリオのパラメータ

- **実体**: 背景スペクトル `R(E) = A·exp(-E/B) + C` の係数 `(A, B, C)` を持つシナリオ表。`"none"`/`"a"` は背景ゼロ（信号のみ）。
- **出所**: [`config.py`](src/quantum_sensor/config.py) の `BACKGROUND_SCENARIOS`（背景の数値が存在する唯一の場所）。

### （参考）クロスチェック用レガシーデータ — 通常の流れでは未使用

`data_loader.py` には次の読み込み関数もあるが、**標準の解析経路では呼ばれない**（単位・露光較正や過去結果との突き合わせ用）:
- `load_ratebin()` — 旧 5 ビン観測カウント `data/Ratebin/Event5...csv`
- `load_legacy_response_matrix()` — 旧露光込み行列 `data/curlyRplotting/...csv`

---

## 2. 途中で計算される量 — 何が、いつ、どこで

計算は **2つのタイミング**に分かれる: ① `DarkMatterQuantumAnalysis(config)` の **コンストラクタ**（順方向モデルの構築）と、② `optimize()` の **呼び出し時**（逆問題の求解）。すべて [`analysis.py`](src/quantum_sensor/analysis.py) が司令塔で、実際の計算は各モジュールへ委譲する。

### フェーズ① — `__init__`（順方向モデル）: オブジェクト生成時に即計算

| 生成される量 | 計算内容 | 計算コード（どこで） |
|---|---|---|
| `self.rm` | 応答行列とグリッドを読み込む | `data_loader.load_response_matrix` |
| `self.eta` | レガシー `eta`（`[1,800] km/s` の 800 点）を**行列の `vmin_mid` グリッドへ線形補間で再標本化**。範囲外は端値に外挿（高速側で `eta→0`） | `model.align_eta` |
| `self.m_phys` | **物理順方向演算子** `M_phys = AL_EXP × matrix`（露光を掛けて生カウントスケールにする） | `model.response_operator`（露光は `model.exposure_factor` 経由で `constants.AL_EXP`） |
| `self.signal` | 信号カウント `signal = M_phys @ eta`（エネルギービンごとの期待イベント数） | `analysis.py`（`m_phys @ eta`） |
| `self.background` | 背景シナリオを `bins.csv` のエネルギー端で各ビンに積分。`A·exp(-E/B)+C` を区間積分し、単位係数 `AMP_SCALE = (eV/keV)·365` を掛ける（旧コードと同じ「/keV/day」由来）。`none`/`a` はゼロ | `backgrounds.background_counts`（→ `_integrate`） |
| `self.observed` | **観測カウント** `observed = signal + background`（逆問題が再現すべきターゲット） | `analysis.py` |

> この時点で `self.result = None`, `self.flux = None`（逆問題は未実行）。

### フェーズ② — `optimize(solver="osqp", ...)`: 明示的に呼んだとき計算

1. **前処理コンディショニング** — `model.condition(m_phys, observed, background)`
   未知数が `O(1)` になるよう、**データから導出した**2つのスケールで問題を無次元化する（旧 `cons1`/`cons2` の置き換え。手調整なし）:
   - `data_scale` = 正のカウントの中央値（各方程式を割って行を `O(1)` に）
   - `flux_scale` = `data_scale / median(M_phys @ 1)`（一様単位ベクトルがデータを出すのに必要なフラックス量 → `x` の自然な大きさ）

   返り値: スケール済み `M_cond`, `data_cond`, `bkg_cond` と、物理単位へ戻す関数 `unscale(u) = flux_scale·u`。
   *フラックスの「形」は `K`(=露光)・スケール因子に依存しない*（`argmin_x |KM·eta − KM·x|` で消えるため）。

2. **逆問題の求解** — `solver ∈ {"osqp","qp","clarabel"}` または `fix` 指定なら `optimizer.run_optimize_qp`、それ以外は `trust-constr`（`optimizer.run_optimize`）。
   - χ²（Neyman）`Σ (data − Bkg − M·x)² / data`（`data>0` のビンのみ）を **QP** として組む（`optimizer._build_qp`）。
   - 制約: **単調非減少**でなく **非増加** `x[i] ≥ x[i+1]`（`ordering_constraint`）＋ **非負** `x ≥ 0`。
   - **ルーティング**（`run_optimize_qp`）:
     - `fix` 指定あり → CLARABEL（固定パラメータ解、`_CLARABELBackend`）
     - 自由解 → OSQP（`_OSQPBackend`）。**劣決定（`n_vmin > 有効ビン数`）**のとき χ²_min=0 の最適面が広がるので、`vertex_select=True`（既定）で **HiGHS 双対シンプレックス LP** に切り替え、`Σx` を最小化する **基底頂点 = 階段状（piecewise-constant）フラックス**を選ぶ（`_staircase_vertex`）。
3. **物理単位へ復元** — `self.flux = unscale(res.x)`。`self.result` に求解メタ（`fun`=χ², `backend`, `nit`, `solve_time`, `staircase` 等）を保持。

> 逆問題は本質的に劣決定（少数のエネルギービン × 多数の `v_min` ビン）。`vertex_select=False` にすると最小ノルム点（入力 `eta` の滑らかな形に近いが綺麗な階段ではない）が返る。いずれにせよ高速側の裾はわずかなエネルギービンでしか拘束されず弱い。

---

## 3. 最終的に得られる量 — 何が出てくるか / どこで作るか

| 出力 | 中身 | 作るコード |
|---|---|---|
| **`flux`（戻り値）** | 物理単位の復元フラックス `x(v_min)`、長さ `n_vmin`。`a.vmin_mid` 上の値。`optimize()` の戻り値かつ `a.flux` | `analysis.optimize` |
| **フラックス CSV** | `vmin_low, vmin_high, vmin_mid, flux` の表。`results/flux_<material>_q<q>_M<mass>_R<nbins>_<eta>_bkg-<background>.csv` | `plotting.save_flux`（`a.save_flux()`） |
| **比較図 PDF** | 入力 `eta(v_min)`（赤線）と復元フラックス（階段）を重ねた図。`results/scenario_bkg_<background>/flux_<stem>.pdf`（x 軸 log, 単位 `cm⁻¹`） | `plotting.plot_flux_comparison`（`a.plot()`） |

これら以外に、`a.result`（求解の診断）、`a.signal` / `a.observed`（順方向カウント）、`a.eta`（整列済み入力）もオブジェクト属性として参照できる。

---

## 4. データフロー全体図

```
[Mathematica 08→10]                         [data/Eta_data]          [constants.py]   [config.py]
 response_matrix/.../                         Eta<eta>M<mass>.csv      AL_EXP, ...      BACKGROUND_SCENARIOS
  ├ matrix.csv (n_ebins × n_vmin)                  │                      │                 │
  ├ vmin.csv  {lo,hi,mid}                          │                      │                 │
  └ bins.csv  {E_lo,E_hi}                          │                      │                 │
        │ load_response_matrix                     │ load_eta             │                 │
        ▼                                          ▼                      ▼                 ▼
   self.rm (ResponseMatrix) ───────────────► align_eta(→vmin_mid) = self.eta                │
        │                                                                                    │
        │ response_operator: m_phys = AL_EXP × matrix                                        │
        ▼                                                                                    │
   self.signal = m_phys @ eta ──────────────►(+)◄──── background_counts(bins) = self.background
                                               │
                                               ▼
                                     self.observed  (順方向の観測カウント)
                                               │  optimize()
                                               ▼
                        condition(): m_cond / data_cond / bkg_cond / unscale   ← データ由来スケール
                                               │
                                               ▼
                        run_optimize_qp: χ² 最小 (単調非増加・非負)
                          劣決定 → HiGHS 頂点 = 階段フラックス
                                               │ unscale
                                               ▼
                                  self.flux  ──►  save_flux (CSV) / plot (PDF)
```

---

## 5. モジュール一覧（責務）

| モジュール | 役割 |
|---|---|
| [`config.py`](src/quantum_sensor/config.py) | `RunConfig`, `BackgroundModel`, `BACKGROUND_SCENARIOS`（実行パラメータと背景数値の唯一の置き場） |
| [`constants.py`](src/quantum_sensor/constants.py) | 物理定数（Mathematica 01/05 から写した単一の真実の源） |
| [`data_loader.py`](src/quantum_sensor/data_loader.py) | 新応答行列フォルダ＋レガシー eta/ratebin の読み込み、`ResponseMatrix` |
| [`model.py`](src/quantum_sensor/model.py) | `eta` 整列、順方向演算子 `m_phys`、データ由来コンディショニング |
| [`backgrounds.py`](src/quantum_sensor/backgrounds.py) | 行列のエネルギービン上での背景カウント |
| [`optimizer.py`](src/quantum_sensor/optimizer.py) | 単調非負 χ²/QP ソルバ（OSQP / CLARABEL / HiGHS 頂点 / trust-constr） |
| [`analysis.py`](src/quantum_sensor/analysis.py) | `DarkMatterQuantumAnalysis(config)` — エントリポイント・司令塔 |
| [`plotting.py`](src/quantum_sensor/plotting.py) | eta vs 復元フラックス図 ＋ CSV 出力 |

---

## 6. 設計方針: マジックナンバーを置かない

- 形・エネルギービン・`v_min` グリッドはすべて行列ファイル（`matrix.csv` / `vmin.csv` / `bins.csv`）から読む。コードには書かない。
- 物理定数は [`constants.py`](src/quantum_sensor/constants.py) に一度だけ（Mathematica 由来）。
- 露光スケール `K`(=`AL_EXP`) は**導出値**であってチューニング値ではない。ソルバのコンディショニング（旧 `cons1`/`cons2`）は `model.condition` で**データから計算**する。復元フラックスの**形は `K` に依存しない**。

---

## 7. 使い方

```python
from quantum_sensor import DarkMatterQuantumAnalysis, RunConfig

a = DarkMatterQuantumAnalysis(RunConfig(material="Al", q="0", mass="3",
                                        nbins=5, eta="Halo", background="none"))
flux = a.optimize(solver="osqp")   # a.vmin_mid 上の復元フラックス
a.plot()                            # 比較図 PDF を保存
a.save_flux()                       # フラックス＋v_min グリッドを CSV 保存
```

最小のエンドツーエンド例は [`examples/run_example.py`](examples/run_example.py)（Al・重媒介子 q0・5 ビン・Halo・信号のみで M1/M2/M3 を回す）。

```
python examples/run_example.py
```

> **スコープ**: TES（Al）。レガシーの `data/`（eta, ratebin）と旧 `curlyRplotting` 行列はクロスチェック用に残してある。R10 / q2 も動く。観測カウントはレガシーの 5 ビン Ratebin ではなく、自己無撞着な順方向モデルから生成する。
