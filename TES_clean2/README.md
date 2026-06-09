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

### 1-B. 速度分布 `eta(v_min)`（ハローモデル、**自然単位系**）

- **実体**: 1 列 CSV。長さは行列の列数 `n_vmin` と一致し、各値は **自然単位系**の `eta(v_min)`（SHM 速度積分 `\[Eta]th`、`ρDM·σe/mχ` のプレファクタ込み）。**行列と同じ自然単位系・同じ `v_min` グリッド**で作られているので、`M_phys @ eta` がそのまま真の期待イベント数になる（補間も単位換算も不要）。
- **出所**: Mathematica の **stage 12** が生成。各行列フォルダの中に書き出される。
  ```
  ../Mathematica/output/TES/response_matrix/M<mass>/<...>/eta_<model>.csv
  例: .../M1/Al_q0M1_R5_v1-800_N1000/eta_Halo.csv
  ```
  - 生成コマンド（**Python 解析の前に 1 回実行が必要**）:
    ```
    wolframscript -file Mathematica/src/TES/12_eta.wl ALL Halo
    ```
    [`12_eta.wl`](../Mathematica/src/TES/12_eta.wl) は各フォルダの `vmin.csv` を読み、その `v_mid` で `\[Eta]th[dmMass][v_mid·kps][v0,ve,vesc]` を評価して `eta_<model>.csv` を出力する（パラメータは `01_setup.wl`、`σe` は `05_parameters.wl`、`\[Eta]th` は `03_functions_response.wl`）。
  - **スコープ**: 現状 **Halo (SHM) のみ**。`.wl` パイプラインには `\[Eta]th ≡ \[Eta]td` で `Bound` の定義が無く、Disk/Bound のレガシー CSV は旧パイプライン由来のため、ここでは再生成しない。
- **読み込むコード**: [`data_loader.py`](src/quantum_sensor/data_loader.py) の `load_natural_eta(rm, model)`。長さが `n_vmin` と合わない／ファイルが無い場合は、上記コマンドを促す明示的なエラーを出す。

> レガシーの物理単位 eta（`data/Eta_data/Eta<eta>M<mass>_Ko.csv`、≈ cm⁻¹、800 サンプル）は `load_eta()` / `model.align_eta()` で読めるが、**自然単位の行列と不整合なのでクロスチェック専用**。主経路では使わない。

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
| `self.eta` | **自然単位の `eta`** を行列フォルダの `eta_<model>.csv` から読む。行列の `v_min` グリッド上で生成済みなので補間・単位換算は無し | `data_loader.load_natural_eta` |
| `self.m_phys` | **物理順方向演算子** `M_phys = AL_EXP × matrix`（露光を掛けて期待カウントスケールにする。自然単位） | `model.response_operator`（露光は `model.exposure_factor` 経由で `constants.AL_EXP`） |
| `self.signal` | 信号カウント `signal = M_phys @ eta`（エネルギービンごとの**真の期待イベント数**。`M_phys` も `eta` も自然単位なので単位整合） | `analysis.py`（`m_phys @ eta`） |
| `self.background` | 背景シナリオを `bins.csv` のエネルギー端で各ビンに積分。`A·exp(-E/B)+C` を区間積分し、単位係数 `AMP_SCALE = (eV/keV)·365` を掛ける（旧コードと同じ「/keV/day」由来）。`none`/`a` はゼロ | `backgrounds.background_counts`（→ `_integrate`） |
| `self.observed` | **観測カウント** `observed = signal + background`（逆問題が再現すべきターゲット） | `analysis.py` |

> この時点で `self.result = None`, `self.flux = None`（逆問題は未実行）。

### フェーズ② — `optimize(solver="osqp", ...)`: 明示的に呼んだとき計算

1. **前処理コンディショニング** — `model.condition(m_phys, observed, background)`
   コンディショニングは2段構え（neutrinoAnalysis 方式）。**データ・背景には一切触れない**ので、最小化される目的関数は**真の Neyman χ²**のまま（`res.fun` が真の χ² → Δχ²=1 等の**統計誤差議論がそのまま成立**。旧 `cons2` のようにデータを掛けると χ² が定数倍され破綻する）。
   - **単一の共通定数 `c`（`config.CONDITION_C`）**: 未知数を `x = c·u` と置き換える全体スケール。`data_cond = data`（そのまま）、`unscale(u) = c·u`。
   - **列正規化（本命のコンディショナ）**: ソルバ内で各列を単位ノルムに正規化する（`optimizer._column_scale`、`D_j = 1/‖M列j‖`）。**`c` は列ノルムで割り消される**ので、復元フラックスは実用域（`c≈1e-31〜1e-25`）で `c` に依存しない（旧 `cons1`/`cons2`・質量別表を完全に置換、手調整不要）。

2. **逆問題の求解（2段階: 列正規化 QP → 裾重み付き頂点選択）** — `solver ∈ {"osqp","qp","clarabel"}` または `fix` 指定なら `optimizer.run_optimize_qp`、それ以外は `trust-constr`（`optimizer.run_optimize`）。
   - χ²（Neyman）`Σ (data − Bkg − M·x)² / data`（`data>0` のビンのみ）を **QP** として組む（`optimizer._build_qp`）。制約: **単調非増加** `x[i] ≥ x[i+1]`＋**非負** `x ≥ 0`。
   - **ルーティング**（`run_optimize_qp`）:
     - `fix` 指定あり → CLARABEL（固定パラメータ解、`_CLARABELBackend`）
     - 自由解（既定）→ `_OSQPBackend`:
       1. **QP を列正規化変数 `x = D⊙z` で OSQP 求解** → 滑らかな内部解（ランプ）と当てはめ値 `μ = M·x` を得る。
       2. `vertex_select=True`（既定）→ **同じ `μ` を再現する頂点（階段状フラックス）を HiGHS シンプレックスで選ぶ**（`_vertex_select`）。目的は **`Σ R^(i/(n-1))·xᵢ` の最小化（幾何的な裾重み、`R=1000`）**で、高 `v_min`（高エネルギー側）を 0 に駆動（`eta→0` と整合）。χ² は `μ` 不変なので変わらない。
3. **物理単位へ復元** — `self.flux = unscale(res.x)`。`self.result` に求解メタ（`fun`=真の χ², `backend`, `nit`, `staircase` 等）を保持。

> **方式は kyphys-creator/neutrinoAnalysis に準拠**: ①列正規化 QP で `μ=M·x` を確定 → ②裾重み付きシンプレックスで階段頂点を選ぶ。χ² 最小は劣決定（χ²=0 が面）なので、OSQP の内部解（滑らかなランプ）ではなく**物理的に意味のある頂点（区分定数フラックス）**を選ぶのが要点。
>
> **コンディショニングは列正規化が本命**で、単一定数 `CONDITION_C` は列ノルムで割り消されるため実用域で結果に効かない（質量別の手調整は不要）。**データを掛けないので真の χ² が保存**され、将来の Δχ² 統計誤差解析にそのまま使える。
>
> **裾の挙動**: 裾重み目的は高 `v_min` を 0 に落とす。応答が高 `v_min` で縮退する 10 MeV/M1 では `eta` が残る帯（~120–565 km/s）も 0 に切られ低速側が過大評価される（rms/eta₀≈0.28）。M2/M3 は `eta` の肩まで追従して裾で 0 に落ちる。裾を残したい場合は別の頂点目的（例: 最上段 `x₀` 最小）に差し替え可能。

---

## 3. 最終的に得られる量 — 何が出てくるか / どこで作るか

| 出力 | 中身 | 作るコード |
|---|---|---|
| **`flux`（戻り値）** | 物理単位の復元フラックス `x(v_min)`、長さ `n_vmin`。`a.vmin_mid` 上の値。`optimize()` の戻り値かつ `a.flux` | `analysis.optimize` |
| **フラックス CSV** | `vmin_low, vmin_high, vmin_mid, flux` の表。`flux` は**自然単位のまま**。`results/flux_<material>_q<q>_M<mass>_R<nbins>_<eta>_bkg-<background>.csv` | `plotting.save_flux`（`a.save_flux()`） |
| **比較図 PDF** | 入力 `eta(v_min)`（赤線）と復元フラックス（階段）を重ねた図。`results/scenario_bkg_<background>/flux_<stem>.pdf`（x 軸 log）。**この図を描くときだけ自然単位 → 物理単位 `cm⁻¹` に変換**（`× CM`、`plotting.ETA_TO_CM_INV`） | `plotting.plot_flux_comparison`（`a.plot()`） |

これら以外に、`a.result`（求解の診断）、`a.signal` / `a.observed`（順方向カウント）、`a.eta`（整列済み入力）もオブジェクト属性として参照できる。

---

## 4. データフロー全体図

```
[Mathematica 08→10]   [Mathematica 12]          [constants.py]   [config.py]
 response_matrix/M*/<name>/                       AL_EXP, ...      BACKGROUND_SCENARIOS
  ├ matrix.csv (n_ebins × n_vmin)                     │                 │
  ├ vmin.csv  {lo,hi,mid}                             │                 │
  ├ bins.csv  {E_lo,E_hi}                             │                 │
  └ eta_<model>.csv  (自然単位, vmin グリッド上)        │                 │
        │ load_response_matrix   │ load_natural_eta    │                 │
        ▼                        ▼                     ▼                 ▼
   self.rm (ResponseMatrix)   self.eta (自然単位)                        │
        │                        │                                       │
        │ response_operator: m_phys = AL_EXP × matrix                    │
        ▼                        │                                       │
   self.signal = m_phys @ eta ──►(+)◄──── background_counts(bins) = self.background
                                               │
                                               ▼
                                     self.observed  (順方向の観測カウント)
                                               │  optimize()
                                               ▼
                        condition(): m_cond=c·M / data そのまま / unscale=c·u   ← 単一定数 CONDITION_C
                                               │
                                               ▼
                        run_optimize_qp: 真の Neyman χ² 最小 (単調非増加・非負)
                          ①列正規化 QP(OSQP) → μ=M·x  ②裾重み付きシンプレックス頂点 → 階段解
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
| [`data_loader.py`](src/quantum_sensor/data_loader.py) | 応答行列＋自然単位 eta（`eta_<model>.csv`）の読み込み、`ResponseMatrix`。レガシー eta/ratebin も保持 |
| [`model.py`](src/quantum_sensor/model.py) | 順方向演算子 `m_phys`、データ由来コンディショニング（`align_eta` はレガシー専用） |
| [`backgrounds.py`](src/quantum_sensor/backgrounds.py) | 行列のエネルギービン上での背景カウント |
| [`optimizer.py`](src/quantum_sensor/optimizer.py) | 単調非負 χ²/QP ソルバ（OSQP / CLARABEL / HiGHS 頂点 / trust-constr） |
| [`analysis.py`](src/quantum_sensor/analysis.py) | `DarkMatterQuantumAnalysis(config)` — エントリポイント・司令塔 |
| [`plotting.py`](src/quantum_sensor/plotting.py) | eta vs 復元フラックス図 ＋ CSV 出力 |

---

## 6. 設計方針

- 形・エネルギービン・`v_min` グリッドはすべて行列ファイル（`matrix.csv` / `vmin.csv` / `bins.csv`）から読む。コードには書かない。
- 物理定数は [`constants.py`](src/quantum_sensor/constants.py) に一度だけ（Mathematica 由来）。
- 露光スケール `K`(=`AL_EXP`) は**導出値**。復元フラックスの**形・値は `K` にもコンディショニングにも依存しない**（`unscale` で打ち消される）。
- **コンディショニングは列正規化（`optimizer._column_scale`、neutrinoAnalysis 方式）が本命**。各列を単位ノルム化するため、単一定数 `config.CONDITION_C` は割り消され実用域で結果に効かない（旧 `cons1`/`cons2`・質量別表を統合・実質不要化）。**データを掛けないので真の Neyman χ² が保存**され、後の Δχ² 統計誤差解析がそのまま成立する。
- **求解は2段階**（`optimizer._OSQPBackend`）: 列正規化 QP で `μ=M·x` を確定 → 裾重み付きシンプレックスで階段頂点を選択。χ² は `μ` 不変なので頂点選択で変わらない。
- **単位系は一貫して自然単位系**。`M_phys`・`AL_EXP`・`eta` がすべて自然単位なので `signal = M_phys @ eta` は真の期待イベント数になり、背景カウントとも整合的に比較できる。`eta` は Mathematica（[`12_eta.wl`](../Mathematica/src/TES/12_eta.wl)）で行列と同じグリッド・同じ単位で生成し、Python は補間も単位換算もせずそのまま使う（単一の真実の源）。
- **物理単位への変換は描画時のみ**。内部・保存 CSV は自然単位を保ち、`plot()` で図を描くときだけ `eta`/`flux`（次元 1/length）を `× CM`（`plotting.ETA_TO_CM_INV`）で `cm⁻¹` に直して表示する。

---

## 7. 使い方

**前提**: 解析に使う行列フォルダに自然単位 eta（`eta_Halo.csv`）が必要。最初に一度だけ Mathematica で生成する:

```
wolframscript -file ../Mathematica/src/TES/12_eta.wl ALL Halo
```

その後:

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

> **スコープ**: TES（Al）。`eta` は当面 **Halo (SHM) のみ**自然単位で生成（[`12_eta.wl`](../Mathematica/src/TES/12_eta.wl)）。Disk/Bound は現行 `.wl` に定義が無いため未対応で、レガシー物理単位 CSV はクロスチェック用に残す。レガシーの `data/`（eta, ratebin）と旧 `curlyRplotting` 行列も同様。R10 / q2 も動く。観測カウントはレガシーの 5 ビン Ratebin ではなく、自己無撞着な順方向モデルから生成する。
