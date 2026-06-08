# QuantumSensor / Mathematica — 応答計算パイプライン (01–06)

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

## ファイル構成と依存関係

各ファイルは番号順に、末尾の番号のファイルが先頭で `Get` して前段を読み込む**直線的な依存チェーン**です。
`06_response_defs.wl` を 1 つロードすれば 01–05 が芋づる式に読み込まれます。

```
01_setup ─▶ 02_functions_math ─▶ 03_functions_response ─▶ 04_material ─▶ 05_parameters ─▶ 06_response_defs
 単位・定数      補間ユーティリティ      運動学・FDM・η            誘電関数・DOS          断面積・分解能・露光      応答核（本体）
```

- **01–06 は `.wl`**（手書きパッケージ、ビルドスクリプト `src/build_wl.py` で TES/MKID 両方を生成）
- **07–11 は `.nb`**（プロット・応答行列・データ・フィット。本 README の対象外）

---

## 各段階が「何を準備しているか」

### 01_setup — 単位系と物理定数
- **自然単位系**を採用し、すべての量を **GeV 建て**で表現します（`GeV = 10^9`）。
  `eV`, `keV`, `cm`, `sec`, `Kel`（ケルビン）, `grams`, `kg` などは GeV 換算の係数として定義。
- **DM ハロー速度**：`v0`（=220 km/s）, `ve`（地球速度の年平均）, `vesc`（脱出速度=544 km/s）。
- **物理定数・DM 密度**：`alpha`（微細構造定数）, `\[Rho]DM`（局所 DM 密度 0.4 GeV/cm³）, `me`（電子質量）。
- 物質密度 `rhoAl`, `rhoTiN`。

> ここで決めた単位換算が以降すべての数式に一貫して効きます。「数値の桁」を読むときは必ずこの単位定義に立ち返ってください。

### 02_functions_math — 補間ユーティリティ
- `LIntpl1D[x, y]` — 並行ベクトルから 1 次（区分線形）の `InterpolatingFunction` を構築。
- `IntplArray[m]` — 共通 x グリッド（`m[[1]]`）に対する複数 y 列を一括補間。04 で `.dat` を関数化するのに使用。
- `Inword` — 図注釈用のテキスト inset（プロット系で使用）。

### 03_functions_response — 運動学・形状因子・速度積分
DM–電子散乱の**運動学**と、レート計算に必要な解析的部品を定義します。
- **運動量移行 q**：エネルギー移行 `Ee` に対する 2 つの分岐 `ql`（左）/ `qr`（右）と、そのヤコビアン `Jacobl/Jacobr`。
- `vmin` — 与えられた `(Ee, q)` に対する DM 最小速度。
- **DM 形状因子 `FDM[q][n] = (α·me / q)^n`** — 媒介子の性質を表す（`n=0` 重い媒介子、`n=2` 軽い媒介子）。`FmedH/FmedL` は別パラメータ化。
- **ハロー速度積分 η**（`\[Eta]th`, `\[Eta]td`）— 切り詰めマクスウェル分布を `vmin` 以上で積分した `η(vmin)`。規格化は `KKf`。

### 04_material — 物質の誘電応答（**TES と MKID で唯一大きく異なる段階**）
散乱レートの「物質側」の入力 `Im[-1/ε(ω, q)]` を準備します。

- **TES（Al, Mermin）**：`input/TES/Al_mermin.dat` から ε₁, ε₂ を読み込み 2 次元補間し、
  `ImepsAlf[w, q]`, `eps1Alf`, `eps2Alf` を構築（定義域 `0≤w≤99.3`, `0≤q≤37289.5` 外は 0）。
  さらにフォノン状態密度 `Al_pDoS.dat → Dw[w]`、フォームファクター表 `Al_Fn.dat → Fnw10GeV1..5` も準備。
- **MKID（TiN, Lindhard）**：データ不要。フェルミ波数 `kFTiN`、プラズマ振動数 `wpTiN`、減衰 `GammaTiN` から
  **解析的 Lindhard 誘電関数 `ImepsLTiN[w, q]`** を直接定義。

> TES は実測ベースの数値テーブル、MKID は解析式という違い。`input/MKID` ディレクトリは存在しませんが、
> TiN は解析式なのでロードに影響しません。

### 05_parameters — 計算パラメータ
- **参照断面積** `\[Sigma]e`（電子 10⁻³⁰ cm²）, `\[Sigma]N`。
- **エネルギー分解能** `TESsig`, `MKIDsig`（FWHM をガウス σ に換算した相対分解能）。
- **露光量**（有効質量 × 時間）`TiNexp`, `Alexp`。

### 06_response_defs — 応答核（パイプラインの本体）
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

## 計算の流れ（物理的に何が起きているか）

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

## ロード方法

```wolfram
(* 06 を Get すれば 01–05 も連鎖ロードされる *)
Get["src/TES/06_response_defs.wl"];      (* または src/MKID/... *)

(* 例：質量 mχ=1 GeV, 媒介子 n=1 での速度積分応答 *)
CRintTES[1 GeV, 1][E1, E2][vmin, vmax, ns]
```

---

## 動作確認（スモークテスト）

01–06 がすべて定義され、数値が健全に返ることを検証するテストを同梱しています。

```bash
cd QuantumSensor/Mathematica
./src/tests/run_tests.sh          # TES と MKID 両方
./src/tests/run_tests.sh TES      # 片方だけ
```

- ロード成否・主要シンボルの定義有無・代表引数での数値スモークの 3 層を検査。
- レポートは `output/<DET>/test/report.txt`（人間用）と `report.m`（機械可読）に出力。
- 現状：**TES 20/20・MKID 18/18 パス**。
```
