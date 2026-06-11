# TODO / 今後やりたいこと

(2026-06-10 時点。会話ベースのメモ — 完了したら消すか ✅ を付ける)

## 統計まわり

- [ ] **背景ペナルティ項(bkg penalty)の実装** — neutrinoAnalysis
  (`1eV/neutrino_analysis_band.py` §3)の移植。背景の規格化が測定誤差を持つ
  状況を扱う層で、背景ありシナリオの band / 有意度を論文の数字にする前に必要。
  - 疑似実験ごとに「測定された背景」を Beta 分布でサンプル
    (`B_varied = f·O^MC`, `f ~ Beta(α,β)`、`0 < B_varied < O^MC` を保証)
  - χ² にペナルティ `Σ (B_fit − B_varied)²/B_varied` を追加、`B_fit ≥ 0` を
    ビンごとのナイサンスとしてフィット(QP のまま解ける形に
    `optimizer._build_qp` を拡張)
  - エッジケース処理(μ≥1 クリップ、σ² 上限、ゼロ背景ビン除外)+統計カウンタ
    +サンプリング健全性チェック図
  - 入口: `RunConfig` か `find_confidence_band` 引数に `bkg_penalty` フラグ
  - 現状の `bkg-none` の結果には影響なし。背景あり解析を始めるときに着手。
- [ ] **背景ありシナリオの本解析** — `results/` は現状ほぼ `bkg-none` のみ。
  `c` / `b` / `b2` / `flat` で point-wise band と判別有意度を回す。その際
  `backgrounds.py` の注記(レガシー由来の背景規格化が新パイプラインの信号と
  整合か要再較正)を先に確認する。
- [ ] **point-wise profile band の拡張** — いま走っているのは fine-binning ×
  Halo × bkg-none の 12 設定。mixture(mix5/mix25)・Bound・粗ビン(R5)への
  拡張、走査点数(`n_indices`)や `n_pseudo_edge` の収束チェック。
- [ ] (任意)toy アンサンブルバンド(`fit_toys`)を point-wise band の
  大域的クロスチェックとして全設定で揃える(3/12 で中断したまま)。

## 物理・パイプライン

- [ ] **応答行列の窓パラメータ α の正当化** — 現状「TES q0 のみ α=0.3、
  他は 0.01」という非対称設定。復元可能な v_min 範囲に直接効くので、
  論文に書ける物理的根拠を整理する(または統一する)。
- [ ] **Bound モデルの扱いの明文化** — M3 のみ非ゼロ(v_esc=11.2 km/s の
  窓に届くのが最重質量だけ)である旨と、Bound = Halo+Bound(上乗せ)の定義。
- [ ] レガシー `.nb`(`10_data`, `11_minimization`)の扱いを決める
  (`.wl` 化するか、obsolete 送りにするか)。

## 論文戦略

- [ ] **第1論文**: 方法+応答関数ライブラリ(ELF 形式、Al-TES / TiN-MKID)+
  最小限の定量主張として **Asimov 判別有意度グリッド**
  (`results/stats/discrimination.csv`、計算済み)を載せる。
  結論部で第2論文の範囲を予告してサラミ批判を予防。
- [ ] **第2論文**: フル統計 — point-wise profile band、背景ペナルティ込みの
  band、背景シナリオ別の判別力マップ。
- [ ] **既存研究との差別化の明示** — 電子散乱のハロー非依存解析は既存
  (Chen–Gelmini–Takhistov arXiv:2105.08101, arXiv:2311.04957;
  sub-GeV の速度分布依存性 arXiv:2001.09156, arXiv:2011.02493;
  フォノン励起 arXiv:2606.04091)。差別化軸は
  (a) meV–eV 級超伝導量子センサー(TES vs MKID)の具体設計での比較、
  (b) Disk/Bound 成分の復元・判別という DM 天文学的な問い。
- [ ] 主要な判別有意度の結果(bkg-none、Asimov):
  mix25 はほぼ全設定で >4σ(最強 TiN q2 M1 の ~50σ、最弱 Al q2 M3 の 2.5σ)、
  mix5 は 10–100 MeV で 2.4–9.8σ / 1 GeV で <1σ、Bound は M3 のみ・桁外れに有意。
  → 背景ありでどこまで残るかが第2論文の核。

## 細かいこと

- [ ] git のコミット作者がホスト名由来(`koichiromacstudio@Mac.lan`)のまま。
  GitHub アカウントに紐付けるなら `git config --global user.name / user.email`。
