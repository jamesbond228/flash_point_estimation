# Marrero-Gani法による沸点推算

SMILESを入力すると、Marrero-Gani法（Marrero & Gani, Fluid Phase Equilib. 183-184, 2001）の第1〜3次グループ寄与から正常沸点 Tb [K] を推算する。

## 構成

| ファイル | 役割 |
|---|---|
| `marrero_gani_tb.py` | 沸点推算（本体。CLI・Python API） |
| `marrero_gani_groups.py` | SMILESから第1〜3次グループへの分解 |
| `marrero_gani_batch.py` | 1化合物の解析ラッパー（エラー・複数成分の処理） |
| `Table_6.csv` | 第1次グループの寄与値（論文 Table 6） |
| `Table_7.csv` | 第2次グループの寄与値（論文 Table 7） |
| `Table_8.csv` | 第3次グループの寄与値（論文 Table 8） |

全ファイルを同じディレクトリに置く。CSVの場所は `--table-dir` で変更できる。

## 環境

- Python 3.9 以上
- RDKit（`pip install rdkit`）

推算スクリプトは pandas を必要としない。

## 使い方

```bash
# SMILESを直接指定（複数可）
python marrero_gani_tb.py "CCCCO" "c1ccccc1"

# グループ内訳を表示
python marrero_gani_tb.py "NCCCCO" --detail

# CSV一括処理（入力列に推算結果の列を追加して出力）
python marrero_gani_tb.py --csv in.csv -o out.csv --smiles-col smiles

# 対話入力（引数なし）
python marrero_gani_tb.py

# 自己テスト
python marrero_gani_tb.py --selftest
```

### オプション

| オプション | 内容 |
|---|---|
| `--detail` | グループ内訳と寄与和を表示 |
| `--csv PATH` / `-o PATH` | CSV一括処理の入力 / 出力 |
| `--smiles-col NAME` | CSVのSMILES列名（既定: `smiles`） |
| `--table-dir DIR` | `Table_6/7/8.csv` のあるディレクトリ |
| `--allow-incomplete` | 第1次で未割当の原子があっても推算する（信頼できない） |
| `--largest-fragment` | 複数成分SMILESは最大成分だけを使う |

### 実行例

```
$ python marrero_gani_tb.py "NCCCCO"
SMILES : NCCCCO
status : ok
Tb(1st) :   448.64 K
Tb(2nd) :   448.64 K
Tb(3rd) :   478.42 K
Tb (final, 3rd order): 478.42 K (205.27 degC)
```

### CSV出力の列

入力の列に、次の列が追加される。

`tb_k`, `tb_c`, `tb_1st_k`, `tb_2nd_k`, `tb_3rd_k`, `tb_status`, `tb_messages`

## Python API

```python
from marrero_gani_tb import estimate_boiling_point, estimate_many

res = estimate_boiling_point("NCCCCO")
res.status        # "ok" / "uncertain" / "error"
res.tb_k          # 最終推算値 [K]（第3近似）。errorのときはNone
res.tb_c          # 最終推算値 [degC]
res.tb_by_level   # {1: ..., 2: ..., 3: ...} 各近似の推算値 [K]
res.sums          # {1: ..., 2: ..., 3: ...} 各近似の寄与和
res.terms         # 使用したグループ (次数, 番号, 名前, 出現回数, 寄与値)
res.messages      # 注意事項・エラー理由

results = estimate_many(["CCO", "c1ccccc1"])  # 入力順のリスト
```

関数は例外を投げず、結果を `status` に反映する。

## 推算式

```
Tb = Tb0 * ln( sum_i N_i*Tb1i + sum_j M_j*Tb2j + sum_k O_k*Tb3k )
Tb0 = 222.543 K  (論文 Table 2)
```

| 近似 | 使うグループ |
|---|---|
| 第1近似 | 第1次のみ |
| 第2近似 | 第1次 + 第2次 |
| 第3近似 | 第1次 + 第2次 + 第3次 |

最終推算値は第3近似。該当する第2・3次グループがなければ下位の近似と同値になる。

## status の意味

| status | 内容 |
|---|---|
| `ok` | 問題なし |
| `uncertain` | 推算は実行したが信頼性に注意。`messages` を確認する |
| `error` | 推算不能。`tb_k` は None |

### uncertain になる場合

- 第2・3次グループに論文表のTb寄与値がなく（`*****`）、その項を無視した
- グループ分解にタイブレークや例外適用の注記がある
- `--allow-incomplete` で未割当原子を許した

### error になる場合

- SMILESが空、または解釈できない
- 第1次で未割当の原子がある（既定）
- 複数成分SMILES（塩・溶媒和物など。既定）
- Tb寄与値が論文表にない第1次グループを使う
- 寄与和が 0 以下で ln を計算できない

## 検証状況

`--selftest` で次を確認する。

- 論文 Appendix B の沸点例題3件（N-phenyl-1,4-benzenediamine、pyrene、4-aminobutanol）の第1・第2・第3近似が、論文値と0.02 K以内で一致する
- 論文表の寄与値から手計算した5化合物（エタノール、1-プロパノール、2-ブタノン、エチルベンゼン、ナフタレン）と一致する
- 不正入力・複数成分が error になり、例外を投げない

実測沸点との精度検証は含まない。

## 制約と注意

- 論文の適用範囲は炭素数3〜60。範囲外（エタノールなど）や小分子では誤差が大きくなる傾向がある。
- 論文表のとおりに計算した結果でも、文献値から数十K離れる化合物がある。実測値つきのデータで精度を確認してから使うこと。
- グループ分解の限界（無置換尿素、ホルムアミド、アシルハライド、チオエステル、ジスルフィド、ピリジンN-オキシド、荷電した環Nなどは未割当になる）は `marrero_gani_groups.py` の判断記録（D1〜D9以降）を参照。
- 論文表にTb寄与値がない第2・3次グループは無視される。該当する場合は `uncertain` と `messages` で分かる。
- CSVは論文Table 6〜8の転記。Appendix Bで使われる寄与値は一致を確認済みだが、全件のPDF照合は別途必要。

## 参考文献

J. Marrero, R. Gani, "Group-contribution based estimation of pure component properties", Fluid Phase Equilibria 183-184 (2001) 183-208.
