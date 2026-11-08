"""Marrero-Gani法の第1次から第3次グループを、化合物データへ一括で適用する。

marrero_gani_groups.py (同じディレクトリに置く) の第1次・第2次・第3次の
解析関数を呼び、化合物ごとの結果と、機械学習の特徴量に使える行列を返す。

主なAPI:
    analyze_compound(smiles)         1化合物の第1-3次の結果 (CompoundAnalysis)
    analyze_many(smiles_list)        複数化合物の結果 (重複SMILESは1回だけ解析)
    compute_group_features(df, ...)  DataFrame -> グループ数の特徴量 (pandas)
    group_catalog()                  全370グループの一覧 (列名との対応表)

コマンドライン:
    python marrero_gani_batch.py input.csv -o output.csv --smiles-col smiles
    python marrero_gani_batch.py --selftest

設計上の注意:
    * 列は第1次182 + 第2次122 + 第3次66 = 370列で固定する。学習用と予測用で
      列の構成がずれないようにするため (only_present=True は固定しない)。
    * 解析できない行 (SMILES不正など) は、グループ数を欠損 (<NA>) にする。
      未割当原子などで結果が信頼できない行は、グループ数は返すが
      mg_status を "incomplete" にする。必ず mg_status を確認すること。
"""

from __future__ import annotations

import argparse
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field
from functools import partial
from pathlib import Path
from typing import TYPE_CHECKING, Dict, Iterable, List, Optional, Sequence, Tuple

from rdkit import Chem, rdBase

import marrero_gani_groups as mg

if TYPE_CHECKING:  # pandas は DataFrame を使う関数でのみ必要
    import pandas as pd

STATUS_OK = "ok"
STATUS_INCOMPLETE = "incomplete"
STATUS_ERROR = "error"
META_COLUMNS = (
    "mg_status", "mg_unassigned", "mg_warnings", "mg_notes", "mg_error",
)
_SEPARATOR = " | "


# ---------------------------------------------------------------------------
# グループの一覧 (列名)
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class GroupColumn:
    """1つのグループと、その特徴量の列名。

    Attributes:
        order: 次数 (1, 2, 3)。
        no: 論文 Table 6-8 のグループ番号。
        name: グループ名。
    """

    order: int
    no: int
    name: str

    @property
    def column_id(self) -> str:
        """機械的に扱いやすい列名 (例: o2_g026)。"""
        return f"o{self.order}_g{self.no:03d}"

    @property
    def column_name(self) -> str:
        """グループ名つきの列名 (例: 2:CHm(NHn)-COOH)。"""
        return f"{self.order}:{self.name}"

    def column(self, style: str) -> str:
        """style は "id" または "name"。"""
        if style == "id":
            return self.column_id
        if style == "name":
            return self.column_name
        raise ValueError(f"column_style must be 'id' or 'name': {style!r}")


def group_catalog() -> List[GroupColumn]:
    """第1-3次の全グループ (370個) を、次数・番号の順に返す。"""
    catalog = [GroupColumn(1, g.no, g.name) for g in mg.GROUP_DEFS]
    catalog += [GroupColumn(2, g.no, g.name) for g in mg.SECOND_ORDER_DEFS]
    catalog += [GroupColumn(3, g.no, g.name) for g in mg.THIRD_ORDER_DEFS]
    return sorted(catalog, key=lambda c: (c.order, c.no))


# ---------------------------------------------------------------------------
# 1化合物の解析
# ---------------------------------------------------------------------------
@dataclass
class CompoundAnalysis:
    """1化合物の第1-3次の解析結果。

    Attributes:
        smiles: 入力のSMILES。
        status: "ok", "incomplete" (結果が信頼できない), "error" (解析不能)。
        first: 第1次のグループ名 -> 出現回数。
        second: 第2次のグループ名 -> 出現回数。
        third: 第3次のグループ名 -> 出現回数。
        unassigned: 第1次で未割当の重原子 (インデックス, 元素記号)。
        warnings: 前処理の警告 (空でなければ結果を信頼できない)。
        notes: 解釈に注意が要る事項 (タイブレーク、フォールバック、例外適用など)。
        error: status が "error" のときの理由。
    """

    smiles: str
    status: str
    first: Dict[str, int] = field(default_factory=dict)
    second: Dict[str, int] = field(default_factory=dict)
    third: Dict[str, int] = field(default_factory=dict)
    unassigned: List[Tuple[int, str]] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)
    error: str = ""

    def counts(self, style: str = "id") -> Dict[str, int]:
        """グループ数を、列名 (style) -> 出現回数の辞書で返す。0の列は含まない。"""
        by_order = {1: self.first, 2: self.second, 3: self.third}
        result: Dict[str, int] = {}
        for col in group_catalog():
            n = by_order[col.order].get(col.name, 0)
            if n:
                result[col.column(style)] = n
        return result


def _error(smiles: object, message: str) -> CompoundAnalysis:
    return CompoundAnalysis(
        smiles="" if smiles is None else str(smiles),
        status=STATUS_ERROR,
        error=message,
    )


def analyze_compound(
    smiles: object,
    paper_aromaticity: bool = True,
    keep_largest_fragment: bool = False,
) -> CompoundAnalysis:
    """1化合物を第1-3次グループへ分解する。例外は投げず status に反映する。

    Args:
        smiles: 化合物のSMILES。文字列以外 (欠損値など) は error になる。
        paper_aromaticity: Trueなら論文流の芳香環定義を適用する。
        keep_largest_fragment: Trueなら、塩や溶媒和物のような複数成分のSMILESで
            重原子数が最大の成分だけを解析する。Falseなら全成分を解析し、
            複数成分であることを注記して status を "incomplete" にする。

    Returns:
        CompoundAnalysis。
    """
    if not isinstance(smiles, str) or not smiles.strip():
        return _error(smiles, "SMILESが空、または文字列ではない")
    text = smiles.strip()
    with rdBase.BlockLogs():  # 不正SMILESのログは error として結果に残す
        mol = Chem.MolFromSmiles(text)
    if mol is None:
        return _error(smiles, "SMILESを解釈できない")

    notes: List[str] = []
    multi_component = False
    fragments = Chem.GetMolFrags(mol, asMols=True)
    if len(fragments) > 1:
        if keep_largest_fragment:
            largest = max(fragments, key=lambda m: m.GetNumHeavyAtoms())
            text = Chem.MolToSmiles(largest)
            notes.append(f"複数成分のため最大成分 {text} だけを解析")
        else:
            multi_component = True
            notes.append("複数成分のSMILES (全成分を解析。純成分の推算には不向き)")

    try:
        first = mg.analyze_smiles(text, paper_aromaticity)
        second = mg.analyze_second_order(text, paper_aromaticity)
        third = mg.analyze_third_order(text, paper_aromaticity)
    except Exception as exc:  # バッチ処理を止めないため、解析器の例外は行の誤りとして扱う
        return _error(smiles, f"解析中の例外: {type(exc).__name__}: {exc}")

    if first.ties:
        notes.append(f"第1次のタイブレーク (D7): {first.ties}")
    if second.ties:
        notes.append(f"第2次で原子集合が一致: {second.ties}")
    notes.extend(third.notes)

    incomplete = bool(first.unassigned or first.warnings or multi_component)
    return CompoundAnalysis(
        smiles=smiles,
        status=STATUS_INCOMPLETE if incomplete else STATUS_OK,
        first=first.groups,
        second=second.groups,
        third=third.groups,
        unassigned=first.unassigned,
        warnings=first.warnings,
        notes=notes,
    )


# ---------------------------------------------------------------------------
# 複数化合物の解析
# ---------------------------------------------------------------------------
def analyze_many(
    smiles_list: Iterable[object],
    paper_aromaticity: bool = True,
    keep_largest_fragment: bool = False,
    n_jobs: int = 1,
) -> List[CompoundAnalysis]:
    """複数化合物を解析する。入力と同じ順序・同じ長さのリストを返す。

    同じSMILESは1回だけ解析する。

    Args:
        smiles_list: SMILESの列。
        paper_aromaticity: analyze_compound を参照。
        keep_largest_fragment: analyze_compound を参照。
        n_jobs: 2以上ならプロセス並列で解析する (大量データ向け)。

    Returns:
        入力順の CompoundAnalysis のリスト。
    """
    items = list(smiles_list)
    # 欠損値などの文字列以外は、キーを None にまとめる
    keys = [s if isinstance(s, str) else None for s in items]
    unique = list(dict.fromkeys(keys))
    func = partial(
        analyze_compound,
        paper_aromaticity=paper_aromaticity,
        keep_largest_fragment=keep_largest_fragment,
    )
    if n_jobs > 1 and len(unique) > 1:
        chunk = max(1, len(unique) // (n_jobs * 8))
        with ProcessPoolExecutor(max_workers=n_jobs) as pool:
            results = list(pool.map(func, unique, chunksize=chunk))
    else:
        results = [func(s) for s in unique]
    by_key = dict(zip(unique, results))
    return [by_key[k] for k in keys]


# ---------------------------------------------------------------------------
# DataFrame (特徴量行列)
# ---------------------------------------------------------------------------
def compute_group_features(
    df: "pd.DataFrame",
    smiles_col: str = "smiles",
    column_style: str = "id",
    only_present: bool = False,
    paper_aromaticity: bool = True,
    keep_largest_fragment: bool = False,
    n_jobs: int = 1,
) -> "pd.DataFrame":
    """DataFrameの各行のSMILESから、第1-3次のグループ数の特徴量を作る。

    Args:
        df: SMILESの列を持つDataFrame。
        smiles_col: SMILESの列名。
        column_style: 列名の形式。"id" (o1_g001 など) か "name" (1:CH3 など)。
        only_present: Trueなら、データ中で1回も現れないグループの列を除く。
            列の構成がデータに依存するため、学習用と予測用で別々に作ると
            列がずれる。モデルに使う場合は False (固定の370列) にすること。
        paper_aromaticity: analyze_compound を参照。
        keep_largest_fragment: analyze_compound を参照。
        n_jobs: analyze_many を参照。

    Returns:
        dfと同じindexのDataFrame。グループ数の列 (Int64。解析不能な行は <NA>) と、
        META_COLUMNS のメタ列 (mg_status など) を持つ。

    Raises:
        KeyError: smiles_col が df にない場合。
    """
    import pandas as pd

    if smiles_col not in df.columns:
        raise KeyError(f"column not found: {smiles_col!r}")
    results = analyze_many(
        df[smiles_col], paper_aromaticity, keep_largest_fragment, n_jobs
    )
    catalog = group_catalog()
    columns = [c.column(column_style) for c in catalog]

    rows = []
    for res in results:
        counts = res.counts(column_style)
        if res.status == STATUS_ERROR:
            rows.append({col: pd.NA for col in columns})
        else:
            rows.append({col: counts.get(col, 0) for col in columns})
    features = pd.DataFrame(rows, index=df.index, columns=columns).astype("Int64")

    meta = pd.DataFrame(
        {
            "mg_status": [r.status for r in results],
            "mg_unassigned": [len(r.unassigned) for r in results],
            "mg_warnings": [_SEPARATOR.join(r.warnings) for r in results],
            "mg_notes": [_SEPARATOR.join(r.notes) for r in results],
            "mg_error": [r.error for r in results],
        },
        index=df.index,
    )
    if only_present:
        features = features.loc[:, (features.fillna(0) != 0).any(axis=0)]
    return pd.concat([features, meta], axis=1)


def catalog_frame() -> "pd.DataFrame":
    """group_catalog() をDataFrame (列名との対応表) にして返す。"""
    import pandas as pd

    return pd.DataFrame(
        [
            {
                "order": c.order,
                "no": c.no,
                "name": c.name,
                "column_id": c.column_id,
                "column_name": c.column_name,
            }
            for c in group_catalog()
        ]
    )


# ---------------------------------------------------------------------------
# 自己テスト
# ---------------------------------------------------------------------------
def _self_test() -> List[str]:
    """ラッパーの動作を検証し、失敗メッセージのリストを返す (空なら成功)。"""
    import pandas as pd

    failures: List[str] = []
    catalog = group_catalog()
    if len(catalog) != 182 + 122 + 66:
        failures.append(f"catalog size {len(catalog)} != 370")
    for style in ("id", "name"):
        cols = [c.column(style) for c in catalog]
        if len(set(cols)) != len(cols):
            failures.append(f"duplicate column names (style={style})")

    smiles = [
        "CC(N)C(=O)O",                                        # l-alanine
        "OC1=C(Cc2c(O)c3ccccc3oc2=O)C(=O)Oc2ccccc12",         # dicoumarol
        "CN1C(=O)C(O)N=C(c2ccccc2F)c2cc(Cl)ccc21",            # flutemazepan
        "c1cc2ccc3cccc4ccc(c1)c2c34",                         # pyrene
        "NCCNCCN",                                            # diethylenetriamine
        "CC(N)C(=O)O",                                        # 重複
        "not a smiles",                                       # 不正
        "",                                                   # 空
        None,                                                 # 欠損
        "CC(=O)[O-].[Na+]",                                   # 複数成分
    ]
    results = analyze_many(smiles)
    if len(results) != len(smiles):
        failures.append("analyze_many length mismatch")
        return failures

    # 個別の解析関数と結果が一致すること
    for smi, res in zip(smiles[:5], results[:5]):
        if res.first != mg.analyze_smiles(smi).groups:
            failures.append(f"first-order mismatch: {smi}")
        if res.second != mg.count_second_order_groups(smi):
            failures.append(f"second-order mismatch: {smi}")
        if res.third != mg.count_third_order_groups(smi):
            failures.append(f"third-order mismatch: {smi}")
    expected_status = [
        STATUS_OK, STATUS_OK, STATUS_OK, STATUS_OK, STATUS_OK, STATUS_OK,
        STATUS_ERROR, STATUS_ERROR, STATUS_ERROR, STATUS_INCOMPLETE,
    ]
    got_status = [r.status for r in results]
    if got_status != expected_status:
        failures.append(f"status {got_status} != {expected_status}")
    if results[0] is not results[5] and results[0].first != results[5].first:
        failures.append("duplicate SMILES gave different results")

    # 最大成分の指定
    largest = analyze_compound("CCO.[Na+].[Cl-]", keep_largest_fragment=True)
    if largest.first != mg.analyze_smiles("CCO").groups:
        failures.append(f"keep_largest_fragment: {largest.first}")
    if largest.status != STATUS_OK or not any("最大成分" in n for n in largest.notes):
        failures.append("keep_largest_fragment: status/note wrong")
    if analyze_compound("CCO.[Na+].[Cl-]").status != STATUS_INCOMPLETE:
        failures.append("multi-component without option should be incomplete")

    # DataFrame
    df = pd.DataFrame({"id": range(len(smiles)), "smiles": smiles},
                      index=[f"r{i}" for i in range(len(smiles))])
    feat = compute_group_features(df)
    if feat.shape != (len(smiles), 370 + len(META_COLUMNS)):
        failures.append(f"feature shape {feat.shape}")
    if list(feat.index) != list(df.index):
        failures.append("feature index mismatch")
    if not feat.loc["r6", "o1_g001":"o3_g066"].isna().all():
        failures.append("error row should be NA")
    if feat.loc["r0", "o2_g026"] != 1 or feat.loc["r0", "o2_g013"] != 0:
        failures.append("l-alanine: second-order columns wrong")
    if not feat.loc["r0"].equals(feat.loc["r5"]):
        failures.append("duplicate rows differ")
    if feat.loc["r9", "mg_status"] != STATUS_INCOMPLETE:
        failures.append("multi-component row should be incomplete")
    named = compute_group_features(df, column_style="name")
    if "2:CHm(NHn)-COOH" not in named.columns:
        failures.append("name-style columns missing")
    present = compute_group_features(df, only_present=True)
    if present.shape[1] >= feat.shape[1]:
        failures.append("only_present did not drop columns")

    # 並列実行と逐次実行が同じ結果になること
    parallel = compute_group_features(df, n_jobs=2)
    if not parallel.equals(feat):
        failures.append("n_jobs=2 differs from n_jobs=1")

    # CSV入出力 (コマンドライン)
    with tempfile.TemporaryDirectory() as tmp:
        src, dst, cat = (Path(tmp) / n for n in ("in.csv", "out.csv", "cat.csv"))
        df.to_csv(src, index=False)
        code = main([str(src), "-o", str(dst), "--catalog", str(cat)])
        out = pd.read_csv(dst)
        if code != 0 or out.shape != (len(smiles), 2 + 370 + len(META_COLUMNS)):
            failures.append(f"cli output shape {out.shape}, code {code}")
        if len(pd.read_csv(cat)) != 370:
            failures.append("catalog csv size")
    return failures


# ---------------------------------------------------------------------------
# コマンドライン
# ---------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    """コマンドライン引数のパーサを作る。"""
    parser = argparse.ArgumentParser(
        description="CSVのSMILES列から、Marrero-Gani法の第1-3次グループ数を求める。"
    )
    parser.add_argument("input", nargs="?", help="入力CSV")
    parser.add_argument("-o", "--output", help="出力CSV (元の列に特徴量を追加)")
    parser.add_argument("--smiles-col", default="smiles", help="SMILES列の名前")
    parser.add_argument("--sep", default=",", help="区切り文字 (TSVなら '\\t')")
    parser.add_argument("--column-style", choices=("id", "name"), default="id",
                        help="列名の形式 (id: o1_g001, name: 1:CH3)")
    parser.add_argument("--only-present", action="store_true",
                        help="データ中に現れないグループの列を除く (列が固定されない)")
    parser.add_argument("--largest-fragment", action="store_true",
                        help="複数成分のSMILESは最大成分だけを解析する")
    parser.add_argument("--no-paper-aromaticity", action="store_true",
                        help="論文流の芳香環定義を使わず、RDKitの判定をそのまま使う")
    parser.add_argument("--n-jobs", type=int, default=1, help="並列プロセス数")
    parser.add_argument("--catalog", help="グループ一覧 (列名との対応表) の出力先CSV")
    parser.add_argument("--selftest", action="store_true", help="自己テストを実行")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    """コマンドラインの入口。終了コードを返す。"""
    args = build_parser().parse_args(argv)
    if args.selftest:
        failures = _self_test()
        for failure in failures:
            print(f"[FAIL] {failure}", file=sys.stderr)
        print(f"selftest: {'OK' if not failures else 'FAILED'}")
        return 1 if failures else 0
    if not args.input or not args.output:
        build_parser().error("input と -o/--output が必要です (または --selftest)")

    import pandas as pd

    sep = "\t" if args.sep == "\\t" else args.sep
    df = pd.read_csv(args.input, sep=sep)
    features = compute_group_features(
        df,
        smiles_col=args.smiles_col,
        column_style=args.column_style,
        only_present=args.only_present,
        paper_aromaticity=not args.no_paper_aromaticity,
        keep_largest_fragment=args.largest_fragment,
        n_jobs=args.n_jobs,
    )
    pd.concat([df, features], axis=1).to_csv(args.output, index=False)
    if args.catalog:
        catalog_frame().to_csv(args.catalog, index=False)
    counts = features["mg_status"].value_counts().to_dict()
    print(f"{len(df)} compounds -> {args.output}  status: {counts}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
