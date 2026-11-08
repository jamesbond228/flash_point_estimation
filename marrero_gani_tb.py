"""Marrero-Gani法 (Marrero & Gani, Fluid Phase Equilib. 183-184, 2001) による
正常沸点 (Tb) の推算。SMILESを入力すると、第1-3次グループ寄与から Tb [K] を返す。

推算式 (論文 Table 1, Eq.(1)):
    exp(Tb / Tb0) = sum_i N_i*Tb1i + sum_j M_j*Tb2j + sum_k O_k*Tb3k
    Tb = Tb0 * ln(右辺)        Tb0 = 222.543 K (論文 Table 2)

    第1近似: 第1次のみ / 第2近似: 第1+2次 / 第3近似: 第1+2+3次。
    論文の最終推算値は第3近似 (該当する第2・3次グループがなければ下位と同値)。

必要なファイル (同じディレクトリ):
    marrero_gani_groups.py   グループ分解 (第1-3次)
    marrero_gani_batch.py    1化合物の解析ラッパー
    Table_6.csv, Table_7.csv, Table_8.csv   論文の寄与値表 (--table-dir で変更可)

使い方:
    python marrero_gani_tb.py "CCO" "c1ccccc1"            SMILESを直接指定
    python marrero_gani_tb.py "CCO" --detail              グループ内訳つき
    python marrero_gani_tb.py --csv in.csv -o out.csv     CSV一括 (--smiles-col)
    python marrero_gani_tb.py                             対話入力
    python marrero_gani_tb.py --selftest                  自己テスト

結果の信頼性 (status):
    ok         問題なし
    uncertain  推算は実行したが信頼性に注意 (寄与値のない第2/3次グループの
               無視、タイブレーク、--allow-incomplete の使用など)。messages を確認。
    error      推算不能 (SMILES不正、未割当原子、Tb寄与値のない第1次グループ等)。
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import marrero_gani_batch as mgb
import marrero_gani_groups as mg

TB0 = 222.543  # 論文 Table 2 の Tb0 [K]
MISSING_MARK = "*****"  # 論文表で寄与値が存在しないことを示す記号
KELVIN_OFFSET = 273.15

# 次数 -> (CSVファイル名, Tb列名, 期待行数)
_TABLES: Dict[int, Tuple[str, str, int]] = {
    1: ("Table_6.csv", "Tb1i", 182),
    2: ("Table_7.csv", "Tb2j", 122),
    3: ("Table_8.csv", "Tb3k", 66),
}

STATUS_OK = "ok"
STATUS_UNCERTAIN = "uncertain"
STATUS_ERROR = "error"


# ---------------------------------------------------------------------------
# 寄与値テーブル
# ---------------------------------------------------------------------------
class TableError(RuntimeError):
    """寄与値CSVの読み込み・検証エラー。"""


def _parse_value(text: str, path: Path, line: int) -> Optional[float]:
    """寄与値の文字列を float にする。'*****' は None (値なし) を返す。"""
    value = text.strip()
    if value == MISSING_MARK:
        return None
    try:
        return float(value)
    except ValueError as exc:
        raise TableError(f"{path.name}:{line}: 数値でない値 {text!r}") from exc


def _load_order_table(
    path: Path, column: str, expected_rows: int
) -> Dict[int, Optional[float]]:
    """1つのCSVから Group_No -> Tb寄与値 (値なしは None) の辞書を作る。"""
    if not path.is_file():
        raise TableError(f"寄与値CSVが見つからない: {path}")
    table: Dict[int, Optional[float]] = {}
    with path.open(encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh)
        fields = reader.fieldnames or []
        for required in ("Group_No", column):
            if required not in fields:
                raise TableError(f"{path.name}: 列 {required!r} がない")
        for row in reader:
            line = reader.line_num
            try:
                no = int(row["Group_No"])
            except (TypeError, ValueError) as exc:
                raise TableError(
                    f"{path.name}:{line}: Group_No が整数でない"
                ) from exc
            if no in table:
                raise TableError(f"{path.name}:{line}: Group_No {no} が重複")
            table[no] = _parse_value(row[column], path, line)
    if len(table) != expected_rows:
        raise TableError(
            f"{path.name}: 行数 {len(table)} != 期待値 {expected_rows}"
        )
    return table


@lru_cache(maxsize=4)
def load_tb_tables(
    table_dir: Optional[str] = None,
) -> Dict[int, Dict[int, Optional[float]]]:
    """第1-3次の Tb 寄与値を読み込む (次数 -> グループ番号 -> 寄与値)。

    Args:
        table_dir: Table_6/7/8.csv のあるディレクトリ。None なら本ファイルと
            同じディレクトリ。

    Raises:
        TableError: ファイル欠落、列欠落、行数不一致、数値不正の場合。
    """
    base = Path(table_dir) if table_dir else Path(__file__).resolve().parent
    return {
        order: _load_order_table(base / name, column, rows)
        for order, (name, column, rows) in _TABLES.items()
    }


@lru_cache(maxsize=1)
def _name_to_no() -> Dict[int, Dict[str, int]]:
    """次数 -> (グループ名 -> グループ番号)。名前の重複は例外にする。"""
    defs = {
        1: mg.GROUP_DEFS,
        2: mg.SECOND_ORDER_DEFS,
        3: mg.THIRD_ORDER_DEFS,
    }
    result: Dict[int, Dict[str, int]] = {}
    for order, group_defs in defs.items():
        mapping = {g.name: g.no for g in group_defs}
        if len(mapping) != len(group_defs):
            raise TableError(f"第{order}次のグループ名が重複している")
        result[order] = mapping
    return result


# ---------------------------------------------------------------------------
# 推算
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class GroupTerm:
    """推算に使ったグループ1種類分の項。

    Attributes:
        order: 次数 (1, 2, 3)。
        no: 論文 Table 6-8 のグループ番号。
        name: グループ名。
        count: 出現回数。
        contribution: Tb寄与値。None は論文表に値がないこと。
    """

    order: int
    no: int
    name: str
    count: int
    contribution: Optional[float]

    @property
    def total(self) -> float:
        """count * contribution (値なしは 0)。"""
        return 0.0 if self.contribution is None else (
            self.count * self.contribution
        )


@dataclass
class BoilingPointResult:
    """沸点推算の結果。

    Attributes:
        smiles: 入力のSMILES。
        status: "ok", "uncertain", "error"。
        tb_k: 最終推算値 [K] (第3近似)。error のとき None。
        tb_by_level: 近似次数 (1, 2, 3) -> 推算値 [K]。
        sums: 近似次数 -> exp(Tb/Tb0) の右辺の累積和。
        terms: 使用したグループの項 (次数、番号順)。
        messages: 注意事項・エラー理由。
    """

    smiles: str
    status: str
    tb_k: Optional[float] = None
    tb_by_level: Dict[int, float] = field(default_factory=dict)
    sums: Dict[int, float] = field(default_factory=dict)
    terms: List[GroupTerm] = field(default_factory=list)
    messages: List[str] = field(default_factory=list)

    @property
    def tb_c(self) -> Optional[float]:
        """最終推算値 [degC]。"""
        return None if self.tb_k is None else self.tb_k - KELVIN_OFFSET


def _error(smiles: object, message: str) -> BoilingPointResult:
    return BoilingPointResult(
        smiles="" if smiles is None else str(smiles),
        status=STATUS_ERROR,
        messages=[message],
    )


def _build_terms(
    analysis: mgb.CompoundAnalysis,
    tables: Dict[int, Dict[int, Optional[float]]],
) -> List[GroupTerm]:
    """グループ分解の結果を、Tb寄与値つきの項に変換する。"""
    mapping = _name_to_no()
    by_order = {1: analysis.first, 2: analysis.second, 3: analysis.third}
    terms: List[GroupTerm] = []
    for order in (1, 2, 3):
        for name, count in by_order[order].items():
            no = mapping[order][name]
            terms.append(
                GroupTerm(order, no, name, count, tables[order][no])
            )
    return sorted(terms, key=lambda t: (t.order, t.no))


def _log_estimate(total: float) -> Optional[float]:
    """Tb = Tb0 * ln(total)。total <= 0 のときは None。"""
    return TB0 * math.log(total) if total > 0 else None


def estimate_boiling_point(
    smiles: object,
    table_dir: Optional[str] = None,
    allow_incomplete: bool = False,
    keep_largest_fragment: bool = False,
    paper_aromaticity: bool = True,
) -> BoilingPointResult:
    """SMILESから正常沸点 [K] を推算する。例外は投げず status に反映する。

    Args:
        smiles: 化合物のSMILES。
        table_dir: 寄与値CSVのディレクトリ。None なら本ファイルと同じ場所。
        allow_incomplete: True なら、第1次グループで未割当の原子があっても
            推算する (status は "uncertain")。未割当原子の寄与は欠落するため、
            結果は信頼できない。False なら error にする。
        keep_largest_fragment: True なら複数成分SMILESの最大成分だけを使う。
        paper_aromaticity: 論文流の芳香環定義を使う (groups.py の P4)。

    Returns:
        BoilingPointResult。
    """
    try:
        tables = load_tb_tables(table_dir)
    except TableError as exc:
        return _error(smiles, f"寄与値テーブルのエラー: {exc}")

    analysis = mgb.analyze_compound(
        smiles,
        paper_aromaticity=paper_aromaticity,
        keep_largest_fragment=keep_largest_fragment,
    )
    if analysis.status == mgb.STATUS_ERROR:
        return _error(smiles, analysis.error)

    messages: List[str] = []
    uncertain = False
    if analysis.status == mgb.STATUS_INCOMPLETE:
        detail = []
        if analysis.unassigned:
            atoms = ", ".join(f"{s}[{i}]" for i, s in analysis.unassigned)
            detail.append(f"第1次で未割当の原子: {atoms}")
        detail.extend(analysis.warnings)
        detail.extend(n for n in analysis.notes if "複数成分" in n)
        text = "; ".join(detail) or "グループ分解が不完全"
        if not allow_incomplete:
            return _error(
                smiles,
                f"{text} (推算不能。--allow-incomplete で欠落を許して実行可)",
            )
        messages.append(f"{text} (欠落を許して推算: 信頼できない)")
        uncertain = True
    messages.extend(analysis.notes)
    uncertain = uncertain or bool(analysis.notes)

    terms = _build_terms(analysis, tables)
    if not terms:
        return _error(smiles, "該当する第1次グループがない")

    # 第1次の寄与値がなければ推算できない。第2・3次は無視して続行する。
    missing_first = [t.name for t in terms if t.order == 1 and
                     t.contribution is None]
    if missing_first:
        return _error(
            smiles,
            "Tb寄与値が論文表にない第1次グループ: " + ", ".join(missing_first),
        )
    for term in terms:
        if term.order > 1 and term.contribution is None:
            messages.append(
                f"第{term.order}次 No.{term.no} {term.name} は論文表に"
                "Tb寄与値がないため無視"
            )
            uncertain = True

    sums: Dict[int, float] = {}
    running = 0.0
    for order in (1, 2, 3):
        running += sum(t.total for t in terms if t.order == order)
        sums[order] = running

    tb_by_level: Dict[int, float] = {}
    for order in (1, 2, 3):
        value = _log_estimate(sums[order])
        if value is None:
            return BoilingPointResult(
                smiles=str(smiles),
                status=STATUS_ERROR,
                sums=sums,
                terms=terms,
                messages=messages + [
                    f"第{order}近似の寄与和が {sums[order]:.4f} <= 0 で "
                    "ln を計算できない (適用範囲外)"
                ],
            )
        tb_by_level[order] = value

    return BoilingPointResult(
        smiles=str(smiles),
        status=STATUS_UNCERTAIN if uncertain else STATUS_OK,
        tb_k=tb_by_level[3],
        tb_by_level=tb_by_level,
        sums=sums,
        terms=terms,
        messages=messages,
    )


def estimate_many(
    smiles_list: Sequence[object], **kwargs: object
) -> List[BoilingPointResult]:
    """複数のSMILESを推算する。入力と同じ順序・長さのリストを返す。"""
    cache: Dict[object, BoilingPointResult] = {}
    results = []
    for smiles in smiles_list:
        key = smiles if isinstance(smiles, str) else None
        if key not in cache:
            cache[key] = estimate_boiling_point(smiles, **kwargs)
        results.append(cache[key])
    return results


# ---------------------------------------------------------------------------
# 表示
# ---------------------------------------------------------------------------
def format_result(result: BoilingPointResult, detail: bool = False) -> str:
    """結果を人が読む形式の文字列にする。"""
    lines = [f"SMILES : {result.smiles}", f"status : {result.status}"]
    if result.tb_k is not None:
        for order, label in ((1, "1st"), (2, "2nd"), (3, "3rd")):
            lines.append(
                f"Tb({label}) : {result.tb_by_level[order]:8.2f} K"
            )
        lines.append(
            f"Tb (final, 3rd order): {result.tb_k:.2f} K "
            f"({result.tb_c:.2f} degC)"
        )
    if detail and result.terms:
        lines.append("groups :")
        for t in result.terms:
            contrib = "n/a" if t.contribution is None else (
                f"{t.contribution:+.4f}"
            )
            lines.append(
                f"  [{t.order}] No.{t.no:<3d} {t.name:<45s} "
                f"x{t.count:<2d} {contrib}"
            )
        lines.append(
            "sums   : " + ", ".join(
                f"L{k}={v:.4f}" for k, v in result.sums.items()
            )
        )
    for message in result.messages:
        lines.append(f"note   : {message}")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# 自己テスト
# ---------------------------------------------------------------------------
# 論文 Appendix B の Tb 例題 (SMILES, 第1/2/3近似の推算値 [K])
APPENDIX_B_TB: List[Tuple[str, str, Tuple[float, float, float]]] = [
    ("Ex.1 N-phenyl-1,4-benzenediamine", "Nc1ccc(Nc2ccccc2)cc1",
     (614.62, 616.03, 623.94)),
    ("Ex.2 pyrene", "c1cc2ccc3cccc4ccc(c1)c2c34",
     (652.43, 652.43, 673.96)),
    ("Ex.3 4-aminobutanol", "NCCCCO", (448.64, 448.64, 478.42)),
]
APPENDIX_TOL = 0.02  # 論文の丸め誤差 [K]

# 論文 Table 6-8 の寄与値から手計算した、グループ分解と寄与値参照の回帰テスト。
# (SMILES, 最終近似の寄与和 sum, 手計算式の説明)
HAND_CASES: List[Tuple[str, float, str]] = [
    ("CCO", 0.8491 + 0.7141 + 2.5670, "CH3 + CH2 + OH"),
    ("CCCO", 0.8491 + 2 * 0.7141 + 2.5670, "CH3 + 2 CH2 + OH"),
    ("CCC(C)=O", 0.8491 + 0.7141 + 3.1178 - 0.0215,
     "CH3 + CH2 + CH3CO + 第2次 CH3COCH2"),
    ("CCc1ccccc1", 0.8491 + 5 * 0.8365 + 1.4925, "CH3 + 5 aCH + aC-CH2"),
    ("c1ccc2ccccc2c1", 8 * 0.8365 + 2 * 1.7324 + 2 * 0.0441,
     "8 aCH + 2 aC(fused) + 第3次 2 AROMFUSED[2]"),
]
HAND_TOL = 0.01  # [K] (表の4桁丸めによる差)


def _self_test() -> List[str]:
    """推算器の動作を検証し、失敗メッセージのリストを返す (空なら成功)。"""
    failures: List[str] = []

    for label, smiles, expected in APPENDIX_B_TB:
        res = estimate_boiling_point(smiles)
        if res.status == STATUS_ERROR or res.tb_k is None:
            failures.append(f"{label}: {res.messages}")
            continue
        for order, want in zip((1, 2, 3), expected):
            got = res.tb_by_level[order]
            if abs(got - want) > APPENDIX_TOL:
                failures.append(
                    f"{label}: level{order} {got:.2f} != {want:.2f}"
                )

    for smiles, total, label in HAND_CASES:
        res = estimate_boiling_point(smiles)
        want = TB0 * math.log(total)
        if res.tb_k is None:
            failures.append(f"hand {smiles}: {res.messages}")
        elif abs(res.tb_k - want) > HAND_TOL:
            failures.append(
                f"hand {smiles} ({label}): {res.tb_k:.2f} != {want:.2f}"
            )

    # 異常系: 例外を投げず error になること
    for bad in ("not a smiles", "", None):
        if estimate_boiling_point(bad).status != STATUS_ERROR:
            failures.append(f"invalid input should be error: {bad!r}")
    # 複数成分は既定で error、最大成分指定なら推算できる
    salt = "CCO.[Na+].[Cl-]"
    if estimate_boiling_point(salt).status != STATUS_ERROR:
        failures.append("multi-component should be error by default")
    kept = estimate_boiling_point(salt, keep_largest_fragment=True)
    ethanol = estimate_boiling_point("CCO")
    if kept.tb_k is None or ethanol.tb_k is None or (
        abs(kept.tb_k - ethanol.tb_k) > 1e-9
    ):
        failures.append("keep_largest_fragment should equal ethanol")
    # 重複SMILESと順序
    many = estimate_many(["CCO", "c1ccccc1", "CCO"])
    if many[0] is not many[2] or len(many) != 3:
        failures.append("estimate_many cache/order")
    return failures


# ---------------------------------------------------------------------------
# コマンドライン
# ---------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    """コマンドライン引数のパーサを作る。"""
    parser = argparse.ArgumentParser(
        description="Marrero-Gani法で、SMILESから正常沸点 [K] を推算する。"
    )
    parser.add_argument("smiles", nargs="*", help="SMILES (複数可)")
    parser.add_argument("--csv", help="入力CSV (SMILES列を持つ)")
    parser.add_argument("-o", "--output", help="CSV入力時の出力先")
    parser.add_argument("--smiles-col", default="smiles",
                        help="CSVのSMILES列名")
    parser.add_argument("--table-dir",
                        help="Table_6/7/8.csv のあるディレクトリ")
    parser.add_argument("--detail", action="store_true",
                        help="グループ内訳を表示")
    parser.add_argument("--allow-incomplete", action="store_true",
                        help="未割当原子があっても推算する (信頼できない)")
    parser.add_argument("--largest-fragment", action="store_true",
                        help="複数成分は最大成分だけを使う")
    parser.add_argument("--selftest", action="store_true",
                        help="自己テストを実行")
    return parser


def _run_csv(args: argparse.Namespace, options: Dict[str, object]) -> int:
    """CSV一括処理。入力列に推算結果の列を追加して出力する。"""
    if not args.output:
        print("--csv には -o/--output が必要です", file=sys.stderr)
        return 2
    with open(args.csv, encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh)
        rows = list(reader)
        fieldnames = list(reader.fieldnames or [])
    if args.smiles_col not in fieldnames:
        print(f"列がない: {args.smiles_col!r}", file=sys.stderr)
        return 2
    results = estimate_many([r[args.smiles_col] for r in rows], **options)
    extra = ["tb_k", "tb_c", "tb_1st_k", "tb_2nd_k", "tb_3rd_k",
             "tb_status", "tb_messages"]
    with open(args.output, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames + extra)
        writer.writeheader()
        for row, res in zip(rows, results):
            lv = res.tb_by_level
            row.update({
                "tb_k": "" if res.tb_k is None else f"{res.tb_k:.2f}",
                "tb_c": "" if res.tb_c is None else f"{res.tb_c:.2f}",
                "tb_1st_k": f"{lv[1]:.2f}" if 1 in lv else "",
                "tb_2nd_k": f"{lv[2]:.2f}" if 2 in lv else "",
                "tb_3rd_k": f"{lv[3]:.2f}" if 3 in lv else "",
                "tb_status": res.status,
                "tb_messages": " | ".join(res.messages),
            })
            writer.writerow(row)
    counts: Dict[str, int] = {}
    for res in results:
        counts[res.status] = counts.get(res.status, 0) + 1
    print(f"{len(rows)} compounds -> {args.output}  status: {counts}",
          file=sys.stderr)
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    """コマンドラインの入口。終了コードを返す。"""
    args = build_parser().parse_args(argv)
    if args.selftest:
        failures = _self_test()
        for failure in failures:
            print(f"[FAIL] {failure}", file=sys.stderr)
        print(f"selftest: {'OK' if not failures else 'FAILED'}")
        return 1 if failures else 0

    options: Dict[str, object] = {
        "table_dir": args.table_dir,
        "allow_incomplete": args.allow_incomplete,
        "keep_largest_fragment": args.largest_fragment,
    }
    if args.csv:
        return _run_csv(args, options)

    smiles_list = list(args.smiles)
    if not smiles_list:
        text = input("SMILES: ").strip()
        smiles_list = [text]
    exit_code = 0
    for smiles in smiles_list:
        result = estimate_boiling_point(smiles, **options)
        print(format_result(result, detail=args.detail))
        print()
        if result.status == STATUS_ERROR:
            exit_code = 1
    return exit_code


if __name__ == "__main__":
    sys.exit(main())
