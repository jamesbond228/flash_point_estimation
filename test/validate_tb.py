"""実測沸点CSV (Book1.csv) に Marrero-Gani法の推算値を追加し、精度を検証する。

入力CSVの1列目は化合物名、'BoilingPont' 列は実測沸点 [degC]。
化合物名 -> SMILES の対応は本ファイルの NAME_TO_SMILES (人手で作成) を使う。
SMILESは自動変換ではないため、必ず出力CSVの smiles 列を確認すること。

使い方:
    python validate_tb.py Book1.csv -o Book1_with_pred.csv

必要なファイル (同じディレクトリ):
    marrero_gani_tb.py, marrero_gani_groups.py, marrero_gani_batch.py,
    Table_6.csv, Table_7.csv, Table_8.csv
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from rdkit import Chem

import marrero_gani_tb as tb

NAME_COL = ""              # 1列目 (見出しが空)
EXP_COL = "BoilingPont"    # 実測沸点 [degC] (原ファイルの綴りのまま)
KELVIN_OFFSET = 273.15

# 化合物名 (CSVの表記) -> (SMILES, 評価から除外する理由。除外しないなら "")
# SMILES が None の行は推算しない。
NAME_TO_SMILES: Dict[str, Tuple[Optional[str], str]] = {
    "Acetaldehyde": ("CC=O", ""),
    "Acetic Acid": ("CC(=O)O", ""),
    "Acetic Anhydride": ("CC(=O)OC(C)=O", ""),
    "Acetone": ("CC(C)=O", ""),
    "Acetonitrile": ("CC#N", ""),
    "Acetophenone": ("CC(=O)c1ccccc1", ""),
    "Acetylacetone": ("CC(=O)CC(C)=O", ""),
    "Allyl Alcohol": ("C=CCO", ""),
    "Ethanolamine": ("NCCO", ""),
    "Aniline": ("Nc1ccccc1", ""),
    "Benzaldehyde": ("O=Cc1ccccc1", ""),
    "Benzene": ("c1ccccc1", ""),
    "Benzyl Alcohol": ("OCc1ccccc1", ""),
    "Benzyl Benzoate": ("O=C(OCc1ccccc1)c1ccccc1", ""),
    "1-Butanol": ("CCCCO", ""),
    "2-Butanol": ("CCC(C)O", ""),
    "Methyl Ethyl Ketone": ("CCC(C)=O", ""),
    "Butyl Acetate": ("CCCCOC(C)=O", ""),
    "tert-Butyl Alcohol": ("CC(C)(C)O", ""),
    "Dibutyl ether": ("CCCCOCCCC", ""),
    "Carbon Disulfide": ("S=C=S", ""),
    "Chloroform": ("ClC(Cl)Cl", ""),
    "Epichlorohydrin": ("ClCC1CO1", ""),
    "o-Cresol": ("Cc1ccccc1O", ""),
    "m-Cresol": ("Cc1cccc(O)c1", ""),
    "p-Cresol": ("Cc1ccc(O)cc1", ""),
    "Cyclohexane": ("C1CCCCC1", ""),
    "Cyclohexanol": ("OC1CCCCC1", ""),
    "1,2-Dichloroethane": ("ClCCCl", ""),
    "Dichloromethane": ("ClCCl", ""),
    "Diethyl Carbonate": ("CCOC(=O)OCC", ""),
    "Diethylene Glycol": ("OCCOCCO", ""),
    "Diethylene Glycol Monobutyl Ether": ("CCCCOCCOCCO", ""),
    "Butyl carbitol acetate": ("CCCCOCCOCCOC(C)=O", ""),
    "Diethylene Glycol Monoethyl Ether": ("CCOCCOCCO", ""),
    "Diethylene Glycol Monoethyl Ether Acetate": ("CCOCCOCCOC(C)=O", ""),
    "Diethylene Glycol Monomethyl Ether": ("COCCOCCO", ""),
    "Diethyl Ether": ("CCOCC", ""),
    "Dimethylacetamide": ("CC(=O)N(C)C", ""),
    "N,N-Dimethylformamide": ("CN(C)C=O", ""),
    "Dimethyl Sulfoxide": ("CS(C)=O", ""),
    "1,4-Dioxane": ("C1COCCO1", ""),
    "Ethanol(95)": ("CCO", "95%水溶液 (混合物)。純エタノールとして推算し評価から除外"),
    "Absolute Ethanol (99.5)": ("CCO", ""),
    "Ethyl Acetate": ("CCOC(C)=O", ""),
    "Ethyl Benzoate": ("CCOC(=O)c1ccccc1", ""),
    "2-Chloroethanol": ("OCCCl", ""),
    "Ethylene Glycol": ("OCCO", ""),
    "1,2-Dimethoxyethane": ("COCCOC", ""),
    "2-Butoxyethanol": ("CCCCOCCO", ""),
    "Ethylene Glycol Monobutyl Ether Acetate": ("CCCCOCCOC(C)=O", ""),
    "2-Ethoxyethanol": ("CCOCCO", ""),
    "Ethylene Glycol Monoethyl Ether Acetate": ("CCOCCOC(C)=O", ""),
    "Ethylene Glycol Monomethyl Ether Acetate": ("COCCOC(C)=O", ""),
    "Ethyl Formate": ("CCOC=O", ""),
    "2-Ethylhexanol": ("CCCCC(CC)CO", ""),
    "Ethyl Acetoacetate": ("CCOC(=O)CC(C)=O", ""),
    "Ethyl Propionate": ("CCOC(=O)CC", ""),
    "Formamide": ("NC=O", ""),
    "Formic Acid 98%(T)以上": ("OC=O", ""),
    "Furfuryl Alcohol": ("OCc1ccco1", ""),
    "glycerine": ("OCC(O)CO", ""),
    "Heptane": ("CCCCCCC", ""),
    "Hexane": ("CCCCCC", ""),
    "1-Hexanol": ("CCCCCCO", ""),
    "Ligroin": (None, "混合物 (石油留分)"),
    "2,6-Lutidine": ("Cc1cccc(C)n1", ""),
    "Methanol": ("CO", ""),
    "2-Methoxyethanol": ("COCCO", ""),
    "Methyl Acetate": ("COC(C)=O", ""),
    "2-Methyl-2-butanol": ("CCC(C)(C)O", ""),
    "3-Methyl-1-butanol": ("CC(C)CCO", ""),
    "Isoamyl Acetate": ("CC(C)CCOC(C)=O", ""),
    "Methyl Propionate": ("CCC(=O)OC", ""),
    "Triethanolamine": ("OCCN(CCO)CCO", ""),
    "Nitrobenzene": ("O=[N+]([O-])c1ccccc1", ""),
    "Nitromethane": ("C[N+](=O)[O-]", ""),
    "n-Octane": ("CCCCCCCC", ""),
    "1-Octanol": ("CCCCCCCCO", ""),
    "2-Octanol": ("CCCCCCC(C)O", ""),
    "Pentane": ("CCCCC", ""),
    "1-Pentanol": ("CCCCCO", ""),
    "3-Pentanol": ("CCC(O)CC", ""),
    "n-Pentyl Acetate": ("CCCCCOC(C)=O", ""),
    "Petroleum Benzine(留分50-80℃)": (None, "混合物 (石油留分)"),
    "Petroleum Benzine(留分30-70℃)": (None, "混合物 (石油留分)"),
    "Petroleum Benzine(留分60-70℃)": (None, "混合物 (石油留分)"),
    "Phenol": ("Oc1ccccc1", ""),
    "1-Propanol": ("CCCO", ""),
    "n-Propyl Acetate": ("CCCOC(C)=O", ""),
    "Propylene Glycol": ("CC(O)CO", ""),
    "Propylene Oxide": ("CC1CO1", ""),
    "n-Propyl Ether": ("CCCOCCC", ""),
    "Pyridine": ("c1ccncc1", ""),
    "1,1,2,2-Tetrachloroethane": ("ClC(Cl)C(Cl)Cl", ""),
    "Tetrachloroethylene": ("ClC(Cl)=C(Cl)Cl", ""),
    "Tetrahydrofuran": ("C1CCOC1", ""),
    "Tetralin": ("c1ccc2c(c1)CCCC2", ""),
    "Toluene": ("Cc1ccccc1", ""),
    "1,1,2-Trichloroethane": ("ClCC(Cl)Cl", ""),
    "1,1,2-TCA": (None, "略称のみの行 (値なし。直前行の別名)"),
    "Trichloroethylene": ("ClC=C(Cl)Cl", ""),
    "Triethylamine": ("CCN(CC)CC", ""),
    "TEA": (None, "略称のみの行 (値なし。直前行の別名)"),
    "Trifluoroacetic Acid": ("OC(=O)C(F)(F)F", ""),
    "TFA": (None, "略称のみの行 (値なし。直前行の別名)"),
    "Xylene": (None, "異性体混合物"),
    "o-Xylene": ("Cc1ccccc1C", ""),
    "m-Xylene": ("Cc1cccc(C)c1", ""),
    "p-Xylene": ("Cc1ccc(C)cc1", ""),
}

_NUM = r"(\d+(?:\.\d+)?)"
_RANGE = re.compile(rf"^\s*{_NUM}\s*[–\-~]\s*{_NUM}\s*$")
_APPROX = re.compile(rf"^\s*約\s*{_NUM}\s*$")
_PLAIN = re.compile(rf"^\s*{_NUM}\s*(?:℃)?\s*$")


def parse_exp(text: str) -> Tuple[Optional[float], str]:
    """実測沸点の文字列を数値 [degC] にする。

    Returns:
        (値, 区分)。区分は "single", "range" (中点を採用), "approx" (約)、
        数値化できなければ (None, "")。
    """
    match = _PLAIN.match(text)
    if match:
        return float(match.group(1)), "single"
    match = _RANGE.match(text)
    if match:
        low, high = float(match.group(1)), float(match.group(2))
        return (low + high) / 2.0, "range"
    match = _APPROX.match(text)
    if match:
        return float(match.group(1)), "approx"
    return None, ""


def n_carbon(smiles: str) -> int:
    """SMILESの炭素原子数。"""
    mol = Chem.MolFromSmiles(smiles)
    return sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "C")


def build_table(rows: List[Dict[str, str]]) -> pd.DataFrame:
    """各行に推算値・誤差・評価用フラグを付けたDataFrameを作る。"""
    records = []
    for row in rows:
        name = (row.get(NAME_COL) or "").strip()
        exp_c, kind = parse_exp(row.get(EXP_COL) or "")
        smiles, reason = NAME_TO_SMILES.get(name, (None, ""))
        rec: Dict[str, object] = {
            "smiles": smiles or "",
            "tb_exp_c": exp_c,
            "tb_exp_kind": kind,
            "tb_exp_k": None if exp_c is None else exp_c + KELVIN_OFFSET,
            "n_carbon": None,
            "tb_pred_k": None, "tb_pred_c": None,
            "tb_1st_k": None, "tb_2nd_k": None,
            "err_k": None, "err_pct": None,
            "tb_status": "", "tb_messages": "",
            "eval_use": 0, "exclude_reason": reason,
        }
        if not name:
            rec["exclude_reason"] = "名前なし (値の継ぎ足し行)"
        elif name not in NAME_TO_SMILES:
            rec["exclude_reason"] = "SMILES未登録"
        elif smiles is None:
            rec["tb_status"] = "skipped"
        else:
            res = tb.estimate_boiling_point(smiles)
            rec["tb_status"] = res.status
            rec["tb_messages"] = " | ".join(res.messages)
            rec["n_carbon"] = n_carbon(smiles)
            if res.tb_k is not None:
                rec["tb_pred_k"] = res.tb_k
                rec["tb_pred_c"] = res.tb_c
                rec["tb_1st_k"] = res.tb_by_level[1]
                rec["tb_2nd_k"] = res.tb_by_level[2]
            if res.tb_k is not None and exp_c is not None:
                rec["err_k"] = res.tb_k - rec["tb_exp_k"]
                rec["err_pct"] = rec["err_k"] / rec["tb_exp_k"] * 100.0
                if not reason:
                    rec["eval_use"] = 1
        records.append(rec)
    return pd.DataFrame(records)


def _metrics(
    df: pd.DataFrame, pred_col: str = "tb_pred_k"
) -> Dict[str, float]:
    """誤差指標 (単位 K)。ARE は実測 [K] に対する相対誤差の平均 (%)。"""
    err = df[pred_col] - df["tb_exp_k"]
    ss_res = float((err ** 2).sum())
    ss_tot = float(((df["tb_exp_k"] - df["tb_exp_k"].mean()) ** 2).sum())
    return {
        "N": int(len(df)),
        "MAE": float(err.abs().mean()),
        "RMSE": float(np.sqrt((err ** 2).mean())),
        "bias": float(err.mean()),
        "ARE%": float((err.abs() / df["tb_exp_k"]).mean() * 100.0),
        "R2": 1.0 - ss_res / ss_tot,
    }


def report(table: pd.DataFrame, name_col: pd.Series) -> str:
    """検証結果の要約文字列を作る。"""
    ev = table[table["eval_use"] == 1].copy()
    ev["name"] = name_col[ev.index]
    lines: List[str] = []

    def line(label: str, df: pd.DataFrame, col: str = "tb_pred_k") -> None:
        if df.empty:
            lines.append(f"{label}: N=0")
            return
        m = _metrics(df, col)
        lines.append(
            f"{label:<34s} N={m['N']:<3d} MAE={m['MAE']:6.2f} K  "
            f"RMSE={m['RMSE']:6.2f} K  bias={m['bias']:+6.2f} K  "
            f"ARE={m['ARE%']:5.2f}%  R2={m['R2']:.4f}"
        )

    lines.append("status 内訳: " + str(
        table["tb_status"].replace("", "none").value_counts().to_dict()))
    lines.append(f"評価対象 (推算成功・実測あり・除外なし): {len(ev)} 件")
    lines.append("")
    line("全評価対象 (第3近似)", ev)
    line("  第1近似", ev, "tb_1st_k")
    line("  第2近似", ev, "tb_2nd_k")
    line("status=ok のみ", ev[ev["tb_status"] == "ok"])
    line("実測が単一値のみ", ev[ev["tb_exp_kind"] == "single"])
    strict = (ev["tb_exp_kind"] == "single") & (ev["tb_status"] == "ok")
    line("実測が単一値 & status=ok", ev[strict])
    line("C>=3 (論文の適用範囲)", ev[ev["n_carbon"] >= 3])
    line("C<3", ev[ev["n_carbon"] < 3])
    lines.append("")
    lines.append("誤差の大きい順 (上位12, 全評価対象):")
    worst = ev.reindex(ev["err_k"].abs().sort_values(ascending=False).index)
    for _, r in worst.head(12).iterrows():
        lines.append(
            f"  {r['name'][:42]:<42s} exp={r['tb_exp_k']:7.2f}  "
            f"pred={r['tb_pred_k']:7.2f}  err={r['err_k']:+7.2f} K "
            f"({r['tb_status']})"
        )
    errors = table[table["tb_status"] == "error"]
    if len(errors):
        lines.append("")
        lines.append("推算不能 (error):")
        for idx, r in errors.iterrows():
            lines.append(f"  {name_col[idx]}: {r['tb_messages']}")
    return "\n".join(lines)


def main(argv: Optional[Sequence[str]] = None) -> int:
    """コマンドラインの入口。"""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("input", help="実測沸点CSV")
    parser.add_argument("-o", "--output", required=True, help="出力CSV")
    args = parser.parse_args(argv)

    with open(args.input, encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh)
        rows = list(reader)
        fields = list(reader.fieldnames or [])
    if EXP_COL not in fields:
        print(f"列がない: {EXP_COL!r}", file=sys.stderr)
        return 2

    table = build_table(rows)
    base = pd.DataFrame(rows, columns=fields)
    out = pd.concat([base, table], axis=1)
    out.to_csv(args.output, index=False, encoding="utf-8-sig")
    print(report(table, base[NAME_COL].fillna("")))
    return 0


if __name__ == "__main__":
    sys.exit(main())
