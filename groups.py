"""Marrero-Gani法 (Marrero & Gani, Fluid Phase Equilib. 183-184, 2001) の
第1次 (Table 6)、第2次 (Table 7)、第3次 (Table 8) の各グループを
SMILESから自動カウントするモジュール。

主なAPI:
    analyze_smiles / count_first_order_groups         第1次 (完全かつ互いに素な分割)
    analyze_second_order / count_second_order_groups  第2次 (重なりを許す)
    analyze_third_order / count_third_order_groups    第3次 (重なりを許す)

設計原則 (論文 Section 2.1 の基準に対応):
    P1 完全かつ互いに素な分割 (基準5): 各原子は高々1つのグループにのみ
       割り当てる。割り当てられない原子は推測で埋めず、未割当として報告する。
    P2 重量規則 + 例外2つ: 同一原子団が複数グループに該当する場合は重い
       グループを優先する (重量 = 構成原子の原子量 + 結合H)。例外は
       芳香族置換基 (aC-R 型) のみを優先度で表現し、アミド・尿素の専用
       グループはパターンの構造条件で表現する。
    P3 環/非環の分離 (基準4): 環内の原子は環状グループだけで記述する。
       非環状の多原子グループは !R を要求し、環内原子を取り込まない。
    P4 芳香族は論文の運用定義に従う: aCH / aC / aN で表せる環 (6員、C/Nのみ、
       環外にヘテロ原子への二重結合なし) だけを芳香環とし、それ以外の
       RDKit芳香環はKekule構造 (通常の環) として扱う。
    P5 Table 6 にないグループは作らない: 寄与値が存在しないため。

SMARTSの規約:
    * アトムマップ番号 (":n") を付けた原子がそのグループの構成原子 (所有原子)。
    * マップ番号のない原子は文脈 (マッチ条件) としてのみ使われ、所有されない。

拡張方法:
    GROUP_DEFS に GroupDef を1行追加するだけでよい。
    SMARTSでは表現できない条件は predicate (Python関数) で与える。

--------------------------------------------------------------------------
判断記録 (Decision log): 論文に記載がなく、実装者が判断した事項
--------------------------------------------------------------------------
D1 [未解決] aC-R 型とアミド専用グループの優先順位
    論文は2つの例外の相互順序を述べていない。本実装は、論文が aC-R を
    「必ず使う」と表現していることから aC-R を最優先にした。副作用として、
    芳香環に結合した CH2 を持つ第三級アミド (例: N-ベンジル-N-メチル
    アセトアミド) では aC-CH2 が CH2 を取り、残る N-C(=O) に専用グループが
    ないため N, C, O が未割当になる (LIMIT_CASES に記録)。
    代替案: アミド専用 > aC-R にすると、この例は完結する一方、ベンジル位の
    CH2 が aC-CH2 にならず「aC-R 必須」の文言から外れる。
D2 [未解決] 完全分割を最大化する探索は行わない
    割り当ては重い候補からの貪欲法である。未割当の最小化を目的とする
    最適化 (厳密被覆) は、重量規則から外れる結果を生みうるため採用しない。
D3 [解釈] 論文に定義文がないグループの判定規則
    * No.16/17: 芳香環に含まれる回数で判定 (芳香環2つ以上 = 16、芳香環1つ
      + 非芳香環 = 17)。
    * No.165-167 (エポキシド): 環炭素のH合計で区別する (3 = C2H3O、
      2 = C2H2O、1以下 = C2O)。C2O の条件は論文の例 (トリメチルオキシラン、
      H=1) から逆算した。無置換エチレンオキシド (H=4) は専用グループが
      なく、CH2 (cyclic) x2 + O (cyclic) になる。
    * No.82/83, 131-137, 164 など: 名称と例示化合物からの推定。
    * CCl4 / CHCl3 / CH2Cl2: 専用グループがなく、部分群 + "Cl except as
      above" に分割される。論文の扱いは確認できない。
D4 [外挿] 芳香環の運用定義 (P4) は論文の例からの推定である。
    根拠: Ex.4 (クマリンのピラノン環は非芳香族)、Table 6 の例 (イミダゾール =
    CH=N (cyclic)、チオフェン = S (cyclic))、aN (ピリジン系)。4-ピラノン、
    ウラシル型 (環外C=O) などは Ex.4 からの外挿であり、論文に明示はない。
    ピリジニウム、N-オキシド (電荷を持つ環N) は対象外で未割当になる。
D5 [解釈] 環内カルボニル/環内Nの分解 (A1)
    環状エステル/アミド/カーボネート/イミドは、環原子単位のグループに分解する
    (Ex.4: ラクトン = CO (cyclic) + O (cyclic)、Ex.5: ラクタム = N (cyclic) +
    CO (cyclic))。環外にカルボニルを持つ環状アミン (N-アシルピロリジン等) は
    該当グループがなく、C=O が未割当になる (LIMIT_CASES に記録)。
D6 [解釈] OH は C または S に結合した水酸基 (A3)
    根拠: Table 7 の aC-SO2-OH。N-OH, O-OH, P-OH は専用グループ (オキシム、
    ペルオキシド、リン酸) が所有し、専用グループのない構造は未割当になる。
D7 [仕様] 同優先度・同重量の競合は "except as above" 系を後にし、さらに
    グループ番号で決定する (定義順に依存しない)。別グループ同士が同優先度・
    同重量で原子を奪い合った場合は GroupAnalysis.ties に記録する。
D8 [仕様] 論文 Appendix B の誤植 (回帰テストで補正済み):
    Ex.5 の本文は ln(25.1421) と記載されているが、469.19 K は ln(24.0959)
    の値である。Ex.6 の本文は 222.543 と記載されているが、303.66 K は
    Tm0 = 147.450 の値である。
D9 [適用範囲外] Table 6 に対応グループがなく、未割当になる構造:
    無置換尿素、ホルムアミド、アシルハライド、チオエステル、ジスルフィド、
    ピリジンN-オキシド、N に CH/C だけが付く第三級アミン (CHN, CN 系なし)。

--- 第2次グループ (Table 7) の判断記録 ---
D10 [仕様] 完全重複の禁止 (論文 2.2 基準4, Fig.2)
    各出現の原子集合 (extent) が別グループの原子集合の真部分集合なら除外する
    (l-アラニン: CHm-COOH を除外し CHm(NHn)-COOH を採用)。部分的な重なりは
    許す。原子集合が完全に一致する別グループは両方残し、ties に記録する。
D11 [解釈] 出現回数の数え方 (identity と extent の分離)
    マップ番号 1-99 を出現の識別子、100以上を原子集合のみの原子とする。
    分岐アルカン (No.1-5) は中心炭素のペアを数え、メチル基の選び方では数えない。
    根拠: 2,2,3,4,4-ペンタメチルペンタン = No.4 が2、N,N-ジメチルシクロ
    ヘキサンアミン = No.87 が1 (表の出現回数と一致することをテストで確認)。
D12 [解釈] 名称に cyc / aC が付かない第2次グループは非環状の炭素にのみ適用する。
    環内 sp3 炭素は CHcyc / Ccyc 系を使う。第1次の環/非環分離 (基準4) からの類推。
D13 [解釈] 範囲表記 (n in 1..2) は H を1個以上要求し、(n in 0..2) は制約なし。
    No.29 は表が (n in 1..2) だが、例示化合物 2-ヒドロキシイソ酪酸は炭素に H が
    ないため矛盾する。例示化合物を優先して制約なしとした。
D14 [解釈] AROMRINGs / PYRIDINEs は縮合していない6員環のみ対象とし、環外に
    重原子が結合した位置の集合を回転・鏡映で正規化して判定する。縮合環は
    対象外 (第3次の領分)。根拠: Ex.5 は AROMRINGs1s2 が1つだけで、縮合した
    ベンゼン環は数えられていない。置換数が定義にない環 (1個、5-6個) は該当なし。
D15 [未解決] No.24 ジエチレントリアミン: 表の出現回数は (1) だが、構造上は
    NH-CH2CH2-NH2 が2か所ある。構造どおり2と数える (テストは2を期待)。
D16 [解釈] No.22, 26, 51, 86, 87 の窒素はアミン窒素 (アミド・ニトロを除く)。
    No.22 は OCR の文字化けから CHm(OH)CHn(NHp) と解釈した (本文 2.2 の例)。
D17 [解釈] No.6-9 の C=C は非環状のみ (環内 C=C は No.69-75 の担当)。
    No.6 は共役 C=C-C=C (非環状)。
D18 [未解決] (CHn=C)cyc 系 (No.69-75) は、C=C を含む環が芳香環と縮合している
    場合は適用しない。根拠は Ex.4 (ジクマロールは「第2次グループなし」だが、
    縮合したピラノン環の C=C に CH2 が結合している) のみで、他の例からの確認は
    できていない。Ex.5 の Ncyc-CH3 / CHcyc-OH は縮合環でも適用される。

--- 第3次グループ (Table 8) の判断記録 ---
D19 [解釈] 完全重複の除外は第3次グループの中だけで行い、第2次とは混ぜない
    (論文 2.3 は第2次と同様の規則とし、回帰も段階ごとに独立)。鎖状グループ
    (No.1-13, 18-21, 33) は末端の組ごとに1回数える。鎖は非環状 sp3 炭素 (分岐可)
    で、m は鎖の炭素数。根拠: glyceryl tridodecanoate = No.13 が1 (1,3位の組のみ m>2)。
D20 [解釈] 「(fused rings)」は接続する2原子が同じ環に含まれる (環結合)、
    「(different rings)」は環外の結合または鎖で接続する。根拠: biphenylene の
    aC-aC が2 (4員環内の結合でも異なる芳香環どうしなら数える)。
D21 [解釈] AROMFUSED[k] の k は、他の環と共有する環原子の数。[4a] は非縮合原子が
    対位 (アントラセン型)、[4p] は隣接 (フェナントレン型)。根拠: ピレン (Ex.2) が
    [3] で2、[4p] で2。置換位置 s は非縮合原子を1から番号付けし、鏡映で同一視する。
    表にない置換パターンは基本グループを使い、notes に記録する。
D22 [原典記載・未検証] Ex.5 は AROMFUSED[2]s3 (係数 2.2589) と記載するが、Table 8 に
    s3 はなく、2,7-ジメチルナフタレンの例 (s2 が2) から s3 は s2 と同一視される。
    表の s2 は 0.3378 で 2.2589 と一致しない。このため 2.2589 は表の値と区別した
    明示係数 (ExplicitCoefficient) として保持し、Ex.5 の第3近似 (435.22 K) の
    再現にだけ使う。グループ構成は本実装の分類 (AROMFUSED[2]s2) と一致する。
D23 [原典内の不整合] No.17 は表が (n in 0..1) だが、例示のテトラリン・インダンは
    CH2 を持つ。範囲 (H0, H1) はそのまま保持し、芳香環と縮合した5-6員炭素環の
    CH2 だけを THIRD_ORDER_EXCEPTIONS で個別管理する例外として数える。
    例外を False にすると表の範囲どおりになる。適用時は notes に記録する。
D24 [解釈] No.40 (P) / No.41 (S) は、P または S に芳香族炭素 (直接または O 経由) が
    つく組を数える。根拠: triphenyl phosphate = 3 (3つのフェニルの組)。
    No.41 は S に O が1個以上つくもの、No.39 (aC-S-aC) は O なしの S。
D25 [解釈] PYRIDINE.FUSED: [2] は N が縮合炭素に隣接 (キノリン型)、[2-iso] は隣接
    しない (イソキノリン型)、[4] は N と C が対位 (アクリジン型)。表に置換体がない
    ため置換基は区別しない。フェナントリジン型などは未分類として notes に記録する。
D26 [実装上の仮定] No.8 の例示 "Butoxypropanol" は構造が特定できないため、範囲 (m>2)
    と整合する 3-ブトキシ-1-プロパノールと仮定する。エーテルの反対側は CHp
    (p in 0..2) なので、メチルエーテルは対象外とする。1-ブトキシ-2-プロパノール
    (m=2) は該当しない。
D27 [原典記載・構造解釈未確定] No.58 のフェナレンは原典の出現回数 3 を保持する。
    1H-フェナレンは CH2 を含む環が非芳香族のため、本実装は2を返す。原典の値は
    変更せず、不一致は THIRD_ORDER_UNRESOLVED として報告する (失敗とは扱わない)。
"""
from __future__ import annotations

import math
import random
import sys
from collections import Counter
from dataclasses import dataclass, field
from itertools import combinations
from typing import Callable, Dict, List, Optional, Sequence, Tuple

from rdkit import Chem

# predicate(mol, 所有原子インデックス) -> bool
Predicate = Callable[[Chem.Mol, Tuple[int, ...]], bool]

H_MASS = 1.008  # 水素の原子量

# 優先度 (大きいほど先に割り当てる)。論文の例外規則のうち、芳香族置換基
# (aC-R 型) のみを優先度で表現する。アミド・尿素の専用グループは、パターンが
# カルボニル炭素に炭素隣接を要求する (尿素は構造的に競合しない) ため、
# 優先度を必要としない (感度実験で優先度の有無により結果が変わらないことを
# 確認済み)。aC-R とアミド専用グループの相互順序は論文に記載がない (D1)。
P_DEFAULT = 0  # 重量規則のみ
P_AROM = 1     # 芳香族置換基 aC-R (論文の例外規則)

# アミン窒素から除外する条件 (アミド・チオアミド・アミジン型N)
_AMINE = "NX3;!$(N[CX3]=[O,S,N])"
# アミド炭素 (尿素・カルバメート炭素は除外、非環状のみ)
_AMIDE_C = "CX3;!R;!$(C(=O)([#7,#8,#16])[#7,#8,#16])"
# ニトロ基 (RDKitの電荷分離表現)
_NITRO_N = "NX3+"
# チオエーテル硫黄 (両隣が炭素の非環状S; ジスルフィドは除外)
_THIOETHER_S = "SX2H0;!R;$(S([#6])[#6])"
# エーテル酸素 (エステル酸素は除外)
_ETHER_O = "OX2H0;!R;!$(O[CX3]=O)"


# ---------------------------------------------------------------------------
# グループ定義
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class GroupDef:
    """第1次グループ定義。

    Attributes:
        no: 論文 Table 6 のグループ番号。
        name: グループ名。
        smarts: アトムマップ付きSMARTS。
        priority: P_AROM なら重量規則より優先 (aC-R 型の例外規則)。
        predicate: SMARTSで表せない追加条件。
        catch_all: "except as above" 系。同優先度・同重量の競合で他グループに
            譲る (論文の「具体的なグループを優先する」趣旨)。
    """

    no: int
    name: str
    smarts: str
    priority: int = P_DEFAULT
    predicate: Optional[Predicate] = None
    catch_all: bool = False


def _ring_counts(mol: Chem.Mol, atom_idx: int) -> Tuple[int, int]:
    """原子を含む (芳香族環の数, 非芳香族環の数) を返す。"""
    n_arom = 0
    n_non = 0
    for ring in mol.GetRingInfo().AtomRings():
        if atom_idx not in ring:
            continue
        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            n_arom += 1
        else:
            n_non += 1
    return n_arom, n_non


def _fused_aromatic(mol: Chem.Mol, owned: Tuple[int, ...]) -> bool:
    """芳香環同士の縮合炭素 (Table 6 No.16)。判定規則は判断記録 D3。"""
    return _ring_counts(mol, owned[0])[0] >= 2


def _fused_nonaromatic(mol: Chem.Mol, owned: Tuple[int, ...]) -> bool:
    """芳香環と非芳香環の縮合炭素 (Table 6 No.17)。"""
    n_arom, n_non = _ring_counts(mol, owned[0])
    return n_arom == 1 and n_non >= 1


def _epoxide_h_count(mol: Chem.Mol, owned: Tuple[int, ...]) -> int:
    """エポキシド環 (C, C, O) の炭素上の水素の合計数 (判断記録 D3)。"""
    return sum(mol.GetAtomWithIdx(i).GetTotalNumHs() for i in owned)


def _epoxide_h3(mol: Chem.Mol, owned: Tuple[int, ...]) -> bool:
    """環炭素のH合計が3 (C2H3O, Table 6 No.165)。"""
    return _epoxide_h_count(mol, owned) == 3


def _epoxide_h2(mol: Chem.Mol, owned: Tuple[int, ...]) -> bool:
    """環炭素のH合計が2 (C2H2O, Table 6 No.166)。"""
    return _epoxide_h_count(mol, owned) == 2


def _epoxide_h01(mol: Chem.Mol, owned: Tuple[int, ...]) -> bool:
    """環炭素のH合計が0または1 (C2O, Table 6 No.167)。"""
    return _epoxide_h_count(mol, owned) <= 1


# 定義順は同重量時のタイブレークに使われる (No.16/17 は No.18 より前)。
GROUP_DEFS: List[GroupDef] = [
    # --- 非環状 脂肪族炭素 ---
    GroupDef(1, "CH3", "[CX4H3;!R:1]"),
    GroupDef(2, "CH2", "[CX4H2;!R:1]"),
    GroupDef(3, "CH", "[CX4H1;!R:1]"),
    GroupDef(4, "C", "[CX4H0;!R:1]"),
    # --- 非環状 二重結合・三重結合 ---
    GroupDef(5, "CH2=CH", "[CX3H2;!R:1]=[CX3H1;!R:2]"),
    GroupDef(6, "CH=CH", "[CX3H1;!R:1]=[CX3H1;!R:2]"),
    GroupDef(7, "CH2=C", "[CX3H2;!R:1]=[CX3H0;!R:2]"),
    GroupDef(8, "CH=C", "[CX3H1;!R:1]=[CX3H0;!R:2]"),
    GroupDef(9, "C=C", "[CX3H0;!R:1]=[CX3H0;!R:2]"),
    GroupDef(13, "CH#C", "[CX2H1;!R:1]#[CX2;!R:2]"),
    GroupDef(14, "C#C", "[CX2H0;!R:1]#[CX2H0;!R:2]"),
    # --- 芳香族炭素・窒素 ---
    GroupDef(15, "aCH", "[cH1:1]"),
    GroupDef(16, "aC fused with aromatic ring", "[c;R:1]",
             predicate=_fused_aromatic),
    GroupDef(17, "aC fused with non-aromatic subring", "[c;R:1]",
             predicate=_fused_nonaromatic),
    GroupDef(18, "aC except as above", "[c;H0:1]",
             catch_all=True),
    GroupDef(19, "aN in aromatic ring", "[nX2;H0:1]"),
    # --- 芳香族置換基 aC-R (例外規則: 重量に関わらず優先) ---
    GroupDef(20, "aC-CH3", "[c:1][CX4H3;!R:2]", priority=P_AROM),
    GroupDef(21, "aC-CH2", "[c:1][CX4H2;!R:2]", priority=P_AROM),
    GroupDef(22, "aC-CH", "[c:1][CX4H1;!R:2]", priority=P_AROM),
    GroupDef(23, "aC-C", "[c:1][CX4H0;!R:2]", priority=P_AROM),
    GroupDef(24, "aC-CH=CH2", "[c:1][CX3H1;!R:2]=[CX3H2;!R:3]", priority=P_AROM),
    GroupDef(25, "aC-CH=CH", "[c:1][CX3H1;!R:2]=[CX3H1;!R:3]", priority=P_AROM),
    GroupDef(26, "aC-C=CH2", "[c:1][CX3H0;!R:2]=[CX3H2;!R:3]", priority=P_AROM),
    GroupDef(27, "aC-C#CH", "[c:1][CX2H0;!R:2]#[CX2H1;!R:3]", priority=P_AROM),
    GroupDef(28, "aC-C#C", "[c:1][CX2H0;!R:2]#[CX2H0;!R:3]", priority=P_AROM),
    GroupDef(30, "aC-OH", "[c:1][OX2H1:2]", priority=P_AROM),
    GroupDef(32, "aC-COOH", "[c:1][CX3:2](=[O:3])[OX2H1:4]", priority=P_AROM),
    GroupDef(37, "aC-CO", "[c:1][CX3;!R:2](=[O:3])[#6]", priority=P_AROM),
    GroupDef(39, "aC-CHO", "[c:1][CX3H1;!R:2]=[O:3]", priority=P_AROM),
    GroupDef(45, "aC-COO", "[c:1][CX3;!R:2](=[O:3])[OX2H0;!R:4][#6]",
             priority=P_AROM),
    GroupDef(46, "aC-OOCH", "[c:1][OX2H0;!R:2][CX3H1;!R:3]=[O:4]",
             priority=P_AROM),
    GroupDef(47, "aC-OOC", "[c:1][OX2H0;!R:2][CX3;!R:3](=[O:4])[#6]",
             priority=P_AROM),
    GroupDef(53, "aC-O", f"[c:1][{_ETHER_O}:2]", priority=P_AROM),
    GroupDef(62, "aC-NH2", "[c:1][NX3H2:2]", priority=P_AROM),
    GroupDef(63, "aC-NH", "[c:1][NX3H1;!R;!$(N[CX3]=[O,S,N]):2]",
             priority=P_AROM),
    GroupDef(64, "aC-N", "[c:1][NX3H0;!R;!$(N[CX3]=[O,S,N]):2]",
             priority=P_AROM),
    GroupDef(71, "aC-CN", "[c:1][CX2:2]#[NX1:3]", priority=P_AROM),
    GroupDef(123, "aC-Cl", "[c:1][Cl:2]", priority=P_AROM),
    GroupDef(124, "aC-F", "[c:1][F:2]", priority=P_AROM),
    GroupDef(125, "aC-I", "[c:1][I:2]", priority=P_AROM),
    GroupDef(126, "aC-Br", "[c:1][Br:2]", priority=P_AROM),
    # --- アルコール・カルボン酸 ---
    GroupDef(29, "OH", "[OX2H1;$(O[#6,#16]):1]"),
    GroupDef(31, "COOH", "[CX3:1](=[O:2])[OX2H1:3]"),
    # --- ケトン・アルデヒド ---
    GroupDef(33, "CH3CO", "[CX4H3;!R:1][CX3;!R:2](=[O:3])[#6]"),
    GroupDef(34, "CH2CO", "[CX4H2;!R:1][CX3;!R:2](=[O:3])[#6]"),
    GroupDef(35, "CHCO", "[CX4H1;!R:1][CX3;!R:2](=[O:3])[#6]"),
    GroupDef(36, "CCO", "[CX4H0;!R:1][CX3;!R:2](=[O:3])[#6]"),
    GroupDef(38, "CHO", "[CX3H1;!R:1](=[O:2])[#6]"),
    # --- エステル ---
    GroupDef(40, "CH3COO",
             "[CX4H3;!R:1][CX3;!R:2](=[O:3])[OX2H0;!R:4][#6]"),
    GroupDef(41, "CH2COO",
             "[CX4H2;!R:1][CX3;!R:2](=[O:3])[OX2H0;!R:4][#6]"),
    GroupDef(42, "CHCOO",
             "[CX4H1;!R:1][CX3;!R:2](=[O:3])[OX2H0;!R:4][#6]"),
    GroupDef(43, "CCOO",
             "[CX4H0;!R:1][CX3;!R:2](=[O:3])[OX2H0;!R:4][#6]"),
    GroupDef(44, "HCOO", "[CX3H1;!R:1](=[O:2])[OX2H0;!R:3][#6]"),
    GroupDef(48, "COO except as above",
             "[CX3;!R:1](=[O:2])[OX2H0;!R:3][#6]",
             catch_all=True),
    # --- エーテル ---
    GroupDef(49, "CH3O", f"[CX4H3;!R:1][{_ETHER_O}:2]"),
    GroupDef(50, "CH2O", f"[CX4H2;!R:1][{_ETHER_O}:2]"),
    GroupDef(51, "CH-O", f"[CX4H1;!R:1][{_ETHER_O}:2]"),
    GroupDef(52, "C-O", f"[CX4H0;!R:1][{_ETHER_O}:2]"),
    # --- アミン ---
    GroupDef(54, "CH2NH2", "[CX4H2;!R:1][NX3H2:2]"),
    GroupDef(55, "CHNH2", "[CX4H1;!R:1][NX3H2:2]"),
    GroupDef(56, "CNH2", "[CX4H0;!R:1][NX3H2:2]"),
    GroupDef(57, "CH3NH", f"[CX4H3;!R:1][{_AMINE};H1;!R:2]"),
    GroupDef(58, "CH2NH", f"[CX4H2;!R:1][{_AMINE};H1;!R:2]"),
    GroupDef(59, "CHNH", f"[CX4H1;!R:1][{_AMINE};H1;!R:2]"),
    GroupDef(60, "CH3N", f"[CX4H3;!R:1][{_AMINE};H0;!R:2]"),
    GroupDef(61, "CH2N", f"[CX4H2;!R:1][{_AMINE};H0;!R:2]"),
    GroupDef(65, "NH2 except as above", "[NX3H2;!$(N[CX3]=[O,S,N]);!$(Nc):1]",
             catch_all=True),
    # --- ニトリル ---
    GroupDef(68, "CH2CN", "[CX4H2;!R:1][CX2:2]#[NX1:3]"),
    GroupDef(69, "CHCN", "[CX4H1;!R:1][CX2:2]#[NX1:3]"),
    GroupDef(70, "CCN", "[CX4H0;!R:1][CX2:2]#[NX1:3]"),
    GroupDef(72, "CN except as above", "[CX2:1]#[NX1:2]",
             catch_all=True),
    # --- 脂肪族ハロゲン ---
    GroupDef(108, "CH2Cl", "[CX4H2;!R:1][Cl:2]"),
    GroupDef(109, "CHCl", "[CX4H1;!R:1][Cl:2]"),
    GroupDef(110, "CCl", "[CX4H0;!R:1][Cl:2]"),
    # --- アミド (脂肪族アシル; 専用グループ) ---
    GroupDef(84, "HCON(CH2)2",
             "[CX3H1;!R:1](=[O:2])[NX3H0;!R:3]([CX4H2;!R:4])[CX4H2;!R:5]"),
    GroupDef(85, "HCONHCH2",
             "[CX3H1;!R:1](=[O:2])[NX3H1;!R:3][CX4H2;!R:4]"),
    GroupDef(86, "CONH2", f"[{_AMIDE_C}:1](=[O:2])([#6])[NX3H2;!R:3]"),
    GroupDef(87, "CONHCH3",
             f"[{_AMIDE_C}:1](=[O:2])([#6])[NX3H1;!R:3][CX4H3;!R:4]"),
    GroupDef(88, "CONHCH2",
             f"[{_AMIDE_C}:1](=[O:2])([#6])[NX3H1;!R:3][CX4H2;!R:4]"),
    GroupDef(89, "CON(CH3)2",
             f"[{_AMIDE_C}:1](=[O:2])([#6])"
             "[NX3H0;!R:3]([CX4H3;!R:4])[CX4H3;!R:5]"),
    GroupDef(90, "CONCH3CH2",
             f"[{_AMIDE_C}:1](=[O:2])([#6])"
             "[NX3H0;!R:3]([CX4H3;!R:4])[CX4H2;!R:5]"),
    GroupDef(91, "CON(CH2)2",
             f"[{_AMIDE_C}:1](=[O:2])([#6])"
             "[NX3H0;!R:3]([CX4H2;!R:4])[CX4H2;!R:5]"),
    GroupDef(92, "CONHCO",
             f"[{_AMIDE_C}:1](=[O:2])([#6])[NX3H1;!R:3]"
             f"[{_AMIDE_C}:4](=[O:5])[#6]"),
    GroupDef(93, "CONCO",
             f"[NX3H0;!R:3]([{_AMIDE_C}:1](=[O:2])[#6])"
             f"([{_AMIDE_C}:4](=[O:5])[#6])"),
    GroupDef(107, "NHCO except as above",
             "[NX3H1;!R:1][CX3;!R;!$(C(=O)([#7,#8,#16])[#7,#8,#16]):2]=[O:3]",
             catch_all=True),
    # --- アミド (芳香族; aC-R 型) ---
    GroupDef(94, "aC-CONH2", "[c:1][CX3;!R:2](=[O:3])[NX3H2;!R:4]",
             priority=P_AROM),
    GroupDef(95, "aC-NH(CO)H", "[c:1][NX3H1;!R:2][CX3H1;!R:3]=[O:4]",
             priority=P_AROM),
    GroupDef(96, "aC-N(CO)H", "[c:1][NX3H0;!R:2][CX3H1;!R:3]=[O:4]",
             priority=P_AROM),
    GroupDef(97, "aC-CONH", "[c:1][CX3;!R:2](=[O:3])[NX3H1;!R:4]",
             priority=P_AROM),
    GroupDef(98, "aC-NHCO",
             "[c:1][NX3H1;!R:2][CX3;!R:3](=[O:4])[#6]", priority=P_AROM),
    GroupDef(99, "aC-NCO (amide)",
             "[c:1][NX3H0;!R:2][CX3;!R:3](=[O:4])[#6]", priority=P_AROM),
    # --- 尿素 ---
    GroupDef(100, "NHCONH",
             "[NX3H1;!R:1][CX3;!R:2](=[O:3])[NX3H1;!R:4]"),
    GroupDef(101, "NH2CONH",
             "[NX3H2;!R:1][CX3;!R:2](=[O:3])[NX3H1;!R:4]"),
    GroupDef(102, "NH2CON",
             "[NX3H2;!R:1][CX3;!R:2](=[O:3])[NX3H0;!R:4]"),
    GroupDef(103, "NHCON",
             "[NX3H1;!R:1][CX3;!R:2](=[O:3])[NX3H0;!R:4]"),
    GroupDef(104, "NCON",
             "[NX3H0;!R:1][CX3;!R:2](=[O:3])[NX3H0;!R:4]"),
    GroupDef(105, "aC-NHCONH2",
             "[c:1][NX3H1;!R:2][CX3;!R:3](=[O:4])[NX3H2;!R:5]",
             priority=P_AROM),
    GroupDef(106, "aC-NHCONH",
             "[c:1][NX3H1;!R:2][CX3;!R:3](=[O:4])[NX3H1;!R:5]",
             priority=P_AROM),
    # --- ニトロ・亜硝酸/硝酸エステル ---
    GroupDef(77, "CH2NO2",
             f"[CX4H2;!R:1][{_NITRO_N}:2](=[O:3])[O-:4]"),
    GroupDef(78, "CHNO2",
             f"[CX4H1;!R:1][{_NITRO_N}:2](=[O:3])[O-:4]"),
    GroupDef(79, "CNO2",
             f"[CX4H0;!R:1][{_NITRO_N}:2](=[O:3])[O-:4]"),
    GroupDef(80, "aC-NO2",
             f"[c:1][{_NITRO_N}:2](=[O:3])[O-:4]", priority=P_AROM),
    GroupDef(81, "NO2 except as above",
             f"[{_NITRO_N}:1](=[O:2])[O-:3]",
             catch_all=True),
    GroupDef(82, "ONO", "[OX2;!R:1][NX2;!R:2]=[OX1:3]"),
    GroupDef(83, "ONO2",
             f"[OX2;!R:1][{_NITRO_N}:2](=[O:3])[O-:4]"),
    # --- 多ハロゲン化物 (Cl) ---
    GroupDef(111, "CHCl2", "[CX4H1;!R:1]([Cl:2])[Cl:3]"),
    GroupDef(112, "CCl2", "[CX4H0;!R:1]([Cl:2])[Cl:3]"),
    GroupDef(113, "CCl3", "[CX4H0;!R:1]([Cl:2])([Cl:3])[Cl:4]"),
    # --- 多ハロゲン化物 (F) ---
    GroupDef(114, "CH2F", "[CX4H2;!R:1][F:2]"),
    GroupDef(115, "CHF", "[CX4H1;!R:1][F:2]"),
    GroupDef(116, "CF", "[CX4H0;!R:1][F:2]"),
    GroupDef(117, "CHF2", "[CX4H1;!R:1]([F:2])[F:3]"),
    GroupDef(118, "CF2", "[CX4H0;!R:1]([F:2])[F:3]"),
    GroupDef(119, "CF3", "[CX4H0;!R:1]([F:2])([F:3])[F:4]"),
    # --- 混合ハロゲン ---
    GroupDef(120, "CCl2F", "[CX4H0;!R:1]([Cl:2])([Cl:3])[F:4]"),
    GroupDef(121, "HCClF", "[CX4H1;!R:1]([Cl:2])[F:3]"),
    GroupDef(122, "CClF2", "[CX4H0;!R:1]([Cl:2])([F:3])[F:4]"),
    # --- ハロゲン (上記以外) ---
    GroupDef(127, "I except as above", "[I:1]",
             catch_all=True),
    GroupDef(128, "Br except as above", "[Br:1]",
             catch_all=True),
    GroupDef(129, "F except as above", "[F:1]",
             catch_all=True),
    GroupDef(130, "Cl except as above", "[Cl:1]",
             catch_all=True),
    # --- 硫黄: チオール・チオエーテル ---
    GroupDef(138, "CH2SH", "[CX4H2;!R:1][SX2H1:2]"),
    GroupDef(139, "CHSH", "[CX4H1;!R:1][SX2H1:2]"),
    GroupDef(140, "CSH", "[CX4H0;!R:1][SX2H1:2]"),
    GroupDef(141, "aC-SH", "[c:1][SX2H1:2]", priority=P_AROM),
    GroupDef(142, "SH except as above", "[SX2H1:1]",
             catch_all=True),
    GroupDef(143, "CH3S", f"[CX4H3;!R:1][{_THIOETHER_S}:2]"),
    GroupDef(144, "CH2S", f"[CX4H2;!R:1][{_THIOETHER_S}:2]"),
    GroupDef(145, "CHS", f"[CX4H1;!R:1][{_THIOETHER_S}:2]"),
    GroupDef(146, "CS", f"[CX4H0;!R:1][{_THIOETHER_S}:2]"),
    GroupDef(147, "aC-S-", f"[c:1][{_THIOETHER_S}:2]", priority=P_AROM),
    # --- 硫黄: スルホキシド・スルホン・硫酸系 ---
    GroupDef(148, "SO", "[SX3;!R:1](=[O:2])([#6])[#6]"),
    GroupDef(149, "SO2", "[SX4;!R:1](=[O:2])(=[O:3])([#6])[#6]"),
    GroupDef(150, "SO3 (sulfite)",
             "[SX3;!R:1](=[O:2])([OX2:3][#6])[OX2:4][#6]"),
    GroupDef(151, "SO3 (sulfonate)",
             "[SX4;!R:1](=[O:2])(=[O:3])([#6])[OX2:4][#6]"),
    GroupDef(152, "SO4 (sulfate)",
             "[SX4;!R:1](=[O:2])(=[O:3])([OX2:4][#6])[OX2:5][#6]"),
    GroupDef(153, "aC-SO", "[c:1][SX3;!R:2](=[O:3])[#6]", priority=P_AROM),
    GroupDef(154, "aC-SO2", "[c:1][SX4;!R:2](=[O:3])=[O:4]",
             priority=P_AROM),
    # --- リン ---
    GroupDef(155, "PH (phosphine)", "[PX3H1;!R:1]([#6])[#6]"),
    GroupDef(156, "P (phosphine)", "[PX3H0;!R:1]([#6])([#6])[#6]"),
    GroupDef(157, "PO3 (phosphite)",
             "[PX3;!R:1]([OX2:2][#6])([OX2:3][#6])[OX2:4][#6]"),
    GroupDef(158, "PHO3 (phosphonate)",
             "[PX4H1;!R:1](=[O:2])([OX2:3][#6])[OX2:4][#6]"),
    GroupDef(159, "PO3 (phosphonate)",
             "[PX4H0;!R:1](=[O:2])([#6])([OX2:3][#6])[OX2:4][#6]"),
    GroupDef(160, "PHO4 (phosphate)",
             "[PX4;!R:1](=[O:2])([OX2H1:3])([OX2:4][#6])[OX2:5][#6]"),
    GroupDef(161, "PO4 (phosphate)",
             "[PX4;!R:1](=[O:2])([OX2:3][#6])([OX2:4][#6])[OX2:5][#6]"),
    GroupDef(162, "aC-PO4",
             "[c:1][OX2:2][PX4;!R:3](=[O:4])([OX2:5][#6])[OX2:6][#6]",
             priority=P_AROM),
    GroupDef(163, "aC-P", "[c:1][PX3;!R:2]", priority=P_AROM),
    # --- 環状 硫黄 ---
    GroupDef(181, "S (cyclic)", "[SX2;R:1]"),
    GroupDef(182, "SO2 (cyclic)", "[SX4;R:1](=[O:2])=[O:3]"),
    # --- アレン ---
    GroupDef(10, "CH2=C=CH", "[CX3H2;!R:1]=[CX2;!R:2]=[CX3H1;!R:3]"),
    GroupDef(11, "CH2=C=C", "[CX3H2;!R:1]=[CX2;!R:2]=[CX3H0;!R:3]"),
    GroupDef(12, "CH=C=CH", "[CX3H1;!R:1]=[CX2;!R:2]=[CX3H1;!R:3]"),
    # --- イミン (非環状) ---
    GroupDef(66, "CH=N", "[CX3H1;!R:1]=[NX2;!R:2]"),
    GroupDef(67, "C=N", "[CX3H0;!R:1]=[NX2;!R:2]"),
    # --- イソシアナート ---
    GroupDef(73, "CH2NCO", "[CX4H2;!R:1][NX2:2]=[CX2:3]=[O:4]"),
    GroupDef(74, "CHNCO", "[CX4H1;!R:1][NX2:2]=[CX2:3]=[O:4]"),
    GroupDef(75, "CNCO", "[CX4H0;!R:1][NX2:2]=[CX2:3]=[O:4]"),
    GroupDef(76, "aC-NCO (isocyanate)", "[c:1][NX2:2]=[CX2:3]=[O:4]", priority=P_AROM),
    # --- オキシム ---
    GroupDef(131, "CHNOH", "[CX3H1;!R:1]=[NX2;!R:2][OX2H1:3]"),
    GroupDef(132, "CNOH", "[CX3H0;!R:1]=[NX2;!R:2][OX2H1:3]"),
    GroupDef(133, "aC-CHNOH", "[c:1][CX3H1;!R:2]=[NX2;!R:3][OX2H1:4]",
             priority=P_AROM),
    # --- グリコールエーテル・ヒドロペルオキシド・カーボネート ---
    GroupDef(134, "OCH2CH2OH",
             f"[{_ETHER_O}:1][CX4H2;!R:2][CX4H2;!R:3][OX2H1:4]"),
    GroupDef(135, "OCHCH2OH",
             f"[{_ETHER_O}:1][CX4H1;!R:2][CX4H2;!R:3][OX2H1:4]"),
    GroupDef(136, "OCH2CHOH",
             f"[{_ETHER_O}:1][CX4H2;!R:2][CX4H1;!R:3][OX2H1:4]"),
    GroupDef(137, "-O-OH", "[OX2;!R:1]([#6])[OX2H1:2]"),
    GroupDef(164, "CO3 (carbonate)",
             "[CX3;!R:1](=[O:2])([OX2;!R:3][#6])[OX2;!R:4][#6]"),
    # --- エポキシド (環炭素のH数で区別) ---
    GroupDef(165, "C2H3O", "[CX4;R:1]1[CX4;R:2][OX2;R:3]1",
             predicate=_epoxide_h3),
    GroupDef(166, "C2H2O", "[CX4;R:1]1[CX4;R:2][OX2;R:3]1",
             predicate=_epoxide_h2),
    GroupDef(167, "C2O", "[CX4;R:1]1[CX4;R:2][OX2;R:3]1",
             predicate=_epoxide_h01),
    # --- 環状イミン ---
    GroupDef(177, "CH=N (cyclic)", "[CX3H1;R:1]=[NX2;R:2]"),
    GroupDef(178, "C=N (cyclic)", "[CX3H0;R:1]=[NX2;R:2]"),
    # --- 環状 (非芳香族) ---
    GroupDef(168, "CH2 (cyclic)", "[CX4H2;R:1]"),
    GroupDef(169, "CH (cyclic)", "[CX4H1;R:1]"),
    GroupDef(170, "C (cyclic)", "[CX4H0;R:1]"),
    GroupDef(171, "CH=CH (cyclic)", "[CX3H1;R:1]=[CX3H1;R:2]"),
    GroupDef(172, "CH=C (cyclic)", "[CX3H1;R:1]=[CX3H0;R:2]"),
    GroupDef(173, "C=C (cyclic)", "[CX3H0;R:1]=[CX3H0;R:2]"),
    GroupDef(174, "CH2=C (cyclic)", "[CX3H2;!R:1]=[CX3H0;R:2]"),
    GroupDef(175, "NH (cyclic)",
             "[NX3H1;R:1]"),
    GroupDef(176, "N (cyclic)",
             "[NX3H0;R:1]"),
    GroupDef(179, "O (cyclic)", "[OX2;R:1]"),
    GroupDef(180, "CO (cyclic)", "[CX3;R:1]=[O:2]"),
]


# ---------------------------------------------------------------------------
# コンパイル済みパターン
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class _Compiled:
    group: GroupDef
    query: Chem.Mol
    owned_query_idx: Tuple[int, ...]


def _compile(defs: Sequence[GroupDef]) -> List[_Compiled]:
    """SMARTSをコンパイルし、所有原子 (マップ番号付き) を検出する。"""
    compiled: List[_Compiled] = []
    names = [g.name for g in defs]
    numbers = [g.no for g in defs]
    if len(set(names)) != len(names):
        raise ValueError("Duplicate group names in GROUP_DEFS")
    if len(set(numbers)) != len(numbers):
        raise ValueError("Duplicate group numbers in GROUP_DEFS")
    for gdef in defs:
        if gdef.catch_all != ("except as above" in gdef.name):
            raise ValueError(f"catch_all flag mismatch: No.{gdef.no}")
        query = Chem.MolFromSmarts(gdef.smarts)
        if query is None:
            raise ValueError(f"Invalid SMARTS: No.{gdef.no} {gdef.smarts}")
        owned = tuple(
            a.GetIdx() for a in query.GetAtoms() if a.GetAtomMapNum() > 0
        )
        if not owned:
            raise ValueError(f"No mapped atoms: No.{gdef.no} {gdef.smarts}")
        compiled.append(_Compiled(gdef, query, owned))
    return compiled


_COMPILED = _compile(GROUP_DEFS)


# ---------------------------------------------------------------------------
# 解析
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class _Candidate:
    priority: int
    weight: float
    atoms: Tuple[int, ...]
    group: GroupDef


@dataclass
class GroupAnalysis:
    """解析結果。

    Attributes:
        groups: グループ名 -> 出現回数。
        assignments: (グループ番号, グループ名, 所有原子インデックス) の一覧。
        unassigned: 未割当の重原子 (インデックス, 元素記号)。空なら完全分割。
        ties: 同優先度・同重量で原子を奪い合った別グループの組
            (採用グループ名, 不採用グループ名)。空でなければ結果が
            タイブレーク規則 (D7) に依存している。
        warnings: 前処理 (芳香族性の再定義) で生じた警告。空でなければ
            結果を信頼できない。
    """

    groups: Dict[str, int]
    assignments: List[Tuple[int, str, Tuple[int, ...]]]
    unassigned: List[Tuple[int, str]]
    ties: List[Tuple[str, str]] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)

    @property
    def is_complete(self) -> bool:
        """未割当原子も警告もなければ True (完全かつ互いに素な分割)。"""
        return not self.unassigned and not self.warnings


def _group_weight(mol: Chem.Mol, atoms: Sequence[int]) -> float:
    """グループ重量 = 構成原子の原子量 + 結合水素の質量。"""
    total = 0.0
    for idx in atoms:
        atom = mol.GetAtomWithIdx(idx)
        total += atom.GetMass() + atom.GetTotalNumHs() * H_MASS
    return total


def _find_candidates(mol: Chem.Mol) -> List[_Candidate]:
    """全パターンの全マッチを候補として列挙する。"""
    candidates: List[_Candidate] = []
    for comp in _COMPILED:
        seen = set()
        matches = mol.GetSubstructMatches(
            comp.query, uniquify=False, maxMatches=100000
        )
        for match in matches:
            owned = tuple(sorted(match[i] for i in comp.owned_query_idx))
            if owned in seen:
                continue
            seen.add(owned)
            if comp.group.predicate and not comp.group.predicate(mol, owned):
                continue
            candidates.append(
                _Candidate(
                    priority=comp.group.priority,
                    weight=_group_weight(mol, owned),
                    atoms=owned,
                    group=comp.group,
                )
            )
    return candidates


def _is_paper_aromatic_ring(mol: Chem.Mol, ring: Sequence[int]) -> bool:
    """論文の運用定義による芳香環か (D4)。

    条件: 6員環、全員がRDKit芳香族、環員がCまたはNのみ、環外にヘテロ原子への
    二重結合がない (2-ピリドン型を除外)。
    """
    if len(ring) != 6:
        return False
    members = set(ring)
    for idx in ring:
        atom = mol.GetAtomWithIdx(idx)
        if not atom.GetIsAromatic() or atom.GetSymbol() not in ("C", "N"):
            return False
        for bond in atom.GetBonds():
            other = bond.GetOtherAtom(atom)
            if (
                other.GetIdx() not in members
                and bond.GetBondType() == Chem.BondType.DOUBLE
                and other.GetSymbol() != "C"
            ):
                return False
    return True


def _apply_paper_aromaticity(
    mol: Chem.Mol,
) -> Tuple[Chem.Mol, Optional[str]]:
    """RDKitの芳香環のうち、論文流の芳香環以外をKekule構造へ戻す (P4)。

    チオフェン、フラン、ピロール、イミダゾール、ピラノン、ピリドン、
    アズレン等は論文では環状グループで記述されるため、非芳香族にする。
    芳香環と共有する原子・結合は芳香族のまま残す
    (-> aC fused with non-aromatic subring)。

    Returns:
        (変換後の分子, 警告)。警告はKekule化や価数検査に失敗した場合のみ。
    """
    ring_info = mol.GetRingInfo()
    arom_rings = [
        r for r in ring_info.AtomRings()
        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in r)
    ]
    keep_rings = [
        set(r) for r in arom_rings if _is_paper_aromatic_ring(mol, r)
    ]
    if len(keep_rings) == len(arom_rings):
        return mol, None

    keep_atoms = set().union(*keep_rings) if keep_rings else set()
    keep_bonds = {
        b.GetIdx()
        for b in mol.GetBonds()
        if any(
            b.GetBeginAtomIdx() in r and b.GetEndAtomIdx() in r
            for r in keep_rings
        )
    }
    was_aromatic = {b.GetIdx() for b in mol.GetBonds() if b.GetIsAromatic()}
    try:
        rw = Chem.RWMol(mol)
        Chem.Kekulize(rw, clearAromaticFlags=False)
    except Exception as exc:  # noqa: BLE001 (RDKitの例外型が多岐にわたる)
        return mol, f"Kekulization failed: {exc}"
    for bond in rw.GetBonds():
        if bond.GetIdx() in keep_bonds:
            bond.SetBondType(Chem.BondType.AROMATIC)
        elif bond.GetIdx() in was_aromatic:
            bond.SetIsAromatic(False)
    for atom in rw.GetAtoms():
        if atom.GetIsAromatic() and atom.GetIdx() not in keep_atoms:
            atom.SetIsAromatic(False)
    new_mol = rw.GetMol()
    warning: Optional[str] = None
    try:
        new_mol.UpdatePropertyCache(strict=True)
    except Exception as exc:  # noqa: BLE001
        warning = f"Valence check failed after aromaticity change: {exc}"
        new_mol.UpdatePropertyCache(strict=False)
    Chem.GetSymmSSSR(new_mol)
    return new_mol, warning


# 同優先度・同重量・同種別 (catch-all か否か) のタイブレーク判定用
def _tie_key(cand: _Candidate) -> Tuple[int, float, bool]:
    return (cand.priority, round(cand.weight, 3), cand.group.catch_all)


def _sort_key(cand: _Candidate) -> Tuple:
    """割り当て順: 優先度 > 重量 > 具体的グループ > グループ番号 > 原子。

    リスト上の定義順には依存しない (D7)。
    """
    return (
        -cand.priority,
        -round(cand.weight, 3),
        cand.group.catch_all,
        cand.group.no,
        cand.atoms,
    )


def analyze_smiles(
    smiles: str, paper_aromaticity: bool = True
) -> GroupAnalysis:
    """SMILESを第1次グループへ分割する。

    候補を (優先度, 重量, 具体性, グループ番号) の順に並べ、既使用原子と
    重ならないものを貪欲に採用する (完全かつ互いに素な分割 + ゴールデン
    ルール)。割り当てアルゴリズムの限界は判断記録 D1, D2 を参照。

    Args:
        smiles: 有機化合物のSMILES。
        paper_aromaticity: Trueなら論文流の芳香環定義 (P4) を適用する。
            Falseの場合はRDKitの芳香族判定をそのまま使う。

    Raises:
        ValueError: SMILESが不正な場合。
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles!r}")
    warnings: List[str] = []
    if paper_aromaticity:
        mol, warning = _apply_paper_aromaticity(mol)
        if warning:
            warnings.append(warning)

    candidates = _find_candidates(mol)
    candidates.sort(key=_sort_key)

    owner: Dict[int, _Candidate] = {}
    counter: Counter = Counter()
    assignments: List[Tuple[int, str, Tuple[int, ...]]] = []
    ties = set()
    for cand in candidates:
        winners = [owner[a] for a in cand.atoms if a in owner]
        if winners:
            for win in winners:
                if (
                    win.group.no != cand.group.no
                    and _tie_key(win) == _tie_key(cand)
                ):
                    ties.add((win.group.name, cand.group.name))
            continue
        for atom_idx in cand.atoms:
            owner[atom_idx] = cand
        counter[cand.group.name] += 1
        assignments.append((cand.group.no, cand.group.name, cand.atoms))

    unassigned = [
        (a.GetIdx(), a.GetSymbol())
        for a in mol.GetAtoms()
        if a.GetIdx() not in owner
    ]
    return GroupAnalysis(
        dict(counter), assignments, unassigned, sorted(ties), warnings
    )


def count_first_order_groups(
    smiles: str, strict: bool = False, paper_aromaticity: bool = True
) -> Dict[str, int]:
    """SMILESから第1次グループ名と出現回数の辞書を返す。

    Args:
        smiles: 有機化合物のSMILES。
        strict: Trueの場合、完全な分割でないとValueErrorを送出する。
        paper_aromaticity: 論文流の芳香環定義を適用する (P4)。

    Returns:
        グループ名 -> 出現回数。
    """
    result = analyze_smiles(smiles, paper_aromaticity)
    if strict and not result.is_complete:
        raise ValueError(
            f"Incomplete split of {smiles!r}: "
            f"unassigned={result.unassigned}, warnings={result.warnings}"
        )
    return result.groups


# ---------------------------------------------------------------------------
# 第2次グループ (Table 7)
# ---------------------------------------------------------------------------
# 第1次との違い (論文 Section 2.2):
#   * 分子全体を覆う必要はなく、第1次で記述しきれない断片だけを表す。
#   * 原子の重複を許す。ただし別のグループの原子集合に完全に含まれるグループは
#     採用しない (完全重複の禁止。例: l-アラニンでは CHm(NHn)-COOH を採用し、
#     CHm-COOH は採用しない)。
#
# SMARTSの規約 (第1次と異なる点に注意):
#   * 1-99 のアトムマップ番号: 出現を数える際の識別子 (identity) になる原子。
#   * 100以上のアトムマップ番号: 原子集合 (extent, 完全重複の判定用) にだけ
#     入る原子。メチル基の選び方などによる二重計上を避けるために使う (D11)。
#   * マップ番号のない原子: 文脈 (マッチ条件) としてのみ使う。
_EXTENT_ONLY = 100

# 二重結合・三重結合・アミン等の共通部品
_AMINE2 = "NX3;+0;!$(N[CX3]=[O,S,N])"      # アミン窒素 (アミド・ニトロを除く)
_KETONE_C = "CX3;!R;$(C(=O)([#6])[#6])"    # ケトンのカルボニル炭素
# 環内 C=C の炭素 2つ (マップ 1, 2) と、環外へ出る単結合 (D12, D18)
_CYC_CC = "[CX3;R;!a:1]=;@[CX3;R;!a:2]-;!@"

Finder = Callable[[Chem.Mol], List[Tuple[frozenset, frozenset]]]


@dataclass(frozen=True)
class SecondOrderDef:
    """第2次グループ定義。

    Attributes:
        no: 論文 Table 7 のグループ番号。
        name: グループ名 (ASCII表記)。
        smarts: アトムマップ付きSMARTS。finder と排他。
        finder: SMARTSで表せないグループ (環の置換位置など) の検出関数。
            分子から (identity, extent) の組のリストを返す。smarts と排他。
        predicate: SMARTSで表せない追加条件。引数は (分子, identity 原子の
            インデックス。マップ番号の昇順)。
    """

    no: int
    name: str
    smarts: str = ""
    finder: Optional[Finder] = None
    predicate: Optional[Predicate] = None


def _cyc_ring_nonaromatic(mol: Chem.Mol, owned: Tuple[int, ...]) -> bool:
    """環内 C=C (owned[0], owned[1]) を含む環に芳香族原子がないこと (D18)。"""
    first, second = owned[0], owned[1]
    for ring in mol.GetRingInfo().AtomRings():
        if first in ring and second in ring:
            if not any(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
                return True
    return False


# --- 環の置換位置による分類 (AROMRINGs, PYRIDINEs) -------------------------
def _canon_benzene(positions: Sequence[int]) -> Tuple[int, ...]:
    """ベンゼン環の置換位置 (0-5) を回転・鏡映で正規化する。"""
    forms = []
    for shift in range(6):
        for sign in (1, -1):
            forms.append(tuple(sorted((sign * p + shift) % 6 for p in positions)))
    return min(forms)


def _canon_pyridine(offsets: Sequence[int]) -> Tuple[int, ...]:
    """ピリジン環の置換位置 (Nからの距離 1-5) を鏡映のみで正規化する。"""
    mirrored = tuple(sorted((6 - p) % 6 for p in offsets))
    return min(tuple(sorted(offsets)), mirrored)


def _make_ring_finder(kind: str, target: Tuple[int, ...]) -> Finder:
    """縮合していない6員芳香環のうち、置換位置が target と一致する環を返す。

    Args:
        kind: "benzene" (環員がすべてC) または "pyridine" (Nを1つだけ含む)。
        target: 正規化済みの置換位置。
    """

    def finder(mol: Chem.Mol) -> List[Tuple[frozenset, frozenset]]:
        info = mol.GetRingInfo()
        found: List[Tuple[frozenset, frozenset]] = []
        for ring in info.AtomRings():
            if len(ring) != 6 or any(info.NumAtomRings(i) > 1 for i in ring):
                continue  # 縮合環は第3次グループの領分 (D14)
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            if not all(a.GetIsAromatic() for a in atoms):
                continue
            if any(a.GetSymbol() not in ("C", "N") for a in atoms):
                continue
            members = set(ring)
            n_pos = [k for k, a in enumerate(atoms) if a.GetSymbol() == "N"]
            subs = [
                k for k, a in enumerate(atoms)
                if any(nb.GetIdx() not in members for nb in a.GetNeighbors())
            ]
            if kind == "benzene":
                if n_pos:
                    continue
                key = _canon_benzene(subs)
            else:
                if len(n_pos) != 1 or n_pos[0] in subs:
                    continue
                key = _canon_pyridine([(k - n_pos[0]) % 6 for k in subs])
            if key == target:
                found.append((frozenset(ring), frozenset(ring)))
        return found

    return finder


# (グループ番号, 1始まりの置換位置)
_BENZENE_RING_GROUPS: List[Tuple[int, Tuple[int, ...]]] = [
    (104, (1, 2)), (105, (1, 3)), (106, (1, 4)),
    (107, (1, 2, 3)), (108, (1, 2, 4)), (109, (1, 3, 5)),
    (110, (1, 2, 3, 4)), (111, (1, 2, 3, 5)), (112, (1, 2, 4, 5)),
]
# (グループ番号, Nを1とした置換位置)
_PYRIDINE_RING_GROUPS: List[Tuple[int, Tuple[int, ...]]] = [
    (113, (2,)), (114, (3,)), (115, (4,)),
    (116, (2, 3)), (117, (2, 4)), (118, (2, 5)), (119, (2, 6)),
    (120, (3, 4)), (121, (3, 5)), (122, (2, 3, 6)),
]


def _ring_defs() -> List[SecondOrderDef]:
    """AROMRINGs (No.104-112) と PYRIDINEs (No.113-122) の定義を作る。"""
    defs: List[SecondOrderDef] = []
    for no, pos in _BENZENE_RING_GROUPS:
        target = _canon_benzene([p - 1 for p in pos])
        name = "AROMRINGs" + "s".join(str(p) for p in pos)
        defs.append(SecondOrderDef(no, name, finder=_make_ring_finder("benzene", target)))
    for no, pos in _PYRIDINE_RING_GROUPS:
        target = _canon_pyridine([p - 1 for p in pos])
        name = "PYRIDINEs" + "s".join(str(p) for p in pos)
        defs.append(SecondOrderDef(no, name, finder=_make_ring_finder("pyridine", target)))
    return defs


# (番号, 名前, SMARTS)。炭素は、名前に cyc / aC が付かないものは非環状 (D12)。
_SECOND_SMARTS: List[Tuple[int, str, str]] = [
    # --- 分岐アルカン (中心炭素が識別子、メチルは原子集合のみ: D11) ---
    (1, "(CH3)2CH", "[CX4H1;!R:1]([CH3:101])[CH3:102]"),
    (2, "(CH3)3C", "[CX4H0;!R:1]([CH3:101])([CH3:102])[CH3:103]"),
    (3, "CH(CH3)CH(CH3)", "[CX4H1;!R:1]([CH3:101])-[CX4H1;!R:2][CH3:102]"),
    (4, "CH(CH3)C(CH3)2",
     "[CX4H1;!R:1]([CH3:101])-[CX4H0;!R:2]([CH3:102])[CH3:103]"),
    (5, "C(CH3)2C(CH3)2",
     "[CX4H0;!R:1]([CH3:101])([CH3:102])-[CX4H0;!R:2]([CH3:103])[CH3:104]"),
    # --- アルケン ---
    (6, "CHn=CHm-CHp=CHk", "[CX3;!R:1]=[CX3;!R:2]-[CX3;!R:3]=[CX3;!R:4]"),
    (7, "CH3-CHm=CHn", "[CH3:1]-[CX3;!R:2]=[CX3;!R:3]"),
    (8, "CH2-CHm=CHn", "[CH2;!R:1]-[CX3;!R:2]=[CX3;!R:3]"),
    (9, "CHp-CHm=CHn", "[CX4;H0,H1;!R:1]-[CX3;!R:2]=[CX3;!R:3]"),
    # --- カルボニル化合物・アルコール ---
    (10, "CHCHO or CCHO", "[CX4;H0,H1;!R:1]-[CX3H1:2]=[O:3]"),
    (11, "CH3COCH2", "[CH3:1][CX3;!R:2](=[O:3])[CH2;!R:4]"),
    (12, "CH3COCH or CH3COC",
     "[CH3:1][CX3;!R:2](=[O:3])[CX4;H0,H1;!R:4]"),
    (13, "CHCOOH or CCOOH",
     "[CX4;H0,H1;!R:1]-[CX3:2](=[O:3])-[OX2H1:4]"),
    (14, "CH3COOCH or CH3COOC",
     "[CH3:1][CX3:2](=[O:3])[OX2:4][CX4;H0,H1;!R:5]"),
    (15, "CO-O-CO", "[CX3:1](=[O:2])-[OX2:3]-[CX3:4]=[O:5]"),
    (16, "CHOH", "[CX4H1;!R:1]-[OX2H1:2]"),
    (17, "COH", "[CX4H0;!R:1]-[OX2H1:2]"),
    (18, "CH3COCHn(OH)", "[CH3:1][CX3;!R:2](=[O:3])[CX4;!R:4]-[OX2H1:5]"),
    (19, "NCCHOH or NCCOH", "[NX1:1]#[CX2:2]-[CX4;H0,H1;!R:3]-[OX2H1:4]"),
    (20, "OH-CHn-COO",
     "[OX2H1:1]-[CX4;!R:2]-[CX3:3](=[O:4])-[OX2H0:5]"),
    (21, "CHm(OH)CHn(OH)",
     "[OX2H1:1]-[CX4;!R:2]-[CX4;!R:3]-[OX2H1:4]"),
    (22, "CHm(OH)CHn(NHp)",
     f"[OX2H1:1]-[CX4;!R:2]-[CX4;!R:3]-[{_AMINE2}:4]"),
    (23, "CHm(NH2)CHn(NH2)",
     "[NX3H2:1]-[CX4;!R:2]-[CX4;!R:3]-[NX3H2:4]"),
    (24, "CHm(NH)CHn(NH2)",
     f"[{_AMINE2};H1:1]-[CX4;H1,H2;!R:2]-[CX4;H1,H2;!R:3]-[NX3H2:4]"),
    (25, "H2NCOCHnCHmCONH2",
     "[NX3H2:1]-[CX3:2](=[O:3])-[CX4;H1,H2;!R:4]-[CX4;H1,H2;!R:5]"
     "-[CX3:6](=[O:7])-[NX3H2:8]"),
    (26, "CHm(NHn)-COOH",
     f"[CX4;!R:1](-[{_AMINE2}:2])-[CX3:3](=[O:4])-[OX2H1:5]"),
    (27, "HOOC-CHn-COOH",
     "[OX2H1:1]-[CX3:2](=[O:3])-[CX4;H1,H2;!R:4]-[CX3:5](=[O:6])-[OX2H1:7]"),
    (28, "HOOC-CHn-CHm-COOH",
     "[OX2H1:1]-[CX3:2](=[O:3])-[CX4;H1,H2;!R:4]-[CX4;H1,H2;!R:5]"
     "-[CX3:6](=[O:7])-[OX2H1:8]"),
    (29, "HO-CHn-COOH", "[OX2H1:1]-[CX4;!R:2]-[CX3:3](=[O:4])-[OX2H1:5]"),
    (30, "NH2-CHn-CHm-COOH",
     "[NX3H2:1]-[CX4;H1,H2;!R:2]-[CX4;H1,H2;!R:3]-[CX3:4](=[O:5])-[OX2H1:6]"),
    (31, "CH3-O-CHn-COOH",
     "[CH3:1]-[OX2:2]-[CX4;H1,H2;!R:3]-[CX3:4](=[O:5])-[OX2H1:6]"),
    (32, "HS-CH-COOH", "[SX2H1:1]-[CX4H1;!R:2]-[CX3:3](=[O:4])-[OX2H1:5]"),
    (33, "HS-CHn-CHm-COOH",
     "[SX2H1:1]-[CX4;H1,H2;!R:2]-[CX4;H1,H2;!R:3]-[CX3:4](=[O:5])-[OX2H1:6]"),
    (34, "NC-CHn-CHm-CN",
     "[NX1:1]#[CX2:2]-[CX4;H1,H2;!R:3]-[CX4;H1,H2;!R:4]-[CX2:5]#[NX1:6]"),
    (35, "OH-CHn-CHm-CN",
     "[OX2H1:1]-[CX4;H1,H2;!R:2]-[CX4;H1,H2;!R:3]-[CX2:4]#[NX1:5]"),
    (36, "HS-CHn-CHm-SH",
     "[SX2H1:1]-[CX4;H1,H2;!R:2]-[CX4;H1,H2;!R:3]-[SX2H1:4]"),
    (37, "COO-CHn-CHm-OOC",
     "[CX3:1](=[O:2])-[OX2:3]-[CX4;H1,H2;!R:4]-[CX4;H1,H2;!R:5]"
     "-[OX2:6]-[CX3:7]=[O:8]"),
    (38, "OOC-CHn-CHm-COO",
     "[OX2H0:1]-[CX3:2](=[O:3])-[CX4;H1,H2;!R:4]-[CX4;H1,H2;!R:5]"
     "-[CX3:6](=[O:7])-[OX2H0:8]"),
    (39, "NC-CHn-COO",
     "[NX1:1]#[CX2:2]-[CX4;H1,H2;!R:3]-[CX3:4](=[O:5])-[OX2H0:6]"),
    (40, "COCHnCOO",
     f"[{_KETONE_C}:1](=[O:2])-[CX4;H1,H2;!R:3]-[CX3:4](=[O:5])-[OX2H0:6]"),
    # --- ビニル化合物 ---
    (41, "CHm-O-CHn=CHp", "[CX4;!R:1]-[OX2:2]-[CX3;!R:3]=[CX3;!R:4]"),
    (42, "CHm=CHn-F", "[CX3;!R:1]=[CX3;!R:2]-[F:3]"),
    (43, "CHm=CHn-Br", "[CX3;!R:1]=[CX3;!R:2]-[Br:3]"),
    (44, "CHm=CHn-I", "[CX3;!R:1]=[CX3;!R:2]-[I:3]"),
    (45, "CHm=CHn-Cl", "[CX3;!R:1]=[CX3;!R:2]-[Cl:3]"),
    (46, "CHm=CHn-CN", "[CX3;!R:1]=[CX3;!R:2]-[CX2:3]#[NX1:4]"),
    (47, "CHn=CHm-COO-CHp",
     "[CX3;!R:1]=[CX3;!R:2]-[CX3:3](=[O:4])-[OX2H0:5]-[CX4;!R:6]"),
    (48, "CHm=CHn-CHO", "[CX3;!R:1]=[CX3;!R:2]-[CX3H1:3]=[O:4]"),
    (49, "CHm=CHn-COOH", "[CX3;!R:1]=[CX3;!R:2]-[CX3:3](=[O:4])-[OX2H1:5]"),
    # --- ベンジル型 (aC-CHn-R) ---
    (50, "aC-CHn-X", "[c:1]-[CX4;H1,H2;!R:2]-[F,Cl,Br,I:3]"),
    (51, "aC-CHn-NHm", f"[c:1]-[CX4;H1,H2;!R:2]-[{_AMINE2}:3]"),
    (52, "aC-CHn-O-", f"[c:1]-[CX4;H1,H2;!R:2]-[{_ETHER_O}:3]"),
    (53, "aC-CHn-OH", "[c:1]-[CX4;H1,H2;!R:2]-[OX2H1:3]"),
    (54, "aC-CHn-CN", "[c:1]-[CX4;H1,H2;!R:2]-[CX2:3]#[NX1:4]"),
    (55, "aC-CHn-CHO", "[c:1]-[CX4;H1,H2;!R:2]-[CX3H1:3]=[O:4]"),
    (56, "aC-CHn-SH", "[c:1]-[CX4;H1,H2;!R:2]-[SX2H1:3]"),
    (57, "aC-CHn-COOH",
     "[c:1]-[CX4;H1,H2;!R:2]-[CX3:3](=[O:4])-[OX2H1:5]"),
    (58, "aC-CHn-CO-", f"[c:1]-[CX4;H1,H2;!R:2]-[{_KETONE_C}:3]=[O:4]"),
    (59, "aC-CHn-S-", f"[c:1]-[CX4;H1,H2;!R:2]-[{_THIOETHER_S}:3]"),
    (60, "aC-CHn-OOC-H",
     "[c:1]-[CX4;H1,H2;!R:2]-[OX2H0:3]-[CX3H1:4]=[O:5]"),
    (61, "aC-CHm-NO2",
     f"[c:1]-[CX4;H1,H2;!R:2]-[{_NITRO_N}:3](=[O:4])[O-:5]"),
    (62, "aC-CHn-CONH2",
     "[c:1]-[CX4;H1,H2;!R:2]-[CX3:3](=[O:4])-[NX3H2:5]"),
    (63, "aC-CHn-OOC",
     "[c:1]-[CX4;H1,H2;!R:2]-[OX2H0:3]-[CX3H0:4](=[O:5])[#6]"),
    (64, "aC-CHn-COO",
     "[c:1]-[CX4;H1,H2;!R:2]-[CX3:3](=[O:4])-[OX2H0:5]"),
    (65, "aC-SO2-OH", "[c:1]-[SX4:2](=[O:3])(=[O:4])-[OX2H1:5]"),
    (66, "aC-CH(CH3)2", "[c:1]-[CX4H1;!R:2]([CH3:3])[CH3:4]"),
    (67, "aC-C(CH3)3", "[c:1]-[CX4H0;!R:2]([CH3:3])([CH3:4])[CH3:5]"),
    (68, "aC-CF3", "[c:1]-[CX4:2]([F:3])([F:4])[F:5]"),
    # --- 環内 C=C への置換基 ((CHn=C)cyc-R) ---
    (69, "(CHn=C)cyc-CHO", f"{_CYC_CC}[CX3H1;!R:3]=[O:4]"),
    (70, "(CHn=C)cyc-COO-CHm",
     f"{_CYC_CC}[CX3;!R:3](=[O:4])-[OX2H0:5]-[CX4;!R:6]"),
    (71, "(CHn=C)cyc-CO-", f"{_CYC_CC}[{_KETONE_C}:3]=[O:4]"),
    (72, "(CHn=C)cyc-CH3", f"{_CYC_CC}[CH3:3]"),
    (73, "(CHn=C)cyc-CH2", f"{_CYC_CC}[CX4H2;!R:3]"),
    (74, "(CHn=C)cyc-CN", f"{_CYC_CC}[CX2:3]#[NX1:4]"),
    (75, "(CHn=C)cyc-Cl", f"{_CYC_CC}[Cl:3]"),
    # --- 環状 CH, C への置換基 ---
    (76, "CHcyc-CH3", "[CX4H1;R:1]-;!@[CH3:2]"),
    (77, "CHcyc-CH2", "[CX4H1;R:1]-;!@[CX4H2;!R:2]"),
    (78, "CHcyc-CH", "[CX4H1;R:1]-;!@[CX4H1;!R:2]"),
    (79, "CHcyc-C", "[CX4H1;R:1]-;!@[CX4H0;!R:2]"),
    (80, "CHcyc-CH=CHn", "[CX4H1;R:1]-;!@[CX3H1;!R:2]=[CX3;H1,H2;!R:3]"),
    (81, "CHcyc-C=CHn", "[CX4H1;R:1]-;!@[CX3H0;!R:2]=[CX3;H1,H2;!R:3]"),
    (82, "CHcyc-Cl", "[CX4H1;R:1]-[Cl:2]"),
    (83, "CHcyc-F", "[CX4H1;R:1]-[F:2]"),
    (84, "CHcyc-OH", "[CX4H1;R:1]-[OX2H1:2]"),
    (85, "CHcyc-NH2", "[CX4H1;R:1]-[NX3H2:2]"),
    (86, "CHcyc-NH-CHn",
     f"[CX4H1;R:1]-[{_AMINE2};H1;!R:2]-[CX4;!R:102]"),
    (87, "CHcyc-N-CHn",
     f"[CX4H1;R:1]-[{_AMINE2};H0;!R:2]-[CX4;!R:102]"),
    (88, "CHcyc-SH", "[CX4H1;R:1]-[SX2H1:2]"),
    (89, "CHcyc-CN", "[CX4H1;R:1]-[CX2:2]#[NX1:3]"),
    (90, "CHcyc-COOH", "[CX4H1;R:1]-[CX3:2](=[O:3])-[OX2H1:4]"),
    (91, "CHcyc-CO", f"[CX4H1;R:1]-[{_KETONE_C}:2]=[O:3]"),
    (92, "CHcyc-NO2", f"[CX4H1;R:1]-[{_NITRO_N}:2](=[O:3])[O-:4]"),
    (93, "CHcyc-S-", f"[CX4H1;R:1]-[{_THIOETHER_S}:2]"),
    (94, "CHcyc-CHO", "[CX4H1;R:1]-[CX3H1;!R:2]=[O:3]"),
    (95, "CHcyc-O-", f"[CX4H1;R:1]-[{_ETHER_O}:2]"),
    (96, "CHcyc-OOCH", "[CX4H1;R:1]-[OX2H0:2]-[CX3H1:3]=[O:4]"),
    (97, "CHcyc-COO", "[CX4H1;R:1]-[CX3;!R:2](=[O:3])-[OX2H0:4]"),
    (98, "CHcyc-OOC", "[CX4H1;R:1]-[OX2H0:2]-[CX3H0:3](=[O:4])[#6]"),
    (99, "Ccyc-CH3", "[CX4H0;R:1]-;!@[CH3:2]"),
    (100, "Ccyc-CH2", "[CX4H0;R:1]-;!@[CX4H2;!R:2]"),
    (101, "Ccyc-OH", "[CX4H0;R:1]-[OX2H1:2]"),
    (102, ">Ncyc-CH3", "[NX3;R:1]-;!@[CH3:2]"),
    (103, ">Ncyc-CH2", "[NX3;R:1]-;!@[CX4H2;!R:2]"),
]

# (CHn=C)cyc 系 (No.69-75) は、芳香環と縮合した環には適用しない (D18)
_CYC_CC_NOS = range(69, 76)

SECOND_ORDER_DEFS: List[SecondOrderDef] = [
    SecondOrderDef(
        no, name, smarts,
        predicate=_cyc_ring_nonaromatic if no in _CYC_CC_NOS else None,
    )
    for no, name, smarts in _SECOND_SMARTS
] + _ring_defs()


@dataclass(frozen=True)
class _SecondCompiled:
    group: SecondOrderDef
    query: Optional[Chem.Mol]
    identity_idx: Tuple[int, ...]  # マップ番号の昇順
    extent_idx: Tuple[int, ...]


def _compile_second(defs: Sequence[SecondOrderDef]) -> List[_SecondCompiled]:
    """第2次グループのSMARTSをコンパイルする。"""
    if len({g.name for g in defs}) != len(defs):
        raise ValueError("Duplicate second-order group names")
    if len({g.no for g in defs}) != len(defs):
        raise ValueError("Duplicate second-order group numbers")
    compiled: List[_SecondCompiled] = []
    for gdef in defs:
        if bool(gdef.smarts) == bool(gdef.finder):
            raise ValueError(f"Specify exactly one of smarts/finder: No.{gdef.no}")
        if gdef.finder:
            compiled.append(_SecondCompiled(gdef, None, (), ()))
            continue
        query = Chem.MolFromSmarts(gdef.smarts)
        if query is None:
            raise ValueError(f"Invalid SMARTS: No.{gdef.no} {gdef.smarts}")
        mapped = sorted(
            (a.GetAtomMapNum(), a.GetIdx())
            for a in query.GetAtoms() if a.GetAtomMapNum() > 0
        )
        if len({m for m, _ in mapped}) != len(mapped):
            raise ValueError(f"Duplicate atom map: No.{gdef.no}")
        identity = tuple(i for m, i in mapped if m < _EXTENT_ONLY)
        if not identity:
            raise ValueError(f"No identity atoms: No.{gdef.no}")
        compiled.append(
            _SecondCompiled(gdef, query, identity, tuple(i for _, i in mapped))
        )
    return compiled


_SECOND_COMPILED = _compile_second(SECOND_ORDER_DEFS)


@dataclass(frozen=True)
class _Occurrence:
    no: int
    name: str
    identity: frozenset
    extent: frozenset


@dataclass
class SecondOrderAnalysis:
    """第2次グループの解析結果。

    Attributes:
        groups: グループ名 -> 出現回数 (完全重複を除いた最終結果)。
        occurrences: (グループ番号, グループ名, identity 原子) の一覧。
        eliminated: 完全重複のため除外された (除外グループ名, 包含グループ名)。
        ties: 原子集合が完全に一致した別グループの組。空でなければ、
            どちらを採るかの規則が論文にないため結果に注意が必要 (D10)。
        warnings: 前処理 (芳香族性の再定義) で生じた警告。
        notes: 表にない置換パターンへのフォールバックや、個別管理している
            例外の適用など、結果の解釈に注意が要る事項 (第3次で使用)。
    """

    groups: Dict[str, int]
    occurrences: List[Tuple[int, str, Tuple[int, ...]]]
    eliminated: List[Tuple[str, str]] = field(default_factory=list)
    ties: List[Tuple[str, str]] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)


def _find_occurrences(
    mol: Chem.Mol, compiled: Sequence[_SecondCompiled]
) -> List[_Occurrence]:
    """全グループの出現を列挙する。同じ identity の出現は原子集合を統合する。"""
    merged: Dict[Tuple[int, frozenset], set] = {}
    for comp in compiled:
        gdef = comp.group
        if comp.query is None:
            for ident, extent in gdef.finder(mol):  # type: ignore[misc]
                merged.setdefault((gdef.no, ident), set()).update(extent)
            continue
        matches = mol.GetSubstructMatches(
            comp.query, uniquify=False, maxMatches=100000
        )
        for match in matches:
            owned = tuple(match[i] for i in comp.identity_idx)
            if gdef.predicate and not gdef.predicate(mol, owned):
                continue
            ident = frozenset(owned)
            extent = {match[i] for i in comp.extent_idx}
            merged.setdefault((gdef.no, ident), set()).update(extent | ident)
    names = {c.group.no: c.group.name for c in compiled}
    return [
        _Occurrence(no, names[no], ident, frozenset(extent))
        for (no, ident), extent in merged.items()
    ]


NotesFn = Callable[[Chem.Mol, List[_Occurrence]], List[str]]


def _analyze_with_overlap(
    smiles: str,
    compiled: Sequence[_SecondCompiled],
    paper_aromaticity: bool = True,
    notes_fn: Optional[NotesFn] = None,
) -> SecondOrderAnalysis:
    """SMILESを第2次/第3次グループへ分解する (重なりを許し、完全重複を除外)。

    別グループの原子集合に真部分集合として含まれる出現は除外する
    (論文 Section 2.2 の完全重複の禁止, D10)。

    Args:
        smiles: 有機化合物のSMILES。
        compiled: コンパイル済みのグループ定義。
        paper_aromaticity: Trueなら論文流の芳香環定義 (P4) を適用する。
        notes_fn: (分子, 採用された出現) から注記を作る関数。

    Raises:
        ValueError: SMILESが不正な場合。
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles!r}")
    warnings: List[str] = []
    if paper_aromaticity:
        mol, warning = _apply_paper_aromaticity(mol)
        if warning:
            warnings.append(warning)

    occs = _find_occurrences(mol, compiled)
    kept: List[_Occurrence] = []
    eliminated = set()
    ties = set()
    for occ in occs:
        bigger = [
            p for p in occs if p.no != occ.no and occ.extent < p.extent
        ]
        if bigger:
            eliminated.add((occ.name, min(bigger, key=lambda p: p.no).name))
            continue
        kept.append(occ)
        for other in occs:
            if other.no > occ.no and other.extent == occ.extent:
                ties.add((occ.name, other.name))

    kept.sort(key=lambda o: (o.no, sorted(o.identity)))
    counter: Counter = Counter(o.name for o in kept)
    return SecondOrderAnalysis(
        groups=dict(counter),
        occurrences=[(o.no, o.name, tuple(sorted(o.identity))) for o in kept],
        eliminated=sorted(eliminated),
        ties=sorted(ties),
        warnings=warnings,
        notes=notes_fn(mol, kept) if notes_fn else [],
    )


def analyze_second_order(
    smiles: str, paper_aromaticity: bool = True
) -> SecondOrderAnalysis:
    """SMILESを第2次グループへ分解する (原子の重複を許す)。"""
    return _analyze_with_overlap(smiles, _SECOND_COMPILED, paper_aromaticity)


def count_second_order_groups(
    smiles: str, paper_aromaticity: bool = True
) -> Dict[str, int]:
    """SMILESから第2次グループ名と出現回数の辞書を返す。"""
    return analyze_second_order(smiles, paper_aromaticity).groups


# ---------------------------------------------------------------------------
# 第3次グループ (Table 8)
# ---------------------------------------------------------------------------
# 第2次と同じ枠組み (SecondOrderDef, identity/extent, 完全重複の除外) を使う。
# 完全重複の除外は第3次グループの中だけで行う (D19)。

# 原典内の不整合などを個別に管理する例外 (D23)。Falseにすると表の範囲どおり。
THIRD_ORDER_EXCEPTIONS: Dict[str, bool] = {
    # No.17 は表が (n in 0..1) だが、例示のテトラリン・インダンは CH2 を持つ。
    "No17_CH2_fused_carbocycle": True,
}

_RINGC = "CX3;R;!a;$(C=;@C)"  # 環内 C=C の炭素 (非芳香族)
_CHAIN_Q = Chem.MolFromSmarts("[CX4;!R]")  # 鎖を構成する炭素 (非環状 sp3)


def _in_fused_carbocycle(mol: Chem.Mol, aro: int, sp3: int) -> bool:
    """aro と sp3 が、芳香環と縮合した5-6員の炭素環に共に含まれるか。"""
    for ring in mol.GetRingInfo().AtomRings():
        if aro not in ring or sp3 not in ring or len(ring) not in (5, 6):
            continue
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if all(a.GetSymbol() == "C" for a in atoms) and not all(
            a.GetIsAromatic() for a in atoms
        ):
            return True
    return False


def _no17_range(mol: Chem.Mol, owned: Tuple[int, ...]) -> bool:
    """No.17 の水素数条件。表は (n in 0..1) で、CH2 は例外として個別管理 (D23)。"""
    n_h = mol.GetAtomWithIdx(owned[1]).GetTotalNumHs()
    if n_h <= 1:
        return True
    return (
        n_h == 2
        and THIRD_ORDER_EXCEPTIONS["No17_CH2_fused_carbocycle"]
        and _in_fused_carbocycle(mol, owned[0], owned[1])
    )


# --- 鎖状グループ (HOOC-(CHn)m-COOH など) ----------------------------------
EndCompiled = Tuple[Chem.Mol, int, Tuple[int, ...]]


def _compile_end(smarts: str) -> EndCompiled:
    """末端官能基のSMARTSをコンパイルする。マップ番号1が鎖に結合する原子。"""
    query = Chem.MolFromSmarts(smarts)
    if query is None:
        raise ValueError(f"Invalid SMARTS: {smarts}")
    mapped = {a.GetAtomMapNum(): a.GetIdx() for a in query.GetAtoms()
              if a.GetAtomMapNum() > 0}
    if 1 not in mapped:
        raise ValueError(f"No anchor (map 1): {smarts}")
    return query, mapped[1], tuple(mapped.values())


def _end_occurrences(
    mol: Chem.Mol, end: EndCompiled
) -> List[Tuple[int, frozenset]]:
    """末端官能基の出現を (結合原子, 原子集合) で返す。"""
    query, anchor_idx, ext_idx = end
    found: Dict[int, set] = {}
    for match in mol.GetSubstructMatches(query, uniquify=False, maxMatches=100000):
        found.setdefault(match[anchor_idx], set()).update(match[i] for i in ext_idx)
    return [(a, frozenset(e)) for a, e in found.items()]


def _chain_atoms(mol: Chem.Mol) -> set:
    return {m[0] for m in mol.GetSubstructMatches(_CHAIN_Q)}


def _chain_path(
    mol: Chem.Mol, chain: set, start: int, goal: int
) -> Optional[List[int]]:
    """鎖の炭素だけを通る start から goal への経路 (非環状なので一意)。"""
    prev: Dict[int, Optional[int]] = {start: None}
    queue = [start]
    while queue:
        cur = queue.pop(0)
        if cur == goal:
            path = []
            node: Optional[int] = cur
            while node is not None:
                path.append(node)
                node = prev[node]
            return path[::-1]
        for nb in mol.GetAtomWithIdx(cur).GetNeighbors():
            j = nb.GetIdx()
            if j in chain and j not in prev:
                prev[j] = cur
                queue.append(j)
    return None


def _chain_starts(mol: Chem.Mol, anchor: int, chain: set) -> List[int]:
    return [n.GetIdx() for n in mol.GetAtomWithIdx(anchor).GetNeighbors()
            if n.GetIdx() in chain]


def _make_chain_finder(end_a: str, end_b: str, min_len: int) -> Finder:
    """末端 A と B が、炭素数 min_len 以上の鎖でつながる組を検出する。

    出現は末端の組ごとに1回数える (D19)。
    """
    comp_a, comp_b = _compile_end(end_a), _compile_end(end_b)

    def finder(mol: Chem.Mol) -> List[Tuple[frozenset, frozenset]]:
        chain = _chain_atoms(mol)
        ends_a = _end_occurrences(mol, comp_a)
        ends_b = ends_a if end_a == end_b else _end_occurrences(mol, comp_b)
        merged: Dict[frozenset, set] = {}
        for anchor_a, ext_a in ends_a:
            starts_a = _chain_starts(mol, anchor_a, chain)
            for anchor_b, ext_b in ends_b:
                if anchor_a == anchor_b:
                    continue
                for sa in starts_a:
                    for sb in _chain_starts(mol, anchor_b, chain):
                        path = _chain_path(mol, chain, sa, sb)
                        if path and len(path) >= min_len:
                            merged.setdefault(
                                frozenset({anchor_a, anchor_b}), set()
                            ).update(ext_a | ext_b | set(path))
        return [(i, frozenset(e)) for i, e in merged.items()]

    return finder


def _make_bridge_finder(end_smarts: str, mid_smarts: str) -> Finder:
    """OH-(CHp)k-X-(CHn)m-OH 型 (X は O, S, N)。k, m は1以上 (No.5-7)。"""
    comp_end, comp_mid = _compile_end(end_smarts), _compile_end(mid_smarts)

    def finder(mol: Chem.Mol) -> List[Tuple[frozenset, frozenset]]:
        chain = _chain_atoms(mol)
        ends = _end_occurrences(mol, comp_end)
        merged: Dict[frozenset, set] = {}
        for mid, ext_mid in _end_occurrences(mol, comp_mid):
            arms = []
            for c in _chain_starts(mol, mid, chain):
                for anchor, ext in ends:
                    for s in _chain_starts(mol, anchor, chain):
                        path = _chain_path(mol, chain, c, s)
                        if path:
                            arms.append((c, anchor, ext | set(path)))
            for (c1, e1, x1), (c2, e2, x2) in combinations(arms, 2):
                if c1 != c2 and e1 != e2:
                    merged.setdefault(frozenset({e1, mid, e2}), set()).update(
                        x1 | x2 | ext_mid
                    )
        return [(i, frozenset(e)) for i, e in merged.items()]

    return finder


def _make_ether_alcohol_finder() -> Finder:
    """CHp-O-(CHn)m-OH (m>2, p in 0..2): エーテルの片側に CH3 以外の炭素 (No.8)。"""
    comp_oh = _compile_end("[OX2H1:1]")

    def finder(mol: Chem.Mol) -> List[Tuple[frozenset, frozenset]]:
        chain = _chain_atoms(mol)
        ends = _end_occurrences(mol, comp_oh)
        merged: Dict[frozenset, set] = {}
        for atom in mol.GetAtoms():
            if (atom.GetSymbol() != "O" or atom.IsInRing() or atom.GetDegree() != 2
                    or atom.GetTotalNumHs() != 0):
                continue
            nbrs = [n.GetIdx() for n in atom.GetNeighbors()]
            if not all(i in chain for i in nbrs):
                continue
            for p_side, start in ((nbrs[0], nbrs[1]), (nbrs[1], nbrs[0])):
                if mol.GetAtomWithIdx(p_side).GetTotalNumHs() > 2:
                    continue  # CHp は p in 0..2。メチルエーテルは対象外
                for anchor, ext in ends:
                    for s in _chain_starts(mol, anchor, chain):
                        path = _chain_path(mol, chain, start, s)
                        if path and len(path) >= 3:
                            merged.setdefault(
                                frozenset({atom.GetIdx(), anchor}), set()
                            ).update(ext | set(path) | {atom.GetIdx(), p_side})
        return [(i, frozenset(e)) for i, e in merged.items()]

    return finder


# --- 芳香族炭素2つをP, Sでつなぐグループ (No.40, 41) -------------------------
def _make_hetero_aryl_finder(symbol: str, need_oxygen: bool) -> Finder:
    """P または S に、芳香族炭素 (直接または O 経由) が2つ以上つくものを検出する。

    出現は芳香族炭素の組ごとに1回数える (D24。リン酸トリフェニル = 3)。
    """

    def finder(mol: Chem.Mol) -> List[Tuple[frozenset, frozenset]]:
        info = mol.GetRingInfo()
        found: List[Tuple[frozenset, frozenset]] = []
        for atom in mol.GetAtoms():
            if atom.GetSymbol() != symbol:
                continue
            oxygens = {n.GetIdx() for n in atom.GetNeighbors() if n.GetSymbol() == "O"}
            if need_oxygen and not oxygens:
                continue
            branches: List[Tuple[int, Tuple[int, ...]]] = []
            for nb in atom.GetNeighbors():
                if nb.GetSymbol() == "C" and nb.GetIsAromatic():
                    branches.append((nb.GetIdx(), ()))
                elif nb.GetSymbol() == "O" and nb.GetDegree() == 2:
                    for nb2 in nb.GetNeighbors():
                        if (nb2.GetIdx() != atom.GetIdx() and nb2.GetSymbol() == "C"
                                and nb2.GetIsAromatic()):
                            branches.append((nb2.GetIdx(), (nb.GetIdx(),)))
            for (c1, via1), (c2, via2) in combinations(branches, 2):
                if any(c1 in r and c2 in r for r in info.AtomRings()):
                    continue  # 同じ環に属する組は「different rings」ではない
                ident = frozenset({atom.GetIdx(), c1, c2})
                found.append((ident, ident | set(via1) | set(via2) | oxygens))
        return found

    return finder


# --- 縮合芳香環 (AROMFUSED, PYRIDINE.FUSED) ---------------------------------
# AROMFUSED[2] の置換パターン (非縮合原子を1-4と番号付け、鏡映 p -> 5-p で同一視)
_AROMFUSED2_NAMES: Dict[Tuple[int, ...], str] = {
    (): "AROMFUSED[2]", (1,): "AROMFUSED[2]s1", (2,): "AROMFUSED[2]s2",
    (2, 3): "AROMFUSED[2]s2s3", (1, 4): "AROMFUSED[2]s1s4",
    (1, 2): "AROMFUSED[2]s1s2", (1, 3): "AROMFUSED[2]s1s3",
}


def _aromatic_fused_rings(mol: Chem.Mol) -> List[Tuple[Tuple[int, ...], str]]:
    """縮合した6員芳香環を (環原子, 種別) で返す。種別は C (炭素のみ) か N (N1個)。"""
    info = mol.GetRingInfo()
    found = []
    for ring in info.AtomRings():
        if len(ring) != 6:
            continue
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        symbols = [a.GetSymbol() for a in atoms]
        if not all(a.GetIsAromatic() for a in atoms):
            continue
        if any(sym not in ("C", "N") for sym in symbols) or symbols.count("N") > 1:
            continue
        if sum(1 for i in ring if info.NumAtomRings(i) > 1) < 2:
            continue
        found.append((tuple(ring), "N" if "N" in symbols else "C"))
    return found


def _classify_fused_ring(
    mol: Chem.Mol, ring: Tuple[int, ...], kind: str
) -> Tuple[Optional[str], Optional[str]]:
    """縮合芳香環をグループ名に分類する。(グループ名 or None, 注記 or None)。

    AROMFUSED[k] の k は、他の環と共有する環原子の数 (D21)。
    """
    info = mol.GetRingInfo()
    members = set(ring)
    fus = [info.NumAtomRings(i) > 1 for i in ring]
    subs = [
        (not fus[t]) and any(
            nb.GetIdx() not in members
            for nb in mol.GetAtomWithIdx(ring[t]).GetNeighbors()
        )
        for t in range(6)
    ]
    k = sum(fus)
    fusion = [t for t in range(6) if fus[t]]
    nonf = [t for t in range(6) if not fus[t]]
    n_sub = sum(subs)

    if kind == "C":
        if k == 2:
            a, b = fusion
            if (b - a) % 6 == 1:
                second = b
            elif (a - b) % 6 == 1:
                second = a
            else:
                return None, None  # 縮合原子が隣接しない
            order = [(second + 1 + j) % 6 for j in range(4)]
            pos = tuple(j + 1 for j, t in enumerate(order) if subs[t])
            canon = min(pos, tuple(sorted(5 - p for p in pos)))
            if canon in _AROMFUSED2_NAMES:
                return _AROMFUSED2_NAMES[canon], None
            return _AROMFUSED2_NAMES[()], (
                f"AROMFUSED[2] の置換位置 {pos} は表になく、基本グループを使用"
            )
        if k == 3:
            note = ("AROMFUSED[3] の置換基は区別しない (表に置換体なし)"
                    if n_sub else None)
            return "AROMFUSED[3]", note
        if k == 4:
            d = (nonf[1] - nonf[0]) % 6
            if d == 3:  # アントラセン型
                names = {0: "AROMFUSED[4a]", 1: "AROMFUSED[4a]s1",
                         2: "AROMFUSED[4a]s1s4"}
                return names[n_sub], None
            if d in (1, 5):  # フェナントレン型
                if n_sub == 2:
                    return "AROMFUSED[4p]s3s4", None
                note = ("AROMFUSED[4p] の置換位置は表になく、基本グループを使用"
                        if n_sub == 1 else None)
                return "AROMFUSED[4p]", note
        return None, None

    # kind == "N" (ピリジン環の縮合)
    n_pos = next(t for t in range(6) if mol.GetAtomWithIdx(ring[t]).GetSymbol() == "N")
    if fus[n_pos]:
        return None, None
    if k == 2:
        a, b = fusion
        if (b - a) % 6 not in (1, 5):
            return None, None
        near = fus[(n_pos - 1) % 6] or fus[(n_pos + 1) % 6]
        return ("PYRIDINE.FUSED[2]" if near else "PYRIDINE.FUSED[2-iso]"), None
    if k == 4 and (nonf[1] - nonf[0]) % 6 == 3:
        return "PYRIDINE.FUSED[4]", None
    return None, None


def _make_fused_finder(target: str) -> Finder:
    """縮合芳香環のうち、分類結果が target のものを検出する。"""

    def finder(mol: Chem.Mol) -> List[Tuple[frozenset, frozenset]]:
        found = []
        for ring, kind in _aromatic_fused_rings(mol):
            if _classify_fused_ring(mol, ring, kind)[0] == target:
                found.append((frozenset(ring), frozenset(ring)))
        return found

    return finder


# --- 定義 ---------------------------------------------------------------------
_COOH_END = "[CX3:1](=[O:2])[OX2H1:3]"
_OH_END = "[OX2H1:1]"
_NH2_END = "[NX3H2:1]"

# (番号, 名前, SMARTS)
_THIRD_SMARTS: List[Tuple[int, str, str]] = [
    (14, "aC-(CHn=CHm)cyc (fused rings)",
     "[c:1]-@[CX3;R;!a:2]=;@[CX3;R;!a:3]"),
    (15, "aC-aC (different rings)", "[c:1]!:[c:2]"),
    (16, "aC-CHn,cyc (different rings)", "[c:1]-!@[CX4;R;H0,H1:2]"),
    (17, "aC-CHn,cyc (fused rings)", "[c:1]-@[CX4;R:2]"),
    (20, "CHcyc-CHcyc (different rings)", "[CX4H1;R:1]-!@[CX4H1;R:2]"),
    (22, "CH multiring", "[CX4H1;!R0;!R1:1]"),
    (23, "C multiring", "[CX4H0;!R0;!R1:1]"),
    (24, "aC-CHm-aC (different rings)", "[c:1]-[CX4;!R:2]-[c:3]"),
    (25, "aC-(CHm=CHn)-aC (different rings)",
     "[c:1]-[CX3;!R:2]=[CX3;!R:3]-[c:4]"),
    (26, "(CHm=C)cyc-CH=CH-(C=CHn)cyc (different rings)",
     f"[{_RINGC}:1]-[CX3H1;!R:2]=[CX3H1;!R:3]-[{_RINGC}:4]"),
    (27, "(CHm=C)cyc-CHp-(C=CHn)cyc (different rings)",
     f"[{_RINGC}:1]-[CX4;!R:2]-[{_RINGC}:3]"),
    (28, "aC-CO-aC (different rings)", "[c:1]-[CX3;!R:2](=[O:3])-[c:4]"),
    (29, "aC-CHm-CO-aC (different rings)",
     "[c:1]-[CX4;!R:2]-[CX3;!R:3](=[O:4])-[c:5]"),
    (30, "aC-CO-(C=CHn)cyc (different rings)",
     f"[c:1]-[CX3;!R:2](=[O:3])-[{_RINGC}:4]"),
    (31, "aC-CO-CO-aC (different rings)",
     "[c:1]-[CX3;!R:2](=[O:3])-[CX3;!R:4](=[O:5])-[c:6]"),
    (32, "aC-CO(cyc) (fused rings)", "[c:1]-@[CX3;R:2]=[O:3]"),
    (34, "aC-CO-CHn,cyc (different rings)",
     "[c:1]-[CX3;!R:2](=[O:3])-[CX4;R;H0,H1:4]"),
    (35, "aC-CO-NHn-aC (different rings)",
     "[c:1]-[CX3;!R:2](=[O:3])-[NX3;H0,H1;!R:4]-[c:5]"),
    (36, "aC-NHn-CONHm-aC (different rings)",
     "[c:1]-[NX3;H0,H1:2]-[CX3:3](=[O:4])-[NX3;H0,H1:5]-[c:6]"),
    (37, "aC-CO-Ncyc (different rings)",
     "[c:1]-[CX3;!R:2](=[O:3])-[NX3;R:4]"),
    (38, "aC-Scyc (fused rings)", "[c:1]-@[SX2;R:2]"),
    (39, "aC-S-aC (different rings)", "[c:1]-[SX2;!R:2]-[c:3]"),
    (42, "aC-NHn,cyc (fused rings)", "[c:1]-@[NX3;R;H0,H1:2]"),
    (43, "aC-NH-aC (different rings)", "[c:1]-[NX3H1;!R:2]-[c:3]"),
    (44, "aC-(C=N)cyc (different rings)",
     "[c:1]-!@[CX3;R;!a:2]=;@[NX2;R:3]"),
    (45, "aC-(N=CHn)cyc (fused rings)",
     "[c:1]-@[NX2;R:2]=;@[CX3;R;!a:3]"),
    (46, "aC-(CHn=N)cyc (fused rings)",
     "[c:1]-@[CX3;R;!a:2]=;@[NX2;R:3]"),
    (47, "aC-O-CHn-aC (different rings)",
     "[c:1]-[OX2;!R:2]-[CX4;!R:3]-[c:4]"),
    (48, "aC-O-aC (different rings)", "[c:1]-[OX2;!R:2]-[c:3]"),
    (49, "aC-CHn-O-CHm-aC (different rings)",
     "[c:1]-[CX4;!R:2]-[OX2;!R:3]-[CX4;!R:4]-[c:5]"),
    (50, "aC-Ocyc (fused rings)", "[c:1]-@[OX2;R:2]"),
]

_FUSED_NAMES: List[Tuple[int, str]] = [
    (51, "AROMFUSED[2]"), (52, "AROMFUSED[2]s1"), (53, "AROMFUSED[2]s2"),
    (54, "AROMFUSED[2]s2s3"), (55, "AROMFUSED[2]s1s4"),
    (56, "AROMFUSED[2]s1s2"), (57, "AROMFUSED[2]s1s3"),
    (58, "AROMFUSED[3]"), (59, "AROMFUSED[4a]"), (60, "AROMFUSED[4a]s1"),
    (61, "AROMFUSED[4a]s1s4"), (62, "AROMFUSED[4p]"),
    (63, "AROMFUSED[4p]s3s4"), (64, "PYRIDINE.FUSED[2]"),
    (65, "PYRIDINE.FUSED[2-iso]"), (66, "PYRIDINE.FUSED[4]"),
]


def _third_finder_defs() -> List[SecondOrderDef]:
    """SMARTSで書けない第3次グループ (鎖状、P/S、縮合環) の定義を作る。"""
    chain_specs = [
        (1, "HOOC-(CHn)m-COOH", _COOH_END, _COOH_END, 3),
        (2, "NHn-(CHn)m-COOH", f"[{_AMINE2}:1]", _COOH_END, 3),
        (3, "NH2-(CHn)m-OH", _NH2_END, _OH_END, 3),
        (4, "OH-(CHn)m-OH", _OH_END, _OH_END, 3),
        (9, "NH2-(CHn)m-NH2", _NH2_END, _NH2_END, 3),
        (10, "NHk-(CHn)m-NH2", f"[{_AMINE2};H0,H1:1]", _NH2_END, 3),
        (11, "SH-(CHn)m-SH", "[SX2H1:1]", "[SX2H1:1]", 3),
        (12, "NC-(CHn)m-CN", "[CX2:1]#[NX1:2]", "[CX2:1]#[NX1:2]", 3),
        (13, "COO-(CHn)m-OOC", "[OX2;!R:1][CX3:2]=[O:3]",
         "[OX2;!R:1][CX3:2]=[O:3]", 3),
        (18, "aC-(CHn)m-aC (different rings)", "[c:1]", "[c:1]", 2),
        (19, "aC-(CHn)m-CHcyc (different rings)", "[c:1]", "[CX4H1;R:1]", 1),
        (21, "CHcyc-(CHn)m-CHcyc (different rings)",
         "[CX4H1;R:1]", "[CX4H1;R:1]", 1),
        (33, "aC-CO-(CHn)m-CO-aC (different rings)",
         "[CX3;!R:1](=[O:2])-[c:3]", "[CX3;!R:1](=[O:2])-[c:3]", 1),
    ]
    defs = [
        SecondOrderDef(no, name, finder=_make_chain_finder(a, b, m))
        for no, name, a, b, m in chain_specs
    ]
    defs += [
        SecondOrderDef(5, "OH-(CHp)k-O-(CHn)m-OH",
                       finder=_make_bridge_finder(_OH_END, "[OX2H0;!R:1]")),
        SecondOrderDef(6, "OH-(CHp)k-S-(CHn)m-OH",
                       finder=_make_bridge_finder(_OH_END, "[SX2H0;!R:1]")),
        SecondOrderDef(7, "OH-(CHp)k-NHx-(CHn)m-OH",
                       finder=_make_bridge_finder(_OH_END, f"[{_AMINE2};!R:1]")),
        SecondOrderDef(8, "CHp-O-(CHn)m-OH", finder=_make_ether_alcohol_finder()),
        SecondOrderDef(40, "aC-PO(n)-aC (different rings)",
                       finder=_make_hetero_aryl_finder("P", False)),
        SecondOrderDef(41, "aC-SOn-aC (different rings)",
                       finder=_make_hetero_aryl_finder("S", True)),
    ]
    defs += [SecondOrderDef(no, name, finder=_make_fused_finder(name))
             for no, name in _FUSED_NAMES]
    return defs


THIRD_ORDER_DEFS: List[SecondOrderDef] = sorted(
    [
        SecondOrderDef(no, name, smarts,
                       predicate=_no17_range if no == 17 else None)
        for no, name, smarts in _THIRD_SMARTS
    ] + _third_finder_defs(),
    key=lambda g: g.no,
)
_THIRD_COMPILED = _compile_second(THIRD_ORDER_DEFS)


def _third_order_notes(mol: Chem.Mol, kept: List[_Occurrence]) -> List[str]:
    """フォールバック、未分類の縮合環、例外の適用を注記する。"""
    notes: List[str] = []
    for ring, kind in _aromatic_fused_rings(mol):
        name, note = _classify_fused_ring(mol, ring, kind)
        if name is None:
            notes.append(f"縮合芳香環 {sorted(ring)} は Table 8 のどれにも該当しない")
        elif note:
            notes.append(f"{name}: {note}")
    for occ in kept:
        if occ.no == 17 and any(
            mol.GetAtomWithIdx(i).GetTotalNumHs() == 2 for i in occ.identity
        ):
            notes.append("No.17: CH2 の例外を適用 (原典の範囲 n in 0..1 と不整合, D23)")
    return sorted(set(notes))


def analyze_third_order(
    smiles: str, paper_aromaticity: bool = True
) -> SecondOrderAnalysis:
    """SMILESを第3次グループへ分解する。

    完全重複の除外は第3次グループの中だけで行う (D19)。結果の解釈に注意が
    要る事項は notes に入る。
    """
    return _analyze_with_overlap(
        smiles, _THIRD_COMPILED, paper_aromaticity, _third_order_notes
    )


def count_third_order_groups(
    smiles: str, paper_aromaticity: bool = True
) -> Dict[str, int]:
    """SMILESから第3次グループ名と出現回数の辞書を返す。"""
    return analyze_third_order(smiles, paper_aromaticity).groups


# ---------------------------------------------------------------------------
# テスト・実行例
# ---------------------------------------------------------------------------
TEST_CASES: List[Tuple[str, str, Dict[str, int]]] = [
    ("ethanol", "CCO", {"CH3": 1, "CH2": 1, "OH": 1}),
    ("toluene", "Cc1ccccc1", {"aC-CH3": 1, "aCH": 5}),
    ("acetone", "CC(=O)C", {"CH3CO": 1, "CH3": 1}),
    ("acetophenone", "CC(=O)c1ccccc1",
     {"aC-CO": 1, "aCH": 5, "CH3": 1}),
    ("3,3-dimethyl-2-butanone", "CC(=O)C(C)(C)C",
     {"CH3CO": 1, "CH3": 3, "C": 1}),
    ("phenol", "Oc1ccccc1", {"aC-OH": 1, "aCH": 5}),
    ("benzoic acid", "OC(=O)c1ccccc1", {"aC-COOH": 1, "aCH": 5}),
    ("acetic acid", "CC(=O)O", {"CH3": 1, "COOH": 1}),
    ("butyl acetate", "CCCCOC(C)=O",
     {"CH3": 1, "CH2": 3, "CH3COO": 1}),
    ("methyl benzoate", "COC(=O)c1ccccc1",
     {"CH3": 1, "aC-COO": 1, "aCH": 5}),
    ("methyl phenyl ether", "COc1ccccc1",
     {"CH3": 1, "aC-O": 1, "aCH": 5}),
    ("cyclohexanol", "OC1CCCCC1",
     {"OH": 1, "CH (cyclic)": 1, "CH2 (cyclic)": 5}),
    ("naphthalene", "c1ccc2ccccc2c1",
     {"aCH": 8, "aC fused with aromatic ring": 2}),
    ("indane", "C1Cc2ccccc2C1",
     {"CH2 (cyclic)": 3, "aCH": 4,
      "aC fused with non-aromatic subring": 2}),
    ("pyrene (paper Ex.2)", "c1cc2ccc3cccc4ccc(c1)c2c34",
     {"aC fused with aromatic ring": 6, "aCH": 10}),
    ("pyridine", "c1ccncc1", {"aN in aromatic ring": 1, "aCH": 5}),
    ("N-phenyl-1,4-benzenediamine (paper Ex.1)", "Nc1ccc(Nc2ccccc2)cc1",
     {"aC-NH2": 1, "aC-NH": 1, "aC except as above": 1, "aCH": 9}),
    ("4-aminobutanol (paper Ex.3)", "NCCCCO",
     {"OH": 1, "CH2NH2": 1, "CH2": 3}),
    ("1,9-nonanediol (paper Ex.6)", "OCCCCCCCCCO",
     {"OH": 2, "CH2": 9}),
    ("styrene", "C=Cc1ccccc1", {"aC-CH=CH2": 1, "aCH": 5}),
    ("diethyl ether", "CCOCC", {"CH3": 2, "CH2": 1, "CH2O": 1}),
    ("propionitrile", "CCC#N", {"CH3": 1, "CH2CN": 1}),
    ("chlorobenzene", "Clc1ccccc1", {"aC-Cl": 1, "aCH": 5}),
    ("acetamide", "CC(N)=O", {"CH3": 1, "CONH2": 1}),
    ("N-methylacetamide", "CC(=O)NC", {"CH3": 1, "CONHCH3": 1}),
    ("N-ethylacetamide", "CCNC(C)=O", {"CH3": 2, "CONHCH2": 1}),
    ("N,N-dimethylacetamide", "CC(=O)N(C)C",
     {"CH3": 1, "CON(CH3)2": 1}),
    ("N-methyl-N-ethylacetamide", "CCN(C)C(C)=O",
     {"CH3": 2, "CONCH3CH2": 1}),
    ("N,N-diethylacetamide", "CCN(CC)C(C)=O",
     {"CH3": 3, "CON(CH2)2": 1}),
    ("diacetamide", "CC(=O)NC(C)=O", {"CH3": 2, "CONHCO": 1}),
    ("methyldiacetamide", "CC(=O)N(C)C(C)=O", {"CH3": 3, "CONCO": 1}),
    ("N,N-diethylformamide", "CCN(CC)C=O",
     {"CH3": 2, "HCON(CH2)2": 1}),
    ("N-ethylformamide", "CCNC=O", {"CH3": 1, "HCONHCH2": 1}),
    ("benzamide", "NC(=O)c1ccccc1", {"aC-CONH2": 1, "aCH": 5}),
    ("N-methylbenzamide", "CNC(=O)c1ccccc1",
     {"aC-CONH": 1, "CH3": 1, "aCH": 5}),
    ("N-phenylformamide", "O=CNc1ccccc1",
     {"aC-NH(CO)H": 1, "aCH": 5}),
    ("N-methyl-N-phenylformamide", "CN(C=O)c1ccccc1",
     {"aC-N(CO)H": 1, "CH3": 1, "aCH": 5}),
    ("N-(2-methylphenyl)acetamide", "CC(=O)Nc1ccccc1C",
     {"aC-NHCO": 1, "CH3": 1, "aC-CH3": 1, "aCH": 4}),
    ("N-methyl-N-phenylacetamide", "CC(=O)N(C)c1ccccc1",
     {"aC-NCO (amide)": 1, "CH3": 2, "aCH": 5}),
    ("methylurea", "CNC(N)=O", {"NH2CONH": 1, "CH3": 1}),
    ("N,N-dimethylurea", "CN(C)C(N)=O", {"NH2CON": 1, "CH3": 2}),
    ("1,3-dimethylurea", "CNC(=O)NC", {"NHCONH": 1, "CH3": 2}),
    ("trimethylurea", "CNC(=O)N(C)C", {"NHCON": 1, "CH3": 3}),
    ("tetramethylurea", "CN(C)C(=O)N(C)C", {"NCON": 1, "CH3": 4}),
    ("phenylurea", "NC(=O)Nc1ccccc1", {"aC-NHCONH2": 1, "aCH": 5}),
    ("N,N'-diphenylurea", "O=C(Nc1ccccc1)Nc1ccccc1",
     {"aC-NHCONH": 1, "aC except as above": 1, "aCH": 10}),
    ("1-nitropropane", "CCC[N+](=O)[O-]",
     {"CH3": 1, "CH2": 1, "CH2NO2": 1}),
    ("1-nitropropane (N(=O)=O)", "CCCN(=O)=O",
     {"CH3": 1, "CH2": 1, "CH2NO2": 1}),
    ("2-nitropropane", "CC(C)[N+](=O)[O-]", {"CH3": 2, "CHNO2": 1}),
    ("2-methyl-2-nitropropane", "CC(C)(C)[N+](=O)[O-]",
     {"CH3": 3, "CNO2": 1}),
    ("nitrobenzene", "[O-][N+](=O)c1ccccc1", {"aC-NO2": 1, "aCH": 5}),
    ("nitrocyclohexane", "[O-][N+](=O)C1CCCCC1",
     {"NO2 except as above": 1, "CH (cyclic)": 1, "CH2 (cyclic)": 5}),
    ("n-butyl nitrate", "CCCCO[N+](=O)[O-]",
     {"CH3": 1, "CH2": 3, "ONO2": 1}),
    ("n-butyl nitrite", "CCCCON=O", {"CH3": 1, "CH2": 3, "ONO": 1}),
    ("1,1-dichloroethane", "CC(Cl)Cl", {"CH3": 1, "CHCl2": 1}),
    ("2,2-dichloropropane", "CC(C)(Cl)Cl", {"CH3": 2, "CCl2": 1}),
    ("1,1,1-trichloroethane", "CC(Cl)(Cl)Cl", {"CH3": 1, "CCl3": 1}),
    ("1,1-difluoroethane", "CC(F)F", {"CH3": 1, "CHF2": 1}),
    ("hexafluoroethane", "FC(F)(F)C(F)(F)F", {"CF3": 2}),
    ("perfluorohexane", "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
     {"CF3": 2, "CF2": 4}),
    ("tetrachloro-1,2-difluoroethane", "ClC(Cl)(F)C(Cl)(Cl)F",
     {"CCl2F": 2}),
    ("1-chloro-1,2,2,2-tetrafluoroethane", "FC(Cl)C(F)(F)F",
     {"HCClF": 1, "CF3": 1}),
    ("iodoethane", "CCI", {"CH3": 1, "CH2": 1, "I except as above": 1}),
    ("ethanethiol", "CCS", {"CH3": 1, "CH2SH": 1}),
    ("2-propanethiol", "CC(C)S", {"CH3": 2, "CHSH": 1}),
    ("2-methyl-2-propanethiol", "CC(C)(C)S", {"CH3": 3, "CSH": 1}),
    ("benzenethiol", "Sc1ccccc1", {"aC-SH": 1, "aCH": 5}),
    ("cyclohexanethiol", "SC1CCCCC1",
     {"SH except as above": 1, "CH (cyclic)": 1, "CH2 (cyclic)": 5}),
    ("dimethyl sulfide", "CSC", {"CH3S": 1, "CH3": 1}),
    ("diethyl sulfide", "CCSCC", {"CH3": 2, "CH2": 1, "CH2S": 1}),
    ("diisopropyl sulfide", "CC(C)SC(C)C",
     {"CH3": 4, "CH": 1, "CHS": 1}),
    ("di-tert-butyl sulfide", "CC(C)(C)SC(C)(C)C",
     {"CH3": 6, "C": 1, "CS": 1}),
    ("phenyl methyl sulfide", "CSc1ccccc1",
     {"aC-S-": 1, "CH3": 1, "aCH": 5}),
    ("dimethyl sulfoxide", "CS(C)=O", {"CH3": 2, "SO": 1}),
    ("phenyl methyl sulfoxide", "CS(=O)c1ccccc1",
     {"aC-SO": 1, "CH3": 1, "aCH": 5}),
    ("dimethyl sulfone", "CS(C)(=O)=O", {"CH3": 2, "SO2": 1}),
    ("diphenyl sulfone", "O=S(=O)(c1ccccc1)c1ccccc1",
     {"aC-SO2": 1, "aC except as above": 1, "aCH": 10}),
    ("dimethyl sulfite", "COS(=O)OC", {"CH3": 2, "SO3 (sulfite)": 1}),
    ("methyl methanesulfonate", "CS(=O)(=O)OC",
     {"CH3": 2, "SO3 (sulfonate)": 1}),
    ("dimethyl sulfate", "COS(=O)(=O)OC", {"CH3": 2, "SO4 (sulfate)": 1}),
    ("dimethylphosphine", "C[PH]C", {"CH3": 2, "PH (phosphine)": 1}),
    ("trimethylphosphine", "CP(C)C", {"CH3": 3, "P (phosphine)": 1}),
    ("triethyl phosphite", "CCOP(OCC)OCC",
     {"CH3": 3, "CH2": 3, "PO3 (phosphite)": 1}),
    ("dimethyl phosphonate", "CO[PH](=O)OC",
     {"CH3": 2, "PHO3 (phosphonate)": 1}),
    ("trimethyl phosphonate", "CP(=O)(OC)OC",
     {"CH3": 3, "PO3 (phosphonate)": 1}),
    ("diethyl phosphate", "CCOP(=O)(O)OCC",
     {"CH3": 2, "CH2": 2, "PHO4 (phosphate)": 1}),
    ("trimethyl phosphate", "COP(=O)(OC)OC",
     {"CH3": 3, "PO4 (phosphate)": 1}),
    ("triphenyl phosphate", "O=P(Oc1ccccc1)(Oc1ccccc1)Oc1ccccc1",
     {"aC-PO4": 1, "aC except as above": 2, "aCH": 15}),
    ("triphenylphosphine", "c1ccc(cc1)P(c1ccccc1)c1ccccc1",
     {"aC-P": 1, "aC except as above": 2, "aCH": 15}),
    ("tetrahydrothiophene", "C1CCSC1", {"S (cyclic)": 1, "CH2 (cyclic)": 4}),
    ("sulfolane", "O=S1(=O)CCCC1", {"SO2 (cyclic)": 1, "CH2 (cyclic)": 4}),
    ("3-sulfolene", "O=S1(=O)CC=CC1",
     {"SO2 (cyclic)": 1, "CH2 (cyclic)": 2, "CH=CH (cyclic)": 1}),
    ("propionaldehyde oxime", "CCC=NO", {"CH3": 1, "CH2": 1, "CHNOH": 1}),
    ("diethyl ketoxime", "CCC(CC)=NO", {"CH3": 2, "CH2": 2, "CNOH": 1}),
    ("benzaldoxime", "ON=Cc1ccccc1", {"aC-CHNOH": 1, "aCH": 5}),
    ("2-ethoxyethanol", "CCOCCO", {"CH3": 1, "CH2": 1, "OCH2CH2OH": 1}),
    ("2-ethoxy-1-propanol", "CCOC(C)CO",
     {"CH3": 2, "CH2": 1, "OCHCH2OH": 1}),
    ("1-methoxy-2-propanol", "COCC(C)O", {"CH3": 2, "OCH2CHOH": 1}),
    ("tert-butyl hydroperoxide", "CC(C)(C)OO",
     {"CH3": 3, "C": 1, "-O-OH": 1}),
    ("diethyl carbonate", "CCOC(=O)OCC",
     {"CH3": 2, "CH2": 2, "CO3 (carbonate)": 1}),
    ("1,2-butadiene", "CC=C=C", {"CH3": 1, "CH2=C=CH": 1}),
    ("3-methyl-1,2-butadiene", "CC(C)=C=C", {"CH3": 2, "CH2=C=C": 1}),
    ("2,3-pentadiene", "CC=C=CC", {"CH3": 2, "CH=C=CH": 1}),
    ("acetaldazine", "CC=NN=CC", {"CH3": 2, "CH=N": 2}),
    ("ketazine", "CC(C)=NN=C(C)C", {"CH3": 4, "C=N": 2}),
    ("ethyl isocyanate", "CCN=C=O", {"CH3": 1, "CH2NCO": 1}),
    ("isopropyl isocyanate", "CC(C)N=C=O", {"CH3": 2, "CHNCO": 1}),
    ("tert-butyl isocyanate", "CC(C)(C)N=C=O", {"CH3": 3, "CNCO": 1}),
    ("phenyl isocyanate", "O=C=Nc1ccccc1", {"aC-NCO (isocyanate)": 1, "aCH": 5}),
    ("ethyl oxirane", "CCC1CO1", {"CH3": 1, "CH2": 1, "C2H3O": 1}),
    ("2,2-dimethyloxirane", "CC1(C)CO1", {"CH3": 2, "C2H2O": 1}),
    ("trimethyloxirane", "CC1OC1(C)C", {"CH3": 3, "C2O": 1}),
    ("1-pyrroline", "C1CC=NC1", {"CH=N (cyclic)": 1, "CH2 (cyclic)": 3}),
    ("2-methyl-1-pyrroline", "CC1=NCCC1",
     {"CH3": 1, "C=N (cyclic)": 1, "CH2 (cyclic)": 3}),
    ("2-methylthiophene", "Cc1cccs1",
     {"CH3": 1, "S (cyclic)": 1, "CH=C (cyclic)": 1,
      "CH=CH (cyclic)": 1}),
    ("furan", "c1ccoc1", {"O (cyclic)": 1, "CH=CH (cyclic)": 2}),
    ("pyrrole", "c1cc[nH]c1", {"NH (cyclic)": 1, "CH=CH (cyclic)": 2}),
    ("imidazole", "c1c[nH]cn1",
     {"NH (cyclic)": 1, "CH=N (cyclic)": 1, "CH=CH (cyclic)": 1}),
    ("2-methylimidazole", "Cc1ncc[nH]1",
     {"CH3": 1, "C=N (cyclic)": 1, "CH=CH (cyclic)": 1,
      "NH (cyclic)": 1}),
    ("1-hexene", "C=CCCCC", {"CH2=CH": 1, "CH2": 3, "CH3": 1}),
    ("2-hexene", "CC=CCCC", {"CH3": 2, "CH=CH": 1, "CH2": 2}),
    ("2-methyl-1-butene", "C=C(C)CC", {"CH2=C": 1, "CH3": 2, "CH2": 1}),
    ("2-methyl-2-butene", "CC=C(C)C", {"CH3": 3, "CH=C": 1}),
    ("2,3-dimethyl-2-butene", "CC(C)=C(C)C", {"CH3": 4, "C=C": 1}),
    ("1-pentyne", "C#CCCC", {"CH#C": 1, "CH2": 2, "CH3": 1}),
    ("3-decyne", "CCC#CCCCCCC", {"CH3": 2, "CH2": 6, "C#C": 1}),
    ("ethylbenzene", "CCc1ccccc1", {"CH3": 1, "aC-CH2": 1, "aCH": 5}),
    ("cumene", "CC(C)c1ccccc1", {"CH3": 2, "aC-CH": 1, "aCH": 5}),
    ("tert-butylbenzene", "CC(C)(C)c1ccccc1",
     {"CH3": 3, "aC-C": 1, "aCH": 5}),
    ("1-propenylbenzene", "CC=Cc1ccccc1",
     {"CH3": 1, "aC-CH=CH": 1, "aCH": 5}),
    ("alpha-methylstyrene", "C=C(C)c1ccccc1",
     {"aC-C=CH2": 1, "CH3": 1, "aCH": 5}),
    ("phenylacetylene", "C#Cc1ccccc1", {"aC-C#CH": 1, "aCH": 5}),
    ("1-phenyl-1-propyne", "CC#Cc1ccccc1",
     {"CH3": 1, "aC-C#C": 1, "aCH": 5}),
    ("benzaldehyde", "O=Cc1ccccc1", {"aC-CHO": 1, "aCH": 5}),
    ("phenyl formate", "O=COc1ccccc1", {"aC-OOCH": 1, "aCH": 5}),
    ("phenyl acetate", "CC(=O)Oc1ccccc1",
     {"CH3": 1, "aC-OOC": 1, "aCH": 5}),
    ("N,N-dimethylaniline", "CN(C)c1ccccc1",
     {"aC-N": 1, "CH3": 2, "aCH": 5}),
    ("benzonitrile", "N#Cc1ccccc1", {"aC-CN": 1, "aCH": 5}),
    ("fluorobenzene", "Fc1ccccc1", {"aC-F": 1, "aCH": 5}),
    ("iodobenzene", "Ic1ccccc1", {"aC-I": 1, "aCH": 5}),
    ("bromobenzene", "Brc1ccccc1", {"aC-Br": 1, "aCH": 5}),
    ("3-pentanone", "CCC(=O)CC", {"CH3": 2, "CH2CO": 1, "CH2": 1}),
    ("2,4-dimethyl-3-pentanone", "CC(C)C(=O)C(C)C",
     {"CH3": 4, "CH": 1, "CHCO": 1}),
    ("2,2,4,4-tetramethyl-3-pentanone", "CC(C)(C)C(=O)C(C)(C)C",
     {"CH3": 6, "C": 1, "CCO": 1}),
    ("hexanal", "CCCCCC=O", {"CH3": 1, "CH2": 4, "CHO": 1}),
    ("methyl butyrate", "CCCC(=O)OC", {"CH3": 2, "CH2": 1, "CH2COO": 1}),
    ("ethyl isobutyrate", "CCOC(=O)C(C)C",
     {"CH3": 3, "CH2": 1, "CHCOO": 1}),
    ("ethyl pivalate", "CCOC(=O)C(C)(C)C",
     {"CH3": 4, "CH2": 1, "CCOO": 1}),
    ("propyl formate", "CCCOC=O", {"CH3": 1, "CH2": 2, "HCOO": 1}),
    ("ethyl acrylate", "C=CC(=O)OCC",
     {"CH3": 1, "CH2": 1, "CH2=CH": 1, "COO except as above": 1}),
    ("methyl butyl ether", "CCCCOC", {"CH3": 1, "CH2": 3, "CH3O": 1}),
    ("di-sec-butyl ether", "CCC(C)OC(C)CC",
     {"CH3": 4, "CH2": 2, "CH-O": 1, "CH": 1}),
    ("di-tert-butyl ether", "CC(C)(C)OC(C)(C)C",
     {"CH3": 6, "C": 1, "C-O": 1}),
    ("sec-butylamine", "CCC(C)N", {"CH3": 2, "CH2": 1, "CHNH2": 1}),
    ("tert-butylamine", "CC(C)(C)N", {"CH3": 3, "CNH2": 1}),
    ("dimethylamine", "CNC", {"CH3": 1, "CH3NH": 1}),
    ("dipropylamine", "CCCNCCC", {"CH3": 2, "CH2": 3, "CH2NH": 1}),
    ("diisopropylamine", "CC(C)NC(C)C", {"CH3": 4, "CH": 1, "CHNH": 1}),
    ("methyldiethylamine", "CCN(C)CC", {"CH3N": 1, "CH3": 2, "CH2": 2}),
    ("triethylamine", "CCN(CC)CC", {"CH3": 3, "CH2": 2, "CH2N": 1}),
    ("cyclohexylamine", "NC1CCCCC1",
     {"NH2 except as above": 1, "CH (cyclic)": 1, "CH2 (cyclic)": 5}),
    ("isobutyronitrile", "CC(C)C#N", {"CH3": 2, "CHCN": 1}),
    ("pivalonitrile", "CC(C)(C)C#N", {"CH3": 3, "CCN": 1}),
    ("acrylonitrile", "C=CC#N", {"CH2=CH": 1, "CN except as above": 1}),
    ("1-chlorobutane", "CCCCCl", {"CH3": 1, "CH2": 2, "CH2Cl": 1}),
    ("2-chloropropane", "CC(C)Cl", {"CH3": 2, "CHCl": 1}),
    ("tert-butyl chloride", "CC(C)(C)Cl", {"CH3": 3, "CCl": 1}),
    ("N-isopropylacetamide", "CC(=O)NC(C)C",
     {"CH3": 3, "CH": 1, "NHCO except as above": 1}),
    ("1-fluorobutane", "CCCCF", {"CH3": 1, "CH2": 2, "CH2F": 1}),
    ("2-fluorobutane", "CCC(C)F", {"CH3": 2, "CH2": 1, "CHF": 1}),
    ("tert-butyl fluoride", "CC(C)(C)F", {"CH3": 3, "CF": 1}),
    ("1,2-dichlorotetrafluoroethane", "FC(F)(Cl)C(F)(F)Cl", {"CClF2": 2}),
    ("bromoethane", "CCBr", {"CH3": 1, "CH2": 1, "Br except as above": 1}),
    ("benzyl fluoride", "FCc1ccccc1",
     {"aC-CH2": 1, "F except as above": 1, "aCH": 5}),
    ("ethyl chloroacetate", "ClCC(=O)OCC",
     {"CH3": 1, "CH2": 1, "CH2COO": 1, "Cl except as above": 1}),
    ("1,1-dimethylcyclohexane", "CC1(C)CCCCC1",
     {"CH3": 2, "C (cyclic)": 1, "CH2 (cyclic)": 5}),
    ("1,2-dimethylcyclopentene", "CC1=C(C)CCC1",
     {"CH3": 2, "C=C (cyclic)": 1, "CH2 (cyclic)": 3}),
    ("methylenecyclohexane", "C=C1CCCCC1",
     {"CH2=C (cyclic)": 1, "CH2 (cyclic)": 5}),
    ("N-methylpyrrolidine", "CN1CCCC1",
     {"CH3": 1, "N (cyclic)": 1, "CH2 (cyclic)": 4}),
    ("cyclobutanone", "O=C1CCC1", {"CO (cyclic)": 1, "CH2 (cyclic)": 3}),
    # --- A1: 環状アミド・エステル・カーボネート・イミドは環原子単位 ---
    ("gamma-butyrolactone", "O=C1CCCO1",
     {"CO (cyclic)": 1, "O (cyclic)": 1, "CH2 (cyclic)": 3}),
    ("ethylene carbonate", "O=C1OCCO1",
     {"CO (cyclic)": 1, "O (cyclic)": 2, "CH2 (cyclic)": 2}),
    ("succinimide", "O=C1CCC(=O)N1",
     {"CO (cyclic)": 2, "NH (cyclic)": 1, "CH2 (cyclic)": 2}),
    ("N-methylpyrrolidone", "CN1CCCC1=O",
     {"CH3": 1, "N (cyclic)": 1, "CO (cyclic)": 1, "CH2 (cyclic)": 3}),
    # --- A2: 論文流の芳香環定義 (6員C/N環のみ芳香族) ---
    ("coumarin", "O=c1ccc2ccccc2o1",
     {"CO (cyclic)": 1, "CH=CH (cyclic)": 1, "O (cyclic)": 1, "aCH": 4,
      "aC fused with non-aromatic subring": 2}),
    ("chromone", "O=c1ccoc2ccccc12",
     {"CO (cyclic)": 1, "CH=CH (cyclic)": 1, "O (cyclic)": 1, "aCH": 4,
      "aC fused with non-aromatic subring": 2}),
    ("4-pyranone", "O=c1ccocc1",
     {"CO (cyclic)": 1, "CH=CH (cyclic)": 2, "O (cyclic)": 1}),
    ("2-pyridone", "O=c1cccc[nH]1",
     {"CO (cyclic)": 1, "CH=CH (cyclic)": 2, "NH (cyclic)": 1}),
    ("quinoline", "c1ccc2ncccc2c1",
     {"aN in aromatic ring": 1, "aCH": 7,
      "aC fused with aromatic ring": 2}),
    ("indole", "c1ccc2[nH]ccc2c1",
     {"CH=CH (cyclic)": 1, "NH (cyclic)": 1, "aCH": 4,
      "aC fused with non-aromatic subring": 2}),
    ("benzofuran", "c1ccc2occc2c1",
     {"CH=CH (cyclic)": 1, "O (cyclic)": 1, "aCH": 4,
      "aC fused with non-aromatic subring": 2}),
    ("carbazole", "c1ccc2c(c1)[nH]c1ccccc12",
     {"NH (cyclic)": 1, "aCH": 8,
      "aC fused with non-aromatic subring": 4}),
    ("adenine", "Nc1ncnc2[nH]cnc12",
     {"aC-NH2": 1, "CH=N (cyclic)": 1, "NH (cyclic)": 1,
      "aN in aromatic ring": 2, "aCH": 1,
      "aC fused with non-aromatic subring": 2}),
    ("caffeine", "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
     {"CO (cyclic)": 2, "CH=N (cyclic)": 1, "C=C (cyclic)": 1,
      "CH3": 3, "N (cyclic)": 3}),
    ("azulene", "c1ccc2cccc2cc1",
     {"CH=CH (cyclic)": 3, "CH=C (cyclic)": 2}),
    # --- A3: OH は C または S に結合 ---
    ("benzenesulfonic acid", "OS(=O)(=O)c1ccccc1",
     {"aC-SO2": 1, "OH": 1, "aCH": 5}),
]


# 適用範囲外 (判断記録 D1, D5, D9): 意図どおり未割当になることを記録する。
# (名前, SMILES, 期待グループ, 未割当原子の元素記号 (ソート済み))
LIMIT_CASES: List[Tuple[str, str, Dict[str, int], List[str]]] = [
    ("N-acetylpyrrolidine (D5)", "CC(=O)N1CCCC1",
     {"CH3": 1, "CH2 (cyclic)": 4, "N (cyclic)": 1}, ["C", "O"]),
    ("N-benzyl-N-methylacetamide (D1)", "CC(=O)N(C)Cc1ccccc1",
     {"aC-CH2": 1, "CH3": 2, "aCH": 5}, ["C", "N", "O"]),
    ("urea (D9)", "NC(N)=O", {}, ["C", "N", "N", "O"]),
    ("formamide (D9)", "NC=O", {}, ["C", "N", "O"]),
    ("triisopropylamine (D9)", "CC(C)N(C(C)C)C(C)C",
     {"CH3": 6, "CH": 3}, ["N"]),
    ("di-tert-butylamine (D9)", "CC(C)(C)NC(C)(C)C",
     {"CH3": 6, "C": 2}, ["N"]),
    ("acetyl chloride (D9)", "CC(Cl)=O",
     {"CH3": 1, "Cl except as above": 1}, ["C", "O"]),
    ("S-methyl thioacetate (D9)", "CC(=O)SC",
     {"CH3S": 1, "CH3": 1}, ["C", "O"]),
    ("dimethyl disulfide (D9)", "CSSC", {"CH3": 2}, ["S", "S"]),
    ("pyridine N-oxide (D4)", "[O-][n+]1ccccc1", {"aCH": 5}, ["N", "O"]),
]


@dataclass(frozen=True)
class AppendixBCase:
    """論文 Appendix B の第1次近似の再現テスト。

    Attributes:
        label: 例題名。
        smiles: 化合物のSMILES。
        t0: 推算式 T = T0 * ln(sum) の定数 (Table 2: Tb0 または Tm0)。
        groups: 論文記載のグループ出現回数 (本実装のグループ名で表記)。
        contributions: 論文記載の寄与値 (Table 6)。
        paper_sum: 論文記載の sum(N_i * C_i)。
        paper_estimate: 論文記載の第1次近似の推算値 [K]。
    """

    label: str
    smiles: str
    t0: float
    groups: Dict[str, int]
    contributions: Dict[str, float]
    paper_sum: float
    paper_estimate: float


TB0 = 222.543  # Table 2
TM0 = 147.450  # Table 2

# 論文本文の誤植 (D8): Ex.5 は ln(25.1421) と記載されているが 469.19 K は
# ln(24.0959) の値、Ex.6 は 222.543 と記載されているが 303.66 K は
# Tm0 = 147.450 の値。ここでは数値が一致する正しい式を使う。
APPENDIX_B_CASES: List[AppendixBCase] = [
    AppendixBCase(
        "Ex.1 N-phenyl-1,4-benzenediamine (Tb)",
        "Nc1ccc(Nc2ccccc2)cc1", TB0,
        {"aC-NH2": 1, "aC-NH": 1, "aC except as above": 1, "aCH": 9},
        {"aC-NH2": 3.8298, "aC-NH": 2.9230, "aC except as above": 1.5468,
         "aCH": 0.8365},
        15.8281, 614.62),
    AppendixBCase(
        "Ex.2 pyrene (Tb)",
        "c1cc2ccc3cccc4ccc(c1)c2c34", TB0,
        {"aC fused with aromatic ring": 6, "aCH": 10},
        {"aC fused with aromatic ring": 1.7324, "aCH": 0.8365},
        18.7593, 652.43),
    AppendixBCase(
        "Ex.3 4-aminobutanol (Tb)",
        "NCCCCO", TB0,
        {"OH": 1, "CH2NH2": 1, "CH2": 3},
        {"OH": 2.5670, "CH2NH2": 2.7987, "CH2": 0.7141},
        7.508, 448.64),
    AppendixBCase(
        "Ex.4 dicoumarol (Tm)",
        "OC1=C(Cc2c(O)c3ccccc3oc2=O)C(=O)Oc2ccccc12", TM0,
        {"OH": 2, "aC fused with non-aromatic subring": 4, "aCH": 8,
         "C=C (cyclic)": 2, "CO (cyclic)": 2, "O (cyclic)": 2, "CH2": 1},
        {"OH": 2.7888, "aC fused with non-aromatic subring": 1.2065,
         "aCH": 0.5860, "C=C (cyclic)": 0.3048, "CO (cyclic)": 3.2119,
         "O (cyclic)": 1.3828, "CH2": 0.2515},
        25.1421, 475.46),
    AppendixBCase(
        "Ex.5 flutemazepan (Tm)",
        "CN1C(=O)C(O)N=C(c2ccccc2F)c2cc(Cl)ccc21", TM0,
        {"OH": 1, "aC-Cl": 1, "aC-F": 1, "CH3": 1, "aC except as above": 1,
         "aC fused with non-aromatic subring": 2, "aCH": 7,
         "CH (cyclic)": 1, "CO (cyclic)": 1, "N (cyclic)": 1,
         "C=N (cyclic)": 1},
        {"OH": 2.7888, "aC-Cl": 1.7134, "aC-F": 0.9782, "CH3": 0.6953,
         "aC except as above": 0.9176,
         "aC fused with non-aromatic subring": 1.2065, "aCH": 0.5860,
         "CH (cyclic)": 0.0335, "CO (cyclic)": 3.2119, "N (cyclic)": 0.6040,
         "C=N (cyclic)": 6.6382},
        24.0959, 469.19),
    AppendixBCase(
        "Ex.6 1,9-nonanediol (Tm)",
        "OCCCCCCCCCO", TM0,
        {"OH": 2, "CH2": 9},
        {"OH": 2.7888, "CH2": 0.2515},
        7.8411, 303.66),
]

SUM_TOL = 1e-3       # sum(N_i C_i) の許容誤差
ESTIMATE_TOL = 0.02  # 推算値 [K] の許容誤差 (論文の丸め誤差)


def check_appendix_b(case: AppendixBCase) -> List[str]:
    """Appendix B の1例を検証し、不一致メッセージのリストを返す。"""
    problems: List[str] = []
    res = analyze_smiles(case.smiles)
    if res.groups != case.groups or not res.is_complete:
        problems.append(
            f"{case.label}: groups {res.groups} != {case.groups} "
            f"(unassigned={res.unassigned}, warnings={res.warnings})"
        )
        return problems
    total = sum(n * case.contributions[g] for g, n in res.groups.items())
    if abs(total - case.paper_sum) > SUM_TOL:
        problems.append(f"{case.label}: sum {total:.4f} != {case.paper_sum}")
    estimate = case.t0 * math.log(total)
    if abs(estimate - case.paper_estimate) > ESTIMATE_TOL:
        problems.append(
            f"{case.label}: estimate {estimate:.2f} != {case.paper_estimate}"
        )
    return problems


def collect_failures() -> List[str]:
    """全テスト (TEST_CASES, LIMIT_CASES, Appendix B) の失敗名を返す。"""
    failures: List[str] = []
    for name, smi, expected in TEST_CASES:
        res = analyze_smiles(smi)
        if res.groups != expected or not res.is_complete:
            failures.append(f"TEST_CASES: {name}")
    for name, smi, expected, unassigned in LIMIT_CASES:
        res = analyze_smiles(smi)
        symbols = sorted(sym for _, sym in res.unassigned)
        if res.groups != expected or symbols != unassigned:
            failures.append(f"LIMIT_CASES: {name}")
    for case in APPENDIX_B_CASES:
        failures.extend(check_appendix_b(case))
    return failures


def run_shuffle_test(n_trials: int = 30, seed: int = 0) -> List[str]:
    """GROUP_DEFS の並びをランダムに入れ替えても結果が変わらないことを検証する。

    結果が定義順に依存すると (D7)、失敗したケース名が返る。
    実行後はモジュールの状態を元に戻す。
    """
    rng = random.Random(seed)
    original = list(GROUP_DEFS)
    failures = set()
    try:
        for _ in range(n_trials):
            shuffled = original[:]
            rng.shuffle(shuffled)
            _COMPILED[:] = _compile(shuffled)
            failures.update(collect_failures())
    finally:
        _COMPILED[:] = _compile(original)
    return sorted(failures)


# ---------------------------------------------------------------------------
# 第2次グループのテスト
# ---------------------------------------------------------------------------
# Table 7 の「例」欄の化合物。(グループ番号, 化合物名, SMILES, 期待出現回数)。
# 期待値は論文記載の出現回数。ただし No.24 は構造に基づく値 (D15)。
SECOND_ORDER_CASES: List[Tuple[int, str, str, int]] = [
    (1, "2-methylpentane", "CC(C)CCC", 1),
    (2, "2,2,4,4-tetramethylpentane", "CC(C)(C)CC(C)(C)C", 2),
    (3, "2,3,4-trimethylpentane", "CC(C)C(C)C(C)C", 2),
    (4, "2,2,3,4,4-pentamethylpentane", "CC(C)(C)C(C)C(C)(C)C", 2),
    (5, "2,2,3,3,4,4-hexamethylpentane", "CC(C)(C)C(C)(C)C(C)(C)C", 2),
    (6, "1,3-butadiene", "C=CC=C", 1),
    (7, "2-methyl-2-butene", "CC=C(C)C", 3),
    (8, "1,4-pentadiene", "C=CCC=C", 2),
    (9, "3-methyl-1-butene", "C=CC(C)C", 1),
    (10, "2-methylbutyraldehyde", "CCC(C)C=O", 1),
    (11, "2-pentanone", "CCCC(C)=O", 1),
    (12, "3-methyl-2-pentanone", "CCC(C)C(C)=O", 1),
    (13, "2-methylbutanoic acid", "CCC(C)C(=O)O", 1),
    (14, "isopropyl acetate", "CC(=O)OC(C)C", 1),
    (15, "propanoic anhydride", "CCC(=O)OC(=O)CC", 1),
    (16, "2-butanol", "CCC(C)O", 1),
    (17, "2-methyl-2-butanol", "CCC(C)(C)O", 1),
    (18, "3-hydroxy-2-butanone", "CC(=O)C(C)O", 1),
    (19, "2-hydroxypropionitrile", "CC(O)C#N", 1),
    (20, "ethyl lactate", "CCOC(=O)C(C)O", 1),
    (21, "ethylene glycol", "OCCO", 1),
    (22, "2-amino-1-butanol", "CCC(N)CO", 1),
    (23, "ethylenediamine", "NCCN", 1),
    (24, "diethylenetriamine", "NCCNCCN", 2),  # 論文は (1) (D15)
    (25, "butanediamide", "NC(=O)CCC(N)=O", 1),
    (26, "l-alanine", "CC(N)C(=O)O", 1),
    (27, "malonic acid", "OC(=O)CC(O)=O", 1),
    (28, "succinic acid", "OC(=O)CCC(O)=O", 1),
    (29, "2-hydroxyisobutyric acid", "CC(C)(O)C(O)=O", 1),
    (30, "beta-alanine", "NCCC(=O)O", 1),
    (31, "methoxyacetic acid", "COCC(=O)O", 1),
    (32, "2-mercaptopropionic acid", "CC(S)C(=O)O", 1),
    (33, "beta-thiolactic acid", "SCCC(=O)O", 1),
    (34, "1,2-dicyanoethane", "N#CCCC#N", 1),
    (35, "3-hydroxypropanenitrile", "OCCC#N", 1),
    (36, "1,2-ethanedithiol", "SCCS", 1),
    (37, "ethylene glycol diacetate", "CC(=O)OCCOC(C)=O", 1),
    (38, "dimethyl succinate", "COC(=O)CCC(=O)OC", 1),
    (39, "methyl cyanoacetate", "COC(=O)CC#N", 1),
    (40, "methyl acetoacetate", "COC(=O)CC(C)=O", 1),
    (41, "ethyl vinyl ether", "CCOC=C", 1),
    (42, "1-fluoro-1-propene", "CC=CF", 1),
    (43, "1-bromo-1-propene", "CC=CBr", 1),
    (44, "1-iodo-1-propene", "CC=CI", 1),
    (45, "1-chloro-2-methylpropene", "CC(C)=CCl", 1),
    (46, "acrylonitrile", "C=CC#N", 1),
    (47, "ethyl acrylate", "C=CC(=O)OCC", 1),
    (48, "propenaldehyde", "C=CC=O", 1),
    (49, "acrylic acid", "C=CC(=O)O", 1),
    (50, "benzyl bromide", "BrCc1ccccc1", 1),
    (51, "benzylamine", "NCc1ccccc1", 1),
    (52, "benzyl ethyl ether", "CCOCc1ccccc1", 1),
    (53, "benzyl alcohol", "OCc1ccccc1", 1),
    (54, "benzyl cyanide", "N#CCc1ccccc1", 1),
    (55, "phenylacetaldehyde", "O=CCc1ccccc1", 1),
    (56, "phenylmethanethiol", "SCc1ccccc1", 1),
    (57, "phenylacetic acid", "OC(=O)Cc1ccccc1", 1),
    (58, "phenylacetone", "CC(=O)Cc1ccccc1", 1),
    (59, "benzyl methyl sulfide", "CSCc1ccccc1", 1),
    (60, "benzyl formate", "O=COCc1ccccc1", 1),
    (61, "phenylnitromethane", "[O-][N+](=O)Cc1ccccc1", 1),
    (62, "phenylacetamide", "NC(=O)Cc1ccccc1", 1),
    (63, "benzyl acetate", "CC(=O)OCc1ccccc1", 1),
    (64, "methyl phenylacetate", "COC(=O)Cc1ccccc1", 1),
    (65, "benzenesulfonic acid", "OS(=O)(=O)c1ccccc1", 1),
    (66, "cumene", "CC(C)c1ccccc1", 1),
    (67, "tert-butylbenzene", "CC(C)(C)c1ccccc1", 1),
    (68, "perfluorotoluene", "FC(F)(F)c1c(F)c(F)c(F)c(F)c1F", 1),
    (69, "furfural", "O=Cc1ccco1", 1),
    (70, "methyl 2-furoate", "COC(=O)c1ccco1", 1),
    (71, "2-acetylfuran", "CC(=O)c1ccco1", 1),
    (72, "1,2-dimethylcyclopentene", "CC1=C(C)CCC1", 2),
    (73, "2-ethylfuran", "CCc1ccco1", 1),
    (74, "3-cyanofuran", "N#Cc1ccoc1", 1),
    (75, "2-chlorofuran", "Clc1ccco1", 1),
    (76, "methylcyclopentane", "CC1CCCC1", 1),
    (77, "ethylcyclohexane", "CCC1CCCCC1", 1),
    (78, "isopropylcyclopentane", "CC(C)C1CCCC1", 1),
    (79, "tert-butylcyclohexane", "CC(C)(C)C1CCCCC1", 1),
    (80, "vinylcyclopentane", "C=CC1CCCC1", 1),
    (81, "limonene", "CC1=CCC(CC1)C(C)=C", 1),
    (82, "chlorocyclopentane", "ClC1CCCC1", 1),
    (83, "fluorocyclohexane", "FC1CCCCC1", 1),
    (84, "cyclohexanol", "OC1CCCCC1", 1),
    (85, "cyclohexylamine", "NC1CCCCC1", 1),
    (86, "N-methylcyclohexylamine", "CNC1CCCCC1", 1),
    (87, "N,N-dimethylcyclohexanamine", "CN(C)C1CCCCC1", 1),
    (88, "cyclohexanethiol", "SC1CCCCC1", 1),
    (89, "cyanocyclopentane", "N#CC1CCCC1", 1),
    (90, "cyclopropanecarboxylic acid", "OC(=O)C1CC1", 1),
    (91, "methyl cyclohexyl ketone", "CC(=O)C1CCCCC1", 1),
    (92, "nitrocyclohexane", "[O-][N+](=O)C1CCCCC1", 1),
    (93, "methyl cyclopentyl sulfide", "CSC1CCCC1", 1),
    (94, "cyclohexanecarboxaldehyde", "O=CC1CCCCC1", 1),
    (95, "methoxycyclohexane", "COC1CCCCC1", 1),
    (96, "cyclohexyl formate", "O=COC1CCCCC1", 1),
    (97, "ethyl cyclobutanecarboxylate", "CCOC(=O)C1CCC1", 1),
    (98, "cyclohexyl acetate", "CC(=O)OC1CCCCC1", 1),
    (99, "1,1-dimethylcyclohexane", "CC1(C)CCCCC1", 2),
    (100, "1-ethyl-1-methylcyclopentane", "CCC1(C)CCCC1", 1),
    (101, "1-methylcyclopentanol", "CC1(O)CCCC1", 1),
    (102, "N-methyl-2-pyrrolidone", "CN1CCCC1=O", 1),
    (103, "N-ethylpyrrole", "CCn1cccc1", 1),
    (104, "2-methylphenol", "Cc1ccccc1O", 1),
    (104, "2-ethyltoluene", "CCc1ccccc1C", 1),
    (105, "3-methylphenol", "Cc1cccc(O)c1", 1),
    (105, "3-ethyltoluene", "CCc1cccc(C)c1", 1),
    (106, "4-methylphenol", "Cc1ccc(O)cc1", 1),
    (106, "4-ethyltoluene", "CCc1ccc(C)cc1", 1),
    (107, "1,2,3-trimethylbenzene", "Cc1cccc(C)c1C", 1),
    (108, "1,2,4-trihydroxybenzene", "Oc1ccc(O)c(O)c1", 1),
    (109, "3,5-diethyltoluene", "CCc1cc(C)cc(CC)c1", 1),
    (110, "3-ethyl-1,2,4-trimethylbenzene", "Cc1c(C)c(CC)c(C)cc1", 1),
    (111, "1,2,3,5-tetramethylbenzene", "Cc1cc(C)c(C)c(C)c1", 1),
    (112, "1,2,4,5-tetramethylbenzene", "Cc1cc(C)c(C)cc1C", 1),
    (113, "2-methylpyridine", "Cc1ccccn1", 1),
    (114, "3-methylpyridine", "Cc1cccnc1", 1),
    (115, "4-methylpyridine", "Cc1ccncc1", 1),
    (116, "2,3-dimethylpyridine", "Cc1cccnc1C", 1),
    (117, "2,4-dimethylpyridine", "Cc1ccnc(C)c1", 1),
    (118, "2,5-dimethylpyridine", "Cc1ccc(C)nc1", 1),
    (119, "2,6-dimethylpyridine", "Cc1cccc(C)n1", 1),
    (120, "3,4-dimethylpyridine", "Cc1ccncc1C", 1),
    (121, "3,5-dimethylpyridine", "Cc1cncc(C)c1", 1),
    (122, "2,3,6-trimethylpyridine", "Cc1ccc(C)nc1C", 1),
]

# 完全重複の禁止と全体の整合を確認するケース (化合物名, SMILES, 期待する全グループ)
SECOND_ORDER_FULL_CASES: List[Tuple[str, str, Dict[str, int]]] = [
    # 論文 Section 2.2 の例: CHm-COOH は CHm(NHn)-COOH に完全に含まれる
    ("l-alanine (paper Fig.2)", "CC(N)C(=O)O", {"CHm(NHn)-COOH": 1}),
    # 論文 Section 2.2 の例: CHOH を併用しない
    ("2-amino-1-butanol (paper 2.2)", "CCC(N)CO", {"CHm(OH)CHn(NHp)": 1}),
    # 論文 Section 2.2 の例: s1s3 や s1s4 を併用しない
    ("3,5-dimethylaniline (paper 2.2)", "Cc1cc(C)cc(N)c1",
     {"AROMRINGs1s3s5": 1}),
    ("ethyl lactate (CHOH embedded)", "CCOC(=O)C(C)O", {"OH-CHn-COO": 1}),
    ("3-hydroxy-2-butanone (CHOH, CH3COCH embedded)", "CC(=O)C(C)O",
     {"CH3COCHn(OH)": 1}),
    ("1-phenylethanol (CHOH embedded)", "CC(O)c1ccccc1",
     {"aC-CHn-OH": 1}),
    ("cumene ((CH3)2CH embedded)", "CC(C)c1ccccc1", {"aC-CH(CH3)2": 1}),
    ("tert-butylbenzene ((CH3)3C embedded)", "CC(C)(C)c1ccccc1",
     {"aC-C(CH3)3": 1}),
    ("biphenyl (no ring group)", "c1ccc(cc1)-c1ccccc1", {}),
    ("naphthalene derivative (fused ring skipped)", "Cc1cccc2ccccc12", {}),
    ("toluene (no second-order group)", "Cc1ccccc1", {}),
    ("ethanol (no second-order group)", "CCO", {}),
]

# 論文 Appendix B で「第2次グループなし」とされている例
NO_SECOND_ORDER_CASES: List[Tuple[str, str]] = [
    ("Ex.2 pyrene", "c1cc2ccc3cccc4ccc(c1)c2c34"),
    ("Ex.3 4-aminobutanol", "NCCCCO"),
    ("Ex.4 dicoumarol", "OC1=C(Cc2c(O)c3ccccc3oc2=O)C(=O)Oc2ccccc12"),
    ("Ex.6 1,9-nonanediol", "OCCCCCCCCCO"),
]


@dataclass(frozen=True)
class SecondOrderAppendixBCase:
    """論文 Appendix B の第2近似の再現テスト。

    Attributes:
        label: 例題名。
        smiles: 化合物のSMILES。
        t0: 推算式 T = T0 * ln(sum) の定数 (Table 2)。
        first_sum: 論文記載の第1次の sum(N_i C_i)。
        groups: 論文記載の第2次グループ出現回数 (本実装のグループ名)。
        contributions: 論文記載の寄与値 (Table 7)。
        paper_estimate: 論文記載の第2次近似の推算値 [K]。
    """

    label: str
    smiles: str
    t0: float
    first_sum: float
    groups: Dict[str, int]
    contributions: Dict[str, float]
    paper_estimate: float


SECOND_ORDER_APPENDIX_B: List[SecondOrderAppendixBCase] = [
    SecondOrderAppendixBCase(
        "Ex.1 N-phenyl-1,4-benzenediamine (Tb)",
        "Nc1ccc(Nc2ccccc2)cc1", TB0, 15.8281,
        {"AROMRINGs1s4": 1}, {"AROMRINGs1s4": 0.1007}, 616.03),
    SecondOrderAppendixBCase(
        "Ex.5 flutemazepan (Tm)",
        "CN1C(=O)C(O)N=C(c2ccccc2F)c2cc(Cl)ccc21", TM0, 24.0959,
        {"AROMRINGs1s2": 1, ">Ncyc-CH3": 1, "CHcyc-OH": 1},
        {"AROMRINGs1s2": -0.6388, ">Ncyc-CH3": -0.0383,
         "CHcyc-OH": 1.3691},
        473.36),
]


def check_second_order_appendix_b(case: SecondOrderAppendixBCase) -> List[str]:
    """Appendix B の第2近似を検証し、不一致メッセージのリストを返す。"""
    res = analyze_second_order(case.smiles)
    if res.groups != case.groups or res.warnings:
        return [f"{case.label}: groups {res.groups} != {case.groups} "
                f"(eliminated={res.eliminated}, warnings={res.warnings})"]
    total = case.first_sum + sum(
        n * case.contributions[g] for g, n in res.groups.items()
    )
    estimate = case.t0 * math.log(total)
    if abs(estimate - case.paper_estimate) > ESTIMATE_TOL:
        return [f"{case.label}: estimate {estimate:.2f} != "
                f"{case.paper_estimate}"]
    return []


def collect_second_order_failures() -> List[str]:
    """第2次グループの全テストの失敗名を返す。"""
    names = {g.no: g.name for g in SECOND_ORDER_DEFS}
    failures: List[str] = []
    for no, label, smi, expected in SECOND_ORDER_CASES:
        got = analyze_second_order(smi).groups.get(names[no], 0)
        if got != expected:
            failures.append(f"No.{no} {names[no]} / {label}: {got} != {expected}")
    for label, smi, expected_groups in SECOND_ORDER_FULL_CASES:
        res = analyze_second_order(smi)
        if res.groups != expected_groups:
            failures.append(f"full: {label}: {res.groups} != {expected_groups}")
    for label, smi in NO_SECOND_ORDER_CASES:
        res = analyze_second_order(smi)
        if res.groups:
            failures.append(f"no-second-order: {label}: {res.groups}")
    for case in SECOND_ORDER_APPENDIX_B:
        failures.extend(check_second_order_appendix_b(case))
    return failures


TABLE7_NUMBERS = range(1, 123)  # 論文 Table 7 の第2次グループ番号


def find_unimplemented_second_order() -> List[int]:
    """Table 7 (No.1-122) のうち SECOND_ORDER_DEFS に未定義の番号を返す。"""
    implemented = {g.no for g in SECOND_ORDER_DEFS}
    return [n for n in TABLE7_NUMBERS if n not in implemented]


def find_untested_second_order() -> List[int]:
    """SECOND_ORDER_CASES に例題がない第2次グループ番号を返す。"""
    tested = {no for no, _, _, _ in SECOND_ORDER_CASES}
    return [g.no for g in SECOND_ORDER_DEFS if g.no not in tested]


def find_second_order_ties() -> Dict[str, List[Tuple[str, str]]]:
    """テスト分子のうち、原子集合が完全一致する別グループが現れたものを返す。"""
    found: Dict[str, List[Tuple[str, str]]] = {}
    smiles_list = [(label, smi) for _, label, smi, _ in SECOND_ORDER_CASES]
    smiles_list += [(label, smi) for label, smi, _ in SECOND_ORDER_FULL_CASES]
    for label, smi in smiles_list:
        ties = analyze_second_order(smi).ties
        if ties:
            found[label] = ties
    return found


def run_second_order_shuffle_test(n_trials: int = 30, seed: int = 0) -> List[str]:
    """SECOND_ORDER_DEFS の並びを入れ替えても結果が変わらないことを検証する。

    実行後はモジュールの状態を元に戻す。
    """
    rng = random.Random(seed)
    original = list(SECOND_ORDER_DEFS)
    failures = set()
    try:
        for _ in range(n_trials):
            shuffled = original[:]
            rng.shuffle(shuffled)
            _SECOND_COMPILED[:] = _compile_second(shuffled)
            failures.update(collect_second_order_failures())
    finally:
        _SECOND_COMPILED[:] = _compile_second(original)
    return sorted(failures)


# ---------------------------------------------------------------------------
# 第3次グループのテスト
# ---------------------------------------------------------------------------
# No.8 の例示 "Butoxypropanol" は構造が特定できない。範囲 (m>2) と整合する
# 3-ブトキシ-1-プロパノールと仮定する (実装上の仮定, D26)。
BUTOXYPROPANOL_ASSUMED_SMILES = "CCCCOCCCO"

# Table 8 の「例」欄の化合物。(グループ番号, 化合物名, SMILES, 期待出現回数)
THIRD_ORDER_CASES: List[Tuple[int, str, str, int]] = [
    (1, "glutaric acid", "OC(=O)CCCC(O)=O", 1),
    (2, "4-aminobutyric acid", "NCCCC(=O)O", 1),
    (3, "4-aminobutanol", "NCCCCO", 1),
    (4, "1,9-nonanediol", "OCCCCCCCCCO", 1),
    (5, "dipropylene glycol", "CC(O)COCC(C)O", 1),
    (6, "thiodiglycol", "OCCSCCO", 1),
    (7, "diethanolamine", "OCCNCCO", 1),
    (8, "butoxypropanol (assumed 3-butoxy-1-propanol)",
     BUTOXYPROPANOL_ASSUMED_SMILES, 1),
    (9, "1,5-diaminopentane", "NCCCCCN", 1),
    (10, "N,N-dimethylpropylenediamine", "CN(C)CCCN", 1),
    (11, "1,5-pentanedithiol", "SCCCCCS", 1),
    (12, "glutaronitrile", "N#CCCCC#N", 1),
    (13, "glyceryl tridodecanoate",
     "CCCCCCCCCCCC(=O)OCC(COC(=O)CCCCCCCCCCC)OC(=O)CCCCCCCCCCC", 1),
    (14, "indene", "C1=Cc2ccccc2C1", 1),
    (14, "acenaphthylene", "C1=Cc2cccc3cccc1c23", 2),
    (15, "biphenyl", "c1ccc(cc1)-c1ccccc1", 1),
    (15, "biphenylene", "c1ccc2c(c1)c1ccccc12", 2),
    (16, "cyclohexylbenzene", "C1CCC(CC1)c1ccccc1", 1),
    (17, "tetralin", "C1CCc2ccccc2C1", 2),
    (17, "indane", "C1Cc2ccccc2C1", 2),
    (18, "bibenzyl", "c1ccc(cc1)CCc1ccccc1", 1),
    (19, "1-cyclopentyl-3-phenylpropane", "C1CCC(C1)CCCc1ccccc1", 1),
    (20, "bicyclohexyl", "C1CCC(CC1)C1CCCCC1", 1),
    (21, "1,2-dicyclohexylethane", "C1CCC(CC1)CCC1CCCCC1", 1),
    (22, "octahydroindene (hexahydroindan)", "C1CCC2CCCC2C1", 2),
    (22, "decalin", "C1CCC2CCCCC2C1", 2),
    (23, "spiropentane", "C1CC12CC2", 1),
    (24, "diphenylmethane", "c1ccc(cc1)Cc1ccccc1", 1),
    (25, "1,2-diphenylethylene", "c1ccc(cc1)C=Cc1ccccc1", 1),
    (26, "1,2-difuranylethene", "C(=Cc1ccco1)c1ccco1", 1),
    (27, "difuranylmethane", "C(c1ccco1)c1ccco1", 1),
    (28, "benzophenone", "O=C(c1ccccc1)c1ccccc1", 1),
    (29, "benzyl phenyl ketone", "O=C(Cc1ccccc1)c1ccccc1", 1),
    (30, "phenyl 2-furyl ketone", "O=C(c1ccccc1)c1ccco1", 1),
    (31, "benzil", "O=C(C(=O)c1ccccc1)c1ccccc1", 1),
    (32, "phenolphthalein", "O=C1OC(c2ccc(O)cc2)(c2ccc(O)cc2)c2ccccc12", 1),
    (33, "1,4-diphenyl-1,4-butanedione", "O=C(CCC(=O)c1ccccc1)c1ccccc1", 1),
    (34, "cyclohexyl phenyl ketone", "O=C(c1ccccc1)C1CCCCC1", 1),
    (35, "benzanilide", "O=C(Nc1ccccc1)c1ccccc1", 1),
    (36, "N,N'-diphenylurea", "O=C(Nc1ccccc1)Nc1ccccc1", 1),
    (37, "N-benzoylpiperidine", "O=C(c1ccccc1)N1CCCCC1", 1),
    (38, "dibenzothiophene", "c1ccc2c(c1)sc1ccccc12", 2),
    (39, "diphenyl sulfide", "c1ccc(cc1)Sc1ccccc1", 1),
    (40, "triphenyl phosphate", "O=P(Oc1ccccc1)(Oc1ccccc1)Oc1ccccc1", 3),
    (41, "diphenyl sulfone", "O=S(=O)(c1ccccc1)c1ccccc1", 1),
    (42, "carbazole", "c1ccc2c(c1)[nH]c1ccccc12", 2),
    (43, "diphenylamine", "c1ccc(cc1)Nc1ccccc1", 1),
    (44, "3-phenylpyrazole", "c1ccc(cc1)-c1cc[nH]n1", 1),
    (45, "benzoxazole", "c1ccc2ocnc2c1", 1),
    (46, "1,2-benzisoxazole", "c1ccc2c(c1)cno2", 1),
    (47, "benzyl phenyl ether", "c1ccc(cc1)COc1ccccc1", 1),
    (48, "diphenyl ether", "c1ccc(cc1)Oc1ccccc1", 1),
    (49, "dibenzyl ether", "c1ccc(cc1)COCc1ccccc1", 1),
    (50, "benzoxazole", "c1ccc2ocnc2c1", 1),
    (51, "naphthalene", "c1ccc2ccccc2c1", 2),
    (52, "1-methylnaphthalene", "Cc1cccc2ccccc12", 1),
    (53, "2,7-dimethylnaphthalene", "Cc1ccc2ccc(C)cc2c1", 2),
    (54, "2,3-dimethylnaphthalene", "Cc1cc2ccccc2cc1C", 1),
    (55, "1,4-dimethylnaphthalene", "Cc1ccc(C)c2ccccc12", 1),
    (56, "1,2-dimethylnaphthalene", "Cc1ccc2ccccc2c1C", 1),
    (57, "1,3-dimethylnaphthalene", "Cc1cc(C)c2ccccc2c1", 1),
    (58, "pyrene", "c1cc2ccc3cccc4ccc(c1)c2c34", 2),
    (59, "anthracene", "c1ccc2cc3ccccc3cc2c1", 1),
    (60, "9-methylanthracene", "Cc1c2ccccc2cc2ccccc12", 1),
    (61, "9,10-dimethylanthracene", "Cc1c2ccccc2c(C)c2ccccc12", 1),
    (62, "phenanthrene", "c1ccc2c(c1)ccc1ccccc12", 1),
    (62, "pyrene", "c1cc2ccc3cccc4ccc(c1)c2c34", 2),
    (63, "9,10-dimethylphenanthrene", "Cc1c(C)c2ccccc2c2ccccc12", 1),
    (64, "quinoline", "c1ccc2ncccc2c1", 1),
    (65, "isoquinoline", "c1ccc2cnccc2c1", 1),
    (66, "acridine", "c1ccc2nc3ccccc3cc2c1", 1),
]

# 原典の出現回数と構造からの結果が一致しない・未確定のもの (D27)。
# (グループ番号, 化合物名, SMILES, 原典の出現回数, 状態)。原典の値を保持し、
# 実装の結果は報告するだけで失敗とは扱わない。
THIRD_ORDER_UNRESOLVED: List[Tuple[int, str, str, int, str]] = [
    (58, "phenalene", "C1C=Cc2cccc3cccc1c23", 3, "原典記載・構造解釈未確定"),
]

# 完全重複の除外と、実装上の仮定の確認 (化合物名, SMILES, 期待する全グループ)
THIRD_ORDER_FULL_CASES: List[Tuple[str, str, Dict[str, int]]] = [
    ("HO(CH2)3O(CH2)3OH (No.8 is contained in No.5)", "OCCCOCCCO",
     {"OH-(CHp)k-O-(CHn)m-OH": 1}),
    ("triethanolamine (3 pairs of arms)", "OCCN(CCO)CCO",
     {"OH-(CHp)k-NHx-(CHn)m-OH": 3}),
    ("glycerol (only the 1,3 pair has m>2)", "OCC(O)CO",
     {"OH-(CHn)m-OH": 1}),
    ("1-butoxy-2-propanol (m=2, not No.8: assumption D26)", "CCCCOCC(C)O", {}),
    ("methoxypropanol (CH3 side is excluded from No.8)", "COCCCO", {}),
    ("succinic acid (m=2 is second-order only)", "OC(=O)CCC(O)=O", {}),
    ("toluene (no third-order group)", "Cc1ccccc1", {}),
    ("ethanol (no third-order group)", "CCO", {}),
]

@dataclass(frozen=True)
class ExplicitCoefficient:
    """表から導けず、論文の例題が直接与える係数 (表の値とは区別して保持する)。

    Attributes:
        paper_label: 論文の例題での表記。
        value: 論文の例題が記載する係数。
        status: 検証状況。
    """

    paper_label: str
    value: float
    status: str


PAPER_EX5_AROMFUSED2_S3 = ExplicitCoefficient(
    "AROMFUSED[2]s3", 2.2589,
    "原典記載・未検証 (Table 8 の AROMFUSED[2]s2 は 0.3378 で一致しない, D22)",
)


@dataclass(frozen=True)
class ThirdOrderAppendixBCase:
    """論文 Appendix B の第3近似の再現テスト。

    Attributes:
        label: 例題名。
        smiles: 化合物のSMILES。
        t0: 推算式 T = T0 * ln(sum) の定数 (Table 2)。
        base_sum: 論文記載の第1次と第2次の合計 sum(N_i C_i) + sum(M_j D_j)。
        groups: 論文記載の第3次グループ出現回数 (本実装のグループ名)。
        contributions: 論文記載の寄与値 (Table 8)。
        paper_estimate: 論文記載の第3次近似の推算値 [K]。
        explicit: 表の値の代わりに明示係数を使うグループ (本実装のグループ名)。
    """

    label: str
    smiles: str
    t0: float
    base_sum: float
    groups: Dict[str, int]
    contributions: Dict[str, float]
    paper_estimate: float
    explicit: Dict[str, ExplicitCoefficient] = field(default_factory=dict)


THIRD_ORDER_APPENDIX_B: List[ThirdOrderAppendixBCase] = [
    ThirdOrderAppendixBCase(
        "Ex.1 N-phenyl-1,4-benzenediamine (Tb)",
        "Nc1ccc(Nc2ccccc2)cc1", TB0, 15.8281 + 0.1007,
        {"aC-NH-aC (different rings)": 1},
        {"aC-NH-aC (different rings)": 0.5768}, 623.94),
    ThirdOrderAppendixBCase(
        "Ex.2 pyrene (Tb)", "c1cc2ccc3cccc4ccc(c1)c2c34", TB0, 18.7593,
        {"AROMFUSED[3]": 2, "AROMFUSED[4p]": 2},
        {"AROMFUSED[3]": 0.0402, "AROMFUSED[4p]": 0.9126}, 673.96),
    ThirdOrderAppendixBCase(
        "Ex.3 4-aminobutanol (Tb)", "NCCCCO", TB0, 7.508,
        {"NH2-(CHn)m-OH": 1}, {"NH2-(CHn)m-OH": 1.0750}, 478.42),
    ThirdOrderAppendixBCase(
        "Ex.4 dicoumarol (Tm)",
        "OC1=C(Cc2c(O)c3ccccc3oc2=O)C(=O)Oc2ccccc12", TM0, 25.1421,
        {"AROMFUSED[2]": 2, "aC-(CHn=CHm)cyc (fused rings)": 2,
         "aC-Ocyc (fused rings)": 2,
         "(CHm=C)cyc-CHp-(C=CHn)cyc (different rings)": 1},
        {"AROMFUSED[2]": 0.2825, "aC-(CHn=CHm)cyc (fused rings)": 0.2479,
         "aC-Ocyc (fused rings)": -0.3545,
         "(CHm=C)cyc-CHp-(C=CHn)cyc (different rings)": 16.8558}, 552.34),
    # Ex.5: 論文は AROMFUSED[2]s3 (2.2589) と記載するが Table 8 に s3 はなく、
    # 本実装は s2 に分類する。2.2589 は明示係数として別に保持する (D22)。
    ThirdOrderAppendixBCase(
        "Ex.5 flutemazepan (Tm)",
        "CN1C(=O)C(O)N=C(c2ccccc2F)c2cc(Cl)ccc21", TM0, 24.0959 + 0.6920,
        {"AROMFUSED[2]s2": 1, "aC-NHn,cyc (fused rings)": 1,
         "aC-(C=N)cyc (different rings)": 1,
         "aC-(CHn=N)cyc (fused rings)": 1},
        {"aC-NHn,cyc (fused rings)": 3.4983,
         "aC-(C=N)cyc (different rings)": -1.3060,
         "aC-(CHn=N)cyc (fused rings)": -10.1007}, 435.22,
        explicit={"AROMFUSED[2]s2": PAPER_EX5_AROMFUSED2_S3}),
    ThirdOrderAppendixBCase(
        "Ex.6 1,9-nonanediol (Tm)", "OCCCCCCCCCO", TM0, 7.8411,
        {"OH-(CHn)m-OH": 1}, {"OH-(CHn)m-OH": 0.6674}, 315.70),
]


def check_third_order_appendix_b(case: ThirdOrderAppendixBCase) -> List[str]:
    """Appendix B の第3近似を検証し、不一致メッセージのリストを返す。"""
    res = analyze_third_order(case.smiles)
    if res.groups != case.groups or res.warnings:
        return [f"{case.label}: groups {res.groups} != {case.groups} "
                f"(eliminated={res.eliminated}, warnings={res.warnings})"]
    total = case.base_sum
    for name, n in res.groups.items():
        coeff = (case.explicit[name].value if name in case.explicit
                 else case.contributions[name])
        total += n * coeff
    estimate = case.t0 * math.log(total)
    if abs(estimate - case.paper_estimate) > ESTIMATE_TOL:
        return [f"{case.label}: estimate {estimate:.2f} != {case.paper_estimate}"]
    return []


def report_third_order_unresolved() -> List[str]:
    """原典の値と実装の結果が異なる未確定ケースを報告用の文字列で返す。"""
    names = {g.no: g.name for g in THIRD_ORDER_DEFS}
    lines = []
    for no, label, smi, paper, status in THIRD_ORDER_UNRESOLVED:
        got = analyze_third_order(smi).groups.get(names[no], 0)
        verdict = "一致" if got == paper else "不一致"
        lines.append(f"No.{no} {names[no]} / {label}: 原典 {paper}, 実装 {got} "
                     f"({verdict}; {status})")
    return lines


def collect_third_order_failures() -> List[str]:
    """第3次グループの全テストの失敗名を返す。"""
    names = {g.no: g.name for g in THIRD_ORDER_DEFS}
    failures: List[str] = []
    for no, label, smi, expected in THIRD_ORDER_CASES:
        got = analyze_third_order(smi).groups.get(names[no], 0)
        if got != expected:
            failures.append(f"No.{no} {names[no]} / {label}: {got} != {expected}")
    for label, smi, expected_groups in THIRD_ORDER_FULL_CASES:
        res = analyze_third_order(smi)
        if res.groups != expected_groups:
            failures.append(f"full: {label}: {res.groups} != {expected_groups}")
    for case in THIRD_ORDER_APPENDIX_B:
        failures.extend(check_third_order_appendix_b(case))
    # No.17 の例外を切ると、表の範囲 (n in 0..1) どおりに CH2 は数えない
    key = "No17_CH2_fused_carbocycle"
    saved = THIRD_ORDER_EXCEPTIONS[key]
    THIRD_ORDER_EXCEPTIONS[key] = False
    try:
        got = analyze_third_order("C1CCc2ccccc2C1").groups.get(names[17], 0)
    finally:
        THIRD_ORDER_EXCEPTIONS[key] = saved
    if got != 0:
        failures.append(f"No.17 exception off: tetralin {got} != 0")
    return failures


TABLE8_NUMBERS = range(1, 67)  # 論文 Table 8 の第3次グループ番号


def find_unimplemented_third_order() -> List[int]:
    """Table 8 (No.1-66) のうち THIRD_ORDER_DEFS に未定義の番号を返す。"""
    implemented = {g.no for g in THIRD_ORDER_DEFS}
    return [n for n in TABLE8_NUMBERS if n not in implemented]


def find_untested_third_order() -> List[int]:
    """THIRD_ORDER_CASES に例題がない第3次グループ番号を返す。"""
    tested = {no for no, _, _, _ in THIRD_ORDER_CASES}
    return [g.no for g in THIRD_ORDER_DEFS if g.no not in tested]


def run_third_order_shuffle_test(n_trials: int = 30, seed: int = 0) -> List[str]:
    """THIRD_ORDER_DEFS の並びを入れ替えても結果が変わらないことを検証する。"""
    rng = random.Random(seed)
    original = list(THIRD_ORDER_DEFS)
    failures = set()
    try:
        for _ in range(n_trials):
            shuffled = original[:]
            rng.shuffle(shuffled)
            _THIRD_COMPILED[:] = _compile_second(shuffled)
            failures.update(collect_third_order_failures())
    finally:
        _THIRD_COMPILED[:] = _compile_second(original)
    return sorted(failures)


TABLE6_NUMBERS = range(1, 183)  # 論文 Table 6 の第1次グループ番号


def find_unimplemented() -> List[int]:
    """Table 6 (No.1-182) のうち GROUP_DEFS に未定義の番号を返す。"""
    implemented = {g.no for g in GROUP_DEFS}
    return [n for n in TABLE6_NUMBERS if n not in implemented]


def find_untested() -> List[Tuple[int, str]]:
    """どのテストケースにも現れないグループ (番号, 名前) を返す。"""
    used = set()
    for _, smi, _ in TEST_CASES:
        used.update(analyze_smiles(smi).groups)
    for case in APPENDIX_B_CASES:
        used.update(analyze_smiles(case.smiles).groups)
    return [(g.no, g.name) for g in GROUP_DEFS if g.name not in used]


def find_ties() -> Dict[str, List[Tuple[str, str]]]:
    """テスト分子のうち、同優先度・同重量の競合が起きたものを返す。"""
    found: Dict[str, List[Tuple[str, str]]] = {}
    smiles_list = [(n, smi) for n, smi, _ in TEST_CASES]
    smiles_list += [(c.label, c.smiles) for c in APPENDIX_B_CASES]
    for name, smi in smiles_list:
        ties = analyze_smiles(smi).ties
        if ties:
            found[name] = ties
    return found


def _print_analysis(smiles: str) -> None:
    res = analyze_smiles(smiles)
    print(f"{smiles}: {res.groups}")
    if res.unassigned:
        print(f"  WARNING unassigned atoms: {res.unassigned}")
    for warning in res.warnings:
        print(f"  WARNING {warning}")
    if res.ties:
        print(f"  NOTE tie-break applied (winner, loser): {res.ties}")
    second = analyze_second_order(smiles)
    print(f"  second-order: {second.groups}")
    if second.eliminated:
        print(f"  NOTE complete overlap removed (removed, kept): "
              f"{second.eliminated}")
    if second.ties:
        print(f"  NOTE identical extent (kept both): {second.ties}")
    third = analyze_third_order(smiles)
    print(f"  third-order: {third.groups}")
    for note in third.notes:
        print(f"  NOTE {note}")
    if third.eliminated:
        print(f"  NOTE complete overlap removed (removed, kept): "
              f"{third.eliminated}")


def main() -> int:
    """引数のSMILESを解析する。引数がなければ内蔵テスト一式を実行する。"""
    args = [a for a in sys.argv[1:] if not a.startswith("-")]
    verbose = "-v" in sys.argv[1:]
    if args:
        for smi in args:
            _print_analysis(smi)
        return 0

    n_fail = 0
    for name, smi, expected in TEST_CASES:
        res = analyze_smiles(smi)
        ok = res.groups == expected and res.is_complete
        n_fail += not ok
        if verbose or not ok:
            print(f"[{'PASS' if ok else 'FAIL'}] {name}: {smi}")
            print(f"       got      : {res.groups}")
        if not ok:
            print(f"       expected : {expected}")
            print(f"       unassigned: {res.unassigned}, "
                  f"warnings: {res.warnings}")
    print(f"unit cases       : {len(TEST_CASES) - n_fail}/{len(TEST_CASES)}")

    n_limit_fail = 0
    for name, smi, expected, unassigned in LIMIT_CASES:
        res = analyze_smiles(smi)
        symbols = sorted(sym for _, sym in res.unassigned)
        ok = res.groups == expected and symbols == unassigned
        n_limit_fail += not ok
        if not ok:
            print(f"[FAIL] limit case {name}: {res.groups} {symbols}")
    print(f"known-limit cases: {len(LIMIT_CASES) - n_limit_fail}/"
          f"{len(LIMIT_CASES)}")

    n_b_fail = 0
    for case in APPENDIX_B_CASES:
        problems = check_appendix_b(case)
        n_b_fail += bool(problems)
        for problem in problems:
            print(f"[FAIL] {problem}")
        if verbose and not problems:
            res = analyze_smiles(case.smiles)
            total = sum(n * case.contributions[g]
                        for g, n in res.groups.items())
            print(f"[PASS] {case.label}: sum={total:.4f} "
                  f"estimate={case.t0 * math.log(total):.2f} K")
    print(f"Appendix B       : {len(APPENDIX_B_CASES) - n_b_fail}/"
          f"{len(APPENDIX_B_CASES)}")

    shuffle_failures = run_shuffle_test()
    for failure in shuffle_failures:
        print(f"[FAIL] order-dependent: {failure}")
    print(f"shuffle test     : {'OK' if not shuffle_failures else 'FAILED'} "
          "(30 random orderings of GROUP_DEFS)")

    second_failures = collect_second_order_failures()
    for failure in second_failures:
        print(f"[FAIL] second-order: {failure}")
    n_second = (len(SECOND_ORDER_CASES) + len(SECOND_ORDER_FULL_CASES)
                + len(NO_SECOND_ORDER_CASES) + len(SECOND_ORDER_APPENDIX_B))
    print(f"second-order     : {n_second - len(second_failures)}/{n_second} "
          "(Table 7 examples, overlap rule, Appendix B)")
    second_shuffle = run_second_order_shuffle_test()
    for failure in second_shuffle:
        print(f"[FAIL] second-order order-dependent: {failure}")
    print(f"second shuffle   : {'OK' if not second_shuffle else 'FAILED'} "
          "(30 random orderings of SECOND_ORDER_DEFS)")
    second_missing = find_unimplemented_second_order()
    print(f"implemented second-order groups: {len(SECOND_ORDER_DEFS)}/"
          f"{len(TABLE7_NUMBERS)}")
    print(f"unimplemented second-order numbers: {second_missing or 'none'}")
    print("second-order groups without a test case: "
          f"{find_untested_second_order() or 'none'}")
    print("second-order identical-extent ties: "
          f"{find_second_order_ties() or 'none'}")

    third_failures = collect_third_order_failures()
    for failure in third_failures:
        print(f"[FAIL] third-order: {failure}")
    n_third = (len(THIRD_ORDER_CASES) + len(THIRD_ORDER_FULL_CASES)
               + len(THIRD_ORDER_APPENDIX_B) + 1)  # +1: No.17 例外の切替
    print(f"third-order      : {n_third - len(third_failures)}/{n_third} "
          "(Table 8 examples, overlap rule, Appendix B)")
    third_shuffle = run_third_order_shuffle_test()
    for failure in third_shuffle:
        print(f"[FAIL] third-order order-dependent: {failure}")
    print(f"third shuffle    : {'OK' if not third_shuffle else 'FAILED'} "
          "(30 random orderings of THIRD_ORDER_DEFS)")
    third_missing = find_unimplemented_third_order()
    print(f"implemented third-order groups: {len(THIRD_ORDER_DEFS)}/"
          f"{len(TABLE8_NUMBERS)}")
    print(f"unimplemented third-order numbers: {third_missing or 'none'}")
    print("third-order groups without a test case: "
          f"{find_untested_third_order() or 'none'}")
    for line in report_third_order_unresolved():
        print(f"[UNRESOLVED] {line}")

    missing = find_unimplemented()
    untested = find_untested()
    ties = find_ties()
    print(f"implemented groups: {len(GROUP_DEFS)}/{len(TABLE6_NUMBERS)}")
    print(f"unimplemented numbers: {missing or 'none'}")
    print(f"groups without a test case: {untested or 'none'}")
    print(f"molecules with tie-break conflicts: {ties or 'none'}")
    failed = (n_fail or n_limit_fail or n_b_fail or shuffle_failures
              or missing or second_failures or second_shuffle
              or second_missing or third_failures or third_shuffle
              or third_missing)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
