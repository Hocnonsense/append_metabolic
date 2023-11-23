# -*- coding: utf-8 -*-
"""
 * @Date: 2023-11-23 17:13:39
 * @LastEditors: hwrn.aou hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-11-23 22:29:47
 * @FilePath: /METABOLIC/append_metabolic/hmm_template.py
 * @Description:
"""
# """

from functools import cached_property
from pathlib import Path
from typing import Generator, Literal, NamedTuple
import pandas as pd

METABOLIC_DIR = Path(__file__).parent.parent


class HmmsearchTbloutLine(NamedTuple):
    """
    https://www.jianshu.com/p/de67da5b84d1

    header: list[str] = [
        "target name",
        "accession",
        "query name",
        "accession",
        *["E-value", "score", "bias"], # full sequence
        *["E-value", "score", "bias"], # best 1 domain
        *["exp", "reg", "clu", "ov", "env", "dom", "rep", "inc"], # domain number estimation
        "description of target",
    ]
    """

    class SeqId(NamedTuple):
        name: str
        accession: str | Literal["-"]

    class Hits(NamedTuple):
        e_value: float
        score: float
        bias: float

    class DomainNumberEstimation(NamedTuple):
        exp: float
        reg: int
        clu: int
        ov: int
        env: int
        dom: int
        rep: int
        inc: int

    target: SeqId
    query: SeqId
    full_sequence: Hits
    best_domain: Hits
    domain_num_est: DomainNumberEstimation
    description: str

    @classmethod
    def from_line(cls, lines: Generator[str, None, None]):
        """
        #                                                                --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----
        # target name         accession  query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
        # ------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------
        M77_2|k141_4450897_13 -          K10944               -            5.9e-59  196.9  21.6   6.7e-59  196.7  21.6   1.0   1   0   0   1   1   1   1 # 7991 # 8642 # 1 # source=Prodigal_v2.6.3;score=47.69;seqid=M77_2|k141_4450897;frame=0;ID=8_13;partial=00;start_type=ATG;rbs_motif=None;rbs_spacer=None;gc_cont=0.429;conf=100.00;cscore=43.55;sscore=4.15;rscore=-0.59;uscore=2.34;tscore=2.40
        M77_2|k141_10148615_7 -          K10944               -            1.2e-58  196.0  20.8   1.3e-58  195.8  20.8   1.0   1   0   0   1   1   1   1 # 4508 # 5159 # -1 # source=Prodigal_v2.6.3;score=46.05;seqid=M77_2|k141_10148615;frame=0;ID=30_7;partial=00;start_type=ATG;rbs_motif=None;rbs_spacer=None;gc_cont=0.418;conf=100.00;cscore=42.09;sscore=3.96;rscore=-0.59;uscore=2.15;tscore=2.40
        #
        # Program:         hmmsearch
        # Version:         3.4 (Aug 2023)
        # Pipeline mode:   SEARCH
        # Query file:      /dssg/home/acct-trench/trench-0/clsxx/Users/hwrn.aou/Templates/METABOLIC/METABOLIC_hmm_db/amoA.hmm
        # Target file:     test/out_compare/total.faa
        # Option settings: hmmsearch --tblout test/out_compare/intermediate_files/Hmmsearch_Outputs/amoA.hmm.total.hmmsearch_result.txt -T 160.90 --cpu 1 /dssg/home/acct-trench/trench-0/clsxx/Users/hwrn.aou/Templates/METABOLIC/METABOLIC_hmm_db/amoA.hmm test/out_compare/total.faa
        # Current dir:     /dssg/home/acct-trench/trench-0/clsxx/Users/hwrn.aou/Templates/METABOLIC
        # Date:            Tue Nov 21 17:34:56 2023
        # [ok]
        """
        for line in (i.strip() for i in lines):
            if line.startswith("#") or not line:
                continue
            values = line.split(maxsplit=18)
            yield cls(
                cls.SeqId(values[0], values[1]),
                cls.SeqId(values[2], values[3]),
                cls.Hits(float(values[4]), float(values[5]), float(values[6])),
                cls.Hits(float(values[7]), float(values[8]), float(values[9])),
                cls.DomainNumberEstimation(
                    float(values[10]),
                    int(values[11]),
                    int(values[12]),
                    int(values[13]),
                    int(values[14]),
                    int(values[15]),
                    int(values[16]),
                    int(values[17]),
                ),
                values[18],
            )


class ThresholdType(NamedTuple):
    threshold: float = 50
    score_type: Literal["full", "domain"] = "full"

    def to_hmmsearch(self):
        if self.score_type == "full":
            return f"-T {self.threshold}"
        return f"--domT {self.threshold}"

    def check(self, hmmsearch_tblout_line: HmmsearchTbloutLine):
        if self.score_type == "full":
            return hmmsearch_tblout_line.full_sequence.score > self.threshold
        else:
            return hmmsearch_tblout_line.best_domain.score > self.threshold

    def __str__(self) -> str:
        return f"{self.threshold}|{self.score_type}"

    @classmethod
    def from_str(cls, threshold_type="50|full"):
        threshold = threshold_type.split("|", 1)[0]
        score_type = threshold_type.split("|", 1)[1]
        return cls(float(threshold), score_type)

    @classmethod
    def from_df_ko(cls, KO: str, df: pd.DataFrame):
        i = df.index[(df["KO"] == KO)][0]
        return cls(float(df.loc[i, "threshold"]), df.loc[i, "score_type"])


class KO2ThresholdTypeDB:
    def __init__(self):
        self._ko2thresholdtype = None

    @staticmethod
    def extract_df(hmm_df, replace="-", default="50|full") -> dict[str, ThresholdType]:
        """
        >>> df
                   KO threshold score_type
        0      K00001    362.77     domain
        1      K00002    443.20       full
        2      K00003    287.10     domain
        3      K00004    365.63     domain
        4      K00005    329.70       full
        ...       ...       ...        ...
        26243  K26928    494.67       full
        26244  K26929    580.07       full
        26245  K26930        50          -
        26246  K26931    413.70       full
        26247  K26932    382.83       full
        """
        threshold_type = ThresholdType.from_str(default)
        df = hmm_df.assign(
            threshold=lambda df: df["threshold"].apply(
                lambda x: threshold_type.threshold if x == replace else x
            ),
            score_type=lambda df: df["score_type"].apply(
                lambda x: threshold_type.score_type if x == replace else x
            ),
        )
        return {
            ko: ThresholdType(t, st)
            for ko, (t, st) in (
                df[["KO", "threshold", "score_type"]]
                .drop_duplicates()
                .set_index("KO")
                .T.to_dict()
            ).items()
        }

    def _extract_df_property(self) -> dict[str, ThresholdType]:
        raise NotImplementedError

    @property
    def ko2thresholdtype(self):
        return self._ko2thresholdtype or self._extract_df_property()


class HmmTemplate(KO2ThresholdTypeDB):
    FILENAME = (
        METABOLIC_DIR / "METABOLIC_template_and_database" / "hmm_table_template.txt"
    )

    def __init__(self) -> None:
        super().__init__()
        self.df = pd.read_csv(self.FILENAME, sep="\t", na_values="N/A")
        self.hmms = (
            self.df[
                ["Hmm file", "Corresponding KO", "Hmm detecting threshold", "#Entry"]
            ]
            .dropna(axis=0)
            .apply(
                lambda s: [
                    {
                        "hmm": i,
                        "KO": j,
                        "threshold": k.split("|", 1)[0],
                        "score_type": k.split("|", 1)[1],
                        "entry": s[3],
                    }
                    for (i, j, k) in zip(
                        s[0].split(", "), s[1].split(", "), s[2].split(", ")
                    )
                ],
                axis=1,
            )
            .explode()
            .apply(pd.Series)
        )

    @cached_property
    def hmm2ko(self):
        return self.hmms.set_index("hmm")["KO"].to_dict()

    def _extract_df_property(self):
        self._ko2thresholdtype = self.extract_df(self.hmms)
        return self._ko2thresholdtype


class HmmTemplate2:
    FILENAME = (
        METABOLIC_DIR / "METABOLIC_template_and_database" / "hmm_table_template_2.txt"
    )

    def __init__(self) -> None:
        self.df = pd.read_csv(self.FILENAME, sep="\t", na_values="N/A")


class KO_LIST(KO2ThresholdTypeDB):
    FILENAME = METABOLIC_DIR / "kofam_database" / "ko_list"

    def __init__(self, ko_list="kofam_database/ko_list") -> None:
        super().__init__()
        self.df = pd.read_csv(ko_list, sep="\t").rename(columns={"knum": "KO"})

    def _extract_df_property(self):
        self._ko2thresholdtype = self.extract_df(self.df)
        return self._ko2thresholdtype


def load_ko2thresholdtype(
    ko_list: str | Path | pd.DataFrame = "kofam_database/ko_list", cache=True
) -> dict[str, ThresholdType]:
    ht = HmmTemplate()
    kofam_db = (
        ko_list if isinstance(ko_list, pd.DataFrame) else KO_LIST(ko_list=ko_list)
    )
    return {
        **KO2ThresholdTypeDB.extract_df(ht.hmms),
        **KO2ThresholdTypeDB.extract_df(kofam_db.df),
    }


def parse_hmmsearch_line(
    lines: Generator[str, None, None],
    ko2thresholdtype: dict[str, ThresholdType],
    default="50|full",
):
    """
    header: list[str] = [
        "target name",
        "accession",
        "query name",
        "accession",
        *["E-value", "score", "bias"], # full sequence
        *["E-value", "score", "bias"], # best 1 domain
        *["exp", "reg", "clu", "ov", "env", "dom", "rep", "inc"], # domain number estimation
        "description of target",
    ]
    """
    default_threshold_type = ThresholdType.from_str(default)
    for htl in HmmsearchTbloutLine.from_line(lines):
        threshold_type = ko2thresholdtype.get(htl.target.name, default_threshold_type)
        yield htl, threshold_type.check(hmmsearch_tblout_line=htl)
