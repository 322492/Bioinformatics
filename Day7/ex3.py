from __future__ import annotations

import argparse
import io
import re
import subprocess
import time
from dataclasses import dataclass
from typing import Iterable
from urllib.parse import urlencode

from Bio.Blast import NCBIXML


HBA_FASTA = """>sp|P69905|HBA_HUMAN Hemoglobin subunit alpha OS=Homo sapiens OX=9606 GN=HBA1 PE=1 SV=2
MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHG
KKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTP
AVHASLDKFLASVSTVLTSKYR
"""

HBB_FASTA = """>NM_000518.5 Homo sapiens hemoglobin subunit beta (HBB), mRNA
ACATTTGCTTCTGACACAACTGTGTTCACTAGCAACCTCAAACAGACACCATGGTGCATCTGACTCCTGA
GGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGC
AGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGATG
CTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGC
TCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACTGTGACAAGCTGCACGTGGAT
CCTGAGAACTTCAGGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGCAAAGAATTCA
CCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAATGCCCTGGCCCACAAGTATCA
CTAAGCTCGCTTTCTTGCTGTCCAATTTCTATTAAAGGTTCCTTTGTTCCCTAAGTCCAACTACTAAACT
GGGGGATATTATGAAGGGCCTTGAGCATCTGGATTCTGCCTAATAAAAAACATTTATTTTCATTGCAA
"""

SPECIES_RE = re.compile(r"\[([^\[\]]+)\]\s*$")
RID_RE = re.compile(r"RID = ([A-Z0-9-]+)")
RTOE_RE = re.compile(r"RTOE = (\d+)")
BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"

@dataclass
class QueryConfig:
    label: str
    program: str
    fasta_text: str


@dataclass
class HitRow:
    target_sequence: str
    species: str
    length: int
    e_value: str
    percent_identity: str


class BlastRequestError(RuntimeError):
    pass


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run BLAST searches for HBA and HBB, then print the top non-human "
            "hits from 10 different species."
        )
    )
    parser.add_argument(
        "--max-species",
        type=int,
        default=10,
        help="Number of distinct non-human species to report per query.",
    )
    return parser.parse_args()


def fasta_sequence(fasta_text: str) -> str:
    lines = [line.strip() for line in fasta_text.strip().splitlines()]
    return "".join(line for line in lines if line and not line.startswith(">"))


def extract_species(hit_description: str) -> str | None:
    match = SPECIES_RE.search(hit_description)
    if match is None:
        return None
    return match.group(1)


def format_e_value(value: float) -> str:
    return f"{value:.2e}".replace(".00e", "e").replace("e-0", "e-").replace("e+0", "e+")


def run_curl(curl_args: list[str]) -> str:
    completed = subprocess.run(
        [
            "curl",
            "--silent",
            "--show-error",
            "--fail",
            "--max-time",
            "60",
            "--retry",
            "3",
            "--retry-delay",
            "2",
            *curl_args,
        ],
        check=True,
        capture_output=True,
        text=True,
        timeout=75,
    )
    return completed.stdout


def submit_blast(program: str, sequence: str) -> tuple[str, int]:
    response = run_curl(
        [
        "--data",
        "CMD=Put",
        "--data",
        f"PROGRAM={program}",
        "--data",
        "DATABASE=nr",
        "--data-urlencode",
        f"QUERY={sequence}",
        "--data",
        "FORMAT_TYPE=XML",
        "--data",
        "FILTER=L",
        "--data",
        "TOOL=bioinformatics_day7_ex3",
        BLAST_URL,
        ]
    )

    rid_match = RID_RE.search(response)
    rtoe_match = RTOE_RE.search(response)
    if rid_match is None or rtoe_match is None:
        raise BlastRequestError("Could not parse RID/RTOE from BLAST submission response.")

    return rid_match.group(1), int(rtoe_match.group(1))


def poll_blast_xml(rid: str, delay_seconds: int) -> str:
    time.sleep(max(delay_seconds, 3))

    params = {
        "CMD": "Get",
        "RID": rid,
        "FORMAT_TYPE": "XML",
        "TOOL": "bioinformatics_day7_ex3",
    }

    for _ in range(30):
        url = f"{BLAST_URL}?{urlencode(params)}"
        response = run_curl([url])

        if "Status=WAITING" in response:
            time.sleep(15)
            continue
        if "Status=FAILED" in response:
            raise BlastRequestError(f"BLAST search {rid} failed on the NCBI server.")
        if "Status=UNKNOWN" in response:
            raise BlastRequestError(f"BLAST search {rid} expired or was not found.")
        if "Status=READY" in response and "ThereAreHits=yes" not in response:
            return ""
        if response.lstrip().startswith("<?xml"):
            return response

        time.sleep(5)

    raise BlastRequestError(f"Timed out while waiting for BLAST search {rid}.")


def fetch_hits(query: QueryConfig, max_species: int) -> list[HitRow]:
    sequence = fasta_sequence(query.fasta_text)
    print(f"Submitting {query.label} using {query.program}...", flush=True)
    rid, rtoe = submit_blast(query.program, sequence)
    print(f"RID {rid} accepted; estimated wait {rtoe}s.", flush=True)
    blast_xml = poll_blast_xml(rid, rtoe)
    if not blast_xml:
        print(f"No hits returned for {query.label}.", flush=True)
        return []

    blast_record = NCBIXML.read(io.StringIO(blast_xml))
    rows: list[HitRow] = []
    seen_species: set[str] = set()

    for alignment in blast_record.alignments:
        species = extract_species(alignment.hit_def)
        if species is None:
            continue
        if species in {"Homo sapiens", "synthetic construct"}:
            continue
        if species in seen_species:
            continue
        if not alignment.hsps:
            continue

        hsp = alignment.hsps[0]
        percent_identity = (hsp.identities / hsp.align_length) * 100
        rows.append(
            HitRow(
                target_sequence=query.label,
                species=species,
                length=alignment.length,
                e_value=format_e_value(hsp.expect),
                percent_identity=f"{percent_identity:.2f}",
            )
        )
        seen_species.add(species)

        if len(rows) >= max_species:
            break

    print(f"Collected {len(rows)} species for {query.label}.", flush=True)
    return rows


def print_table(rows: Iterable[HitRow]) -> None:
    rows = list(rows)
    headers = ["target_sequence", "species", "length", "E-value", "%identity"]
    data = [
        [
            row.target_sequence,
            row.species,
            str(row.length),
            row.e_value,
            row.percent_identity,
        ]
        for row in rows
    ]

    widths = [len(header) for header in headers]
    for row in data:
        for index, value in enumerate(row):
            widths[index] = max(widths[index], len(value))

    def render(row: list[str]) -> str:
        return " | ".join(value.ljust(widths[index]) for index, value in enumerate(row))

    separator = "-+-".join("-" * width for width in widths)

    print(render(headers))
    print(separator)
    for row in data:
        print(render(row))


def main() -> None:
    args = parse_args()
    queries = [
        QueryConfig("HBA_HUMAN_P69905", "blastp", HBA_FASTA),
        QueryConfig("HBB_HUMAN_NM_000518.5", "blastx", HBB_FASTA),
    ]

    all_rows: list[HitRow] = []
    for query in queries:
        all_rows.extend(fetch_hits(query, args.max_species))

    print_table(all_rows)


if __name__ == "__main__":
    main()
