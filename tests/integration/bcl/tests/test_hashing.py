from collections import Counter
from pathlib import Path
import csv
import gzip


def read_tsv(filename: Path):
    data = []
    with open(filename, 'r', newline='', encoding='utf-8') as f:
        tsv_reader = csv.reader(f, delimiter='\t')
        next(tsv_reader)
        for row in tsv_reader:
            data.append(row)
    return data


def test_hashing_metrics():
    hashing_metrics = Path(__file__).parent.parent / "output" / "bcl_one_run" / "hashing" / "bcl_one_run_hashing_metrics.tsv"
    assert hashing_metrics.exists()
    data = read_tsv(hashing_metrics)
    assert len(data) == 1675
    # check top three
    assert data[0] == ['bcl_one_run', 'zfish-hash', '40uM_P7_F12', 'CGCGGCCATA_CGGCGAGACC_TACTTGCGTG_GCGTTGGAGC', 'C06_D03_LIG87_P01-A02', '10', '6']
    assert data[1] == ['bcl_one_run', 'zfish-hash', '20uM_P7_G10', 'AATATTACTT_GAACTGCATC_CCAAGGATTG_TCCAATCATC', 'C07_B03_LIG78_P01-F03', '10', '3']
    assert data[2] == ['bcl_one_run', 'zfish-hash', '40uM_P7_B12', 'ATTGGCAGAT_GTTGACGGCC_GCTCTTAGTG_TCCGGCCTCG', 'C08_A03_LIG68_P01-D04', '8', '5']


def test_hashing_umi_counts_match_hashing_metrics():
    hashing_dir = Path(__file__).parent.parent / "output" / "bcl_one_run" / "hashing"
    hashing_umi_counts = hashing_dir / "bcl_one_run_hashing_umi_counts.tsv.gz"
    assert hashing_umi_counts.exists()
    with gzip.open(hashing_umi_counts, 'rt', newline='', encoding='utf-8') as f:
        tsv_reader = csv.reader(f, delimiter='\t')
        assert next(tsv_reader) == ['experiment_name', 'sample_name', 'hashing_name', 'cell_barcode', 'umi', 'reads', 'reads_rescued']
        umi_rows = list(tsv_reader)

    # Per cell / hashing barcode, the reads per UMI must sum to count and the no. of UMIs must equal n_umi.
    reads = Counter()
    n_umi = Counter()
    for row in umi_rows:
        assert int(row[6]) <= int(row[5])
        reads[(row[1], row[2], row[3])] += int(row[5])
        n_umi[(row[1], row[2], row[3])] += 1
    metrics = read_tsv(hashing_dir / "bcl_one_run_hashing_metrics.tsv")
    assert {(row[1], row[2], row[3]): int(row[5]) for row in metrics} == dict(reads)
    assert {(row[1], row[2], row[3]): int(row[6]) for row in metrics} == dict(n_umi)
