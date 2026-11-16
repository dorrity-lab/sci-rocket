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
    hashing_metrics = Path(__file__).parent.parent / "output" / "fastq_from_aviti_one_run" / "hashing" / "fastq_from_aviti_one_run_hashing_metrics.tsv"
    assert hashing_metrics.exists()
    data = read_tsv(hashing_metrics)
    assert len(data) == 176
    # check top five
    assert data[0] == ['fastq_from_aviti_one_run', 'ZAe-14hpf-28', '10uM_P7_D7', 'TTGATTGGCG_GGCCGTCAAC_AATTCTAGGT_TTAATTGAAT', 'D02_A03_LIG114_P02-B05', '22', '22']
    assert data[1] == ['fastq_from_aviti_one_run', 'ZAe-14hpf-28', '10uM_P7_A8', 'CTAGCCAGCC_CTTAATCTTG_CGTCTTCCTG_TATCATGATC', 'D05_G03_LIG185_P02-D04', '17', '17']
    assert data[2] == ['fastq_from_aviti_one_run', 'ZAe-14hpf-28', '10uM_P7_A7', 'TTGGCAAGCC_GGTCTCGCCG_ATGGTTGGTG_ACTATAGGTT', 'D01_D03_LIG72_P02-E05', '15', '15']
    assert data[3] == ['fastq_from_aviti_one_run', 'ZAe-14hpf-28', '10uM_P7_B7', 'GCATATGAGC_CTCCTGGACC_TGATGCGATG_GTTACGCAAG', 'D03_F03_LIG93_P02-F05', '15', '15']
    assert data[4] == ['fastq_from_aviti_one_run', 'ZAe-14hpf-28', '10uM_P7_C7', 'TGCGACCTCT_GGTCTCGCCG_TGGATTCTAT_GAGCATATGG', 'D06_D03_LIG12_P02-D05', '14', '14']


def test_hashing_umi_counts_match_hashing_metrics():
    hashing_dir = Path(__file__).parent.parent / "output" / "fastq_from_aviti_one_run" / "hashing"
    hashing_umi_counts = hashing_dir / "fastq_from_aviti_one_run_hashing_umi_counts.tsv.gz"
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
    metrics = read_tsv(hashing_dir / "fastq_from_aviti_one_run_hashing_metrics.tsv")
    assert {(row[1], row[2], row[3]): int(row[5]) for row in metrics} == dict(reads)
    assert {(row[1], row[2], row[3]): int(row[6]) for row in metrics} == dict(n_umi)
