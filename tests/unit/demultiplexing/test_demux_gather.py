from collections import Counter

from workflow.scripts.demultiplexing.demux_gather import combine_pickle


def test_combine_pickle_keeps_max_r2_read_length():
    combined = {
        "experiment_name": "exp",
        "version": "0.0.0",
        "max_r2_read_length": 84,
    }
    incoming = {
        "experiment_name": "exp",
        "version": "0.0.0",
        "max_r2_read_length": 101,
    }

    merged = combine_pickle(incoming, combined)
    assert merged["max_r2_read_length"] == 101

    incoming_smaller = {
        "experiment_name": "exp",
        "version": "0.0.0",
        "max_r2_read_length": 90,
    }
    merged = combine_pickle(incoming_smaller, merged)
    assert merged["max_r2_read_length"] == 101


def test_combine_pickle_sets_max_r2_read_length_if_missing():
    combined = {
        "experiment_name": "exp",
        "version": "0.0.0",
    }
    incoming = {
        "experiment_name": "exp",
        "version": "0.0.0",
        "max_r2_read_length": 92,
    }

    merged = combine_pickle(incoming, combined)
    assert merged["max_r2_read_length"] == 92


def test_combine_pickle_sums_hashing_reads_per_umi():
    def chunk(umi, umi_rescued):
        return {
            "experiment_name": "exp",
            "version": "0.0.0",
            "hashing": {
                "sample1": {
                    "hash1": {
                        "n_correct": sum(umi.values()) - sum(umi_rescued.values()),
                        "n_corrected": sum(umi_rescued.values()),
                        "n_correct_upstream": 0,
                        "counts": {
                            "CELL1": {
                                "umi": Counter(umi),
                                "umi_rescued": Counter(umi_rescued),
                                "count": sum(umi.values()),
                                "cell_barcode_label": "A01_A01_LIG1_P01-A01",
                            }
                        },
                    }
                }
            },
        }

    combined = chunk({"AAAAAAAA": 5, "CCCCCCCC": 1}, {"CCCCCCCC": 1})
    merged = combine_pickle(chunk({"AAAAAAAA": 3, "GGGGGGGG": 2}, {}), combined)

    cell = merged["hashing"]["sample1"]["hash1"]["counts"]["CELL1"]
    assert cell["umi"] == Counter({"AAAAAAAA": 8, "GGGGGGGG": 2, "CCCCCCCC": 1})
    assert cell["umi_rescued"] == Counter({"CCCCCCCC": 1})
    assert cell["count"] == 11
    assert len(cell["umi"]) == 3
