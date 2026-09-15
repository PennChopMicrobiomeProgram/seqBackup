import io
from pathlib import Path

import pytest

from seqBackupLib.nanopore import (
    ONT_MACHINE_TYPES,
    NanoporeDir,
    NanoporeFastq,
    extract_ont_device_code,
    natural_sort_key,
)

# Every supported ONT machine type must have a fixture exercising it, mirroring
# the mechanism in test_illumina.py.
ont_machine_fixtures = {
    "ONT-MN": "minion_dir",
    "ONT-P2I": "p2i_dir",
}


def test_all_ont_machine_types_are_tested():
    missing = set(ONT_MACHINE_TYPES.values()) - set(ont_machine_fixtures)
    assert not missing, f"Missing nanopore fixtures for machine types: {missing}"


@pytest.mark.parametrize("machine_type", ont_machine_fixtures.keys())
def test_nanopore_dir(machine_type, request):
    run_dir = request.getfixturevalue(ont_machine_fixtures[machine_type])

    nd = NanoporeDir(run_dir.name)

    assert nd.machine_type == machine_type
    assert nd.folder_info["flowcell_id"] in run_dir.name
    assert nd.build_archive_dir() == run_dir.name


def test_extract_ont_device_code():
    assert extract_ont_device_code("MN47822") == "MN"
    assert extract_ont_device_code("P2I-00513-B") == "P2I"
    assert extract_ont_device_code("GXB01234-X1") == "GXB"


def test_nanopore_dir_parses_minion():
    nd = NanoporeDir("20260401_1506_MN47822_FBE92725_8b1f29fd")
    assert nd.folder_info == {
        "date": "2026-04-01",
        "time": "1506",
        "position": "MN47822",
        "device_code": "MN",
        "flowcell_id": "FBE92725",
        "short_run_id": "8b1f29fd",
        "run_alias": "",
    }
    assert nd.machine_type == "ONT-MN"


def test_nanopore_dir_parses_p2i_with_alias():
    nd = NanoporeDir("20260408_1652_P2I-00513-B_PBK70557_2fe9fcef_Promethion_Training")
    assert nd.folder_info["position"] == "P2I-00513-B"
    assert nd.folder_info["device_code"] == "P2I"
    assert nd.folder_info["flowcell_id"] == "PBK70557"
    assert nd.folder_info["run_alias"] == "Promethion_Training"
    assert nd.machine_type == "ONT-P2I"


@pytest.mark.parametrize(
    "bad_name",
    [
        "not_a_run",
        "20260401_1506_XX12345_FBE92725_8b1f29fd",  # unknown device
        "2026041_1506_MN47822_FBE92725_8b1f29fd",  # bad date
        "20260401_156_MN47822_FBE92725_8b1f29fd",  # bad time
    ],
)
def test_nanopore_dir_rejects_bad_names(bad_name):
    with pytest.raises(ValueError):
        NanoporeDir(bad_name)


def test_natural_sort_key_orders_numerically():
    names = ["r_10.fastq.gz", "r_2.fastq.gz", "r_1.fastq.gz"]
    ordered = sorted((Path(n) for n in names), key=natural_sort_key)
    assert [p.name for p in ordered] == [
        "r_1.fastq.gz",
        "r_2.fastq.gz",
        "r_10.fastq.gz",
    ]


def _handle(header: str) -> io.StringIO:
    return io.StringIO(f"@{header}\nACGT\n+\nIIII\n")


def test_nanopore_fastq_header_parsing_and_checks():
    header = (
        "0a1b2c3d-1111 runid=abc read=5 ch=9 flow_cell_id=FBE92725 barcode=barcode01"
    )
    nf = NanoporeFastq(
        _handle(header), expected_flowcell="FBE92725", expected_barcode="barcode01"
    )
    assert nf.info["read_id"] == "0a1b2c3d-1111"
    assert nf.info["flow_cell_id"] == "FBE92725"
    ok, problems = nf.check_fp_vs_content()
    assert ok and problems == {}


def test_nanopore_fastq_detects_mismatch():
    header = "readx runid=abc flow_cell_id=WRONG barcode=barcode02"
    nf = NanoporeFastq(
        _handle(header), expected_flowcell="FBE92725", expected_barcode="barcode01"
    )
    ok, problems = nf.check_fp_vs_content()
    assert not ok
    assert problems["flow_cell_id"] == ("WRONG", "FBE92725")
    assert problems["barcode"] == ("barcode02", "barcode01")


def test_nanopore_fastq_tolerates_missing_keys():
    nf = NanoporeFastq(_handle("readx"), expected_flowcell="FBE92725")
    ok, problems = nf.check_fp_vs_content()
    assert ok and problems == {}


def test_nanopore_fastq_rejects_non_header():
    with pytest.raises(ValueError):
        NanoporeFastq(io.StringIO("ACGT\n"))
