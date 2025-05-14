import gzip
import pytest
from pathlib import Path
from seqBackupLib.backup import DEFAULT_MIN_FILE_SIZE
from seqBackupLib.illumina import IlluminaDir, IlluminaFastq, MACHINE_TYPES


machine_fixtures = {
    "A": "novaseq_dir",
    "D": "hiseq_dir",
    "LH": "novaseqx_dir",
    "M": "miseq_dir",
    "NB": "miniseq_dir",
    "VH": "nextseq_dir",
}


@pytest.mark.parametrize("machine_type", MACHINE_TYPES.keys())
def test_illumina_fastq(machine_type, request):
    fixture_name = machine_fixtures.get(machine_type)
    if not fixture_name:
        raise ValueError(
            f"All supported machine types must be tested. Missing: {machine_type}"
        )

    fp = request.getfixturevalue(fixture_name)

    with gzip.open(fp / "Undetermined_S0_L001_R1_001.fastq.gz", "rt") as f:
        r1 = IlluminaFastq(f)

    print("FASTQ info: ", r1.fastq_info, "\nFolder info: ", r1.folder_info)
    assert r1.machine_type == MACHINE_TYPES[machine_type]
    assert r1.check_fp_vs_content()[0], r1.check_fp_vs_content()
    assert not r1.check_file_size(DEFAULT_MIN_FILE_SIZE)
    assert r1.check_file_size(100)
    assert r1.check_index_read_exists()


@pytest.mark.parametrize("machine_type", MACHINE_TYPES.keys())
def test_illumina_dir(machine_type, request):
    fixture_name = machine_fixtures.get(machine_type)
    if not fixture_name:
        raise ValueError(
            f"All supported machine types must be tested. Missing: {machine_type}"
        )

    fp = request.getfixturevalue(fixture_name)

    d = IlluminaDir(fp.name)
