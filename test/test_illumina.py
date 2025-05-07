import gzip
import pytest
from pathlib import Path
from seqBackupLib.illumina import IlluminaFastq


def setup_illumina_dir(fp: Path, r1: str, r1_lines: list[str]) -> Path:
    fp.mkdir(parents=True, exist_ok=True)

    r1 = fp / r1
    with gzip.open(r1, "wt") as f:
        f.writelines(r1_lines)

    return fp


@pytest.fixture
def novaseq_dir(tmp_path) -> Path:
    return setup_illumina_dir(
        tmp_path / "250101_A12345_0001_A1234",
        "Undetermined_S0_L001_R1_001.fastq.gz",
        [
            "@A12345:1:1234:1:1101:1078:1091 R1:Y:0:ATTACTCG\n",
            "ACGT\n",
            "+\n",
            "IIII\n",
        ],
    )


@pytest.fixture
def hiseq_dir(tmp_path) -> Path:
    return setup_illumina_dir(
        tmp_path / "250101_D12345_0001_1234",
        "Undetermined_S0_L001_R1_001.fastq.gz",
        [
            "@D12345:1:1234:1:1101:1078:1091 R1:Y:0:ATTACTCG\n",
            "ACGT\n",
            "+\n",
            "IIII\n",
        ],
    )


@pytest.fixture
def novaseqx_dir(tmp_path) -> Path:
    return setup_illumina_dir(
        tmp_path / "250101_LH12345_0001_A1234",
        "Undetermined_S0_L001_R1_001.fastq.gz",
        [
            "@LH12345:1:1234:1:1101:1078:1091 R1:Y:0:ATTACTCG\n",
            "ACGT\n",
            "+\n",
            "IIII\n",
        ],
    )


@pytest.fixture
def miseq_dir(tmp_path) -> Path:
    return setup_illumina_dir(
        tmp_path / "250101_M12345_0001_1234",
        "Undetermined_S0_L001_R1_001.fastq.gz",
        [
            "@M12345:1:1234:1:1101:1078:1091 R1:Y:0:ATTACTCG\n",
            "ACGT\n",
            "+\n",
            "IIII\n",
        ],
    )


@pytest.fixture
def miniseq_dir(tmp_path) -> Path:
    return setup_illumina_dir(
        tmp_path / "250101_N12345_0001_1234",
        "Undetermined_S0_L001_R1_001.fastq.gz",
        [
            "@N12345:1:1234:1:1101:1078:1091 R1:Y:0:ATTACTCG\n",
            "ACGT\n",
            "+\n",
            "IIII\n",
        ],
    )


@pytest.fixture
def nextseq_dir(tmp_path) -> Path:
    return setup_illumina_dir(
        tmp_path / "250101_V12345_0001_1234",
        "Undetermined_S0_L001_R1_001.fastq.gz",
        [
            "@V12345:1:1234:1:1101:1078:1091 R1:Y:0:ATTACTCG\n",
            "ACGT\n",
            "+\n",
            "IIII\n",
        ],
    )


machine_fixtures = {
    "A": "novaseq_dir",
    "D": "hiseq_dir",
    "LH": "novaseqx_dir",
    "M": "miseq_dir",
    "N": "miniseq_dir",
    "V": "nextseq_dir",
}


@pytest.mark.parametrize("machine_type", IlluminaFastq.MACHINE_TYPES.keys())
def test_illumina_fastq(machine_type, request):
    fixture_name = machine_fixtures.get(machine_type)
    if not fixture_name:
        raise ValueError(f"All supported machine types must be tested. Missing: {machine_type}")

    fp = request.getfixturevalue(fixture_name)

    with gzip.open(fp / "Undetermined_S0_L001_R1_001.fastq.gz", "rt") as f:
        r1 = IlluminaFastq(f)