import gzip
import importlib
from urllib.error import URLError

import pytest

from seqBackupLib.backup import DEFAULT_MIN_FILE_SIZE


machine_fixtures = {
    "A": "novaseq_dir",
    "D": "hiseq_dir",
    "LH": "novaseqx_dir",
    "M": "miseq_dir",
    "NB": "miniseq_dir",
    "VH": "nextseq_dir",
    "SH": "sh_dir",
}


@pytest.fixture
def illumina_module(monkeypatch):
    tsv_rows = ["instrument_code\tmachine_type"]
    fallback = {
        "VH": "Illumina-NextSeq",
        "D": "Illumina-HiSeq",
        "M": "Illumina-MiSeq",
        "A": "Illumina-NovaSeq",
        "NB": "Illumina-MiniSeq",
        "LH": "Illumina-NovaSeqX",
        "SH": "Illumina-MiSeq",
    }
    tsv_rows.extend(f"{code}\t{machine}" for code, machine in fallback.items())
    tsv = "\n".join(tsv_rows) + "\n"

    class FakeResponse:
        def __init__(self, data: str):
            self._data = data

        def read(self):
            return self._data.encode("utf-8")

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

    monkeypatch.setattr(
        "urllib.request.urlopen", lambda *args, **kwargs: FakeResponse(tsv)
    )
    import seqBackupLib.illumina as illumina

    return importlib.reload(illumina)


@pytest.mark.parametrize("machine_type", machine_fixtures.keys())
def test_illumina_fastq(machine_type, request, illumina_module):
    fixture_name = machine_fixtures.get(machine_type)
    if not fixture_name:
        raise ValueError(
            f"All supported machine types must be tested. Missing: {machine_type}"
        )

    fp = request.getfixturevalue(fixture_name)

    with gzip.open(fp / "Undetermined_S0_L001_R1_001.fastq.gz", "rt") as f:
        r1 = illumina_module.IlluminaFastq(f)

    print("FASTQ info: ", r1.fastq_info, "\nFolder info: ", r1.folder_info)
    assert r1.machine_type == illumina_module.MACHINE_TYPES[machine_type]
    assert r1.check_fp_vs_content()[0], r1.check_fp_vs_content()
    assert not r1.check_file_size(DEFAULT_MIN_FILE_SIZE)
    assert r1.check_file_size(100)
    assert r1.check_index_read_exists()


@pytest.mark.parametrize("machine_type", machine_fixtures.keys())
def test_illumina_dir(machine_type, request, illumina_module):
    fixture_name = machine_fixtures.get(machine_type)
    if not fixture_name:
        raise ValueError(
            f"All supported machine types must be tested. Missing: {machine_type}"
        )

    fp = request.getfixturevalue(fixture_name)

    d = illumina_module.IlluminaDir(fp.name)


def test_illumina_fastq_without_lane(novaseq_dir, illumina_module):
    original = novaseq_dir / "Undetermined_S0_L001_R1_001.fastq.gz"
    renamed = novaseq_dir / "Undetermined_S0_R1_001.fastq.gz"
    original.rename(renamed)
    with gzip.open(renamed, "rt") as f:
        r1 = illumina_module.IlluminaFastq(f)
    assert r1.check_fp_vs_content()[0]
    assert r1.build_archive_dir().endswith("L001")


def test_load_machine_types_from_tsv(monkeypatch):
    tsv = "instrument_code\tmachine_type\nZZ\tIllumina-Test\n"

    class FakeResponse:
        def __init__(self, data: str):
            self._data = data

        def read(self):
            return self._data.encode("utf-8")

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

    monkeypatch.setattr(
        "urllib.request.urlopen", lambda *args, **kwargs: FakeResponse(tsv)
    )
    import seqBackupLib.illumina as illumina

    illumina = importlib.reload(illumina)
    assert illumina.MACHINE_TYPES["ZZ"] == "Illumina-Test"


def test_load_machine_types_fallback_warning(monkeypatch):
    def raise_url_error(*args, **kwargs):
        raise URLError("network down")

    monkeypatch.setattr("urllib.request.urlopen", raise_url_error)
    import seqBackupLib.illumina as illumina

    with pytest.warns(RuntimeWarning, match="Falling back to bundled machine types"):
        illumina = importlib.reload(illumina)
    assert illumina.MACHINE_TYPES == illumina.MACHINE_TYPES_FALLBACK
