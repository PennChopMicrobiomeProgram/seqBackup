import pytest
from pathlib import Path
from seqBackupLib.backup import backup_fastq, build_fp_to_archive, return_md5


def test_build_fp_to_archive():
    archive = build_fp_to_archive(
        Path("Undetermined_S0_L001_R1_001.fastq.gz"), True, "1"
    )
    assert archive == [
        Path("Undetermined_S0_L001_R1_001.fastq.gz"),
        Path("Undetermined_S0_L001_R2_001.fastq.gz"),
        Path("Undetermined_S0_L001_I1_001.fastq.gz"),
        Path("Undetermined_S0_L001_I2_001.fastq.gz"),
    ]

    archive = build_fp_to_archive(
        Path("Undetermined_S0_L001_R1_001.fastq.gz"), False, "1"
    )
    assert archive == [
        Path("Undetermined_S0_L001_R1_001.fastq.gz"),
        Path("Undetermined_S0_L001_R2_001.fastq.gz"),
    ]

    archive = build_fp_to_archive(
        Path("Undetermined_S0_L002_R1_001.fastq.gz"), True, "2"
    )
    assert archive == [
        Path("Undetermined_S0_L002_R1_001.fastq.gz"),
        Path("Undetermined_S0_L002_R2_001.fastq.gz"),
        Path("Undetermined_S0_L002_I1_001.fastq.gz"),
        Path("Undetermined_S0_L002_I2_001.fastq.gz"),
    ]

    with pytest.raises(IOError):
        build_fp_to_archive(Path("Undetermined_S0_L001_R2_001.fastq.gz"), True, "1")


def test_return_md5(tmp_path):
    test_file = tmp_path / "test.txt"
    with open(test_file, "w") as f:
        f.write("Hello, World!")

    md5_hash = return_md5(test_file)
    assert md5_hash == "65a8e27d8879283831b664bd8b7f0ad4"  # MD5 hash of "Hello, World!"


def test_backup_fastq(tmp_path, full_miseq_dir):
    raw = tmp_path / "raw_reads"
    raw.mkdir(parents=True, exist_ok=True)
    sample_sheet_fp = full_miseq_dir / "sample_sheet.csv"

    backup_fastq(
        full_miseq_dir / "Undetermined_S0_L001_R1_001.fastq.gz",
        raw,
        sample_sheet_fp,
        True,
        100,
    )
    backup_fastq(
        full_miseq_dir / "Undetermined_S0_L002_R1_001.fastq.gz",
        raw,
        sample_sheet_fp,
        True,
        100,
    )

    with pytest.raises(FileNotFoundError):
        backup_fastq(
            full_miseq_dir / "Undetermined_S0_L003_R1_001.fastq.gz",
            raw,
            sample_sheet_fp,
            True,
            100,
        )
