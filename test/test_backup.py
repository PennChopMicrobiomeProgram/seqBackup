import pytest
from pathlib import Path
from seqBackupLib.backup import (
    backup_fastq,
    build_fp_to_archive,
    return_md5,
    main,
)


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

    archive = build_fp_to_archive(Path("Undetermined_S0_R1_001.fastq.gz"), True, "1")
    assert archive == [
        Path("Undetermined_S0_R1_001.fastq.gz"),
        Path("Undetermined_S0_R2_001.fastq.gz"),
        Path("Undetermined_S0_I1_001.fastq.gz"),
        Path("Undetermined_S0_I2_001.fastq.gz"),
    ]


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


def test_backup_fastq_without_lane(tmp_path, full_miseq_dir):
    raw = tmp_path / "raw_reads"
    raw.mkdir(parents=True, exist_ok=True)
    sample_sheet_fp = full_miseq_dir / "sample_sheet.csv"

    for lab in ["R1", "R2", "I1", "I2"]:
        (full_miseq_dir / f"Undetermined_S0_L001_{lab}_001.fastq.gz").rename(
            full_miseq_dir / f"Undetermined_S0_{lab}_001.fastq.gz"
        )

    backup_fastq(
        full_miseq_dir / "Undetermined_S0_R1_001.fastq.gz",
        raw,
        sample_sheet_fp,
        True,
        100,
    )

    out_dir = raw / "250407_M03543_0443_000000000-DTHBL_L001"
    assert (out_dir / "Undetermined_S0_L001_R1_001.fastq.gz").is_file()
    assert (out_dir / "Undetermined_S0_L001_R2_001.fastq.gz").is_file()


def test_main_returns_archive_path(tmp_path, full_miseq_dir):
    raw = tmp_path / "raw_reads"
    raw.mkdir(parents=True, exist_ok=True)
    sample_sheet_fp = full_miseq_dir / "sample_sheet.csv"

    out_dir = main(
        [
            "--forward-reads",
            str(full_miseq_dir / "Undetermined_S0_L001_R1_001.fastq.gz"),
            "--destination-dir",
            str(raw),
            "--sample-sheet",
            str(sample_sheet_fp),
            "--min-file-size",
            "100",
        ]
    )

    expected_dir = raw / "250407_M03543_0443_000000000-DTHBL_L001"
    assert out_dir == expected_dir
    assert expected_dir.is_dir()
