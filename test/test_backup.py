import pytest
from pathlib import Path
import gzip
from seqBackupLib.backup import (
    backup_fastq,
    build_fp_to_archive,
    cli,
    main,
    return_md5,
)


def _write_fastq(fp: Path, header: str) -> None:
    sequence = "N" * 10
    content = f"{header}\n{sequence}\n+\n{'#' * len(sequence)}\n"
    with gzip.open(fp, "wt") as handle:
        handle.write(content)


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
    # main() is also called directly as a library function -- auto_bfx's
    # archive task calls seqBackupLib.backup.main(args) and returns the
    # result as the archive location -- so main() itself must keep
    # returning the Path. See test_cli_* below for the console-script entry
    # point, which wraps main() instead of replacing it.
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


def test_cli_prints_archive_path_and_returns_zero(capsys, tmp_path, full_miseq_dir):
    # cli() is what [project.scripts] points backup_illumina at. The
    # installed console-script wraps it as sys.exit(cli()), and sys.exit()
    # treats any non-None, non-int argument as an error -- printed to
    # stderr, exit code 1 -- so cli() must return 0/None rather than the
    # Path main() returns, or every successful run would report as a
    # failure.
    raw = tmp_path / "raw_reads"
    raw.mkdir(parents=True, exist_ok=True)
    sample_sheet_fp = full_miseq_dir / "sample_sheet.csv"
    expected_dir = raw / "250407_M03543_0443_000000000-DTHBL_L001"

    result = cli(
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

    assert result == 0
    assert capsys.readouterr().out.strip() == str(expected_dir)
    assert expected_dir.is_dir()


def test_allow_check_failures_continues_archive(tmp_path):
    run_dir = tmp_path / "240101_M01234_0001_ABCDEFGX"
    run_dir.mkdir(parents=True, exist_ok=True)
    sample_sheet_fp = run_dir / "sample_sheet.csv"
    sample_sheet_fp.write_text(
        "[Header]\nIEMFileVersion,4\n[Data]\nSample_ID,Sample_Name\nS1,S1\n"
    )

    header = "@M01234:1:ZZZZZZ:1:1101:10000:10000 1:N:0:ATCACG"
    for name in [
        "Undetermined_S0_L001_R1_001.fastq.gz",
        "Undetermined_S0_L001_R2_001.fastq.gz",
        "Undetermined_S0_L001_I1_001.fastq.gz",
        "Undetermined_S0_L001_I2_001.fastq.gz",
    ]:
        _write_fastq(run_dir / name, header)

    raw = tmp_path / "raw_reads"
    raw.mkdir(parents=True, exist_ok=True)

    with pytest.raises(ValueError, match="header information don't match"):
        backup_fastq(
            run_dir / "Undetermined_S0_L001_R1_001.fastq.gz",
            raw,
            sample_sheet_fp,
            True,
            1,
        )

    with pytest.warns(UserWarning, match="header information don't match"):
        out_dir = backup_fastq(
            run_dir / "Undetermined_S0_L001_R1_001.fastq.gz",
            raw,
            sample_sheet_fp,
            True,
            1,
            allow_check_failures=True,
        )

    assert out_dir.is_dir()
    md5_fp = out_dir / f"{out_dir.name}.md5"
    assert md5_fp.is_file()
