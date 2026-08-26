import gzip
import os
import shutil
import stat
from pathlib import Path

import pytest

import seqBackupLib.backup_nanopore as bn
from seqBackupLib.backup_nanopore import (
    DEFAULT_DESTINATION_DIR,
    backup_nanopore,
    find_fastq_groups,
    main,
    return_md5,
)


def _count_records(fp) -> int:
    with gzip.open(fp, "rt") as handle:
        return sum(1 for _ in handle) // 4


def test_find_fastq_groups_barcoded(minion_dir):
    groups = find_fastq_groups(minion_dir / "fastq_pass")
    assert set(groups) == {"barcode01", "barcode02", "unclassified"}
    # chunks come back in natural (numeric) order, not lexical (_10 before _5)
    assert [p.name for p in groups["barcode01"]] == [
        "FBE92725_pass_barcode01_8b1f29fd_117e08fc_0.fastq.gz",
        "FBE92725_pass_barcode01_8b1f29fd_117e08fc_5.fastq.gz",
        "FBE92725_pass_barcode01_8b1f29fd_117e08fc_10.fastq.gz",
    ]


def test_find_fastq_groups_non_barcoded(non_barcoded_nanopore_dir):
    groups = find_fastq_groups(non_barcoded_nanopore_dir / "fastq_pass")
    assert set(groups) == {None}


def test_backup_nanopore_archives_reads_and_reports(
    tmp_path, minion_dir, nanopore_sample_sheet
):
    dest = tmp_path / "archive"

    out_dir = backup_nanopore(minion_dir, dest, nanopore_sample_sheet, min_total_size=1)

    assert out_dir == dest / minion_dir.name
    assert out_dir.is_dir()

    # the fastq_pass/<barcode>/ layout is preserved, one concatenated file per
    # barcode folder (each holding all 3 single-read chunks)
    for barcode in ("barcode01", "barcode02", "unclassified"):
        fq = out_dir / "fastq_pass" / barcode / f"{barcode}.fastq.gz"
        assert fq.is_file()
        assert _count_records(fq) == 3

    # reports kept
    assert (out_dir / "final_summary_FBE92725_x.txt").is_file()
    assert (out_dir / "report_FBE92725_x.html").is_file()
    assert (out_dir / "report_FBE92725_x.json").is_file()
    assert (out_dir / "report_FBE92725_x.md").is_file()
    assert (out_dir / "sample_sheet_FBE92725_x.csv").is_file()
    assert (out_dir / "throughput_FBE92725_x.csv").is_file()
    assert (out_dir / "barcode_alignment_FBE92725_x.tsv").is_file()

    # the supplied (official) sample sheet is archived
    assert (out_dir / nanopore_sample_sheet.name).is_file()

    # noise left behind
    assert not (out_dir / ".DS_Store").exists()
    assert not (out_dir / ".Rhistory").exists()
    assert not (out_dir / "sequencing_summary_FBE92725_x.txt").exists()

    # md5 manifest covers every concatenated fastq, keyed by archive-relative path
    md5_fp = out_dir / f"{minion_dir.name}.md5"
    assert md5_fp.is_file()
    entries = dict(line.split("\t") for line in md5_fp.read_text().splitlines() if line)
    assert set(entries) == {
        "fastq_pass/barcode01/barcode01.fastq.gz",
        "fastq_pass/barcode02/barcode02.fastq.gz",
        "fastq_pass/unclassified/unclassified.fastq.gz",
    }
    bc01 = out_dir / "fastq_pass" / "barcode01" / "barcode01.fastq.gz"
    assert entries["fastq_pass/barcode01/barcode01.fastq.gz"].strip() == return_md5(
        bc01
    )

    # archived fastqs are read-only
    assert not (stat.S_IMODE(os.stat(bc01).st_mode) & stat.S_IWUSR)


def test_backup_nanopore_non_barcoded(
    tmp_path, non_barcoded_nanopore_dir, nanopore_sample_sheet
):
    dest = tmp_path / "archive"
    out_dir = backup_nanopore(
        non_barcoded_nanopore_dir, dest, nanopore_sample_sheet, min_total_size=1
    )
    fq = out_dir / "fastq_pass" / "FBE92725.fastq.gz"
    assert fq.is_file()
    assert _count_records(fq) == 3


def test_backup_nanopore_missing_sample_sheet(tmp_path, minion_dir):
    with pytest.raises(IOError, match="Sample sheet does not exist"):
        backup_nanopore(
            minion_dir,
            tmp_path / "archive",
            tmp_path / "nope.csv",
            min_total_size=1,
        )
    # bailed out before creating anything
    assert not (tmp_path / "archive").exists()


def test_backup_nanopore_dummy_sample_sheet_ok(tmp_path, minion_dir):
    # transfer runs with no metadata use a placeholder file
    dummy = tmp_path / "dummy.tsv"
    dummy.write_text("no metadata; transfer run\n")

    out_dir = backup_nanopore(minion_dir, tmp_path / "archive", dummy, min_total_size=1)
    assert (out_dir / "dummy.tsv").is_file()


def test_backup_nanopore_refuses_existing_archive(
    tmp_path, minion_dir, nanopore_sample_sheet
):
    dest = tmp_path / "archive"
    backup_nanopore(minion_dir, dest, nanopore_sample_sheet, min_total_size=1)
    with pytest.raises(FileExistsError):
        backup_nanopore(minion_dir, dest, nanopore_sample_sheet, min_total_size=1)


def test_backup_nanopore_size_check(tmp_path, minion_dir, nanopore_sample_sheet):
    dest = tmp_path / "archive"

    with pytest.raises(ValueError, match="below the minimum"):
        backup_nanopore(minion_dir, dest, nanopore_sample_sheet, min_total_size=10**12)

    # nothing left behind after the failure
    assert not (dest / minion_dir.name).exists()
    assert list(dest.iterdir()) == []

    with pytest.warns(UserWarning, match="below the minimum"):
        out_dir = backup_nanopore(
            minion_dir,
            dest,
            nanopore_sample_sheet,
            min_total_size=10**12,
            allow_check_failures=True,
        )
    assert out_dir.is_dir()


def test_backup_nanopore_header_mismatch(tmp_path, minion_dir, nanopore_sample_sheet):
    # rewrite barcode01 chunks with the wrong flow cell id
    bc = minion_dir / "fastq_pass" / "barcode01"
    for chunk in bc.glob("*.fastq.gz"):
        with gzip.open(chunk, "wt") as handle:
            handle.write("@readx flow_cell_id=WRONG barcode=barcode01\nACGT\n+\nIIII\n")

    dest = tmp_path / "archive"
    with pytest.raises(ValueError, match="header info does not match"):
        backup_nanopore(minion_dir, dest, nanopore_sample_sheet, min_total_size=1)

    with pytest.warns(UserWarning, match="header info does not match"):
        backup_nanopore(
            minion_dir,
            dest,
            nanopore_sample_sheet,
            min_total_size=1,
            allow_check_failures=True,
        )


def test_backup_nanopore_missing_fastq_pass(
    tmp_path, minion_dir, nanopore_sample_sheet
):
    shutil.rmtree(minion_dir / "fastq_pass")
    with pytest.raises(IOError):
        backup_nanopore(
            minion_dir, tmp_path / "archive", nanopore_sample_sheet, min_total_size=1
        )


def test_main_returns_archive_path(tmp_path, p2i_dir, nanopore_sample_sheet):
    dest = tmp_path / "archive"
    out_dir = main(
        [
            "--run-dir",
            str(p2i_dir),
            "--destination-dir",
            str(dest),
            "--sample-sheet",
            str(nanopore_sample_sheet),
            "--min-total-size",
            "1",
        ]
    )
    assert out_dir == dest / p2i_dir.name
    assert out_dir.is_dir()
    assert (out_dir / "fastq_pass" / "barcode01" / "barcode01.fastq.gz").is_file()
    assert (out_dir / nanopore_sample_sheet.name).is_file()


def test_sample_sheet_is_required(tmp_path, minion_dir):
    with pytest.raises(SystemExit):
        main(["--run-dir", str(minion_dir), "--destination-dir", str(tmp_path)])


def test_destination_dir_defaults_to_raw_data(monkeypatch, tmp_path, minion_dir):
    assert DEFAULT_DESTINATION_DIR == Path("/mnt/isilon/microbiome/raw_data")

    seen = {}

    def fake_backup(run_dir, dest_dir, *args, **kwargs):
        seen["dest"] = dest_dir
        return tmp_path

    monkeypatch.setattr(bn, "backup_nanopore", fake_backup)
    main(["--run-dir", str(minion_dir), "--sample-sheet", str(tmp_path / "s.csv")])
    assert seen["dest"] == DEFAULT_DESTINATION_DIR
