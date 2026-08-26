from __future__ import annotations

import gzip
import pytest
from pathlib import Path


def _write_ont_fastq(fp: Path, headers: list[str]) -> None:
    sequence = "ACGTACGTAC"
    quality = "IIIIIIIIII"
    with gzip.open(fp, "wt") as handle:
        for header in headers:
            handle.write(f"@{header}\n{sequence}\n+\n{quality}\n")


def setup_nanopore_dir(
    base: Path,
    run_name: str,
    flowcell_id: str,
    barcodes: list[str],
    n_chunks: int = 3,
) -> Path:
    """Build a realistic MinKNOW run folder.

    ``barcodes`` is the list of ``fastq_pass`` subdirectory names to create; pass
    an empty list for a non-multiplexed run (chunks land directly in
    ``fastq_pass/``).
    """
    run_dir = base / run_name
    fastq_pass = run_dir / "fastq_pass"
    fastq_pass.mkdir(parents=True, exist_ok=True)

    def _chunks_into(target: Path, label: str, barcode: str | None) -> None:
        for i in range(n_chunks):
            fields = [
                f"{label}-read-{i}-0000-0000-000000000000",
                "runid=8b1f29fd-6233-456a-9c6c-d8d72cecc734",
                "ch=5",
                "start_time=2026-04-01T15:14:56.068015-04:00",
                f"flow_cell_id={flowcell_id}",
                "protocol_group_id=CHOPMC-611_NimaGen_04012026",
                "sample_id=",
            ]
            if barcode is not None:
                fields.append(f"barcode={barcode}")
                fields.append(f"barcode_alias={barcode}")
            fields.append(
                "basecall_model_version_id=dna_r10.4.1_e8.2_400bps_hac@v4.3.0"
            )
            # chunk indices deliberately out of lexical order (_10 sorts after _2)
            _write_ont_fastq(
                target
                / f"{flowcell_id}_pass_{label}_8b1f29fd_117e08fc_{i * 5}.fastq.gz",
                [" ".join(fields)],
            )

    if barcodes:
        for barcode in barcodes:
            sub = fastq_pass / barcode
            sub.mkdir(parents=True, exist_ok=True)
            _chunks_into(sub, barcode, barcode)
    else:
        _chunks_into(fastq_pass, "all", None)

    # Top-level run reports/metadata that should be archived.
    (run_dir / f"final_summary_{flowcell_id}_x.txt").write_text(
        f"instrument=MN47822\nflow_cell_id={flowcell_id}\n"
    )
    (run_dir / f"report_{flowcell_id}_x.html").write_text("<html>ONT report</html>")
    (run_dir / f"report_{flowcell_id}_x.json").write_text('{"run": "ok"}')
    (run_dir / f"report_{flowcell_id}_x.md").write_text("# ONT report\n")
    (run_dir / f"sample_sheet_{flowcell_id}_x.csv").write_text(
        f"flow_cell_id\n{flowcell_id}\n"
    )
    (run_dir / f"throughput_{flowcell_id}_x.csv").write_text("minute,reads\n0,0\n")
    (run_dir / f"barcode_alignment_{flowcell_id}_x.tsv").write_text("barcode\talias\n")

    # Noise that must NOT be archived.
    (run_dir / ".DS_Store").write_text("junk")
    (run_dir / ".Rhistory").write_text("")
    (run_dir / f"sequencing_summary_{flowcell_id}_x.txt").write_text("huge summary\n")

    return run_dir


@pytest.fixture
def minion_dir(tmp_path) -> Path:
    return setup_nanopore_dir(
        tmp_path,
        "20260401_1506_MN47822_FBE92725_8b1f29fd",
        "FBE92725",
        ["barcode01", "barcode02", "unclassified"],
    )


@pytest.fixture
def p2i_dir(tmp_path) -> Path:
    return setup_nanopore_dir(
        tmp_path,
        "20260408_1652_P2I-00513-B_PBK70557_2fe9fcef_Promethion_Training",
        "PBK70557",
        ["barcode01", "barcode02", "unclassified"],
    )


@pytest.fixture
def non_barcoded_nanopore_dir(tmp_path) -> Path:
    return setup_nanopore_dir(
        tmp_path,
        "20260401_1506_MN47822_FBE92725_8b1f29fd",
        "FBE92725",
        [],
    )


def setup_illumina_dir(fp: Path, r1: str, r1_lines: list[str]) -> Path:
    fp.mkdir(parents=True, exist_ok=True)

    r1_fp = fp / r1
    with gzip.open(r1_fp, "wt") as f:
        f.writelines(r1_lines)

    (fp / r1.replace("R1", "R2")).touch()
    (fp / r1.replace("R1", "I1")).touch()
    (fp / r1.replace("R1", "I2")).touch()
    (fp / "sample_sheet.csv").touch()

    return fp


@pytest.fixture
def novaseq_dir(tmp_path) -> Path:
    return setup_illumina_dir(
        tmp_path / "250218_A00901_1295_BHTKCGDRX5",
        "Undetermined_S0_L001_R1_001.fastq.gz",
        [
            "@A00901:1295:HTKCGDRX5:1:2101:1054:1000 1:N:0:NAGTGTTAGG+CGGAACTAGC\n",
            "GTAAAAAGCTAGATTTTCGCGATTTACCAGACGAACTANTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n",
            "+\n",
            "FFFFF,FFF:FFFFFF,FFFFFFFFFFFFFFFFFFFFF#,###############################################################################################################\n",
            "@A00901:1295:HTKCGDRX5:1:2101:1090:1000 1:N:0:AATTCTTGGA+AAGTTGACAA\n",
            "GCTGCAATATGCGCCAACAAAACCGGTGGATAAAAAGGTTTCGTAATATAGTCATCNCNGNCNTNTNCNANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNATTTCTCCGC\n",
            "+\n",
            ",FFFFFFFF:FFFFFFFFFFFF::FFFFFFFFF,FFFFFFF:F::FFFFFFFFFF:#F#F#F#F#:#F#F#######################################################################FFFFFFFFFF\n",
        ],
    )


@pytest.fixture
def hiseq_dir(tmp_path) -> Path:
    return setup_illumina_dir(
        tmp_path / "201118_D00728_0139_ACD5C3ANXX",
        "Undetermined_S0_L001_R1_001.fastq.gz",
        [
            "@D00728:139:CD5C3ANXX:1:1101:1228:2123 1:N:0:ATCTCAGG+CCTAGAGT\n",
            "NTGCGCAGGGGGACCTGCACCGGCATCCCCTGTACCGGCGGGGCGCTCAGGCTGAATGCGCCGTCCTGCATCAGTACCGACTCCGGCTCGATGGCTTTATCCTGTCTCTTATACACATCTCCGAGC\n",
            "+\n",
            "#:<>AE<EGGGGCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGEGGGGGGGGGGGGGGEGGGGDEGGGGGGGGGG;EGGGGGGEDDEGGGGGGGGGGGGGGBEGGDDDCDC\n",
            "@D00728:139:CD5C3ANXX:1:1101:1176:2132 1:N:0:CTCTCTAC+TTCTAGCT\n",
            "NATCTCGACCTTTTGGGTCAGTCCGTTTTTCGATTTATAAACCGCCGGCAATACACGGCTCTTGACATCCGATAATCCCGGAGACGCATGGACCGACAATGTCCCAATCCTGTCTCTTATACACAT\n",
            "+\n",
            "#<<>BGGGGGGGGGDDG=CGGGGF=FGGGGGGGGGGGGGGGGF<EG@EEGGCGGGGGGGGDGGGGEGGGGGGGDGGGGGBGGGGGGGGBGGGGEGGGGG77CEEE=GGG.DGDEGEGGB/8/CDB@\n",
        ],
    )


@pytest.fixture
def novaseqx_dir(tmp_path) -> Path:
    return setup_illumina_dir(
        tmp_path / "20250429_LH00732_0028_A22YJWWLT3",
        "Undetermined_S0_L001_R1_001.fastq.gz",
        [
            "@LH00732:28:22YJWWLT3:1:1101:1213:1080 1:N:0:CCTCCGTCCA+CACCGATGTG\n",
            "ACGT\n",
            "+\n",
            "IIII\n",
        ],
    )


@pytest.fixture
def miseq_dir(tmp_path) -> Path:
    return setup_illumina_dir(
        tmp_path / "250407_M03543_0443_000000000-DTHBL",
        "Undetermined_S0_L001_R1_001.fastq.gz",
        [
            "@M03543:443:000000000-DTHBL:1:1101:16223:1348 1:N:0:TTTTTTTTTTTT+TTCTTTTTCCTT\n",
            "TCTTCCCTCTTTCTTCTTTCTTCCTCCCTTCCCTTCTCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n",
            "+\n",
            ">>>>A1C1BB1B3A333BB33B311100BA000BBBE122B110A//AA/>//>///>///<<<//<<--<---:------9---99-999-9---9-99-99---9----9-9-999-9-9-9-----999--9--9--99-9/99/9/--99---999-9-999--9--999-9-9-9-99-9--9-999-99-999999999@>-9-99--99---999--999@999-9999@>---9------9-9\n",
            "@M03543:443:000000000-DTHBL:1:1101:15497:1351 1:N:0:TTTTTTTTTTTT+TTCTTTTTCCTC\n",
            "TCTTCCCTCTTTCTTCTTTCTTCCTCCCTTCCCTTCTCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCTCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n",
            "+\n",
            ">>>>A1C1BB1B3B333BB33B311100BB000BBCD122B110A//AA/>//>///>////<<//<<--<-:-:-9--9-9---9--999--///9//9/99--------9-9-9999999-9-----999-----9--99-9/99/9/--99-9--99-9-999--9--9-9--9--9-99-9----999-99-9999999999@9--99--99---99---999@999--999=>---9------9-9\n",
        ],
    )


@pytest.fixture
def miniseq_dir(tmp_path) -> Path:
    return setup_illumina_dir(
        tmp_path / "210612_NB551353_0107_AHWJFCAFX2",
        "Undetermined_S0_L001_R1_001.fastq.gz",
        [
            "@NB551353:107:HWJFCAFX2:1:11101:1486:1048 1:N:0:TAATTAGCGT+NNNTTAACCA\n",
            "GAAATNGACCGCCTCAATGAGGTTGCCAAGAATTTAAATGAATCTCTCATCGATCTCCAAGAACTTGGAAAGTA\n",
            "+\n",
            "AAAAA#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEE\n",
            "@NB551353:107:HWJFCAFX2:1:11101:6713:1048 1:N:0:GAAGACTAGA+NNNTTCTAGT\n",
            "GGCTANATCTGAGGACAAGAGGGCAAAAGTTACTAGTGCTATGCAGACAATGCTTTTCACTATGCTTAGAAAGT\n",
            "+\n",
            "AAAAA#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE\n",
        ],
    )


@pytest.fixture
def nextseq_dir(tmp_path) -> Path:
    return setup_illumina_dir(
        tmp_path / "250101_VH12345_0022_222C2NYNX",
        "Undetermined_S0_L001_R1_001.fastq.gz",
        [
            "@VH12345:22:222C2NYNX:1:1101:18286:1000 1:N:0:GGCACTAAGG+GTTGACCTGA\n",
            "ACGT\n",
            "+\n",
            "IIII\n",
        ],
    )


@pytest.fixture
def sh_dir(tmp_path) -> Path:
    return setup_illumina_dir(
        tmp_path / "20250610_SH00024_0006_ASC2107697-SC3",
        "Undetermined_S0_L001_R1_001.fastq.gz",
        [
            "@SH00024:6:SC2107697-SC3:1:1101:18286:1000 1:N:0:ACGT+ACGT\n",
            "ACGT\n",
            "+\n",
            "IIII\n",
        ],
    )


@pytest.fixture
def full_miseq_dir(tmp_path) -> Path:
    # Lane 1
    setup_illumina_dir(
        tmp_path / "250407_M03543_0443_000000000-DTHBL",
        "Undetermined_S0_L001_R1_001.fastq.gz",
        [
            "@M03543:443:000000000-DTHBL:1:1101:16223:1348 1:N:0:TTTTTTTTTTTT+TTCTTTTTCCTT\n",
            "TCTTCCCTCTTTCTTCTTTCTTCCTCCCTTCCCTTCTCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n",
            "+\n",
            ">>>>A1C1BB1B3A333BB33B311100BA000BBBE122B110A//AA/>//>///>///<<<//<<--<---:------9---99-999-9---9-99-99---9----9-9-999-9-9-9-----999--9--9--99-9/99/9/--99---999-9-999--9--999-9-9-9-99-9--9-999-99-999999999@>-9-99--99---999--999@999-9999@>---9------9-9\n",
            "@M03543:443:000000000-DTHBL:1:1101:15497:1351 1:N:0:TTTTTTTTTTTT+TTCTTTTTCCTC\n",
            "TCTTCCCTCTTTCTTCTTTCTTCCTCCCTTCCCTTCTCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCTCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n",
            "+\n",
            ">>>>A1C1BB1B3B333BB33B311100BB000BBCD122B110A//AA/>//>///>////<<//<<--<-:-:-9--9-9---9--999--///9//9/99--------9-9-9999999-9-----999-----9--99-9/99/9/--99-9--99-9-999--9--9-9--9--9-99-9----999-99-9999999999@9--99--99---99---999@999--999=>---9------9-9\n",
        ],
    )
    setup_illumina_dir(
        tmp_path / "250407_M03543_0443_000000000-DTHBL",
        "Undetermined_S0_L001_R2_001.fastq.gz",
        [
            "@M03543:443:000000000-DTHBL:1:1101:16223:1348 2:N:0:TTTTTTTTTTTT+TTCTTTTTCCTT\n",
            "TCTTCCCTCTTTCTTCTTTCTTCCTCCCTTCCCTTCTCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n",
            "+\n",
            ">>>>A1C1BB1B3A333BB33B311100BA000BBBE122B110A//AA/>//>///>///<<<//<<--<---:------9---99-999-9---9-99-99---9----9-9-999-9-9-9-----999--9--9--99-9/99/9/--99---999-9-999--9--999-9-9-9-99-9--9-999-99-999999999@>-9-99--99---999--999@999-9999@>---9------9-9\n",
            "@M03543:443:000000000-DTHBL:1:1101:15497:1351 2:N:0:TTTTTTTTTTTT+TTCTTTTTCCTC\n",
            "TCTTCCCTCTTTCTTCTTTCTTCCTCCCTTCCCTTCTCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCTCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n",
            "+\n",
            ">>>>A1C1BB1B3B333BB33B311100BB000BBCD122B110A//AA/>//>///>////<<//<<--<-:-:-9--9-9---9--999--///9//9/99--------9-9-9999999-9-----999-----9--99-9/99/9/--99-9--99-9-999--9--9-9--9--9-99-9----999-99-9999999999@9--99--99---99---999@999--999=>---9------9-9\n",
        ],
    )
    setup_illumina_dir(
        tmp_path / "250407_M03543_0443_000000000-DTHBL",
        "Undetermined_S0_L001_I1_001.fastq.gz",
        [
            "@M03543:443:000000000-DTHBL:1:2106:17605:1940 1:N:0:TTTTTTTTTTTT+TCTTTCCCTACA\n",
            "TTTTTTTTTTTT\n",
            "+\n",
            "111>111>0000\n",
            "@M03543:443:000000000-DTHBL:1:2106:14807:1943 1:N:0:TTTTTTTTTTTT+TCTTTCCCTACA\n",
            "TTTTTTTTTTTT\n",
            "+\n",
            "1>1111>100>0\n",
        ],
    )
    setup_illumina_dir(
        tmp_path / "250407_M03543_0443_000000000-DTHBL",
        "Undetermined_S0_L001_I2_001.fastq.gz",
        [
            "@M03543:443:000000000-DTHBL:1:2106:17605:1940 2:N:0:TTTTTTTTTTTT+TCTTTCCCTACA\n",
            "TTTTTTTTTTTT\n",
            "+\n",
            "111>111>0000\n",
            "@M03543:443:000000000-DTHBL:1:2106:14807:1943 2:N:0:TTTTTTTTTTTT+TCTTTCCCTACA\n",
            "TTTTTTTTTTTT\n",
            "+\n",
            "1>1111>100>0\n",
        ],
    )

    # Lane 2
    setup_illumina_dir(
        tmp_path / "250407_M03543_0443_000000000-DTHBL",
        "Undetermined_S0_L002_R1_001.fastq.gz",
        [
            "@M03543:443:000000000-DTHBL:2:1101:16223:1348 1:N:0:TTTTTTTTTTTT+TTCTTTTTCCTT\n",
            "TCTTCCCTCTTTCTTCTTTCTTCCTCCCTTCCCTTCTCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n",
            "+\n",
            ">>>>A1C1BB1B3A333BB33B311100BA000BBBE122B110A//AA/>//>///>///<<<//<<--<---:------9---99-999-9---9-99-99---9----9-9-999-9-9-9-----999--9--9--99-9/99/9/--99---999-9-999--9--999-9-9-9-99-9--9-999-99-999999999@>-9-99--99---999--999@999-9999@>---9------9-9\n",
            "@M03543:443:000000000-DTHBL:2:1101:15497:1351 1:N:0:TTTTTTTTTTTT+TTCTTTTTCCTC\n",
            "TCTTCCCTCTTTCTTCTTTCTTCCTCCCTTCCCTTCTCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCTCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n",
            "+\n",
            ">>>>A1C1BB1B3B333BB33B311100BB000BBCD122B110A//AA/>//>///>////<<//<<--<-:-:-9--9-9---9--999--///9//9/99--------9-9-9999999-9-----999-----9--99-9/99/9/--99-9--99-9-999--9--9-9--9--9-99-9----999-99-9999999999@9--99--99---99---999@999--999=>---9------9-9\n",
        ],
    )
    setup_illumina_dir(
        tmp_path / "250407_M03543_0443_000000000-DTHBL",
        "Undetermined_S0_L002_R2_001.fastq.gz",
        [
            "@M03543:443:000000000-DTHBL:2:1101:16223:1348 2:N:0:TTTTTTTTTTTT+TTCTTTTTCCTT\n",
            "TCTTCCCTCTTTCTTCTTTCTTCCTCCCTTCCCTTCTCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n",
            "+\n",
            ">>>>A1C1BB1B3A333BB33B311100BA000BBBE122B110A//AA/>//>///>///<<<//<<--<---:------9---99-999-9---9-99-99---9----9-9-999-9-9-9-----999--9--9--99-9/99/9/--99---999-9-999--9--999-9-9-9-99-9--9-999-99-999999999@>-9-99--99---999--999@999-9999@>---9------9-9\n",
            "@M03543:443:000000000-DTHBL:2:1101:15497:1351 2:N:0:TTTTTTTTTTTT+TTCTTTTTCCTC\n",
            "TCTTCCCTCTTTCTTCTTTCTTCCTCCCTTCCCTTCTCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCTCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n",
            "+\n",
            ">>>>A1C1BB1B3B333BB33B311100BB000BBCD122B110A//AA/>//>///>////<<//<<--<-:-:-9--9-9---9--999--///9//9/99--------9-9-9999999-9-----999-----9--99-9/99/9/--99-9--99-9-999--9--9-9--9--9-99-9----999-99-9999999999@9--99--99---99---999@999--999=>---9------9-9\n",
        ],
    )
    setup_illumina_dir(
        tmp_path / "250407_M03543_0443_000000000-DTHBL",
        "Undetermined_S0_L002_I1_001.fastq.gz",
        [
            "@M03543:443:000000000-DTHBL:2:2106:17605:1940 1:N:0:TTTTTTTTTTTT+TCTTTCCCTACA\n",
            "TTTTTTTTTTTT\n",
            "+\n",
            "111>111>0000\n",
            "@M03543:443:000000000-DTHBL:2:2106:14807:1943 1:N:0:TTTTTTTTTTTT+TCTTTCCCTACA\n",
            "TTTTTTTTTTTT\n",
            "+\n",
            "1>1111>100>0\n",
        ],
    )
    return setup_illumina_dir(
        tmp_path / "250407_M03543_0443_000000000-DTHBL",
        "Undetermined_S0_L002_I2_001.fastq.gz",
        [
            "@M03543:443:000000000-DTHBL:2:2106:17605:1940 2:N:0:TTTTTTTTTTTT+TCTTTCCCTACA\n",
            "TTTTTTTTTTTT\n",
            "+\n",
            "111>111>0000\n",
            "@M03543:443:000000000-DTHBL:2:2106:14807:1943 2:N:0:TTTTTTTTTTTT+TCTTTCCCTACA\n",
            "TTTTTTTTTTTT\n",
            "+\n",
            "1>1111>100>0\n",
        ],
    )
