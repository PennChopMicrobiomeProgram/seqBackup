# seqBackup

Logic for parsing Illumina headers and folders as well as for archiving reads.
Also handles Oxford Nanopore (MinKNOW) run folders.

## Backing up a run

Illumina (paired-end, demultiplexed later from the archived `Undetermined` reads):

```
backup_illumina --forward-reads .../Undetermined_S0_L001_R1_001.fastq.gz \
  --destination-dir /path/to/archive --sample-sheet .../sample_sheet.csv
```

Nanopore (single-end, already basecalled by MinKNOW). Point it at the MinKNOW
run folder (`<date>_<time>_<position>_<flowcell>_<runid>[...]`). Every
`.fastq.gz` chunk under `fastq_pass/` (across all barcode subdirectories, if
any) is concatenated into one `<flowcell>.fastq.gz` in the archive -- the same
convention as Illumina, which archives the undemultiplexed `Undetermined`
reads and splits later. Demultiplexing is not lost: each pass read's fastq
header still carries `barcode=barcodeNN`/`barcode_alias=`, which is enough to
split back out with a header-based script; re-demultiplexing from raw sequence
(e.g. to change stringency) needs `dorado demux`. The `final_summary_*.txt`,
`report_*.html`, and the top-level `*.csv` / `*.json` / `*.tsv` / `*.md` run
reports are copied into the archive alongside a `.md5` manifest.

```
backup_nanopore --run-dir .../20260401_1506_MN47822_FBE92725_8b1f29fd \
  --sample-sheet .../CHOPMC-611_NimaGen_04012026.tsv
```

`--sample-sheet` is required so the run's metadata is always recorded (for a
transfer run with no metadata, point it at a dummy placeholder file).
`--destination-dir` defaults to `/mnt/isilon/microbiome/raw_data`. `--min-file-size`
(default 500MB) is the minimum size of the concatenated fastq.

Both commands take `--allow-check-failures` to archive despite failed validation
checks (a warning is emitted instead of an error).

## Dev

```
git clone https://github.com/PennChopMicrobiomeProgram/seqBackup.git
cd seqBackup/
python -m venv env
source env/bin/activate
pip install -e .
pip install black pytest
```

Before commits, make sure everything is well formatted and working:

```
black .
pytest test/
git commit ...
```

### Adding a new machine type

To add a new Illumina machine type, add the new machine code to the `MACHINE_TYPES` map in `seqBackuplib/illumina.py`. In some cases, you may have to add machine specific parsing in `_parse_header` or `_parse_folder`. In `test/test_illumina.py`, we have a mechanism for requiring tests for each supported machine type. Add the new machine type to the `machine_fixtures` map and then create the fixture that it points to in `test/conftest.py`. Follow the pattern laid out by other fixtures and try to make the test data as realistic as possible.

To add a new Nanopore device, add its position-token prefix to the `ONT_MACHINE_TYPES` map in `seqBackupLib/nanopore.py` (e.g. `MN47822` -> `MN`, `P2I-00513-B` -> `P2I`; see `extract_ont_device_code`). `test/test_nanopore.py` has the same "every machine type must be tested" mechanism: add the device to `ont_machine_fixtures` and create the matching run-folder fixture in `test/conftest.py`.

### Incorporating new version

This software is the "source of truth" for Illumina file handling logic. Other software in our ecosystem depend on this logic including the sample registry and the automation pipeline. When you update this software you will have to then update the installed versions wherever it is deployed as a dependency. We don't bother with official GitHub releases and instead just point directly at the `master` branch, so usually it is a matter of running `pip install git+https://github.com/PennChopMicrobiomeProgram/seqBackup.git@master` from the host machine.