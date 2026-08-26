# seqBackup

Logic for parsing Illumina headers and folders as well as for archiving reads.
Also handles Oxford Nanopore (MinKNOW) run folders.

## Backing up a run

Illumina (paired-end, demultiplexed later from the archived `Undetermined` reads):

```
backup_illumina --forward-reads .../Undetermined_S0_L001_R1_001.fastq.gz \
  --destination-dir /path/to/archive --sample-sheet .../sample_sheet.csv
```

Nanopore (single-end; the reads are already demultiplexed by MinKNOW). Point it
at the MinKNOW run folder (`<date>_<time>_<position>_<flowcell>_<runid>[...]`).
Each `fastq_pass/<barcode>/` subdirectory's chunk files are concatenated into a
single `fastq_pass/<barcode>/<barcode>.fastq.gz` in the archive, preserving the
per-barcode folder layout the ONT tools expect (a non-multiplexed run yields one
`fastq_pass/<flowcell>.fastq.gz`). The `final_summary_*.txt`, `report_*.html`, and
the top-level `*.csv` / `*.json` / `*.tsv` / `*.md` run reports are copied into the
archive alongside a `.md5` manifest.

```
backup_nanopore --run-dir .../20260401_1506_MN47822_FBE92725_8b1f29fd
```

`--destination-dir` defaults to `/mnt/isilon/microbiome/raw_data`.

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