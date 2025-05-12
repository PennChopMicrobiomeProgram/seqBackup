# seqBackup

Logic for parsing Illumina headers and folders as well as for archiving reads.

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

To add a new machine type, add the new machine code to the `MACHINE_TYPES` map in `seqBackuplib/illumina.py`. In some cases, you may have to add machine specific parsing in `_parse_header` or `_parse_folder`. In `test/test_illumina.py`, we have a mechanism for requiring tests for each supported machine type. Add the new machine type to the `machine_fixtures` map and then create the fixture that it points to in `test/conftest.py`. Follow the pattern laid out by other fixtures and try to make the test data as realistic as possible.

### Incorporting new version

This software is the "source of truth" for Illumina file handling logic. Other software in our ecosystem depend on this logic including the sample registry and the automation pipeline. When you update this software you will have to then update the installed versions wherever it is deployed as a dependency. We don't bother with official GitHub releases and instead just point directly at the `master` branch, so usually it is a matter of running `pip install git+https://github.com/PennChopMicrobiomeProgram/seqBackup.git@master` from the host machine.