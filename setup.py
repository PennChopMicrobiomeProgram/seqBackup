#!/usr/bin/env python

from distutils.core import setup

# Get version number from package
exec(open("seqBackupLib/version.py").read())

setup(
    name="seqBackup",
    version=__version__,
    description="Set of rules to organize our fastq storage on the server.",
    author="Ceylan Tanes",
    author_email="ctanes@gmail.com",
    url="https://github.com/PennChopMicrobiomeProgram",
    packages=["seqBackupLib"],
    scripts=["scripts/backup_illumina.py"],  # ,
    # install_requires=["pandas", "biopython"]
)
