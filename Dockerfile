# syntax=docker/dockerfile:1
FROM python:3.11-slim

ENV PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1

WORKDIR /app

# Install the seqBackup package and its dependencies
COPY pyproject.toml README.md ./
COPY seqBackupLib ./seqBackupLib
RUN pip install --upgrade pip \
    && pip install .

ENTRYPOINT ["backup_illumina"]
CMD ["--help"]
