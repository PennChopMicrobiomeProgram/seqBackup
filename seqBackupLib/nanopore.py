from __future__ import annotations

import re
from io import TextIOWrapper
from pathlib import Path

# Oxford Nanopore device codes, keyed by the alphabetic prefix of the position
# token in a MinKNOW run folder name (e.g. ``MN47822`` -> ``MN``,
# ``P2I-00513-B`` -> ``P2I``).  Extend this map when a new device shows up; see
# the "Adding a new machine type" section of the README.
ONT_MACHINE_TYPES = {
    "MN": "ONT-MN",  # MinION Mk1B / Mk1C
    "P2I": "ONT-P2I",  # PromethION 2 Integrated
}

_DIGITS = "0123456789"


def extract_ont_device_code(position: str) -> str:
    """Return the device code embedded in a MinKNOW position token.

    ``MN47822`` -> ``MN`` (trailing serial number stripped)
    ``P2I-00513-B`` -> ``P2I`` (only the first hyphen-delimited segment is used)
    """
    return position.split("-", 1)[0].rstrip(_DIGITS)


def natural_sort_key(path: Path) -> list:
    """Sort key that orders ``reads_2.fastq.gz`` before ``reads_10.fastq.gz``."""
    return [
        int(token) if token.isdigit() else token.lower()
        for token in re.split(r"(\d+)", path.name)
    ]


class NanoporeDir:
    """Parses an Oxford Nanopore (MinKNOW) run folder name.

    Expected shape::

        <YYYYMMDD>_<HHMM>_<position>_<flow_cell_id>_<short_run_id>[_<run_alias>]

    e.g. ``20260401_1506_MN47822_FBE92725_8b1f29fd`` or
    ``20260408_1652_P2I-00513-B_PBK70557_2fe9fcef_Promethion_Training``.
    """

    def __init__(self, run_name: str):
        self.run_name = run_name
        self.folder_info = self._parse_folder()

    def _parse_folder(self) -> dict[str, str]:
        parts = self.run_name.split("_")
        if len(parts) < 5:
            raise ValueError(f"Not a nanopore run folder name: {self.run_name}")

        date, time_, position, flowcell_id, short_run_id = parts[:5]
        run_alias = "_".join(parts[5:])

        if not (len(date) == 8 and date.isdigit()):
            raise ValueError(f"Invalid date in run name: {date}")
        formatted_date = f"{date[0:4]}-{date[4:6]}-{date[6:8]}"

        if not (len(time_) == 4 and time_.isdigit()):
            raise ValueError(f"Invalid time in run name: {time_}")

        device_code = extract_ont_device_code(position)
        if device_code not in ONT_MACHINE_TYPES:
            raise ValueError(f"Unknown nanopore device in run name: {position}")
        self.machine_type = ONT_MACHINE_TYPES[device_code]

        return {
            "date": formatted_date,
            "time": time_,
            "position": position,
            "device_code": device_code,
            "flowcell_id": flowcell_id,
            "short_run_id": short_run_id,
            "run_alias": run_alias,
        }

    def build_archive_dir(self) -> str:
        # Nanopore has no lane concept; archive under the full run folder name.
        return self.run_name


class NanoporeFastq:
    """Parses the first header line of an ONT fastq to sanity check a run.

    ONT/dorado headers are ``@<read_id> key=value key=value ...`` where the
    interesting keys are ``flow_cell_id`` and (for multiplexed runs)
    ``barcode``.  Older basecallers may omit these; missing keys are treated as
    "nothing to check" rather than a failure.
    """

    def __init__(
        self,
        f: TextIOWrapper,
        expected_flowcell: str | None = None,
        expected_barcode: str | None = None,
    ):
        self.file = f
        self.expected_flowcell = expected_flowcell
        self.expected_barcode = expected_barcode
        self.info = self._parse_header()

    def _parse_header(self) -> dict[str, str]:
        try:
            line = next(self.file).strip()
        except StopIteration:
            raise ValueError("FASTQ file is empty")
        if not line.startswith("@"):
            raise ValueError("Not a FASTQ header line")

        tokens = line[1:].split()
        info = {"read_id": tokens[0] if tokens else ""}
        for token in tokens[1:]:
            if "=" in token:
                key, value = token.split("=", 1)
                info[key] = value
        return info

    def check_fp_vs_content(self) -> tuple[bool, dict[str, tuple[str, str]]]:
        problems: dict[str, tuple[str, str]] = {}

        flowcell = self.info.get("flow_cell_id")
        if self.expected_flowcell and flowcell and flowcell != self.expected_flowcell:
            problems["flow_cell_id"] = (flowcell, self.expected_flowcell)

        barcode = self.info.get("barcode")
        if self.expected_barcode and barcode and barcode != self.expected_barcode:
            problems["barcode"] = (barcode, self.expected_barcode)

        return (not problems, problems)
