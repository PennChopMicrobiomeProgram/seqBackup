import csv
import re
import warnings
from collections.abc import Mapping
from typing import Optional
from io import TextIOWrapper
from pathlib import Path
from urllib.error import URLError
from urllib.request import urlopen


MACHINE_TYPES_FALLBACK = {
    "VH": "Illumina-NextSeq",
    "D": "Illumina-HiSeq",
    "M": "Illumina-MiSeq",
    "A": "Illumina-NovaSeq",
    "NB": "Illumina-MiniSeq",
    "LH": "Illumina-NovaSeqX",
    "SH": "Illumina-MiSeq",
}  # Fallback mapping if machine_types.tsv is unavailable.
MACHINE_TYPES_URL = (
    "https://raw.githubusercontent.com/PennChopMicrobiomeProgram/"
    "SampleRegistry/master/sample_registry/data/machine_types.tsv"
)
_machine_types_cache: Optional[dict[str, str]] = None


class MachineTypesMapping(Mapping):
    def _mapping(self) -> dict[str, str]:
        return load_machine_types()

    def __getitem__(self, key: str) -> str:
        return self._mapping()[key]

    def __iter__(self):
        return iter(self._mapping())

    def __len__(self) -> int:
        return len(self._mapping())


def load_machine_types() -> dict[str, str]:
    global _machine_types_cache
    if _machine_types_cache is not None:
        return _machine_types_cache

    try:
        with urlopen(MACHINE_TYPES_URL, timeout=10) as response:
            content = response.read().decode("utf-8")
    except (URLError, TimeoutError) as exc:
        warnings.warn(
            f"Falling back to bundled machine types; unable to load {MACHINE_TYPES_URL}: {exc}",
            RuntimeWarning,
        )
        _machine_types_cache = MACHINE_TYPES_FALLBACK
        return _machine_types_cache

    reader = csv.reader(content.splitlines(), delimiter="\t")
    rows = [row for row in reader if row]
    if not rows:
        warnings.warn(
            "Falling back to bundled machine types; received empty machine_types.tsv.",
            RuntimeWarning,
        )
        _machine_types_cache = MACHINE_TYPES_FALLBACK
        return _machine_types_cache

    if rows[0][0].lower() in {"instrument_code", "code"}:
        rows = rows[1:]

    mapping: dict[str, str] = {}
    for row in rows:
        if len(row) < 2:
            continue
        code = row[0].strip()
        machine_type = row[1].strip()
        if code and machine_type:
            mapping[code] = machine_type

    if not mapping:
        warnings.warn(
            "Falling back to bundled machine types; no valid rows in machine_types.tsv.",
            RuntimeWarning,
        )
        _machine_types_cache = MACHINE_TYPES_FALLBACK
        return _machine_types_cache

    _machine_types_cache = mapping
    return _machine_types_cache


MACHINE_TYPES = MachineTypesMapping()


def extract_instrument_code(instrument: str) -> str:
    return "".join(filter(lambda x: not x.isdigit(), instrument))


class IlluminaDir:
    def __init__(self, run_name: str):
        self.run_name = run_name
        self.folder_info = self._parse_folder()

    def _parse_folder(self) -> dict[str, str]:
        # Extract directory name info
        parts = self.run_name.split("_")

        date = parts[0]
        if len(date) == 8:
            date = f"{date[0:4]}-{date[4:6]}-{date[6:8]}"
        elif len(date) == 6:
            date = f"20{date[0:2]}-{date[2:4]}-{date[4:6]}"
        else:
            raise ValueError(f"Invalid date format in run name: {date}")

        instrument = parts[1]
        instrument_code = extract_instrument_code(instrument)
        if instrument_code not in MACHINE_TYPES:
            raise ValueError(f"Invalid instrument code in run name: {instrument}")
        self.machine_type = MACHINE_TYPES[instrument_code]

        run_number = parts[2]
        if not run_number.isdigit():
            raise ValueError(f"Invalid run number in run name: {run_number}")

        flowcell_id = parts[3]

        if len(parts) > 4:
            raise ValueError(f"Unexpected extra parts in run name: {parts[4:]}")

        vals1 = {
            "date": date,
            "instrument": instrument,
            "run_number": str(int(run_number)),
            "flowcell_id": flowcell_id,
        }

        if instrument_code in {"D", "A", "NB", "LH", "SH"}:
            vals1["flowcell_id"] = vals1["flowcell_id"][1:]

        return vals1


class IlluminaFastq:
    def __init__(self, f: TextIOWrapper):
        self.file = f
        self.fastq_info = self._parse_header()
        self.folder_info = IlluminaDir(self.run_name).folder_info
        self.folder_info.update(self._parse_fastq_file())

    def __str__(self):
        return "_".join(
            [
                self.fastq_info["instrument"],
                self.fastq_info["run_number"],
                self.fastq_info["flowcell_id"],
                self.fastq_info["lane"],
            ]
        )

    def is_same_run(self, other: "IlluminaFastq") -> bool:
        keys = ["run_number", "instrument", "flowcell_id"]
        return all(self.fastq_info[k] == other.fastq_info[k] for k in keys)

    def _parse_header(self) -> dict[str, str]:
        line = next(self.file).strip()
        if not line.startswith("@"):
            raise ValueError("Not a FASTQ header line")
        # Remove first character, @
        line = line[1:]
        word1, _, word2 = line.partition(" ")

        keys1 = ("instrument", "run_number", "flowcell_id", "lane")
        vals1 = dict((k, v) for k, v in zip(keys1, word1.split(":")))

        keys2 = ("read", "is_filtered", "control_number", "index_reads")
        vals2 = dict((k, v) for k, v in zip(keys2, word2.split(":")))

        vals1.update(vals2)
        return vals1

    def _parse_fastq_file(self) -> dict[str, str]:
        # Extract file name info
        filename = self.filepath.name
        if matches := re.match(
            "Undetermined_S0_L00([1-8])_([RI])([12])_001.fastq.gz", filename
        ):
            keys2 = ("lane", "read_or_index", "read")
            return dict(zip(keys2, matches.groups()))

        if matches := re.match("Undetermined_S0_([RI])([12])_001.fastq.gz", filename):
            return {
                "lane": self.fastq_info["lane"],
                "read_or_index": matches.group(1),
                "read": matches.group(2),
            }

        raise ValueError(f"Unexpected FASTQ file name: {filename}")

    @property
    def lane(self) -> str:
        return self.fastq_info["lane"]

    @property
    def filepath(self) -> Path:
        return Path(self.file.name)

    @property
    def machine_type(self) -> str:
        return MACHINE_TYPES[extract_instrument_code(self.fastq_info["instrument"])]

    @property
    def run_name(self) -> str:
        for part in self.filepath.parts:
            segments = part.split("_")
            if (
                len(segments) >= 4
                and segments[0].isdigit()
                and extract_instrument_code(segments[1]) in MACHINE_TYPES
                and segments[2].isdigit()
            ):
                return part
        raise ValueError(f"Run name not found in path: {self.filepath}")

    def build_archive_dir(self) -> str:
        return "_".join([self.run_name, "L{:0>3}".format(self.lane)])

    def check_fp_vs_content(self, verbose: bool = False) -> list[bool]:
        run_check = self.fastq_info["run_number"] == self.folder_info["run_number"]
        instrument_check = (
            self.fastq_info["instrument"] == self.folder_info["instrument"]
        )
        flowcell_check = (
            self.fastq_info["flowcell_id"] == self.folder_info["flowcell_id"]
        )
        lane_check = self.lane == self.folder_info["lane"]
        read_check = self.fastq_info["read"] == self.folder_info["read"]

        if verbose:
            (
                print(
                    "Fastq run number: ",
                    self.fastq_info["run_number"],
                    "Folder run number: ",
                    self.folder_info["run_number"],
                )
                if not run_check
                else None
            )
            (
                print(
                    "Fastq instrument: ",
                    self.fastq_info["instrument"],
                    "Folder instrument: ",
                    self.folder_info["instrument"],
                )
                if not instrument_check
                else None
            )
            (
                print(
                    "Fastq flowcell id: ",
                    self.fastq_info["flowcell_id"],
                    "Folder flowcell id: ",
                    self.folder_info["flowcell_id"],
                )
                if not flowcell_check
                else None
            )
            (
                print(
                    "Fastq lane: ", self.lane, "Folder lane: ", self.folder_info["lane"]
                )
                if not lane_check
                else None
            )
            (
                print(
                    "Fastq read: ",
                    self.fastq_info["read"],
                    "Folder read: ",
                    self.folder_info["read"],
                )
                if not read_check
                else None
            )

        return [
            run_check
            and instrument_check
            and flowcell_check
            and lane_check
            and read_check,
            run_check,
            instrument_check,
            flowcell_check,
            lane_check,
            read_check,
            self.fastq_info["flowcell_id"],
            self.folder_info["flowcell_id"],
        ]

    def check_file_size(self, min_file_size) -> bool:
        return self.filepath.stat().st_size > min_file_size

    def check_index_read_exists(self) -> bool:
        return len(self.fastq_info["index_reads"]) > 2
