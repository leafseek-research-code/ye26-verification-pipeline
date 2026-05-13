#!/usr/bin/env python3

from pathlib import Path

PARAMETER_FILES = {
    "01_jessica.txt": "\n".join([
        "7",
        "612",
        "6",
        "0 20 42 85 107 127",
        "14 22 270 270 22 14",
        "1.2599210498948732",
        "0.5 0.5",
    ]),
    "02_sonetto.txt": "\n".join([
        "7",
        "360",
        "6",
        "0 21 42 85 106 127",
        "8 22 150 150 22 8",
        "1.2599210498948732",
        "0.38 0.38",
    ]),
    "03_regulus.txt": "\n".join([
        "7",
        "360",
        "8",
        "0 42 43 45 77 85 86 127",
        "8 36 20 20 160 64 44 8",
        "1.2599210498948732",
        "0.336 0.414",
    ]),
    "03_regulus_transposed.txt": "\n".join([
        "7",
        "360",
        "8",
        "0 41 42 50 82 84 85 127",
        "8 44 64 160 20 20 36 8",
        "1.2599210498948732",
        "0.414 0.336",
    ]),
}


def main() -> None:
    out_dir = Path("02_circuit_parameters")
    out_dir.mkdir(parents=True, exist_ok=True)

    for filename, content in PARAMETER_FILES.items():
        (out_dir / filename).write_text(content, encoding="utf-8")


if __name__ == "__main__":
    main()

