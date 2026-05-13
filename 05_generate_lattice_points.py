#!/usr/bin/env python3

from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

SCALE = 1_000_000

# Each triple is (start, stop, step) in millionths, inclusive at both ends.
LATTICE_SEGMENTS: Dict[str, Sequence[Tuple[int, int, int]]] = {
    "01_jessica.txt": (
        (100_000, 300_000, 100_000),
        (310_000, 340_000, 10_000),
        (341_000, 360_000, 1_000),
        (360_100, 370_800, 100),
        (370_830, 372_000, 30),
        (372_100, 380_000, 100),
        (381_000, 400_000, 1_000),
        (410_000, 500_000, 10_000),
        (600_000, 900_000, 100_000),
    ),
    "02_sonetto.txt": (
        (10_000, 200_000, 10_000),
        (201_000, 500_000, 1_000),
        (510_000, 990_000, 10_000),
    ),
    "03_regulus.txt": (
        (10_000, 200_000, 10_000),
        (201_000, 320_000, 1_000),
        (330_000, 990_000, 10_000),
    ),
    "03_regulus_transposed.txt": (
        (10_000, 320_000, 10_000),
        (321_000, 440_000, 1_000),
        (450_000, 990_000, 10_000),
    ),
}


def format_point(millionths: int) -> str:
    return f"{millionths / SCALE:.6f}".rstrip("0").rstrip(".")


def iter_segment(start: int, stop: int, step: int) -> Iterable[int]:
    return range(start, stop + 1, step)


def build_lattice(filename: str) -> List[str]:
    try:
        segments = LATTICE_SEGMENTS[filename]
    except KeyError as exc:
        raise ValueError(f"unknown lattice family: {filename}") from exc

    scaled_points: List[int] = []
    for start, stop, step in segments:
        scaled_points.extend(iter_segment(start, stop, step))

    for left, right in zip(scaled_points, scaled_points[1:]):
        if right <= left:
            raise ValueError(
                f"lattice for {filename} is not strictly increasing: {left} then {right}"
            )

    return [format_point(value) for value in scaled_points]


def main() -> None:
    out_dir = Path("06_lattice_points")
    out_dir.mkdir(parents=True, exist_ok=True)

    for filename in LATTICE_SEGMENTS:
        (out_dir / filename).write_text("\n".join(build_lattice(filename)), encoding="utf-8")


if __name__ == "__main__":
    main()
