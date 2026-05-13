#!/usr/bin/env python3

import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List


@dataclass
class Parameters:
    K: int
    H: int
    block_ids: List[int]
    multiplicities: List[int]
    Z: float
    alpha: float
    beta: float

    def block_of_state(self, h: int) -> int:
        total = 0
        for block_id, multiplicity in zip(self.block_ids, self.multiplicities):
            total += multiplicity
            if h < total:
                return block_id
        raise ValueError(f"state {h} is outside the encoded partition")


@dataclass
class Term:
    row_support: int
    col_support: int
    delta: int
    symbol: str
    step_index: int
    raw_value: float
    half_distance: float
    rounded_target: float


def read_parameters(path: Path) -> Parameters:
    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines) != 7:
        raise ValueError(f"{path} should contain exactly 7 lines")

    K = int(lines[0])
    H = int(lines[1])
    L = int(lines[2])
    block_ids = [int(x) for x in lines[3].split()]
    multiplicities = [int(x) for x in lines[4].split()]
    alpha, beta = (float(x) for x in lines[6].split())

    if len(block_ids) != L or len(multiplicities) != L:
        raise ValueError(f"{path} has inconsistent block data: expected {L} entries")
    if sum(multiplicities) != H:
        raise ValueError(f"{path} has multiplicities summing to {sum(multiplicities)} instead of H={H}")

    return Parameters(
        K=K,
        H=H,
        block_ids=block_ids,
        multiplicities=multiplicities,
        Z=float(lines[5]),
        alpha=alpha,
        beta=beta,
    )


def support_weight(bitset: int, subset_weights: List[float]) -> float:
    total = 0.0
    for subset in range(len(subset_weights)):
        if (bitset >> subset) & 1:
            total += subset_weights[subset]
    return total


def round_nearest_integer(x: float) -> int:
    if x >= 0.0:
        return int(math.floor(x + 0.5))
    return int(math.ceil(x - 0.5))


def block_word(block_id: int, K: int) -> str:
    bits = format(block_id, f"0{K}b")
    return "".join("R" if bit == "0" else "C" for bit in bits)


def build_sergeev_terms(
    K: int,
    block_id: int,
    Z: float,
    alpha: float,
    beta: float,
) -> List[Term]:
    n = 1 << K
    subsets_by_size: List[List[int]] = [[] for _ in range(K + 1)]
    row_weights = [0.0] * n
    col_weights = [0.0] * n

    for mask in range(n):
        weight = mask.bit_count()
        subsets_by_size[weight].append(mask)
        row_weights[mask] = (1.0 - alpha) ** (K - weight) * alpha ** weight
        col_weights[mask] = (1.0 - beta) ** (K - weight) * beta ** weight

    terms: List[Term] = []
    word = block_word(block_id, K) + "R"
    current_r = 0
    current_c = 0

    for step_index, symbol in enumerate(word):
        if symbol == "R":
            for subset in subsets_by_size[current_r]:
                row_support = 1 << subset
                col_support = 0
                for other in range(n):
                    if (other & subset) == 0 and other.bit_count() >= current_c:
                        col_support |= 1 << other

                numerator = support_weight(col_support, col_weights)
                denominator = support_weight(row_support, row_weights)
                raw_value = math.log(numerator / denominator, Z)
                rounded = round_nearest_integer(raw_value)
                half_target = math.floor(raw_value) + 0.5
                half_distance = abs(raw_value - half_target)

                terms.append(
                    Term(
                        row_support=row_support,
                        col_support=col_support,
                        delta=rounded,
                        symbol=symbol,
                        step_index=step_index,
                        raw_value=raw_value,
                        half_distance=half_distance,
                        rounded_target=half_target,
                    )
                )

            current_r += 1
        else:
            for subset in subsets_by_size[current_c]:
                col_support = 1 << subset
                row_support = 0
                for other in range(n):
                    if (other & subset) == 0 and other.bit_count() >= current_r:
                        row_support |= 1 << other

                numerator = support_weight(col_support, col_weights)
                denominator = support_weight(row_support, row_weights)
                raw_value = math.log(numerator / denominator, Z)
                rounded = round_nearest_integer(raw_value)
                half_target = math.floor(raw_value) + 0.5
                half_distance = abs(raw_value - half_target)

                terms.append(
                    Term(
                        row_support=row_support,
                        col_support=col_support,
                        delta=rounded,
                        symbol=symbol,
                        step_index=step_index,
                        raw_value=raw_value,
                        half_distance=half_distance,
                        rounded_target=half_target,
                    )
                )

            current_c += 1

    return terms


def clamp_h(value: int, H: int) -> int:
    if value < 0:
        return 0
    if value >= H:
        return H - 1
    return value


def process_file(path: Path) -> None:
    params = read_parameters(path)
    out_dir = Path("04_transition_matrices")
    out_dir.mkdir(parents=True, exist_ok=True)

    block_cache: Dict[int, List[Term]] = {}
    worst_info = None

    for block_id in params.block_ids:
        terms = build_sergeev_terms(
            K=params.K,
            block_id=block_id,
            Z=params.Z,
            alpha=params.alpha,
            beta=params.beta,
        )
        block_cache[block_id] = terms
        for term_index, term in enumerate(terms):
            record = (
                term.half_distance,
                block_id,
                term_index,
                term.symbol,
                term.step_index,
                term.raw_value,
                term.rounded_target,
                term.delta,
            )
            if worst_info is None or record[0] < worst_info[0]:
                worst_info = record

    matrices = [defaultdict(int) for _ in range(params.K + 1)]

    canonical_rows = [(1 << i) - 1 for i in range(params.K + 1)]
    for h_old in range(params.H):
        block_id = params.block_of_state(h_old)
        terms = block_cache[block_id]
        for i, canonical_row in enumerate(canonical_rows):
            for term in terms:
                if (term.row_support >> canonical_row) & 1:
                    h_new = clamp_h(h_old + term.delta, params.H)
                    matrices[i][(h_old, h_new)] += 1

    entries = []
    for i in range(params.K + 1):
        for (h_old, h_new), value in sorted(matrices[i].items()):
            entries.append((i, h_old, h_new, value))

    output_lines = [f"{params.K} {params.H} {len(entries)}"]
    output_lines.extend(f"{i} {h_old} {h_new} {value}" for i, h_old, h_new, value in entries)

    out_path = out_dir / path.name
    out_path.write_text("\n".join(output_lines), encoding="utf-8")

    if worst_info is not None:
        distance, block_id, term_index, symbol, step_index, raw_value, half_target, rounded = worst_info
        print(
            f"[{path.name}] worst rounding: block={block_id}, term={term_index}, "
            f"symbol={symbol}, step={step_index}, raw={raw_value:.12f}, "
            f"nearest_half={half_target:.12f}, distance={distance:.12e}, rounded={rounded}"
        )


def main() -> None:
    in_dir = Path("02_circuit_parameters")
    for path in sorted(in_dir.glob("*.txt")):
        process_file(path)


if __name__ == "__main__":
    main()

