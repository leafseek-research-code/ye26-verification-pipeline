#!/usr/bin/env python3

import bisect
import math
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import matplotlib.colors as mcolors
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np


K = 7
EPS = 1e-11
TIE_TOL = 1e-10
CERTIFY_THRESHOLD = 0.3198989999
VIOLATION_THRESHOLD = 0.319899
NUM_THREADS = 16
EFFECTIVE_LIPSCHITZ_CAP = 13.0

AX_BBOX = [0.08, 0.10, 0.80, 0.80]
COLORBAR_BBOX = [0.90, 0.10, 0.03, 0.80]
STRATEGYBAR_BBOX = [0.90, 0.10, 0.09, 0.80]
ALL_FAMILY_FILES = {
    "Jessica": "01_jessica.txt",
    "Sonetto": "02_sonetto.txt",
    "Regulus": "03_regulus.txt",
    "Regulus^T": "03_regulus_transposed.txt",
}
ENTROPY_ANCHOR = 0.01

palette_final = {
    "0+1": "#F8F4EA",
    "0+2": "#F1E8D8",
    "1+2": "#E6DAC0",
    "2": "#F4BE18",
    "5": "#E69A00",
    "2+5": "#69B2D8",
    "2+6": "#2D6C9B",
    "2+7": "#143454",
    "0+4": "#F7C8C1",
    "4+7": "#F2A08F",
    "1+7": "#E97866",
    "0+7": "#D85A4B",
    "5+7": "#B83A31",
    "1+5": "#8F241D",
    "1+3": "#E1E9BD",
    "3+6": "#C8DF8B",
    "0+6": "#A9CE59",
    "1+6": "#7FB64A",
    "5+6": "#569632",
    "0+5": "#2F7428",
}
STRATEGY_LABELS = [
    "2",
    "5",
    "0+1",
    "0+2",
    "0+4",
    "0+5",
    "0+6",
    "0+7",
    "1+2",
    "1+3",
    "1+5",
    "1+6",
    "1+7",
    "2+5",
    "2+6",
    "2+7",
    "3+6",
    "4+7",
    "5+6",
    "5+7",
]
COLORBAR_LABELS = list(palette_final)
STRATEGY_RGB = np.array([mcolors.to_rgb(palette_final[label]) for label in STRATEGY_LABELS], dtype=float)


@dataclass(frozen=True)
class Strategy:
    label: str
    indices: Tuple[int, ...]

    @property
    def is_one_strategy(self) -> bool:
        return len(self.indices) == 1


STRATEGIES: List[Strategy] = []
for label in STRATEGY_LABELS:
    parts = tuple(int(part) for part in label.split("+"))
    STRATEGIES.append(Strategy(label=label, indices=parts))


def read_lattice(path: Path):
    return [float(line.strip()) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def read_certificates(path: Path):
    lines = path.read_text(encoding="utf-8").splitlines()
    index = 0
    records = []

    while index < len(lines):
        p = float(lines[index].strip())
        index += 1

        alpha = []
        for _ in range(K + 1):
            alpha.append([float(x) for x in lines[index].split()])
            index += 1

        beta = []
        for _ in range(K + 1):
            beta.append([float(x) for x in lines[index].split()])
            index += 1

        records.append((p, np.array(alpha, dtype=float), np.array(beta, dtype=float)))

    return records


def bernstein_coefficients(alpha: np.ndarray) -> np.ndarray:
    coeffs = np.zeros(2 * K + 1, dtype=float)
    for t in range(2 * K + 1):
        total = 0.0
        for i in range(K + 1):
            for j in range(K + 1):
                if i + j == t:
                    total += math.comb(K, i) * math.comb(K, j) * alpha[i, j]
        coeffs[t] = total / math.comb(2 * K, t)
    return coeffs


def evaluate_bernstein(coeffs: np.ndarray, p: float) -> float:
    total = 0.0
    one_minus_p = 1.0 - p
    for t, value in enumerate(coeffs):
        total += value * math.comb(2 * K, t) * (p ** t) * (one_minus_p ** (2 * K - t))
    return total


def binomial_entropy_base2(p: float) -> float:
    total = 0.0
    one_minus_p = 1.0 - p
    for i in range(K + 1):
        probability = math.comb(K, i) * (p ** i) * (one_minus_p ** (K - i))
        if probability > 0.0:
            total -= probability * math.log2(probability)
    return total


_ENTROPY_ANCHOR_VALUE = binomial_entropy_base2(ENTROPY_ANCHOR)
_ENTROPY_LOWER_SLOPE = _ENTROPY_ANCHOR_VALUE / ENTROPY_ANCHOR


def entropy_lower_envelope(p: float) -> float:
    if p <= ENTROPY_ANCHOR:
        return _ENTROPY_LOWER_SLOPE * p
    if p >= 1.0 - ENTROPY_ANCHOR:
        return _ENTROPY_LOWER_SLOPE * (1.0 - p)
    return binomial_entropy_base2(p)


class Envelope:
    def __init__(self, lattice, coeffs):
        self.lattice = lattice
        self.coeffs = coeffs
        self.raw_lipschitz = max(
            14.0 * np.max(np.abs(coeff[1:] - coeff[:-1]))
            for coeff in coeffs
        )
        self.entropy_lipschitz = _ENTROPY_LOWER_SLOPE
        self.lipschitz = (self.raw_lipschitz + self.entropy_lipschitz) / K

    def raw_curve_value(self, index: int, p: float) -> float:
        return evaluate_bernstein(self.coeffs[index], p)

    def corrected_curve_value(self, index: int, p: float) -> float:
        return (self.raw_curve_value(index, p) - entropy_lower_envelope(p)) / K

    def __call__(self, p: float) -> float:
        if p <= self.lattice[0]:
            return self.corrected_curve_value(0, p)
        if p >= self.lattice[-1]:
            return self.corrected_curve_value(len(self.lattice) - 1, p)

        left = bisect.bisect_right(self.lattice, p) - 1
        return min(
            self.corrected_curve_value(left, p),
            self.corrected_curve_value(left + 1, p),
        )


def build_envelope(filename: str) -> Envelope:
    lattice = read_lattice(Path("06_lattice_points") / filename)
    records = read_certificates(Path("08_dual_certificates") / filename)
    coeffs = [bernstein_coefficients(alpha) for _, alpha, _ in records]
    return Envelope(lattice, coeffs)


def effective_verification_lipschitz(all_envelopes) -> float:
    family_lipschitz = [all_envelopes[family_name].lipschitz for family_name in ALL_FAMILY_FILES]
    if all(value <= EFFECTIVE_LIPSCHITZ_CAP for value in family_lipschitz):
        return EFFECTIVE_LIPSCHITZ_CAP
    return max(1.0, *family_lipschitz)


def candidate_points(p: float, q: float, envelopes):
    return [
        (0.0, 1.0 - q),
        (1.0 - p, 0.0),
        (p, q),
        (p * (1.0 - p), 0.5 * (1.0 - q * q)),
        (0.5 * (1.0 - p * p), q * (1.0 - q)),
        (envelopes["Sonetto"](p), envelopes["Sonetto"](q)),
        (envelopes["Regulus"](p), envelopes["Regulus^T"](q)),
        (envelopes["Regulus^T"](p), envelopes["Regulus"](q)),
    ]


def strategy_value(strategy: Strategy, points) -> float:
    if strategy.is_one_strategy:
        x, y = points[strategy.indices[0]]
        return max(x, y)

    i, j = strategy.indices
    x_i, y_i = points[i]
    x_j, y_j = points[j]
    d_i = x_i - y_i
    d_j = x_j - y_j

    if (d_i < -EPS and d_j > EPS) or (d_i > EPS and d_j < -EPS):
        t_star = d_i / (d_i - d_j)
        return (1.0 - t_star) * x_i + t_star * x_j

    return math.inf


def delta_value_and_strategy(p: float, q: float, envelopes):
    points = candidate_points(p, q, envelopes)

    best_value = math.inf
    best_index = -1

    for index, strategy in enumerate(STRATEGIES):
        candidate = strategy_value(strategy, points)
        if candidate < best_value - TIE_TOL:
            best_value = candidate
            best_index = index

    return best_value, best_index


def split_points_evenly(points, num_chunks):
    n = len(points)
    num_chunks = max(1, min(num_chunks, n))
    return [
        points[(worker * n) // num_chunks:((worker + 1) * n) // num_chunks]
        for worker in range(num_chunks)
        if (worker * n) // num_chunks < ((worker + 1) * n) // num_chunks
    ]


def evaluate_certification_chunk(points_chunk, envelopes, radius):
    certified = True
    max_upper = -math.inf
    max_info = None
    max_sample_value = -math.inf
    max_sample_point = None
    max_sample_strategy = None
    violation = None

    for p, q in points_chunk:
        value, strategy_index = delta_value_and_strategy(p, q, envelopes)

        if value > max_sample_value:
            max_sample_value = value
            max_sample_point = (p, q)
            max_sample_strategy = strategy_index

        if value > VIOLATION_THRESHOLD:
            violation = (p, q, value, strategy_index)
            break

        upper = value + radius
        if upper > max_upper:
            max_upper = upper
            max_info = (p, q, value, strategy_index)

        if upper > CERTIFY_THRESHOLD:
            certified = False

    return {
        "certified": certified,
        "max_upper": max_upper,
        "max_info": max_info,
        "max_sample_value": max_sample_value,
        "max_sample_point": max_sample_point,
        "max_sample_strategy": max_sample_strategy,
        "violation": violation,
    }


def certify_rectangle(p0, p1, q0, q1, depth, envelopes, l_star, global_state):
    global_state["recursive_calls"] += 1
    global_state["max_depth"] = max(global_state["max_depth"], depth)

    p_grid = np.linspace(p0, p1, 101)
    q_grid = np.linspace(q0, q1, 101)
    radius = l_star * ((p1 - p0) / 200.0 + (q1 - q0) / 200.0)

    points = [(float(p), float(q)) for p in p_grid for q in q_grid]
    global_state["grid_points_checked"] += len(points)

    num_workers = max(1, min(NUM_THREADS, len(points)))
    point_chunks = split_points_evenly(points, num_workers)

    executor = global_state.get("executor")
    if executor is not None and len(point_chunks) > 1:
        futures = [
            executor.submit(evaluate_certification_chunk, chunk, envelopes, radius)
            for chunk in point_chunks
        ]
        chunk_results = [future.result() for future in futures]
    else:
        chunk_results = [
            evaluate_certification_chunk(chunk, envelopes, radius)
            for chunk in point_chunks
        ]

    certified = True
    max_upper = -math.inf
    max_info = None

    for result in chunk_results:
        if result["max_sample_value"] > global_state["max_sample_value"]:
            global_state["max_sample_value"] = result["max_sample_value"]
            global_state["max_sample_point"] = result["max_sample_point"]
            global_state["max_sample_strategy"] = result["max_sample_strategy"]

        if result["violation"] is not None:
            p, q, value, strategy_index = result["violation"]
            raise RuntimeError(
                "Explicit sampled violation found in "
                f"R=[{p0:.12f}, {p1:.12f}] x [{q0:.12f}, {q1:.12f}], "
                f"point=({p:.12f}, {q:.12f}), value={value:.12f}, "
                f"strategy_index={strategy_index}, label={STRATEGY_LABELS[strategy_index]}"
            )

        if result["max_upper"] > max_upper:
            max_upper = result["max_upper"]
            max_info = result["max_info"]

        if not result["certified"]:
            certified = False

    if certified:
        return

    if depth % 2 == 0:
        mid = 0.5 * (p0 + p1)
        if mid == p0 or mid == p1:
            p, q, value, strategy_index = max_info
            raise RuntimeError(
                "Certification did not terminate on rectangle "
                f"depth={depth}, bounds=[{p0:.12f}, {p1:.12f}] x [{q0:.12f}, {q1:.12f}], "
                f"max_upper={max_upper:.12f}, sample_point=({p:.12f}, {q:.12f}), "
                f"sample_value={value:.12f}, strategy_index={strategy_index}, "
                f"label={STRATEGY_LABELS[strategy_index]}"
            )
        certify_rectangle(p0, mid, q0, q1, depth + 1, envelopes, l_star, global_state)
        certify_rectangle(mid, p1, q0, q1, depth + 1, envelopes, l_star, global_state)
    else:
        mid = 0.5 * (q0 + q1)
        if mid == q0 or mid == q1:
            p, q, value, strategy_index = max_info
            raise RuntimeError(
                "Certification did not terminate on rectangle "
                f"depth={depth}, bounds=[{p0:.12f}, {p1:.12f}] x [{q0:.12f}, {q1:.12f}], "
                f"max_upper={max_upper:.12f}, sample_point=({p:.12f}, {q:.12f}), "
                f"sample_value={value:.12f}, strategy_index={strategy_index}, "
                f"label={STRATEGY_LABELS[strategy_index]}"
            )
        certify_rectangle(p0, p1, q0, mid, depth + 1, envelopes, l_star, global_state)
        certify_rectangle(p0, p1, mid, q1, depth + 1, envelopes, l_star, global_state)


def compute_landscape(p_values, q_values, envelopes):
    best = np.empty((len(p_values), len(q_values)), dtype=np.float32)
    idx = np.empty((len(p_values), len(q_values)), dtype=np.int16)

    for i, p in enumerate(p_values):
        for j, q in enumerate(q_values):
            value, strategy_index = delta_value_and_strategy(float(p), float(q), envelopes)
            best[i, j] = value
            idx[i, j] = strategy_index
        if (i + 1) % 50 == 0 or i + 1 == len(p_values):
            print(f"heatmap row {i + 1}/{len(p_values)}")

    return best, idx


def make_value_cmap():
    return mcolors.LinearSegmentedColormap.from_list(
        "value_landscape",
        ["#03051a", "#0b3d91", "#006d9c", "#24b3b3", "#9fffcf"],
    )


def plot_value_heatmap(array, extent, title, ticks, output_path):
    cmap = make_value_cmap()
    norm = mcolors.Normalize(vmin=0.31, vmax=0.32)

    fig = plt.figure(figsize=(10, 10), dpi=600)
    ax = fig.add_axes(AX_BBOX)
    image = ax.imshow(
        array.T,
        origin="lower",
        extent=extent,
        interpolation="nearest",
        cmap=cmap,
        norm=norm,
        aspect="equal",
    )

    ax.set_title(title)
    ax.set_xlabel("p")
    ax.set_ylabel("q")
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)

    cax = fig.add_axes(COLORBAR_BBOX)
    colorbar = fig.colorbar(image, cax=cax, ticks=[0.31 + 0.001 * k for k in range(11)])
    colorbar.set_label("δ (cost)")

    fig.savefig(output_path)
    plt.close(fig)


def draw_strategy_colorbar(fig):
    ax = fig.add_axes(STRATEGYBAR_BBOX)
    ax.axis("off")
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)

    rows = len(COLORBAR_LABELS)
    margin_y = 0.02
    usable_height = 1.0 - 2.0 * margin_y
    row_height = usable_height / rows
    swatch_height = row_height * 0.75

    for index, label in enumerate(COLORBAR_LABELS):
        y_center = 1.0 - margin_y - (index + 0.5) * row_height
        y0 = y_center - 0.5 * swatch_height
        ax.add_patch(
            patches.Rectangle(
                (0.04, y0),
                0.28,
                swatch_height,
                facecolor=palette_final[label],
                edgecolor="black",
                linewidth=0.5,
            )
        )
        ax.text(0.38, y_center, label, va="center", ha="left", fontsize=9)


def plot_strategy_heatmap(idx_array, extent, title, ticks, output_path):
    rgb = STRATEGY_RGB[idx_array.T]

    fig = plt.figure(figsize=(10, 10), dpi=600)
    ax = fig.add_axes(AX_BBOX)
    ax.imshow(
        rgb,
        origin="lower",
        extent=extent,
        interpolation="nearest",
        aspect="equal",
    )
    ax.set_title(title)
    ax.set_xlabel("p")
    ax.set_ylabel("q")
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)

    draw_strategy_colorbar(fig)

    fig.savefig(output_path)
    plt.close(fig)


def main() -> None:
    start_time = time.perf_counter()
    print("Starting Verification Algorithm V...")
    print(f"Target Bound: {CERTIFY_THRESHOLD:.10f}")
    print(f"Using Threads: {NUM_THREADS}")

    all_envelopes = {
        family_name: build_envelope(filename)
        for family_name, filename in ALL_FAMILY_FILES.items()
    }
    envelopes = {name: all_envelopes[name] for name in ["Sonetto", "Regulus", "Regulus^T"]}
    print("Certificates loaded successfully.")
    print("------------------------------------------------------------")

    print(f"Entropy lower-envelope slope S_E={_ENTROPY_LOWER_SLOPE:.12f}")
    for family_name in ALL_FAMILY_FILES:
        print(f"L_{family_name}={all_envelopes[family_name].lipschitz:.12f}")
    l_star = effective_verification_lipschitz(all_envelopes)
    print(f"L_*={l_star:.12f}")

    global_state = {
        "max_sample_value": -math.inf,
        "max_sample_point": None,
        "max_sample_strategy": None,
        "executor": None,
        "recursive_calls": 0,
        "grid_points_checked": 0,
        "max_depth": 0,
    }

    if NUM_THREADS > 1:
        with ThreadPoolExecutor(max_workers=NUM_THREADS) as executor:
            global_state["executor"] = executor
            certify_rectangle(0.0, 1.0, 0.0, 1.0, 0, envelopes, l_star, global_state)
    else:
        certify_rectangle(0.0, 1.0, 0.0, 1.0, 0, envelopes, l_star, global_state)

    point = global_state["max_sample_point"]
    strategy_index = global_state["max_sample_strategy"]

    elapsed = time.perf_counter() - start_time
    print("------------------------------------------------------------")
    print("Verification Complete.")
    print(f"Time Elapsed: {elapsed:.2f} seconds.")
    print(f"Total Recursive Calls: {global_state['recursive_calls']}")
    print(f"Total Grid Points Checked: {global_state['grid_points_checked']}")
    print(f"Max Recursion Depth: {global_state['max_depth']}")
    print(f"Max Delta_C Found: {global_state['max_sample_value']:.14f}")
    print(f"Location of Max: (p={point[0]:.14f}, q={point[1]:.14f})")
    print(
        f"Winning Strategy: index={strategy_index}, label={STRATEGY_LABELS[strategy_index]}"
    )
    print("RESULT: CERTIFIED (Max Found <= Target)")

    out_dir = Path("11_figures")
    out_dir.mkdir(parents=True, exist_ok=True)

    n = 1600
    full_p = (np.arange(n, dtype=float) + 0.5) / 1600.0
    full_q = (np.arange(n, dtype=float) + 0.5) / 1600.0
    zoom_p = 0.25 + (np.arange(n, dtype=float) + 0.5) * 0.15 / 1600.0
    zoom_q = 0.25 + (np.arange(n, dtype=float) + 0.5) * 0.15 / 1600.0

    best_full, idx_full = compute_landscape(full_p, full_q, envelopes)
    best_zoom, idx_zoom = compute_landscape(zoom_p, zoom_q, envelopes)

    plot_value_heatmap(
        best_full,
        extent=[0.0, 1.0, 0.0, 1.0],
        title="Cost landscape: δ(p,q) on [0,1] × [0,1]",
        ticks=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
        output_path=out_dir / "02_value_heatmap_full.png",
    )
    plot_value_heatmap(
        best_zoom,
        extent=[0.25, 0.4, 0.25, 0.4],
        title="Cost landscape: δ(p,q) on [0.25,0.4] × [0.25,0.4]",
        ticks=[0.25, 0.28, 0.31, 0.34, 0.37, 0.40],
        output_path=out_dir / "03_value_heatmap_zoom.png",
    )

    plot_strategy_heatmap(
        idx_full,
        extent=[0.0, 1.0, 0.0, 1.0],
        title="Strategy landscape: best strategy on [0,1] × [0,1]",
        ticks=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
        output_path=out_dir / "04_strategy_heatmap_full.png",
    )
    plot_strategy_heatmap(
        idx_zoom,
        extent=[0.25, 0.4, 0.25, 0.4],
        title="Strategy landscape: best strategy on [0.25,0.4] × [0.25,0.4]",
        ticks=[0.25, 0.28, 0.31, 0.34, 0.37, 0.40],
        output_path=out_dir / "05_strategy_heatmap_zoom.png",
    )


if __name__ == "__main__":
    main()
