#!/usr/bin/env python3

import bisect
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


K = 7
FAMILIES = [
    ("Jessica", "01_jessica.txt"),
    ("Sonetto", "02_sonetto.txt"),
    ("Regulus", "03_regulus.txt"),
    ("Regulus^T", "03_regulus_transposed.txt"),
]
COLORS = {
    "Jessica": "#D7E0E1",
    "Jessica excess minus one": "#C4CED0",
    "Target minus one": "#B4ABDA",
    "Sonetto": "#D9782B",
    "Regulus": "#D94B3A",
    "Regulus^T": "#9E2F2A",
}
ENTROPY_ANCHOR = 0.01
VERIFICATION_THRESHOLD = 1.244843
EFFECTIVE_LIPSCHITZ_CAP = 13.0


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


def binary_entropy_base2(p: float) -> float:
    if p <= 0.0 or p >= 1.0:
        return 0.0
    return -p * math.log2(p) - (1.0 - p) * math.log2(1.0 - p)


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


def build_envelope(family_filename: str) -> Envelope:
    lattice = read_lattice(Path("06_lattice_points") / family_filename)
    records = read_certificates(Path("08_dual_certificates") / family_filename)
    coeffs = [bernstein_coefficients(alpha) for _, alpha, _ in records]
    return Envelope(lattice, coeffs)


def effective_verification_lipschitz(envelopes) -> float:
    family_lipschitz = [envelopes[family_name].lipschitz for family_name, _ in FAMILIES]
    if all(value <= EFFECTIVE_LIPSCHITZ_CAP for value in family_lipschitz):
        return EFFECTIVE_LIPSCHITZ_CAP
    return max(family_lipschitz)


def verify_size_bound(envelope: Envelope, verification_lipschitz: float) -> None:
    mesh = [0.0] + envelope.lattice + [1.0]

    for left, right in zip(mesh[:-1], mesh[1:]):
        step = (right - left) / 1000.0
        for index in range(1000):
            a = left + index * step
            b = left + (index + 1) * step

            ell_a = envelope(a)
            ell_b = envelope(b)
            entropy_max = binary_entropy_base2(min(max(0.5, a), b))
            upper = ell_b + verification_lipschitz * (b - a) + entropy_max

            if upper > VERIFICATION_THRESHOLD:
                raise RuntimeError(
                    "Jessica size-bound verification failed on "
                    f"[{a:.12f}, {b:.12f}] with upper={upper:.12f}, "
                    f"ell(a)={ell_a:.12f}, ell(b)={ell_b:.12f}, "
                    f"h_max={entropy_max:.12f}, L={verification_lipschitz:.12f}"
                )


def plot_curves(envelopes) -> None:
    out_dir = Path("11_figures")
    out_dir.mkdir(parents=True, exist_ok=True)

    target_minus_one = VERIFICATION_THRESHOLD - 1.0
    target_minus_one_label = f"{target_minus_one:.6f}"

    xs = np.linspace(0.0, 1.0, 10001)
    y_jessica = np.array([envelopes["Jessica"](x) for x in xs], dtype=float)
    y_excess = np.maximum(
        -0.1,
        y_jessica + np.array([binary_entropy_base2(x) for x in xs], dtype=float) - 1.0,
    )
    y_sonetto = np.array([envelopes["Sonetto"](x) for x in xs], dtype=float)
    y_regulus = np.array([envelopes["Regulus"](x) for x in xs], dtype=float)
    y_regulus_t = np.array([envelopes["Regulus^T"](x) for x in xs], dtype=float)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(xs, y_jessica, color=COLORS["Jessica"], linewidth=2.0, label="ell_Jessica(p)")
    ax.plot(xs, y_excess, color=COLORS["Jessica excess minus one"], linewidth=2.0, label="max(-0.1, ell_Jessica(p) + h(p) - 1)")
    ax.axhline(
        target_minus_one,
        color=COLORS["Target minus one"],
        linewidth=2.0,
        label=target_minus_one_label,
    )
    ax.plot(xs, y_sonetto, color=COLORS["Sonetto"], linewidth=2.0, label="ell_Sonetto(p)")
    ax.plot(xs, y_regulus, color=COLORS["Regulus"], linewidth=2.0, label="ell_Regulus(p)")
    ax.plot(xs, y_regulus_t, color=COLORS["Regulus^T"], linewidth=2.0, label="ell_Regulus^T(p)")

    ax.set_xlabel("p")
    ax.set_ylabel("value")
    ax.set_xlim(0.0, 1.0)
    ax.grid(True, alpha=0.3)
    ax.legend()

    fig.savefig(out_dir / "01_certificate_curves.png", dpi=300)
    plt.close(fig)


def main() -> None:
    envelopes = {
        family_name: build_envelope(filename)
        for family_name, filename in FAMILIES
    }

    l_star = effective_verification_lipschitz(envelopes)

    print(f"Entropy lower-envelope slope S_E={_ENTROPY_LOWER_SLOPE:.12f}")
    for family_name, _ in FAMILIES:
        print(f"L_{family_name}={envelopes[family_name].lipschitz:.12f}")
    print(f"L_*={l_star:.12f}")
    print(f"Target Bound: {VERIFICATION_THRESHOLD:.6f}")

    verify_size_bound(envelopes["Jessica"], l_star)
    print("Jessica verification succeeded.")

    plot_curves(envelopes)


if __name__ == "__main__":
    main()
