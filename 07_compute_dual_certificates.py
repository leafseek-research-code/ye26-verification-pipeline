#!/usr/bin/env python3

import math
import sys
from pathlib import Path

import cvxpy as cp
import numpy as np


def read_lattice(path: Path):
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    return [(line, float(line)) for line in lines]


def read_transition_matrices(path: Path):
    lines = path.read_text(encoding="utf-8").splitlines()
    K, H, _ = map(int, lines[0].split())
    matrices = [[[] for _ in range(H)] for _ in range(K + 1)]

    for line in lines[1:]:
        i, u, v, value = map(int, line.split())
        matrices[i][u].append((v, value))

    return K, H, matrices


def build_problem(K: int, H: int, matrices):
    alpha = cp.Variable((K + 1, K + 1))
    beta = cp.Variable((K + 1, H))
    weight_matrix = cp.Parameter((K + 1, K + 1), nonneg=True)

    constraints = []
    for i in range(K + 1):
        for u in range(H):
            terms = []
            for j in range(K + 1):
                for v, value in matrices[j][u]:
                    terms.append(math.log(value) - alpha[i, j] - beta[i, u] + beta[j, v])
            if terms:
                constraints.append(cp.log_sum_exp(cp.hstack(terms)) <= 0.0)

    middle_left = H // 2 - 1
    middle_right = H // 2
    for i in range(K + 1):
        constraints.append(beta[i, middle_left] + beta[i, middle_right] == 0.0)

    objective = cp.Minimize(cp.sum(cp.multiply(weight_matrix, alpha)))
    problem = cp.Problem(objective, constraints)
    return problem, alpha, beta, weight_matrix


def verify_constraints(alpha_base2: np.ndarray, beta_base2: np.ndarray, matrices):
    K = alpha_base2.shape[0] - 1
    H = beta_base2.shape[1]
    largest_sum = 0.0

    for i in range(K + 1):
        for u in range(H):
            total = 0.0
            for j in range(K + 1):
                for v, value in matrices[j][u]:
                    total += value * (2.0 ** (-alpha_base2[i, j] - beta_base2[i, u] + beta_base2[j, v]))
            if total > largest_sum:
                largest_sum = total

    return largest_sum


def write_certificate_file(path: Path, p_texts, alpha_blocks, beta_blocks) -> None:
    lines = []
    for block_index, p_text in enumerate(p_texts):
        lines.append(p_text)
        alpha = alpha_blocks[block_index]
        beta = beta_blocks[block_index]

        for row in alpha:
            lines.append(" ".join(f"{value:14.10f}" for value in row))
        for row in beta:
            lines.append(" ".join(f"{value:14.10f}" for value in row))

    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: python 07_compute_dual_certificates.py <family_name_without_extension>")

    family = sys.argv[1]
    lattice_path = Path("06_lattice_points") / f"{family}.txt"
    matrix_path = Path("04_transition_matrices") / f"{family}.txt"
    out_path = Path("08_dual_certificates") / f"{family}.txt"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    lattice = read_lattice(lattice_path)
    K, H, matrices = read_transition_matrices(matrix_path)
    problem, alpha_var, beta_var, weight_matrix = build_problem(K, H, matrices)

    ln2 = math.log(2.0)
    p_texts = []
    alpha_blocks = []
    beta_blocks = []

    for p_text, p in lattice:
        binomial = np.array(
            [math.comb(K, i) * (p ** i) * ((1.0 - p) ** (K - i)) for i in range(K + 1)],
            dtype=float,
        )
        weight_matrix.value = np.outer(binomial, binomial)

        problem.solve(solver=cp.MOSEK)
        if problem.status not in {cp.OPTIMAL, cp.OPTIMAL_INACCURATE}:
            raise RuntimeError(f"MOSEK failed for p={p_text}: status={problem.status}")

        alpha_base2 = np.asarray(alpha_var.value, dtype=float) / ln2
        beta_base2 = np.asarray(beta_var.value, dtype=float) / ln2

        s0 = verify_constraints(alpha_base2, beta_base2, matrices)
        shift = math.log2(s0) + 1e-8
        alpha_base2 = alpha_base2 + shift

        p_texts.append(p_text)
        alpha_blocks.append(alpha_base2)
        beta_blocks.append(beta_base2)

        print(
            f"[{family}] p={p_text} objective={problem.value / ln2:.10f} "
            f"S0={s0:.12f} alpha_shift={shift:.12e}"
        )

    write_certificate_file(out_path, p_texts, alpha_blocks, beta_blocks)


if __name__ == "__main__":
    main()
