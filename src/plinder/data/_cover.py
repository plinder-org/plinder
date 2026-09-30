"""Compiled linear-time operations for disk-backed directed covers."""

import numpy as np
from numba import njit
from numpy.typing import NDArray


@njit(nogil=True)
def outgoing_counts(queries: NDArray[np.uint32], node_count: int) -> NDArray[np.int64]:
    counts = np.zeros(node_count, dtype=np.int64)
    for query in queries:
        counts[query] += 1
    return counts


@njit(nogil=True)
def transpose_edges(
    offsets: NDArray[np.int64],
    queries: NDArray[np.uint32],
    scores: NDArray[np.uint8],
    outgoing_offsets: NDArray[np.int64],
    targets: NDArray[np.uint32],
    outgoing_scores: NDArray[np.uint8],
) -> None:
    positions = outgoing_offsets[:-1].copy()
    for target in range(len(offsets) - 1):
        for edge in range(offsets[target], offsets[target + 1]):
            query = queries[edge]
            position = positions[query]
            targets[position] = target
            outgoing_scores[position] = scores[edge]
            positions[query] += 1


@njit(nogil=True)
def decrement_gains(
    covered: NDArray[np.int64],
    offsets: NDArray[np.int64],
    targets: NDArray[np.uint32],
    scores: NDArray[np.uint8],
    primary: int,
    fallback: int,
    primary_gains: NDArray[np.int64],
    fallback_gains: NDArray[np.int64],
) -> None:
    for query in covered:
        primary_gains[query] -= 1
        fallback_gains[query] -= 1
        for edge in range(offsets[query], offsets[query + 1]):
            target = targets[edge]
            if target != query:
                score = scores[edge]
                if score >= primary:
                    primary_gains[target] -= 1
                if score >= fallback:
                    fallback_gains[target] -= 1
