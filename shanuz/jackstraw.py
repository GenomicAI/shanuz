from __future__ import annotations

from typing import Optional

import numpy as np


class JackStrawData:
    """Stores JackStraw permutation test results for a DimReduc.

    Mirrors R's JackStraw / JackStrawData from jackstraw.R.
    """

    __slots__ = (
        "empirical_p_values",
        "fake_reduction_scores",
        "overall_p_values",
        "score",
        "method",
    )

    def __init__(
        self,
        empirical_p_values: Optional[np.ndarray] = None,
        fake_reduction_scores: Optional[np.ndarray] = None,
        overall_p_values: Optional[np.ndarray] = None,
        score: Optional[np.ndarray] = None,
        method: Optional[str] = None,
    ) -> None:
        self.empirical_p_values = empirical_p_values
        self.fake_reduction_scores = fake_reduction_scores
        self.overall_p_values = overall_p_values
        self.score = score
        self.method = method

    def is_empty(self) -> bool:
        return self.empirical_p_values is None

    def __repr__(self) -> str:
        if self.is_empty():
            return "JackStrawData(empty)"
        shape = self.empirical_p_values.shape if self.empirical_p_values is not None else "?"
        return f"JackStrawData(p_values shape={shape}, method={self.method!r})"
