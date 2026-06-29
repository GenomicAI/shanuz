from __future__ import annotations

from typing import TYPE_CHECKING, Optional, Union

import numpy as np

if TYPE_CHECKING:
    from scipy.sparse import spmatrix

# Mirrors R's AnyMatrix (dgCMatrix | matrix | any other matrix type)
AnyMatrix = Union[np.ndarray, "spmatrix"]

OptionalStr = Optional[str]
OptionalDict = Optional[dict]
OptionalList = Optional[list]
