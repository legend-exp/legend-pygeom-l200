"""Efficiency for the KETEK SiPMs. This is not part of legendoptics.

.. [Wiest2011] F. Wiest. “SiPM Developments at KETEK.” (Feb. 17, 2011).
    https://indico.cern.ch/event/117424/contributions/1329246/attachments/56776/81752/CERN_SiPM-Status-Ketek_17-Feb-2011.pdf
"""

from __future__ import annotations

import logging

import numpy as np
import pint
from legendoptics import store
from pint import Quantity

log = logging.getLogger(__name__)
u = pint.get_application_registry()


@store.register_pluggable
def ketek_sipm_efficiency() -> tuple[Quantity, Quantity]:
    """Detection efficiency for the KETEK SiPM.

    [Wiest2011]_
    """
    λ = np.array([100, 280, 310, 350, 400, 435, 505, 525, 595, 670][::-1]) * u.nm
    eff = np.array([0.0, 0.19, 0.30, 0.32, 0.33, 0.32, 0.27, 0.19, 0.12, 0.07][::-1])
    return λ, eff
