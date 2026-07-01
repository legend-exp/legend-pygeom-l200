from __future__ import annotations

import matplotlib as mpl

mpl.use("Agg")

import numpy as np


def test_public_hpge_mass_comparison():
    # the public geometry must yield a non-empty, finite comparison without a metadata git checkout.
    from pygeoml200 import plot_hpge_mass_comparison

    ax = plot_hpge_mass_comparison({"public_geom": True})

    finite = []
    for line in ax.get_lines():
        y = np.asarray(line.get_ydata(), dtype=float)
        finite.extend(y[np.isfinite(y)])

    # at least one plotted detector point (besides the y=0 reference line) with a finite value.
    assert any(v != 0 for v in finite)
    assert all(np.isfinite(v) for v in finite)


def test_public_hpge_mass_comparison_keyword():
    # the public_geometry keyword forces the testdata path, matching core.construct.
    from pygeoml200 import plot_hpge_mass_comparison

    ax = plot_hpge_mass_comparison(public_geometry=True)
    assert len(ax.get_lines()) > 0
