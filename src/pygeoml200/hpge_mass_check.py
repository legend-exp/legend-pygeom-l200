from __future__ import annotations

import contextlib
import copy
import logging
from typing import TYPE_CHECKING

import numpy as np
from git import GitCommandError
from legendmeta import LegendMetadata
from pygeomhpges import make_hpge
from pygeomtools.utils import load_dict_from_config

from .core import DEFAULT_METADATA_TIMESTAMP
from .metadata import PublicMetadataProxy

if TYPE_CHECKING:
    from matplotlib.axes import Axes

log = logging.getLogger(__name__)


def plot_hpge_mass_comparison(
    config: dict | None = None,
    *,
    allow_cylindrical_asymmetry: bool = True,
    public_geometry: bool = False,
    ax: Axes | None = None,
) -> Axes:
    """Plot the relative mass difference ``(MC - measured) / measured`` for each HPGe in the geometry.

    The mass of each detector, built with :func:`pygeomhpges.make_hpge`, is compared to the measured mass
    from the metadata and grouped by manufacturer. The axes drawn into are returned.

    Parameters
    ----------
    config
        the geometry configuration dictionary accepted by :func:`pygeoml200.core.construct`. If it carries a
        truthy ``"public_geom"`` flag (or ``public_geometry`` is set), the public testdata metadata is used
        instead of a real LEGEND metadata checkout.
    allow_cylindrical_asymmetry
        passed to :func:`pygeomhpges.make_hpge`; if ``False``, build detectors with the cylindrically
        symmetric base class, ignoring non-symmetric features.
    public_geometry
        force the use of the public testdata metadata, overriding the ``config["public_geom"]``
        auto-detection. Mirrors the ``public_geometry`` argument of :func:`pygeoml200.core.construct`.
    ax
        a :class:`matplotlib.axes.Axes` to draw into. If ``None``, a new figure and axes are created.
    """
    import matplotlib.pyplot as plt

    config = config if config is not None else {}

    lmeta = None
    public = public_geometry or bool(config.get("public_geom", False))
    if not public:
        with contextlib.suppress(GitCommandError):
            lmeta = LegendMetadata(lazy=True)
        if lmeta is None:
            msg = "real LEGEND metadata is required for the HPGe mass comparison but could not be loaded"
            raise RuntimeError(msg)
    if lmeta is None:
        dummy_geom = PublicMetadataProxy()

    timestamp = config.get("metadata_timestamp", DEFAULT_METADATA_TIMESTAMP)
    channelmap = load_dict_from_config(
        config,
        "channelmap",
        lambda: lmeta.channelmap(timestamp) if lmeta is not None else dummy_geom.chmap,
    )
    diodes = lmeta.hardware.detectors.germanium.diodes if lmeta is not None else dummy_geom.diodes

    # only consider the germanium detectors actually built into the geometry, i.e. the geds channels.
    geds = channelmap.map("system", unique=False).get("geds", {})
    names = sorted(ch.name for ch in geds.values())

    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    group_color = {"Mirion": colors[0], "Ortec": colors[1]}

    # group -> {detector name: relative mass difference in percent}
    group_values: dict[str, dict[str, float]] = {}

    for name in names:
        meta = copy.deepcopy(diodes[name])

        measured_mass = meta.production.get("mass_in_g")
        if measured_mass is None:
            log.warning("%s has no measured mass in metadata, skipping", name)
            continue

        # some detectors miss the enrichment value: fall back to a dummy value. The metadata stores it
        # either as a plain number or as a {"val": ..., "unc": ...} mapping.
        enrichment = meta.production.get("enrichment")
        enrichment_val = enrichment.val if hasattr(enrichment, "val") else enrichment
        if enrichment_val is None:
            log.warning("%s has no enrichment in metadata, setting to dummy value 0.92", name)
            meta.production.enrichment = 0.92

        hpge = make_hpge(meta, None, allow_cylindrical_asymmetry)
        sim_mass = hpge.mass.to("g").m

        rel_diff = 100 * (sim_mass - measured_mass) / measured_mass

        group = meta.production.manufacturer
        group_values.setdefault(group, {})[name] = rel_diff

    if ax is None:
        _, ax = plt.subplots(figsize=(16, 3))
    assert ax is not None

    for group, values in sorted(group_values.items()):
        # keep one entry per detector name (NaN elsewhere) so all series share the same, ordered x axis.
        y = [values.get(name, np.nan) for name in names]
        ax.plot(names, y, "o", color=group_color.get(group, colors[2]), markersize=4, label=group)

    ax.axhline(0, color="gray", linewidth=1.5, zorder=0)
    ax.grid(visible=True)
    ax.set_title("difference between simulated and measured HPGe mass")
    ax.set_ylabel("(MC - measured) / measured [%]")
    ax.tick_params(axis="x", labelrotation=90)
    ax.legend(ncol=3)
    return ax
