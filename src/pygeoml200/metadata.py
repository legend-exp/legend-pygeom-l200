from __future__ import annotations

import copy
from importlib import resources

from dbetto import AttrsDict, TextDB


class PublicMetadataProxy:
    """Provides proxies to transparently replace legend hardware metadata with sample data."""

    def __init__(self):
        dummy = TextDB(resources.files("pygeoml200") / "configs" / "dummy_geom")

        self.chmap = dummy.channelmap
        self.diodes = _DiodeProxy(dummy)
        self.fibers = _FiberProxy()

    def update_special_metadata(self, special_metadata: AttrsDict) -> AttrsDict:
        # the string is shorter because of missing special detectors.
        special_metadata = copy.deepcopy(special_metadata)
        special_metadata.hpge_string[7].minishroud_delta_length_in_mm = -200

        # mark as readonly to match real loaded metadata.
        special_metadata.__readonly__ = True
        return special_metadata


class _DiodeProxy:
    def __init__(self, dummy_detectors: TextDB):
        self.dummy_detectors = dummy_detectors

    def __getitem__(self, det_name: str) -> AttrsDict:
        # create the detector from the matching dummy metadata.
        det = self.dummy_detectors[det_name[0] + "99000A"]
        m = copy.deepcopy(det)
        m.name = det_name

        # also test the code paths with no enrichment set.
        if det_name[0] == "P":
            m.production.enrichment.val = None

        # mark as readonly to match real loaded metadata.
        m.__readonly__ = True
        return m


class _FiberProxy:
    def __getitem__(self, det_name: str) -> AttrsDict:
        m = {
            "name": det_name,
            "type": "inner" if det_name.startswith("IB") else "outer",
            "geometry": {"tpb": {"thickness_in_nm": 1000}},
        }
        # mark as readonly to match real loaded metadata.
        return AttrsDict(m, readonly=True)
