from __future__ import annotations

import copy
from importlib import resources

from dbetto import AttrsDict, TextDB


class PublicMetadataProxy:
    """Provides proxies to transparently replace legend hardware metadata with sample data."""

    def __init__(self):
        dummy = TextDB(resources.files("l200geom") / "configs" / "dummy_geom")

        self.chmap = dummy.channelmap
        self.diodes = _DiodeProxy(dummy)
        self.fibers = _FiberProxy()

    def update_special_metadata(self, special_metadata) -> None:
        # the string is shorter because of missing special detectors.
        special_metadata.hpge_string["7"].minishroud_delta_length_in_mm = -200


class _DiodeProxy:
    def __init__(self, dummy_detectors: TextDB):
        self.dummy_detectors = dummy_detectors

    def __getitem__(self, det_name: str) -> AttrsDict:
        det = self.dummy_detectors[det_name[0] + "99000A"]
        m = copy.copy(det)
        m.name = det_name
        return m


class _FiberProxy:
    def __getitem__(self, det_name: str) -> AttrsDict:
        m = {
            "name": det_name,
            "type": "inner" if det_name.startswith("IB") else "outer",
            "geometry": {"tpb": {"thickness_in_nm": 1000}},
        }
        return AttrsDict(m)
