from __future__ import annotations

import logging
from pathlib import Path

import pyg4ometry.geant4 as g4
import vtk
from pyg4ometry import visualisation
from pyg4ometry.gdml.Defines import Auxiliary

log = logging.getLogger(__name__)


def visualize(registry: g4.Registry) -> None:
    v = visualisation.VtkViewerColouredNew()
    v.addLogicalVolume(registry.worldVolume)

    _load_color_auxvals_recursive(registry.worldVolume)
    registry.worldVolume.pygeom_color_rgba = False  # hide the wireframe of the world.
    _color_recursive(registry.worldVolume, v)

    # v.addClipper([0, 0, 0], [1, 0, 0], bClipperCloseCuts=False)

    v.buildPipelinesAppend()
    v.addAxes(length=5000)
    v.axes[0].SetVisibility(False)  # hide axes by default.

    # override the interactor style.
    v.interactorStyle = _KeyboardInteractor(v.ren, v.iren, v)
    v.interactorStyle.SetDefaultRenderer(v.ren)
    v.iren.SetInteractorStyle(v.interactorStyle)

    # set some defaults
    _set_camera(v, up=(1, 0, 0), pos=(0, 0, +20000))

    v.view()


class _KeyboardInteractor(vtk.vtkInteractorStyleTrackballCamera):
    def __init__(self, renderer, iren, vtkviewer):
        self.AddObserver("KeyPressEvent", self.keypress)

        self.ren = renderer
        self.iren = iren
        self.vtkviewer = vtkviewer

    def keypress(self, obj, event):
        # predefined: "e"xit

        if self.iren.GetKeySym() == "a":  # toggle "a"xes
            ax = self.vtkviewer.axes[0]
            ax.SetVisibility(not ax.GetVisibility())

            self.ren.GetRenderWindow().Render()

        if self.iren.GetKeySym() == "u":  # "u"p
            _set_camera(self, up=(0, 0, 1), pos=(-20000, 0, 0))

        if self.iren.GetKeySym() == "b":  # "b"ottom
            _set_camera(self, up=(1, 0, 0), pos=(0, 0, +20000))

        if self.iren.GetKeySym() == "F1":
            _set_camera(self, up=(0.55, 0, 0.82), pos=(-14000, 0, 8000))

        if self.iren.GetKeySym() == "s":  # "s"ave
            _export_png(self.vtkviewer)


def _set_camera(v, up, pos):
    v.ren.GetActiveCamera().SetViewUp(*up)
    v.ren.GetActiveCamera().SetPosition(*pos)
    v.ren.ResetCameraClippingRange()
    v.ren.GetRenderWindow().Render()


def _export_png(v, fileName="scene.png"):
    ifil = vtk.vtkWindowToImageFilter()
    ifil.SetInput(v.renWin)
    ifil.ReadFrontBufferOff()
    ifil.Update()

    png = vtk.vtkPNGWriter()
    png.SetFileName("./" + fileName)
    png.SetInputConnection(ifil.GetOutputPort())
    png.Write()


def _color_recursive(lv: g4.LogicalVolume, viewer: visualisation.ViewerBase) -> None:
    if hasattr(lv, "pygeom_color_rgba"):
        for vis in viewer.instanceVisOptions[lv.name]:
            if lv.pygeom_color_rgba is False:
                vis.alpha = 0
                vis.visible = False
            else:
                vis.colour = lv.pygeom_color_rgba[0:3]
                vis.alpha = lv.pygeom_color_rgba[3]
                vis.visible = vis.alpha > 0

    for pv in lv.daughterVolumes:
        if pv.type == "placement":
            _color_recursive(pv.logicalVolume, viewer)


def _color_macro_recursive(lv: g4.LogicalVolume, macro_lines: dict) -> None:
    if hasattr(lv, "pygeom_color_rgba") and lv.name not in macro_lines:
        mac = f"/vis/geometry/set/forceSolid {lv.name}\n"
        if lv.pygeom_color_rgba is False or lv.pygeom_color_rgba[3] == 0:
            mac += f"/vis/geometry/set/visibility {lv.name} -1 false\n"
        else:
            rgba = " ".join([str(c) for c in lv.pygeom_color_rgba])
            mac += f"/vis/geometry/set/colour {lv.name} 0 {rgba}\n"
        macro_lines[lv.name] = mac

    for pv in lv.daughterVolumes:
        if pv.type == "placement":
            _color_macro_recursive(pv.logicalVolume, macro_lines)


def generate_color_macro(registry: g4.Registry, filename: str) -> None:
    """Create a Geant4 macro file containing the defined visualization attributes."""
    macro_lines = {registry.worldVolume: None}
    _color_macro_recursive(registry.worldVolume, macro_lines)
    macro_contents = "".join([m for m in macro_lines.values() if m is not None])

    with Path.open(filename, "w", encoding="utf-8") as f:
        f.write(macro_contents)


def write_color_auxvals(registry: g4.Registry) -> None:
    """Append an auxiliary structure, storing the visualization color information."""
    written_lvs = set()

    def _append_color_recursive(lv: g4.LogicalVolume) -> None:
        if hasattr(lv, "pygeom_color_rgba") and lv.name not in written_lvs:
            if lv.pygeom_color_rgba is False or lv.pygeom_color_rgba[3] == 0:
                rgba = "-1"
            else:
                rgba = ",".join([str(c) for c in lv.pygeom_color_rgba])
            lv.addAuxiliaryInfo(Auxiliary("rmg_color", rgba, registry, addRegistry=False))
            written_lvs.add(lv.name)

        for pv in lv.daughterVolumes:
            if pv.type == "placement":
                _append_color_recursive(pv.logicalVolume)

    written_lvs.add(registry.worldVolume.name)  # do not store world vis args.
    _append_color_recursive(registry.worldVolume)


def _load_color_auxvals_recursive(lv: g4.LogicalVolume) -> None:
    auxvals = list(filter(lambda aux: aux.auxtype == "rmg_color", lv.auxiliary))
    assert len(auxvals) <= 1
    if len(auxvals) == 1 and not hasattr(lv, "pygeom_color_rgba"):
        rgba = auxvals[0].auxvalue
        if rgba == "-1":
            lv.pygeom_color_rgba = False
        else:
            lv.pygeom_color_rgba = list(map(float, rgba.split(",")))

    for pv in lv.daughterVolumes:
        if pv.type == "placement":
            _load_color_auxvals_recursive(pv.logicalVolume)
