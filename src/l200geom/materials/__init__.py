"""Subpackage to provide all implemented materials and their (optical) material properties."""

from __future__ import annotations

import pyg4ometry.geant4 as g4


class OpticalMaterialRegistry:
    _elements = {}
    _elements_cb = {}

    def __init__(self, reg: g4.Registry):
        self.g4_registry = reg
        self.lar_temperature = 88.8
        self._define_elements()

    def get_element(self, symbol: str) -> g4.Element:
        if (symbol in self._elements_cb) and (symbol not in self._elements):
            self._elements[symbol] = (self._elements_cb[symbol])()
        return self._elements[symbol]

    def _add_element(
        self, name: str, symbol: str, Z: int, A: float  # noqa: N803
    ) -> None:
        """Lazily define an element on the current registry."""
        assert symbol not in self._elements_cb
        self._elements_cb[symbol] = lambda: g4.ElementSimple(
            name=name, symbol=symbol, Z=Z, A=A, registry=self.g4_registry
        )

    def _define_elements(self) -> None:
        """Lazily define all used elements."""
        self._add_element(name="Hydrogen", symbol="H", Z=1, A=1.00794)
        self._add_element(name="Carbon", symbol="C", Z=6, A=12.011)
        self._add_element(name="Nitrogen", symbol="N", Z=7, A=14.01)
        self._add_element(name="Oxygen", symbol="O", Z=8, A=16.00)
        self._add_element(name="Fluorine", symbol="F", Z=9, A=19.00)
        self._add_element(name="Silicon", symbol="Si", Z=14, A=28.09)
        self._add_element(name="argon", symbol="Ar", Z=18, A=39.95)
        self._add_element(name="Chromium", symbol="Cr", Z=24, A=51.9961)
        self._add_element(name="Manganese", symbol="Mn", Z=25, A=54.93805)
        self._add_element(name="Iron", symbol="Fe", Z=26, A=55.845)
        self._add_element(name="Indium", symbol="In", Z=49, A=114.82)
        self._add_element(name="Cobalt", symbol="Co", Z=27, A=58.9332)
        self._add_element(name="Nickel", symbol="Ni", Z=28, A=58.6934)
        self._add_element(name="Copper", symbol="Cu", Z=29, A=63.55)

    @property
    def liquidargon(self) -> g4.Material:
        """LEGEND liquid argon."""
        if hasattr(self, "_liquidargon"):
            return self._liquidargon

        self._liquidargon = g4.Material(
            name="LiquidArgon",
            density=1.390,  # g/cm3
            number_of_components=1,
            state="liquid",
            temperature=self.lar_temperature,  # K
            pressure=1.0 * 1e5,  # pascal
            registry=self.g4_registry,
        )
        self._liquidargon.add_element_natoms(self.get_element("Ar"), natoms=1)

        return self._liquidargon

    @property
    def metal_steel(self) -> g4.Material:
        """Stainless steel of the GERDA cryostat."""
        if hasattr(self, "_metal_steel"):
            return self._metal_steel

        self._metal_steel = g4.Material(
            name="metal_steel",
            density=7.9,
            number_of_components=5,
            registry=self.g4_registry,
        )
        self._metal_steel.add_element_massfraction(
            self.get_element("Si"), massfraction=0.01
        )
        self._metal_steel.add_element_massfraction(
            self.get_element("Cr"), massfraction=0.20
        )
        self._metal_steel.add_element_massfraction(
            self.get_element("Mn"), massfraction=0.02
        )
        self._metal_steel.add_element_massfraction(
            self.get_element("Fe"), massfraction=0.67
        )
        self._metal_steel.add_element_massfraction(
            self.get_element("Ni"), massfraction=0.10
        )

        return self._metal_steel
