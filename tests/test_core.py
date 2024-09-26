from __future__ import annotations

import numpy as np


def test_import():
    import l200geom  # noqa: F401


def test_construct():
    from l200geom import core

    core.construct(use_detailed_fiber_model=False)
    core.construct(use_detailed_fiber_model=True)


def test_material_store():
    # replacing material properties is _not_ a core functionality of this package, but
    # we have to make sure that replaced material properties from the optics package are
    # propagated correctly to the generated GDML files.

    from legendoptics import store
    from legendoptics.fibers import fiber_core_refractive_index

    from l200geom import core

    # test that replaced material properties are reflected in the GDML.
    fiber_core_refractive_index.replace_implementation(lambda: 1234)
    reg = core.construct(use_detailed_fiber_model=False)
    rindex = reg.defineDict["ps_fibers_RINDEX"].eval()
    assert np.all(rindex[:, 1] == [1234, 1234])

    # test that after the reset, the created GDML contains the original values again.
    store.reset_all_to_original()
    reg = core.construct(use_detailed_fiber_model=False)
    rindex = reg.defineDict["ps_fibers_RINDEX"].eval()
    assert np.all(rindex[:, 1] == [1.6, 1.6])
