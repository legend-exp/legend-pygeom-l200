# Changing optical properties

It might be useful for some studies to deviate from the default values provided
for optical material properties.

More information on the available properties can be found on the
{mod}`legendoptics` documentation. Especially important is the functionality of
the {mod}`legendoptics.store` which allows users to easily swap some of the
default implementation with their own.

As an example, consider this simple file `lar-ly.py` that can be specified on
the command line: `legend-pygeom-l200 --pygeom-optics-plugin lar-ly.py`:

```python3
from __future__ import annotations

import pint
from legendoptics.lar import lar_scintillation_params

u = pint.get_application_registry()
lar_scintillation_params.replace_implementation(
    lambda _: lar_scintillation_params.original_impl()(flat_top_yield=40 / u.keV)
)
```

This sets the flat-top light yield (as described in
{meth}`legendoptics.lar.lar_scintillation_params`) to 40 photons/keV.
