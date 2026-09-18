# shadowpy

A cleaner Shadow interface in Python.

`shadowpy` wraps [Shadow3](https://github.com/oasys-kit/shadow3)'s Python API with
small, user-friendly classes for building and tracing synchrotron beamlines.

## Features

- `BeamLine` – orchestrates a sequence of optical elements and caches the beam
  after every element.
- `Source` / `BendingMagnet` – source wrappers.
- `ToroidalMirror`, `Screen`, `Slits` – optical element wrappers with validation
  against the specification file.
- **Reference frames** – elements are positioned relative to their neighbours
  and can be translated/tilted in the lab frame (`tx/ty/tz`, `rx/ry/rz`).
- **Automatic retracing** – reading an `.image` after an element parameter
  changed triggers a trace from the first stale element.
- `Screen.parallel_caustic` / `simple_caustic` – beam size scans around the focus.
- `si_beamlines.CAXSim` – a ready-made simulation of the CARCARÁ-X beamline.

## Installation

```bash
pip install shadowpy
```

Installing from a checkout:

```bash
pip install .
```

Requires Python 3.9+.

### Dependencies

- `shadow3` (provides the `Shadow` module)
- `numpy`
- `attrs`

## Quick start

```python
from shadowpy.si_beamlines import CAXSim

beamline = CAXSim(total_rays=100_000)

print(beamline.mirror.rx)   # mirror tilt, mrad
beamline.mirror.rx = 0.1    # change it...

image = beamline.dvf_B1.image   # ...and read the screen (retraces automatically)
```

Build a custom beamline:

```python
import Shadow
from shadowpy.beamline import BeamLine
from shadowpy.sources import BendingMagnet
from shadowpy.optical_elements import ToroidalMirror, Slits
from shadowpy.utils import ReferenceFrame, rotation_matrix

beam = Shadow.Beam()
mirror = ToroidalMirror("M1", "path/to/mirror.txt")

beamline = BeamLine(
    name="my-beamline",
    beam=beam,
    source=BendingMagnet("B1", "path/to/source.txt"),
    optical_elements=[mirror, Slits("A1", "path/to/slits.txt")],
    beamline_frame=ReferenceFrame(orientation=rotation_matrix("x", 90)),
)
```

## Layout

- `shadowpy/` – package source
- `shadowpy/beamline_config/` – beamline specification files (shipped as
  package data and resolved with `importlib.resources`)

## Notes

Tracing a beamline writes Shadow3's runtime files (e.g. `star.00`, `SRANG`,
`effic.00`) to the current working directory. These are safe to delete and are
covered by the repository's `.gitignore`.