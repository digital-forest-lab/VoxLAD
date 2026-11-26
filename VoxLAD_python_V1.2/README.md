# VoxLAD for Python

Leaf area density (m² of leaf per m³) in cubic voxels from classified ground-based
(terrestrial) lidar scans, after

> Béland, M., Widlowski, J.-L., Fournier, R.A. (2014). A model for deriving voxel-level
> tree leaf area density estimates from ground-based LiDAR. *Environmental Modelling &
> Software* 51, 184–189. doi:10.1016/j.envsoft.2013.09.034

This is the Python version of the MATLAB program `voxlad_v7.m`. It gives the same
numbers as the MATLAB program when asked to (`compat = "matlab"`), corrects a number of
defects by default, and accepts any terrestrial scanner that samples on a regular
angular pattern. The two example scans (1.8 million returns) are processed in about three
seconds.

New users should start with the quick start guide, `docs/QUICKSTART.pdf` (four pages); the
user manual, `docs/MANUAL.pdf`, has the details.

## Install

Python 3.9 or newer.

```
pip install .            # from this folder; installs numpy
pip install numba        # recommended: counting 10 to 30 times faster
pip install laspy lazrs  # only for LAS/LAZ input
```

## Try it

```
cd examples/shea_tree
voxlad inspect shea_tree.toml    # describe the two scans and check their settings
voxlad run shea_tree.toml        # about 3 s; writes results/shea_tree_*.csv and a summary,
                                 # then opens the result in 3-D in the web browser
```

The 3-D viewer is a single HTML file written next to the tables (`shea_tree_viewer.html`):
turn, zoom, colour the voxels by leaf area density, percent explored, quality flag or scan,
filter and cut them, and point at a voxel to read its values. `voxlad view shea_tree.toml`
opens it again; `--no-view` runs without it.

The instrument does not have to be level: for panoramic scanners the axis the head turned
about is found from the scan lines in the file, and used for the pulse count (`axis = "auto"`,
the default). The same analysis checks the declared scanner position and angular step against
the returns, and says so in the summary when they do not agree.

Start a configuration file of your own with `voxlad template my_plot.toml`.

From Python:

```python
import voxlad
res = voxlad.run("examples/shea_tree/shea_tree.toml")
print(res.total_leaf_area())     # m2
lad = res.grid("lad")            # 3-D array, NaN where no value
voxlad.view(res)                 # the 3-D viewer, in the web browser
```

## What is where

| Path | Content |
| --- | --- |
| `src/voxlad/` | the package |
| `examples/shea_tree/` | two ILRIS-3D scans of one tree and their configuration file |
| `tests/` | `pytest` suite, including a simulated scan of a medium of known leaf area density |
| `validation/` | comparison with the original MATLAB code run in GNU Octave, and other checks |
| `docs/` | quick start guide (`QUICKSTART.md`, `QUICKSTART.html`, `QUICKSTART.pdf`) and user manual (`MANUAL.md`, `MANUAL.html`, `MANUAL.pdf`) |

Run the tests with `pip install pytest numba` then `pytest`.

## Input in one paragraph

One file per scan, each return with x, y, z and a class (leaf, wood or noise), all scans
in one coordinate frame with z up. For every scan the program needs the scanner position,
the angular step, the scan pattern (`spherical` for panoramic scanners, `uniform` for
window scanners such as the ILRIS-3D) and, if the scan does not cover the full sphere, its
limits. Files in the original VoxLAD `.asc` layout, plain text/CSV, LAS/LAZ and NumPy
arrays are read. The manual has the details.
