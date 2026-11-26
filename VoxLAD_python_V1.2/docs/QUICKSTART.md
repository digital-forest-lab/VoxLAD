---
title: "VoxLAD for Python: quick start"
subtitle: "Version 1.2.0"
---

VoxLAD turns classified ground lidar scans into leaf area density (LAD, in m² of leaf per m³) for each voxel of a grid. This guide holds what is needed to run it, for version 1.2.0; everything else is in the user manual (`docs/MANUAL.pdf` in the package).

# What you need

You need Python 3.9 or newer, one file per scan position and a few facts about each scan.

| Item | What it is |
| --- | --- |
| Scan files | One file per scan position, with x, y, z and a class code for each return. Text, CSV, LAS or LAZ, and NumPy files are read. |
| Classes | Each return is marked as leaf, wood or noise before the run; VoxLAD does not classify. |
| Common frame | All scans in one coordinate frame, in metres, with z up. |
| Scanner position | x, y, z of the scanner in that frame, for each scan. |
| Scan pattern | `spherical` for panoramic scanners (RIEGL VZ, Leica, FARO), `uniform` for window scanners (Optech ILRIS-3D). |
| Angular step | The angle between consecutive pulses, as set at acquisition. This is the setting with the largest effect on the result. |
| Limits of the scan | The vertical field of view, and the left and right edges if the scan is not a full turn. |
| Voxel size and leaf angles | For instance 0.3 m, and `spherical` when the leaf angle distribution is not known. |

The files must hold the returns as the scanner recorded them:

1. Keep every return, including wood, ground and noise. Do not crop, thin or filter the files.
2. Keep one return per pulse: the first one, for instruments that record several.
3. Keep the order in which the scanner software wrote the returns. VoxLAD uses it to find out whether the instrument was level.

# Install

The four commands below, typed in a terminal, install the program and check it.

```bash
cd VoxLAD_python          # the unzipped folder, which holds pyproject.toml
pip install .             # installs the voxlad command
pip install numba         # recommended: the run is ten to thirty times faster
voxlad --version          # prints: voxlad 1.2.0
```

For LAS or LAZ files, also run `pip install laspy lazrs`.

# Try the example

The package holds two scans of one tree and their configuration file; running them takes a few seconds.

```bash
cd examples/shea_tree
voxlad inspect shea_tree.toml     # describes each scan and checks its settings
voxlad run shea_tree.toml         # computes, writes the results and opens the 3-D viewer
```

The run ends with a summary that holds these lines:

```text
voxels evaluated      4402 (voxels holding a leaf return)
  0 ok                4356 voxels   leaf area     37.188 m2
  1 low_sampling        46 voxels   leaf area      0.408 m2
Total leaf area       37.595 m2  (flags ok + low_sampling)
```

The summary also gives one warning per scan about the angular step. It is expected with these two files: the returns in them are spaced 3% wider than the step written in their headers.

# Run your own scans

A run is described by one text file. Three steps lead to a result.

1. Put the scan files in a folder and save the text below there as `my_plot.toml`, with your own values. `voxlad template my_plot.toml` writes a longer version in which every setting is explained.
2. Check the scans with `voxlad inspect my_plot.toml`. It lists, for each scan, the returns per class, the angles covered next to the limits declared and the tilt found for the instrument, and marks with `<-- check` what does not fit.
3. Run the model with `voxlad run my_plot.toml`.

```toml
[grid]
voxel_size = 0.3                # side of the voxels, in metres

[model]
leaf_angle = "spherical"        # or erectophile, planophile, plagiophile, uniform

[classes]                       # the class codes used in your files
leaf = [1]
wood = [2]
noise = [3]

[scanner]                       # the same for all scans
pattern = "spherical"           # panoramic scanner; "uniform" for a window scanner
angular_step_deg = 0.04         # angle between consecutive pulses
zenith_range_deg = [30, 130]    # vertical field of view; 0 is straight up, 90 the horizon

[[scan]]                        # one block per scan position
file = "scan1.txt"              # columns: x y z class
position = [0.0, 0.0, 1.5]      # scanner x y z

[[scan]]
file = "scan2.txt"
position = [12.0, 0.0, 1.5]

[output]
prefix = "results/my_plot"      # start of the names of the result files
```

| Your case | What to set |
| --- | --- |
| The scan is not a full turn | `left_edge = [x, y]` and `right_edge = [x, y]` in its `[[scan]]` block: any point on the left and on the right limit of the scan, as seen from the scanner |
| Window scanner | `pattern = "uniform"`, with both edges and the `zenith_range_deg` of the window |
| Scans with different settings | Move `pattern`, `angular_step_deg` or `zenith_range_deg` into each `[[scan]]` block |
| Columns in another order, or header lines | `columns = ["x", "y", "z", "class"]` in the order of the file, and `skip_rows` |
| Only part of the scene is of interest | `bounds = [xmin, ymin, zmin, xmax, ymax, zmax]` or `z_range = [zmin, zmax]` under `[grid]`; the files stay whole |

# Read the results

Each run writes four files, whose names start with the output prefix, and opens the viewer in the web browser.

| File | Content |
| --- | --- |
| `_voxels.csv` | One line per voxel: its centre `x`, `y`, `z`, `lad` (m² m⁻³), `leaf_area` (m²), `flag`, `best_scan` (the scan used) and `percent_explored` (percent of the pulses aimed at the voxel that reached it) |
| `_summary.txt` | Settings, voxels per flag, total leaf area, leaf area by height, and the warnings of the run |
| `_scans.csv` | One line per voxel and scan, with every intermediate count |
| `_viewer.html` | The voxels in three dimensions: drag to turn, use the wheel to zoom, point at a voxel to read its values. `voxlad view my_plot.toml` opens it again |

Every voxel carries a quality flag, to be read before its value.

| `flag` | Name | Meaning |
| --- | --- | --- |
| 0 | ok | Value given: at least 15% of the pulses aimed at the voxel reached it |
| 1 | low_sampling | Value given, but fewer pulses reached the voxel: less reliable |
| 2 | unresolved | Pulses reached the voxel, but too few went through it for a solution: no value |
| 3 | occluded | No pulse reached the voxel: no value |
| 4 | not_seen | Outside the field of view of every scan, or too close to the scanner: no value |

A voxel without a value holds `nan`, which is not zero. The total leaf area of the summary adds the voxels with flags 0 and 1.

# Before using the numbers

Three checks take a few minutes and catch most problems.

1. **Warnings.** The summary ends with the warnings of the run. The one on the angular step matters most: the number of pulses varies with the inverse square of the step, and in the example a difference of 3% in the step changes the total leaf area by about 20%. Compare the step with the settings of the instrument and correct `angular_step_deg` if needed. A warning on the scanner position means that the scan lines meet more than 10 cm from the position given: check `position`.
2. **Flags.** Many voxels flagged low_sampling, unresolved or occluded mean that the canopy hides too much of itself: add scan positions or use larger voxels. Colouring the viewer by quality flag shows where.
3. **Voxel size.** The result depends on it: the total leaf area of the example goes from 43.5 m² with voxels of 0.1 m to 32.7 m² with voxels of 1 m. Give the voxel size with any result.

For each scan of a panoramic scanner, the `head axis` line of the summary gives the tilt found for the instrument, and the run corrects for it. When that line says that the tilt could not be found, the scan is treated as level. If the instrument was not level, give the registration matrix of the scan with `matrix_file` (see the manual).

# Going further

The user manual covers what this guide leaves out: the description of the scanner with a table of common instruments, every setting of the configuration file, the quality checks, the instrument that was not level, the use from Python and the comparison with the MATLAB version.

Work that uses VoxLAD should cite the article that describes the model: Béland, M., Widlowski, J.-L., Fournier, R.A. (2014). A model for deriving voxel-level tree leaf area density estimates from ground-based LiDAR. *Environmental Modelling & Software* 51, 184–189. [doi:10.1016/j.envsoft.2013.09.034](https://doi.org/10.1016/j.envsoft.2013.09.034)

The software is free for non-commercial use, for academic and research purposes; see `LICENSE.txt` in the package.
