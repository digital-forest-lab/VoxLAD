# Validation

Everything quoted in the manual and in the code review can be reproduced from this folder.
Run the scripts from here; they use the example scans in `../examples/shea_tree`.

| Script | What it checks | Result (November 2026, Python 3.13, NumPy 2.5, Numba 0.68) |
| --- | --- | --- |
| `compare_octave.py all` | Python in MATLAB-compatible mode against the original `voxlad_v7.m` run in GNU Octave 8.4 | 5,760 voxel-scan records: all pulse counts identical; A_proj, d_bar, G agree to 14 digits; 2,880 result lines identical |
| `engine_checks.py` | Default (traversal) engine against a brute-force test and against the MATLAB rules | Identical to brute force (400 voxel-scans); identical to the MATLAB rules for the 8,256 voxel-scans where the scanner is outside the voxel's extent; the MATLAB rules lose 0.9% of blocked pulses in the other 560 |
| `path_length_checks.py` | Eq. 7, the MATLAB constant and the closed form against a Monte Carlo of chord lengths | Eq. 7 within 1.4% of truth for K <= 1.5 at any voxel size; MATLAB constant right at 0.3 m only |
| `example_effects.py` | Default against MATLAB mode on the example tree; voxel size and path-length treatment | Total leaf area 37.38 m2 (MATLAB mode) and 37.67 m2 (default) on the same grid; 3.8% of voxels differ by more than 5% |
| `single_vs_double.py` | Single-precision storage of v7 against double | Counts differ by 1 in 5 of 1,656 voxel-scans; 834 voxels evaluated instead of 840 |
| `benchmark.py` | Timings | See the manual, section Speed and memory |
| `tilt_checks.py` | Instrument not level, on simulated scans: the head axis found from the scan lines, and LAD of a box of leaves seen by a leaning scanner | Axis found to better than 0.001 degree for tilts of 0 to 90 degrees; mean LAD within 0.7% of the level scan with the axis found, up to 19% off with the axis taken as vertical |
| `real_scans/vz400_checks.py` | The same on real scans of a RIEGL VZ-400, against the pulses the instrument fired | See "Real scans" below |
| `../tests` (`pytest`) | 39 tests, including simulated scans of a medium of known LAD and the whole pipeline in MATLAB mode against the Octave run (same 2,016 voxels; 19 of 4,032 records differ by one return lying exactly on a voxel face) | All pass |

## The Octave reference

`octave/reference/` holds the outputs of the Octave runs:

- `thin400_*`, `thin100_*`: the two scans thinned 1:400 and 1:100 (angular step multiplied by
  the square root of the factor), whole grid, double-precision storage;
- `thin400_single_*`: the same 1:400 run with the single-precision storage of v7 as written;
- `full_*`: the full scans, 24 voxels listed in `subset_voxels.txt` (the six with most leaf
  returns and 18 drawn at random).

`*_dbg.txt` has one line per voxel and scan: `x_min y_min z_min scan A_proj Re(G) Im(G) d_bar
Nt Nb Ni_leaf Ni_wood Ni_noise theta`. `*_out_1.txt` and `*_out_2.txt` are the result files of
the original script, including its `here` / `there` debug lines.

To regenerate them (several hours; Octave runs the per-return loop at about 50 microseconds
per return, i.e. 40 s per voxel and scan on the full scans):

```
python octave/make_reference.py /path/to/voxlad_v7.m ../examples/shea_tree /tmp/voxlad_ref
cd /tmp/voxlad_ref && bash run_all.sh
```

`make_reference.py` applies to the original script only the edits listed at the top of the
file: file paths, number of scans, one array size that otherwise stops v7 with an index
error, a debug file, optionally double-precision storage and a restriction to some voxels.
For the full scans it fills the workspace with a fast loader and uses the script's own
`skip_read` switch, because Octave needs many minutes to read 900,000 returns with the
original loop.

The same comparison can be made against a MATLAB run: write the files with the scripts that
`make_reference.py` produces, copy them into `octave/reference/` under the names above and
run `compare_octave.py`.

## Real scans (instrument not level)

`real_scans/vz400_checks.py` was run in November 2026 on scans of a RIEGL VZ-400 that are not
distributed with VoxLAD: an upright scan and a scan on a 90 degree tilt mount taken at Pasoh
(30 April 2018, 22.5 and 3.1 million pulses), and a tilt-mount scan taken at Harvard Forest
(15 August 2017). The raw `.rxp` files were read pulse by pulse with RIEGL's RiVLib, so that
every pulse fired is known, with or without a return. `real_scans/README.md` says how to
prepare the input.

Upright scan, point cloud in the scanner's own frame (head axis = z exactly), written in
millimetres and turned by known angles:

| Case | Axis accepted | Error of the axis (degrees) |
| --- | --- | --- |
| All returns, cloud turned by 0, 0.3, 1.7, 3, 10, 45 and 90 degrees | yes | 0.0002 to 0.0003 |
| Sector of 30 degrees | yes | 0.0033 |
| Sector of 10 degrees, elevation 10 to 35 degrees (one tree) | yes | 0.0022 |
| The same, returns 8 to 25 m away only | yes | 0.0011 |
| Sector of 4 degrees, elevation 20 to 30 degrees | no (standard error 0.11) | 0.50 |
| Box of 6 m by 6 m, 10 m from the scanner | yes | 0.0022 |
| File rounded to centimetres | yes | 0.0002 |
| Every fourth return kept | yes | 0.0004 |
| Declared position 5 cm off (found: 5.0 cm) | yes | 0.0014 |
| Declared position 30 cm off, 28 cm of it across the axis (found: 28.0 cm) | yes | 0.0078 |
| Declared position 20 cm too high (along the axis: not found) | yes | 0.031 |
| Returns in random order | no (no scan lines) | not estimated |

The two tilt-mount scans, in their own frame and turned by 0 to 90 degrees, gave errors of
0.0020 and 0.0021 degrees. The scan lines of this instrument are cones of 1.04 degrees.

Angular steps measured from the returns: 0.04009 degrees between lines and 0.03997 along them,
against 0.04014 and 0.03999 from all the pulses of the file (the three facets of the mirror
give spacings between lines of 0.046, 0.046 and 0.027 degrees, whose mean is what counts).

N_T of the model against the pulses fired at 150 voxels of 0.5 m, 6 to 25 m from the scanner
(about 1,000 to 35,000 pulses per voxel), in frames turned by known angles:

| Frame turned by | Axis taken as vertical: rms, worst | True axis: rms, worst |
| --- | --- | --- |
| 0 degrees | 0.62%, -2.4% | 0.62%, -2.4% |
| 1 degree | 0.83%, -2.3% | 0.46%, +1.4% |
| 3 degrees | 2.43%, +7.5% | 0.30%, +1.2% |
| 10 degrees | 7.70%, +28.9% | 0.22%, +1.0% |

The tilt-mount scan of Pasoh, placed in the frame of the upright scan with a registration
matrix obtained independently (surface matching of the two clouds): the axis found from its
returns is 0.0058 degrees from the z axis of that matrix, and N_T is within 0.19% (rms) of the
pulses fired with either axis, whereas the vertical axis gives values several times too large.

These scans come from a single instrument; other instruments remain to be tested.
