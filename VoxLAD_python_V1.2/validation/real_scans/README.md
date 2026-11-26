# Checks on real scans of a RIEGL VZ-400

`vz400_checks.py` compares the tilt handling of VoxLAD with what the instrument did, using every
pulse of the raw scan file. It needs data that cannot be distributed: the scans themselves and
RIEGL's RiVLib library (available to owners of an instrument from the RIEGL member area).

## Preparing the input

1. Build a reader that writes every pulse of an `.rxp` file with RiVLib: direction, time and
   number of returns of each pulse (`<name>_shots.bin`, 33 bytes per pulse: three float64 for the
   direction, one float64 for the time, one uint8), and pulse index, return index, range,
   amplitude, reflectance and deviation of each return (`<name>_echoes.bin`, 25 bytes per return:
   uint32, uint8, float64, three float32).
2. Make it also write the beam origin of each pulse, as three float32 per pulse, to
   `<name>_origins.bin`. The beam of a VZ-400 leaves the mirror up to 10 cm from the origin of the
   scanner's frame, and an exported point is `origin + range x direction`.
3. Run it on an upright scan, and optionally on a scan taken on a tilt mount for which the
   registration on the upright scan is known (a text file with the 4 x 4 matrix).

## Running

```
python vz400_checks.py FOLDER UPRIGHT_NAME [TILTED_NAME MATRIX_FILE]
```

The results obtained in November 2026 are in `../README.md`, section "Real scans".
