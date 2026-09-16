# Bayesian Stokes inversion code

---

## Overview

An MPI-parallel Bayesian inversion program for polarized
spectral profiles. It uses the DREAM Markov-chain Monte Carlo algorithm and a
Milne-Eddington forward model. The original method is described in
[MCMC Inversion of Stokes Profiles](https://ui.adsabs.harvard.edu/abs/2019ApJ...875..127L). The GEMC simulation is removed from the current version.

The current implementation was refactored and optimized from the earlier
codebase with assistance from OpenAI Codex. Scientific validation was
performed by comparing the revised implementation with the original program.

The current version supports FITS image cubes, single-profile DAT input,
multiple spectral lines and wavelength regions, Cartesian or spherical
magnetic parameters, per-pixel noise, asymmetric credible errors, and optional
sample output. Noise may be derived from Stokes I or supplied as a shared or
per-pixel wavelength-dependent array.

---

## Build

MPI and CFITSIO are required. GSL is optional.

```sh
./configure
make
```

Use `./configure --enable-gsl` to build with GSL random-number routines. Run
`./configure --help` for the remaining compiler options.

---

## Run

```sh
cd example
mpirun -np 4 ../MCMCINV input.inv
```

The program is MPI-only, but a serial-sized run is available through
`mpirun -np 1`. FITS input can distribute profiles among several islands; DAT
input contains one profile and therefore uses one island. Input and output
paths are resolved from the working directory.

See [example/input.inv](example/input.inv) for the current input keywords.
Results are written to the configured FITS output. When sample output is
enabled, samples are stored in the binary `SMPL` sample format.
[tools/analyze_samples.ipynb](tools/analyze_samples.ipynb) provides plotting
examples, while [tools/inversion_results.py](tools/inversion_results.py)
contains the file readers.

---

## Observation input formats

Select the format with `profile_format = FITS` or `profile_format = DAT`.
Wavelengths must be finite, strictly increasing, and contain at least five
samples.

### FITS input

`data_path` identifies the observation FITS file. Its primary HDU contains the
Stokes profiles with FITS-axis order:

```text
[nw, 4, nx, ny]
```

The three-dimensional form `[nw, 4, nx]` is also accepted and is interpreted
as `ny = 1`. The second axis is ordered as Stokes I, Q, U, and V. Profile data
may use 16-bit integer, 32-bit floating-point, or 64-bit floating-point FITS
images.

`wavelength_path` identifies a separate FITS file whose primary HDU is a
one-dimensional 32-bit or 64-bit floating-point array of length `nw`. The
wavelength values and the `lines` input values are expressed in Angstrom.

For example:

```ini
profile_format = FITS
data_path = ./profiles.fits
wavelength_path = ./wavelength.fits
```

### DAT input

DAT input contains exactly one spectral profile. Each nonempty data row must
contain exactly five finite values:

```text
wavelength  I  Q  U  V
```

Lines beginning with `#` are ignored. DAT mode reads the wavelength directly
from `data_path`, uses one island, and currently requires `noise_mode = INTENSITY`.

Example:

```text
# wavelength        I             Q             U             V
6302.2431       1.0280e6       1.2310e3      -8.4200e2       2.1050e3
6302.2531       1.0275e6       1.1980e3      -8.0100e2       2.2440e3
```

---

## Noise input

All supplied noise values represent the standard deviation `sigma`.
The likelihood code converts them internally to `1 / sigma^2`.
Every noise value must be finite and strictly positive.

### `noise_mode = INTENSITY`

No noise HDU is required. `noise_level` supplies four fractional noise levels
for I, Q, U, and V. For each profile, the corresponding standard deviation is:

```text
sigma_stokes = noise_level_stokes * max(I)
```

For example:

```ini
noise_mode = INTENSITY
noise_level = 1e-3, 1e-3, 1e-3, 1e-3
```

### `noise_mode = GLOBAL`

The second HDU of `data_path` contains one wavelength-dependent noise array
shared by every image pixel. Its FITS-axis order must be:

```text
[nw, 4]
```

The four rows correspond to I, Q, U, and V. The array is read once by world
rank 0 and broadcast to the other MPI ranks. `noise_level` is ignored.

### `noise_mode = PIXEL`

The second HDU of `data_path` contains a separate wavelength-dependent noise
array for every image pixel. It must match the profile geometry:

```text
[nw, 4, nx, ny]
```

For observations with `ny = 1`, `[nw, 4, nx]` is accepted. The program reads
the corresponding `4 * nw` values whenever it reads a pixel. `noise_level` is
ignored.

---

## License

- This project is licensed under the MIT License.
