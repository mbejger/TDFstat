genseg
======

Generates time segments of any length.
Part of the TDFstat pipeline by Polgraw - Polish Virgo group.

This code assembles STS (short time series) chunks stored in a HDF5 file
(produced by extract_band) into time segments of desired length.
In addition:
- it applies science (analysis ready) mask
   (Tukey window is used to smooth edges of continuous regions)
- removes outliers
- writes ephemeris for the detectorr (requires lalsuite)

Compilation:
```
conda activate lal
make genseg-hdf
```

Run:
```
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
./genseg-hdf <config.g2d>
```
