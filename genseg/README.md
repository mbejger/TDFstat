genseg
======

Generates time segments of any length.
Part of the TDFstat pipeline by Polgraw - Polish Virgo group.

This code assembles short time chunks stored in HDF5 files
(produced by extract_band) into time segments of desired length.
In addition:
- it applies science (analysis ready) mask
   (Tukey window is used to smooth edges of continuous regions)
- removes outliers
- writes ephemeris for the detectors (requires lalpulsar)

Compilation:
conda activate lal
make genseg-hdf

Run:
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
./genseg-hdf <config.g2d>
