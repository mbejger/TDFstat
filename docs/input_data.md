---
title: Input data generation
---

# Input data generation

## General workflow

The input data for TDFstat pipeline are generated from SFDB files created by the Rome Group using [PSS library](https://git.ligo.org/pia.astone/create_sfdb) (LVK login required). Those files contain wide-band (e.g. 2048 Hz) Short-time Fourier Transforms. Two codes are used to create long, time domain series suitable for TDFstat pipeline:

* `extract_band_hdf_td` - extracts narrow band Short-Time Fourier Transforms from the SFDB, performs inverse FFT and writes Short Time Series (STS) into a HDF5 file. This code is based on the original `extract_band` code from PSS which only extracts narrow band SFTs and writes them to text files.
* `genseg` - combines STSes into longer time segments, applies regions mask (e.g. regions with Analysis Ready flag) and cleans the data (outliers removal). In addition the code generates ephemeris files for coresponding time segments.

The general scheme of data generation is illustrated in the picture.

[![General scheme](./img/tdfstat_data_gen-v2.png){width=50%}](./img/tdfstat_data_gen-v2.png)

### Band and segment numbering

For practical reasons we assingn integer numbers to our frequency bands according to this general formula:

$$fpo = 10 + (1 - 2^{-ov})\cdot band\cdot B \quad \mathrm{[Hz]}.$$

where:

* $fpo$ is starting/reference frequency of the band
* $band$ is integer band number
* $B$ is bandwidth
* $2^{-ov}$ is band overlap (this form assures alignement of bands with Fourier bins in the SFDB).

In the names of files and directories the $band$ is formatted using `%04d` specifier, e.g. `0072`.

Time segments (or simply segments) have length of integer multiple of 1 day and are assigned subsequent integer numbers starting from 1. In the names of files and directories the segment number is formatted using `%03d` specifier, e.g. `001`.


### Example

Typical All-sky search data generation workflow with numbers:

1. We use SFDB files with 2048 Hz bandwidth which contain STFT of length 1024 s overlapped in time by half (512 s).
2. `extract_band` was used to extract narrow bands with $B=0.25 Hz$. For this $B$ we have sampling interval $dt = 1/(2 B) = 2 s$. Since STFT length is 1024 s we have 512 samples/STFT but we use only middle half of the inverse FFT which results in 256 samples per time chunk.
3. `genseg` was used to assemble time segments of length nod=6 days. The whole run was about 1 yr long so we've got 60 segments. Number of samples in the single segment is given by formula:

$$N = round(nod*C\_SIDDAY/dt)$$

where $C\_SIDDAY$ is the duration of sidereal day in seconds.

----
## extract_band

This program converts Short Fourier Transformation series to time series.
Written by Pia Astone (INFN, Physics Department of University of Rome "La Sapienza") a part of the PSS package used to create SFDB.

### Prerequisites

C compiler & standard C libraries (`math.h`). Links to the PSS library (created by Pia Astone).

### Example

```
% extract_band < input_file
```
where `input_file` is an ASCII file containing the following rows:

* Maximal number of SFT
* The name of the output file
* The list of SFT files
* The frequency band in Hz
* The width of frequency band in Hz

e.g.,
```
100000
J0034+1612_2010-10-10.out
J0034+1612_2010-10-10.list
718.2480
1
```

### Output
```
HDF5 object                  | data type
----------------------------------------------------
/root
├── attr: format_version     | int
├── attr: site               | str[3]
├── attr: fpo                | double
├── attr: bandwidth          | float
├── attr: df                 | double
├── attr: dtype              | str[4]
├── attr: sft_overlap        | int
├── attr: nsamples           | int
├── attr: scaling_factor     | double
├── attr: subsampling_factor | double
├── attr: last_ichunk        | int
├── dataset "ichunk"
│   ├── data                 | float data[nsamples]
│   ├── attr: ichunk         | int
│   ├── attr: gps_sec        | int
│   ├── attr: gps_nsec       | int
│   ├── attr: sft_mjdtime    | double
│   ├── attr: sfdb_file      | str[]
│   ├── attr: ctime          | str[]
│   └── attr: nfft           | int
└── dataset "ichunk+1"
    ├── data
    └── ...
```

### old version
```
% Beginning freq- Band- Samples in one stretch- Subsampling factor- inter (overlapping, 2 if data were overlapped)- Frequency step- Scaling factor- ***The data are real and imag of the FFT
% 908.152344 0.250000 256 8192.000000 2  0.0009766 1.000000e-20
% FFT number in the file; Beginning mjd days; Gps s; Gps ns;
% 100 55099.5879745370 937922816 0
 4.59662571e+02  2.27630825e+01
-3.50387007e+02 -2.20005558e+02
 3.57587904e+02  1.01217077e+02
 1.74400486e+02  2.62086552e+02
 2.21804800e+02 -5.20278366e+02
-3.87826732e+02 -1.55758978e+02
```

----
## `genseg`

*TODO* For the implementation see [here](https://github.com/mbejger/polgraw-allsky/tree/master/gen2day).

## Input data structure

The time series input data is divided into time segments of typically a few days length and consists - for each detector - of the input time series data, the ephemerides and the grid-generating matrix file (defining the parameter space of the search).

A single run requires 2 data files for each detector `DD`, stored in `data_dir/nnn/DD` subdirectory, where `DD` is currently either `H1` (Hanford), `L1` (Livingston) or `V1` (Virgo Cascina):

   * `xdat_nnn_bbbb.bin` - time-domain narrow-band data sequence, sampled at  half second. `nnn` is the number of time frame, `bbbb` is the number of frequency band (see below),
   * `DetSSB.bin`  -  location  of  the  detector  w.r.t. the Solar
   System Barycenter (SSB), in Cartesian coordinates, sampled at `dt` sampling rate (array of size `2N`),
   * Last two records in this file  are the angle `phir`, determining the  position of Earth in  its diurnal motion, and the obliquity of  the ecliptic `epsm`, both calculated for the first sample of the data.

Third file is the sky positions-frequency-spindown grid file in linear coordinates (common for all the detectors), stored in `data_dir/nnn` in case of the network search (one grid file is used by all the detectors) or in each detector directory separately (in case of single-detector searches):

   * `grid.bin` - generator matrix of an optimal grid of templates (defining the parameter space; see [here](./grid_generation) for details).

An example for two LIGO detectors H1 and L1, and data frame segments $nnn=001-008$ with pure Gaussian noise 2-day time segments with sampling time equal to 2s for a fiducial narrow band number $bbbb=1234$ (`xdatc_nnn_1234.bin`) coresponding the the band frequency $fpo=308.859375$ is [available here](https://polgraw.camk.edu.pl/H1L1_2d_0.25.tar.gz).

A typical directory structure is as follows:

```bash
001
├── grid.bin
├── H1
│   ├── DetSSB.bin
│   ├── grid.bin
│   ├── starting_date
│   └── xdatc_001_1234.bin
└── L1
    ├── DetSSB.bin
    ├── grid.bin
    ├── starting_date
    └── xdatc_001_1234.bin
```

Beginning of each time frame is saved in the `nnn/DD/starting_date` file, e.g.,
```
% cat 2d_0.25/001/H1/starting_date
1.1260846080e+09
```
Frames `nnn` are labelled with three-digit consecutive number. For the `O1` data, the bandwidth is `0.25 Hz` ($dt = 2s$). For a given $dt$, the reference band frequency `fpo` is defined as

$$fpo = 10 + (1 - 2^{-5})\cdot bbbb\cdot \frac{1}{2dt}\ \mathrm{[Hz]}.$$

Neighboring bands overlap by $2^{-5}/(2dt)\ \mathrm{Hz}$. `O1` data in the frequency range $10-2000\ \mathrm{Hz}$ contains $8220$ narrow `0.25 Hz` bands. With the $dt = 2s$ sampling time, the total number of data points in time segments of 2 sideral day long is `N=86164`. For lower frequencies (`10-475 Hz`, see [documents and publications](./articles)) 6 day length segments are used (`N=258492` double-precision numbers).



## Gaussian input data (for tests)

The directory `search/network/src-cpu` contains a standalone `gauss-xdat` code to generate time series drawn from the Gaussian distribution.

### Prerequisites

C compiler and standard libraries (`math.h`, `sys/time.h` for `gettimeofday`). The code depends on the [GNU Scientific Library (GSL)](http://www.gnu.org/software/gsl/) random number generation (`gsl/gsl_rng.h`) and random number distributions (`gsl/gsl_randist.h`; using the Marsaglia-Tsang ziggurat implementation).

### Compilation

```
% gcc gauss-xdat.c -o gauss-xdat -lm -lgsl -lgslcblas
```

### Example

The program takes input values from the command line:
```
% ./gauss-xdat N amplitude sigma output-file
```
e.g.,

```
% ./gauss-xdat 86164 1 1 ../../../testdata/2d_0.25/001/H1/xdatc_001_1234.bin
```
The output is a binary file containing `N` double-precision numbers.
