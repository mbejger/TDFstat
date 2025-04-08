---
title: Input data generation
---

# Input data generation

## General scheme

The input data for TDFstat pipeline are generated from SFDB files created by the Rome Group using [PSS library](https://git.ligo.org/pia.astone/create_sfdb) (login required). Those files contain wide-band (e.g. 2048 Hz) Short-time Fourier Transforms. Two codes are used to create long, time domain series suitable for TDFstat pipeline:

* `extract_band_hdf_td` - extracts narrow band Short-Time Fourier Transforms from the SFDB, performs inverse FFT and writes Short Time Series (STS) into a HDF5 file. This code is based on the original `extract_band` code from PSS which only extracts narrow band SFTs and writes them to text files.
* `genseg` - combines STS-es into longer time segments, applies regions mask (e.g. regions with Analysis Ready flag) and cleans the data (outliers removal). In addition the code generates ephemeris files for coresponding time segments.

The general scheme of data generation is illustrated in the picture.

[![General scheme](./img/tdfstat_data_gen-v2.png){width=50%}](./img/tdfstat_data_gen-v2.png)


### Example

Typical All-sky search data generation workflow with O3 numbers:

1. We used SFDB files with 2048 Hz bandwidth which contain STFT of length 1024 s
   overlapped in time by half (512 s).
2. `extract_band` was used to extract narrow bands with $B=0.25\ Hz$. The
   resulting sampling interval was $dt = 1/(2\ B) = 2\ s$. 
   Since STFT length is 1024 s we have 512 samples/STFT, but we use only middle
   half of the inverse FFT which makes 256 samples per time chunk.
3. `genseg` was used to assemble short time chunks into segments of length 6
   sidereal days which contain N samples:

$$N = round(nod\cdot C\_SIDDAY/dt) = 258492$$
  
   where $C\_SIDDAY$ is the duration of sidereal day in seconds.
   The whole O3 run time span was ~1 yr so we've got 60 segments.


### Band and segment numbering

For practical reasons we assingn integer numbers to our **frequency bands**
according to this general formula:

$$fpo = 10 + (1 - 2^{-ov})\cdot bbbb\cdot B \quad \mathrm{[Hz]}.$$

where:

* $fpo$ - starting/reference frequency of the band,
* $bbbb$ - integer band number,
* $B$ - bandwidth,
* $2^{-ov}$ - band overlap; $oi$ is natural number; this form assures alignement of bands with Fourier bins in the SFDB.

In the names of files and directories, the **band number**, $bbbb$ is formatted using `%04d` specifier, e.g. `0072`.

The **time segment** (or simply 'segment') has length $nod$ which is INTEGER
number of days.  
Segments are assigned subsequent natural numbers starting from 1.  
In the names of files and directories the segment number ($nnn$) is formatted using `%03d` specifier, e.g. `001`.


----
## extract_band code

This program converts Short Fourier Transformation series to time series.
It was written by Pia Astone (INFN, Physics Department of University of Rome "La
Sapienza") a part of the PSS package used to create SFDB. The
`extract_band_hdf_td` code is its substantially modified version by Paweł
Ciecieląg, used in TDFstat since run O4.
The major differences comprise: inverse FFT to get short time series (STS) and
output to HDF5 file.

### Compilation

Source: currently part of a [private gitlab
repo](https://gitlab.camk.edu.pl/polgraw-cw/scripts)  
Prerequisities: C compiler, selected files from the PSS library ([PSS
library](https://git.ligo.org/pia.astone/create_sfdb)), FFTW3 & HDF5 libraries.  

To compile, in the `extract_band/pss_sfdb` directory type
```
make extract_band_hdf_td
```

### Running

`extract_band_hdf_td` reads following input parameters from stdin:
```
string[170] out_file: output file name
string[170] in_file: input file name (file with a list of all SFDB files)
float fpo: reference frequency
float B: badwidth
```
e.g. create file `0072_H1.in`
```
sfft_0072_H1.h5
../../sfdb_files_H1.list
26.2
0.25
```
and pass it as input:
```
./extract_band_hdf_td < 0072_H1.in
```

Important remarks

- The file names can be specified with absolute or relative paths.
- The code is designed to work in an *incremental mode* - if you add more files
  to `in_file` and `out_file` exists, the new data chunks will be **appended**
  to it. For chunks already present in the HDF file only attributes gps_sec, nfft and
  sfdb_name are verified, no actual data is read.
- the code is very generic - it does not know which band number is
  processed (but it can be encoded in the name of the input file).


### Helper script

To simplify mass data generation, e.g. for all-sky search, we provide an example
of bash script: `eb2hdf.sh`.
It uses `extract_band` is to create STS data for all bands and all detectors
following certain convention. Before using the script, first edit some parameters in 
section marked "EDIT HERE":

- eb - full path to the `extract_band_hdf_td` executable
- B and ov - bandwidth and overlap as defined above
- ilist - path to file containing list of SFDB files to be processed; default: `./sfdb_${det}.list`
- odir - path to the output directory (will be created); default: `sts_B${B}_ov${ov}/${det}`
- fpo - check the definition of `fpo( <band_number> )` function

To generate HDF file for detector `<det_name>` and for bands between 
`<start_band>` and `<end_band>` as follows run:
```
./eb2hdf.sh <det_name> <start_band> [<end_band>]
```
Currently det_name can be H1, L1 or V1. If `<end_band>` is omitted then only
single band is generated.  
The ilist file should contain one file per line, to create use: `ls -1 /somedirectory/*.SFDB09`  
The input file for `extract_band` will be saved in:
`sts_B${B}_ov${ov}/${det}/${B4}_${det}.in`  
The standard output from `extract_band` will be saved in: `sts_B${B}_ov${ov}/${det}/eb2hdf-${B4}_${det}.out`  


### Output HDF5 file format

To display contents of the HDF5 file use:
```
h5dump -q creation_order <file_name>
```
*Our HDF files have creation_order tracking and indexing flags enabled for easier
sorting of chucks (typically the chunks will be created in an increasing ichunk
order).*

The internal file structure:
```
HDF5 object                  | data type
----------------------------------------------------
/                            | root group
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
├── dataset "ichunk"         | dataset
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

Please refer to the `genseg` code for example how to read this HDF5 file in C.

#### Old/original extract_band

<details>
<summary>Click to expand</summary>
Running:
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
Output to a text file:
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
</details>

----
## genseg code

This code combines STS (short time series) chunks stored in a HDF5 file
(produced by `extract_band`) into time segments of any length. In addition:

* it applies science (analysis ready) mask (Tukey window is used to smooth edges
  of each continuous region)
* removes outliers
* writes ephemeris files (DetSSB.bin) for the detector (requires lalsuite)


### Compilation:

[Source code](https://github.com/Polgraw/TDFstat/tree/main/genseg)  
Prerequisities: C compiler, HDF5 library, lalsuite (if -DUSE_LAL is used in
Makefile, this is required to generate ephemeris)

lalsuite can be installed from conda-forge repository using [miniforge installer](https://conda-forge.org/download/):
```
mamba create -n lal lalsuite
mamba activate lal
```
To compile type `make genseg-hdf`.  
The old genseg version (used before O4) is preserved in subdirectory `old`.

### Running:

`genseg-hdf` requires a configuration file in the INI format. An example file is
provided in the source directory
[H1_0072_6d.g2d](https://github.com/Polgraw/TDFstat/blob/main/genseg/H1_0072_6d.g2d).
All configuration options are explained in the comments in this file.

* To generate list of Analysis Ready segments this python script can be used:

```python
from gwpy.segments import DataQualityFlag, SegmentList
import sys

GPSstart = int(sys.argv[1]) # 1238166018 (1st April 2019, 15 UTC) 
GPSend   = int(sys.argv[2]) # 1253314818 
det_channel_name = sys.argv[3] # e.g. 'H1:DMT-ANALYSIS_READY'

dqf = DataQualityFlag.query(det_channel_name, GPSstart, GPSend)
seg_science = SegmentList(dqf.active)
print(seg_science)
```
```bash
python get_sci_segment_list.py 1238166018 1253314818 H1:DMT-ANALYSIS_READY > C00_H1_gps_sci_segments
python get_sci_segment_list.py 1238166018 1253314818 L1:DMT-ANALYSIS_READY > C00_L1_gps_sci_segments
python get_sci_segment_list.py 1238166018 1253314818 V1:ITF_SCIENCE:1 > C00_V1_gps_sci_segments
```

* To calculate ephemeris the code should be run under proper conda environment
  and two efemerid files must be specified in the config file. They can be
  extracted from lalpulsar to te current directory this way:
  
```bash
cp $CONDA_PREFIX/share/lalpulsar/sun00-40-DE405.dat.gz .
cp $CONDA_PREFIX/share/lalpulsar/earth00-40-DE405.dat.gz .
gunzip sun00-40-DE405.dat.gz earth00-40-DE405.dat.gz
```

To run genseg type `./genseg-hdf <config file>`. 


## TDFstat input data structure

The TDFstat pipeline (actually the mein `search` code) assumes that input data
have fixed directory and fileneme structure.
The input data are divided into time segments of typically a few days length and consists - for each detector - of the input time series data, the ephemerides and the grid-generating matrix file (defining the parameter space of the search).

A single `search` run requires 2 data files for each detector `DD` and segment
`nnn`, stored in `data_dir/nnn/DD` subdirectory, where `DD` is currently either
`H1` (Hanford), `L1` (Livingston) or `V1` (Virgo Cascina): 

* `xdat_nnn_bbbb.bin` - time-domain narrow-band data sequence (`bbbb` is the number of frequency band),
* `DetSSB.bin`  -  location  of  the  detector  w.r.t. the Solar System
  Barycenter (SSB), in Cartesian coordinates, sampled at `dt` sampling rate
  (array of size `2N`).  
  The last two records in this file are the angle `phir`, determining the  position
  of Earth in its diurnal motion, and the obliquity of the ecliptic `epsm`,
  both calculated for the first sample of the data. 

Third file is the sky positions-frequency-spindown grid file in linear coordinates (common for all the detectors), stored in `data_dir/nnn` in case of the network search (one grid file is used by all the detectors) or in each detector directory separately (in case of single-detector searches):

   * `grid.bin` - generator matrix of an optimal grid of templates (defining the parameter space; see [here](/grid_generation) for details).

A typical directory structure is as follows:

```
xdat_O3_C01/
├── 001/
	├── grid.bin
	├── H1/
	│   ├── DetSSB.bin
	│   ├── grid.bin
	│   ├── starting_date
	│   └── xdat_001_1234.bin
	└── L1/
		├── DetSSB.bin
		├── grid.bin
		├── starting_date
		└── xdat_001_1234.bin
```

Beginning of each time frame is saved in the `nnn/DD/starting_date` file, e.g.,
```
% cat 2d_0.25/001/H1/starting_date
1.1260846080e+09
```

With this structure the root directory `xdat_O3_C01` is passed to the search code
as parameter `indir`.

An example for two LIGO detectors H1 and L1, and data frame segments $nnn=001-008$ with pure Gaussian noise 2-day time segments with sampling time equal to 2s for a fiducial narrow band number $bbbb=1234$ (`xdatc_nnn_1234.bin`) coresponding the the band frequency $fpo=308.859375$ is [available here](https://polgraw.camk.edu.pl/H1L1_2d_0.25.tar.gz).


## Gaussian data generator

`genseg` directory contains additional code `gauss-xdat-mask.c` to generate time series drawn from the Gaussian distribution.

### Compilation

Prerequisites: C compiler , GSL library
```
% gcc gauss--xdat-mask.c -o gauss-xdat-mask -lm -lgsl -lgslcblas
```

### Running

The program takes input values from the command line:
```
% ./gauss-xdat-mask N amplitude sigma output-file [xdat-template]
```
where N is the length of the time series (in samples), amplitude and sigma are
Gauss function parameters, output-file is the name of the output file and
xdat-template is an optional parameter: name of the existing xdat file from which
only gaps (zeros) will be copied to the output (e.g. the real xdat file -
usefull for testing).

Example:
```
% ./gauss-xdat-mask 86164 1 1 xdat_001_1234.bin
```
The output is a binary file, `xdat_001_1234.bin` containing  86164 float-precision numbers.
