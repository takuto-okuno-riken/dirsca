[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)

# DirSCA and FFM toolbox
Directional SCA (Seed-based connecvitivy analysis) and Functional flat mapping toolbox for MATLAB

## Introduction
TODO:

<b>Command line tools</b>

| name | description |
|:---|:---|
| dirsca | Calculate and plot MTESS for a group of multivariate time-series data. |
| dummy | Generate a group surrogate model (VAR, PCVAR, VARDNN surrogate) and (multivariate time-series) group surrogate data.|

## Requirements: Software
* MATLAB R2019b or later
* Fuzzy Logic Toolbox ver2.6 or later
* Econometrics Toolbox ver5.3 or later
* Parallel Computing Toolbox ver7.1 or later
* [VARDNN Toolbox](https://github.com/takuto-okuno-riken/vardnn)

Please download the [VARDNN Toolbox](https://github.com/takuto-okuno-riken/vardnn) and "Add Path" in the MATLAB before using GSDGM and MTESS Toolbox.


## Installation
1. Download this Toolbox and [VARDNN Toolbox](https://github.com/takuto-okuno-riken/vardnn) zip files.
2. Extract zip files under your working directory <work_path>.
3. Run the MATLAB software, and "Add Path" extracted directories (i.e. <work_path>/vardnn-master and <work_path>/dirsca-master).
4. Move to <work_path>/dirsca-master directory and run the following demos.

## Command Line Tools Demos
<b>Demo 1</b><br>
The first demo shows the calculation of MTESS among time-series data and figure output.<br>
(Copy and paste this command line. Demo data is included in GSDGM and MTESS Toolbox.)
~~~
>> dirsca --showinsig --showmat --showsig --showprop --shownode data/cx-8x500-demo-surrogate.mat 
...
output mat file : results/cx-8x500-demo-surrogate_mtess.mat
~~~
