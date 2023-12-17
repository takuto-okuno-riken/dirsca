[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)

# DirSCA and FFM toolbox
Directional SCA (Seed-based connectivity analysis) and Functional flat mapping toolbox for MATLAB

## Introduction
TODO:

<b>Command line tools</b>

| name | description |
|:---|:---|
| dirsca | Calculate directional and non-directional Seed-based connectivity analysis. |
| flatmap | Plot functional flat mapping. |
| plotsca | Plot SCA result with background image.|

## Requirements: Software
* MATLAB R2019b or later
* Image Processing Toolbox ver11.0 or later
* Parallel Computing Toolbox ver7.1 or later
* [VARDNN Toolbox](https://github.com/takuto-okuno-riken/vardnn)

Please download the [VARDNN Toolbox](https://github.com/takuto-okuno-riken/vardnn) and "Add Path" in the MATLAB before using DirSCA and FFM Toolbox.


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
