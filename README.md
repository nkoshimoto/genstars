[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-yellow.svg)](http://www.gnu.org/licenses/gpl-3.0)



# genstars
`genstars` is a stellar population synthesis tool for our Galaxy that generates stars using the Galactic model developed by [Koshimoto, Baba & Bennett (2021)](https://ui.adsabs.harvard.edu/abs/2021ApJ...917...78K/abstract).  
A difference from other population synthesis tools is that the Galactic model used is optimized for the bulge direction.  
Currently the simulation is only available inside a box with -9.5 < l < 9.5 and -10.0 < b < 4.5 because it uses the VVV extinction map combining the ones by [Gonzalez et al. (2012)](https://ui.adsabs.harvard.edu/abs/2012A%26A...543A..13G/abstract) and [Surot et al. (2020)](https://ui.adsabs.harvard.edu/abs/2020A%26A...644A.140S/abstract) to calculate extinction values. 

Although the Galactic model was originally developed for the Galactic bulge region in |b| > 1 deg. `genstars` uses a newly developed version of the model that includes the nuclear stellar disk (NSD) model based on [Sormani et al. (2022)](https://ui.adsabs.harvard.edu/abs/2022MNRAS.512.1857S/abstract).
A paper describing the new model is currently in preparation, and will be linked here as long as it becomes available.  

A simulator for microlensing events, [`genulens`](https://github.com/nkoshimoto/genulens), is also available.

The copyright of an included supplementary code, "option.c", belongs to Ian A. Bond and Takahiro Sumi.


## Before installation
Please ensure that GSL (GNU Scientific Library) is installed in your environment.
If you do not have GSL, please download it from the link provided on the [GSL page](https://www.gnu.org/software/gsl/), and install it following README or INSTALL in the downloaded directory.

You can check where you have GSL by
```
gsl-config --libs
```
If the command gsl-config does not work in the terminal, it probably means that the GSL lib is not installed, or unknown to the OS.


## Installation
If you have `git`, you can download the package by
``` 
git clone https://github.com/nkoshimoto/genstars.git
```
This is recommended because that way you can track any future updates with `git`.

If you do not have `git`, you can simply download the package by clicking the green button "Code" on the upper right in [the repository page](https://github.com/nkoshimoto/genstars), selecting "Download ZIP", and then unzipping it.

You need a C compiler to `make` using Makefile to compile the program genstars.c in the downloaded directory.
The default compiler is `clang`, which is available in macOS.
Please replace the first uncommented line in Makefile with `gcc` or any other C compiler if you prefer.
If you are not sure if you have `clang` or `gcc`, you can check it by
```
which clang (or gcc)
```
If your terminal returns a full path for it, then you have the compiler. 

You might also need to change the paths for GSL specified in Makefile.
There are two lines to specify paths for GSL in the file;
> INCLUDE = -I/opt/local/include
> LINK = -L/opt/local/lib

Please replace the path in INCLUDE with the path returned by
```
gsl-config --cflags
```
and replace the path in LINK with the path returned by
``` 
gsl-config --libs 
```



After making sure that you specify your C compiler and the paths for GSL in Makefile, you can compile the program by
```
make
```

If
```
./genstars
```
returns output lines that starts from
> \#   Output of "./genstars "

and ends with
> \# (n\_BD n\_MS n\_WD n\_NS n\_BH)/n\_all= (  78317 144302  25174   1085    503 ) / 249381 = ( 0.314046 0.578641 0.100946 0.004351 0.002017 )


you are ready to use `genstars`. Note that the exact numbers of the end line might depend on your environment because the calculation uses random numbers.
Please make sure that the input\_files/ directory is in the same directory as where you run `genstars`.


## Usage

Some example runs and their explanations are presented in [genstars\_samples.ipynb](https://github.com/nkoshimoto/genstars/blob/main/genstars_samples.ipynb)

