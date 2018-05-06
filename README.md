# CharFunToolR: The Characteristic Functions Package
R repository of characteristic functions and tools for their combinations and numerical inversion.

For the R version of the package see the CharFun package development available at

- https://github.com/gajdosandrej/CharFunToolR current status of the MATLAB toolbox (not an identical clone) see the CharFunTool development available at

- https://github.com/witkovsky/CharFunTool

About
=====

The Characteristic Functions Package (CharFunToolR) consists of a set of algorithms for evaluating selected characteristic funcions
and algorithms for numerical inversion of the (combined and/or compound) characteristic functions, used to evaluate the probability density function (PDF) and the cumulative distribution function (CDF).
                                                                              
The package includes inversion algorithm, including those based on simple trapezoidal rule for computing the integrals defined by the Gil-Pelaez formulae, and/or by using the FFT algorithm for computing the Fourier transform integrals.
                                                                       
Installation and requirements
============================= 

CharFunToolR was developed with R version 3.4.2 (2017-09-28).

To install, you can either clone the directory with Git or download a .zip file or install package from R GitHub repository.

## Option 1: Download .zip file

Download a .zip of CharFunToolR from

- https://github.com/gajdosandrej/CharFunToolR/releases

After unzipping, you will need to open CharFunToolR.Rproj.

## Option 2: Clone with Git

To clone the CharFunToolR repository, first navigate in a terminal to where you want the repository cloned, then type
```
git clone https://github.com/gajdosandrej/CharFunToolR.git
```
and you will need to open CharFunToolR.Rproj.

## Option 3: Install package from source archive file

You can download package from 

- https://github.com/gajdosandrej/CharFunToolR/releases

in install packages in RStudio you chose Package Archive File (.tar.gz).

## Option 4: Install package from GitHub

Just run the following command in your RStudio console:
```
devtools::install_github("gajdosandrej/CharFunToolR")
```

Getting started
===============

We recommend taking a look at the Examples collection. 

To get a taste of what computing with CharFunToolR is like, type
```
   cf <- function(t) exp(-t^2/2)  # the standard normal characteristic function (CF)
   result <- cf2DistGP(cf)   # Invert the CF to get the CDF and PDF   
```


License
=======

See `LICENSE` for CharFunToolR licensing information.
