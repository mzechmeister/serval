# SERVAL (SpEctrum Radial Velocity AnaLyser)
calculate radial velocities from stellar spectra

The concept of SERVAL is described in http://adsabs.harvard.edu/abs/2017A%26A...609A..12Z.

Currently, SERVAL can process data from CARM_VIS, CARM_NIR, HARPS, and HARPN.

## Install instruction

Requirements:
- python 2.7 + numpy, scipy, pyfits
- gnuplot

Setup the path:
```bash
export SERVALHOME=~/mzechmeister
export SERVAL=$SERVALHOME/serval/
export PYTHONPATH=$PYTHONPATH:$SERVALHOME/python/
```
You can include these lines into your `~/.bashrc`.
Above is for bash; for tcsh:
```tcsh
setenv SERVALHOME ~
...
```

Download SERVAL and required tools:
```bash
mkdir $SERVALHOME
cd $SERVALHOME
git clone https://github.com/mzechmeister/serval.git
git clone https://github.com/mzechmeister/python.git
```

Make main files executable:
```bash
chmod u+x $SERVAL/src/serval.py
chmod u+x $SERVAL/src/read_spec.py
```

A few c programs come precompiled. Therefore, it is likely not necessary to compile but in case ...
```bash
cd $SERVAL/src/
gcc -c  -Wall -O2 -ansi -pedantic -fPIC polyregression.c; gcc -o polyregression.so -shared polyregression.o
gcc -c  -Wall -O2 -ansi -pedantic -fPIC cbspline.c; gcc -o cbspline.so -shared cbspline.o
cd $SERVAL/src/BarCor
gfortran bary.f -o bary.e
```

A first try to check whether there are any conflicts. It should list all available options:
```bash
$SERVAL/src/serval.py --help
```

If you have a ~/bin folder, a useful shortcut is:
```bash
ln -s $SERVAL/src/serval.py ~/bin/serval
```
Otherwise, an alias can be create and included in `~/.bashrc`.
```bash
alias serval=$SERVAL/src/serval.py
```
and you can run it as
```bash
serval --help
```

### Install problems and solution

An alternative to pyfits might be (in ds9.py and read_spec.py):
```python
import astropy.io.fits as pyfits
```

## SERVAL in action

A basic example is:
```bash
serval gj3917 /path/to/e2ds/ -inst HARPS -targ gj3917 -vref auto
```
More explanation is in prepapration...
