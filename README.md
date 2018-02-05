# SERVAL (SpEctrum Radial Velocity AnaLyser)
calculate radial velocities from stellar spectra

The concept of SERVAL is described in http://adsabs.harvard.edu/abs/2017A%26A...609A..12Z.

Currently, SERVAL can process data from CARM_VIS, CARM_NIR, HARPS, and HARPN.

## Recent changes
* request to simbad to get RADE, PM, and PL via -targ
* barycentric correction with Wright & Eastman (2014) is now default (requires https://github.com/shbhuk/barycorrpy)

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

Install barycorrpy
```bash
pip install barycorrpy
```
See also https://github.com/shbhuk/barycorrpy/wiki/1.-Installation for other possibilities.

A few c programs come precompiled. Therefore, it is likely not necessary to compile but in case (e.g. Mac OS) ...
```bash
cd $SERVAL/src/
gcc -c  -Wall -O2 -ansi -pedantic -fPIC polyregression.c; gcc -o polyregression.so -shared polyregression.o
gcc -c  -Wall -O2 -ansi -pedantic -fPIC cbspline.c; gcc -o cbspline.so -shared cbspline.o
f2py -c -m spl_int spl_int.f
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
