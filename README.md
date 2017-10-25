# SERVAL (SpEctrum Radial Velocity AnaLyser)
calculate radial velocities from stellar spectra

The concept of SERVAL is described in https://arxiv.org/abs/.

## Install instruction

Requirements:
- python 2.7 + numpy, scipy, pyfits
- gnuplot

Setup the path:
```
export SERVALHOME=~/mzechmeister
export SERVAL=$SERVALHOME/serval/
export PYTHONPATH=$PYTHONPATH:$SERVALHOME/python/
```
You can include these lines into your `~/.bashrc`.
Above is for bash; for tcsh:
```
setenv SERVALHOME ~
...
```

Download SERVAL and required tools
```
mkdir $SERVALHOME
cd $SERVALHOME
git clone https://github.com/mzechmeister/serval.git
git clone https://github.com/mzechmeister/python.git
```

Make main files executable
```
chmod u+x $SERVAL/src/serval.py
chmod u+x $SERVAL/src/read_spec.py
mv $SERVAL/zoom.gnu ~/
```

A few c programs come precompiled. Likely not nessecary
```
cd $SERVAL/src/
gcc -c  -Wall -O2 -ansi -pedantic -fPIC polyregression.c; gcc -o polyregression.so -shared polyregression.o
gcc -c  -Wall -O2 -ansi -pedantic -fPIC psplinelib.c; gcc -o psplinelib.so -shared psplinelib.o
```

A first try to check whether there are any conflict. It should list all available options:
```
$SERVAL/src/serval.py --help
```

If you have a ~/bin folder, a useful shortcut is:
```
ln -s $SERVAL/src/serval.py ~/bin/serval
```
Otherwise, an alias can be create and included in `~/.bashrc`.
```
alias serval=$SERVAL/src/serval.py
```
and you can run it as
```
serval --help
```

### Install problems and solution

An alternative to pyfits might be (in ds9.py and read_spec.py):
```
import astropy.io.fits as pyfits
```

## SERVAL in action

A basic example is:
```
serval gj3917 /path/to/e2ds/ -inst HARPS -targrade 06:00:03.495 +02:42:23.67 -targpm 311.1 -42.4 -vref auto
```
Usage of star.cat and more explanation is in prepapration...
