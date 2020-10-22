# SERVAL (SpEctrum Radial Velocity AnaLyser)
calculate radial velocities from stellar spectra

The concept of SERVAL is described in http://adsabs.harvard.edu/abs/2017A&A...609A..12Z [[pdf](https://www.aanda.org/articles/aa/pdf/2018/01/aa31483-17.pdf)].

Currently, SERVAL can process data from CARM_VIS, CARM_NIR, ELODIE, ESPRESSO, HARPS, HARPN, and SOPHIE.

## Install instruction

Requirements:
- python (2.7 or 3.6+) + numpy, scipy, astropy
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

Install barycorrpy:
```bash
pip install --user barycorrpy
```
See also https://github.com/shbhuk/barycorrpy/wiki/1.-Installation for other possibilities.

A few c programs come precompiled. But probably it is necessary to compile (e.g. Mac OS) ...
```bash
cd $SERVAL/src/
gcc -c  -Wall -O2 -ansi -pedantic -fPIC polyregression.c; gcc -o polyregression.so -shared polyregression.o
gcc -c  -Wall -O2 -ansi -pedantic -fPIC cbspline.c; gcc -o cbspline.so -shared cbspline.o
python -m numpy.f2py -c -m spl_int spl_int.f
cd $SERVAL/src/BarCor
gfortran bary.f -o bary.e
```

A first try to check whether there are any conflicts. It should list all available options:
```bash
$SERVAL/src/serval.py --help
```

If you have a `~/bin` folder, a useful shortcut is:
```bash
ln -s $SERVAL/src/serval.py ~/bin/serval
ln -s $SERVAL/src/srv.py ~/bin/srv
```
Otherwise, an alias can be created and included in `~/.bashrc`.
```bash
alias serval=$SERVAL/src/serval.py
```
and you can run it as
```bash
serval --help
```

## SERVAL in action

A basic example is:
```bash
mkdir data
(cd data; git clone https://github.com/mzechmeister/HARPS.git)
serval gj699 data/HARPS/gj699/ -inst HARPS -targ gj699
```

`-targ` requests the coordinates from simbad (otherwise RA and DEC from fits header is used).

After serval has finished, you can inspect the results with `srv.py`, for instance
```bash
srv gj699 -rv -x
```

### Tip

You may want to include the following lines in your `~/.gnuplot`:
```
set colors classic
load "~/mzechmeister/python/zoom.gnu"
```
In gnuplot 5 this uses the old color scheme from gnuplot 4.

And [zoom.gnu](https://github.com/mzechmeister/python/blob/master/zoom.gnu) gives you some additional features, like pan and zoom with keyboard and arrows keys. In particular, very useful for the `look` and `lookt` options to explore the spectra.


Further tips are giving in the [wiki](../../wiki/).
