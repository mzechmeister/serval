# SERVAL (SpEctrum Radial Velocity AnaLyser)
calculate radial velocities from stellar spectra

## Install instruction
excecute + include to ~/.bashrc
```
export SERVALHOME=~/
export SERVAL=$SERVALHOME/serval/
export PYTHONPATH=$PYTHONPATH:$SERVALHOME/python/zechmeister/
```

above is for bash
for tcsh:
```
setenv SERVALHOME ~
...
```

download serval and required tools
```
https://github.com/mzechmeister/serval.git $SERVALHOME
https://github.com/mzechmeister/python.git $SERVALHOME/python/zechmeister
```

make main files executable
```
chmod u+x $SERVAL/src/serval.py
chmod u+x $SERVAL/src/read_spec.py
mv $SERVALHOME/zoom.gnu ~/
```

likely not nessecary
```
gcc -c  -Wall -O2 -ansi -pedantic -fPIC polyregression.c; gcc -o polyregression.so -shared polyregression.o
gcc -c  -Wall -O2 -ansi -pedantic -fPIC psplinelib.c; gcc -o psplinelib.so -shared psplinelib.o
```

for running the program see
```
./src/serval.py --help
```

a helpful short cut
```
ln -s $SERVAL/src/serval.py ~/bin/serval
```

and you can run it as
```
serval --help
```

an alternative to pyfits might be (in ds9.py and read_spec.py):
```
import astropy.io.fits as pyfits
```

