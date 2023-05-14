# pyFortranCWT
A Python wrapper of the [FortranCWT](https://github.com/WangYun1995/FortranCWT) created with f2py

## Requirements
- FFTW3
- Python3
- Numpy
- gfortran

## Usage
- ## Linux:
```
$:git clone  https://github.com/WangYun1995/pyFortranCWT.git
$:cd pyFortranCWT/fftcwt
$:python -m numpy.f2py -c Parameters.f95 Tools.f95 FFTCWT.f95 -I/.../fftw3/include -L/.../fftw3/lib/libfftw3 -lfftw3 --link-fftw3 --fcompiler=gnu95 --compiler=mingw32 --opt='-O3' -m fftcwt
$:cd pyFortranCWT/v97cwt
$:python -m numpy.f2py -c Parameters.f95 Tools.f95 V97CWT.f95 --fcompiler=gnu95 --compiler=mingw32 --opt='-O3' -m v97cwt
$:cd pyFortranCWT/m02cwt
$:python -m numpy.f2py -c Parameters.f95 Tools.f95 M02CWT.f95 --fcompiler=gnu95 --compiler=mingw32 --opt='-O3' -m m02cwt
$:cd pyFortranCWT/a19cwt
$:python -m numpy.f2py -c Parameters.f95 Tools.f95 A19CWT.f95 --fcompiler=gnu95 --compiler=mingw32 --opt='-O3' -m a19cwt
```
Then add ```export PYTHONPATH="${PYTHONPATH}:/.../.../pyFortranCWT/``` to your ```~/.bashrc``` and run the following commond

```
$:source ~/.bashrc
```
