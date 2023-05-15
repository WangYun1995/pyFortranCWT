# pyFortranCWT
**Python wrappers of the [FortranCWT](https://github.com/WangYun1995/FortranCWT) codes created with f2py**

If you use our codes or your research is related to our paper, please kindly cite the following paper, in which the fast CWT algorithms are reviewed and compared.

``` bib
@ARTICLE{Wang2023,
       author = {{Wang}, Yun and {He}, Ping},
        title = "{Comparisons between fast algorithms for the continuous wavelet transform and applications in cosmology: the one-dimensional case}",
      journal = {arXiv e-prints},
     keywords = {Astrophysics - Instrumentation and Methods for Astrophysics, Physics - Computational Physics},
         year = 2023,
        month = feb,
          eid = {arXiv:2302.03909},
        pages = {arXiv:2302.03909},
          doi = {10.48550/arXiv.2302.03909},
archivePrefix = {arXiv},
       eprint = {2302.03909},
 primaryClass = {astro-ph.IM},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2023arXiv230203909W},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```

## Requirements
- FFTW3
- Python3
- Numpy
- GNU Fortran 95 compiler
- MinGW (for Windows)

## Usage
- ## Linux:
```
$:git clone  https://github.com/WangYun1995/pyFortranCWT.git
$:cd pyFortranCWT/fftcwt
$:python -m numpy.f2py -c Parameters.f95 Tools.f95 FFTCWT.f95 -I/.../fftw3/include -L/.../fftw3/lib -lfftw3 --link-fftw3 --fcompiler=gnu95 --compiler=unix --opt='-O3' -m fftcwt
$:cd ..
$:cd v97cwt
$:python -m numpy.f2py -c Parameters.f95 Tools.f95 V97CWT.f95 --fcompiler=gnu95 --compiler=unix --opt='-O3' -m v97cwt
$:cd ..
$:cd m02cwt
$:python -m numpy.f2py -c Parameters.f95 Tools.f95 M02CWT.f95 --fcompiler=gnu95 --compiler=unix --opt='-O3' -m m02cwt
$:cd ..
$:cd a19cwt
$:python -m numpy.f2py -c Parameters.f95 Tools.f95 A19CWT.f95 --fcompiler=gnu95 --compiler=unix --opt='-O3' -m a19cwt
```
The above commands will yield ".so" files, which can be imported in Python. Then add ```export PYTHONPATH="${PYTHONPATH}:/.../.../pyFortranCWT/``` to your ```~/.bashrc``` and run the following commond

```
$:source ~/.bashrc
```

- ## Windows:
```
> git clone  https://github.com/WangYun1995/pyFortranCWT.git
> cd pyFortranCWT\fftcwt
> python -m numpy.f2py -c Parameters.f95 Tools.f95 FFTCWT.f95 -I"\...\fftw3\include" -L"\...\fftw3\lib" -lfftw3 --link-fftw3 --fcompiler=gnu95 --compiler=mingw32 --opt='-O3' -m fftcwt
> cd ..
> cd v97cwt
> python -m numpy.f2py -c Parameters.f95 Tools.f95 V97CWT.f95 --fcompiler=gnu95 --compiler=mingw32 --opt='-O3' -m v97cwt
> cd ..
> cd m02cwt
> python -m numpy.f2py -c Parameters.f95 Tools.f95 M02CWT.f95 --fcompiler=gnu95 --compiler=mingw32 --opt='-O3' -m m02cwt
> cd ..
> cd a19cwt
> python -m numpy.f2py -c Parameters.f95 Tools.f95 A19CWT.f95 --fcompiler=gnu95 --compiler=mingw32 --opt='-O3' -m a19cwt
```
The above commands will yield ".pyd" files, which can be imported in Python. Then add ```...\...\pyFortranCWT``` to the PYTHONPATH. (see [How to add to the PYTHONPATH in Windows](https://stackoverflow.com/questions/3701646/how-to-add-to-the-pythonpath-in-windows-so-it-finds-my-modules-packages))

## Examples
see [pyFortran_examples.ipynb](https://github.com/WangYun1995/pyFortranCWT/blob/main/pyFortranCWT_examples.ipynb) for examples.
