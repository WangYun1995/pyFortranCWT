# pyFortranCWT
**Python wrappers of the [FortranCWT](https://github.com/WangYun1995/FortranCWT) codes created with f2py**

If you use our codes or your research is related to our paper, please kindly cite the following paper, in which the fast CWT algorithms are reviewed and compared.

``` bib
@article{Wang2023,
    author = {{Wang}, Yun and {He}, Ping},
    title = "{Comparisons between fast algorithms for the continuous wavelet transform and applications in cosmology: The one-dimensional case}",
    journal = {RAS Techniques and Instruments},
    year = {2023},
    month = {06},
    doi = {10.1093/rasti/rzad020},
    url = {https://doi.org/10.1093/rasti/rzad020},
    note = {rzad020},
    eprint = {https://academic.oup.com/rasti/advance-article-pdf/doi/10.1093/rasti/rzad020/50599483/rzad020.pdf},
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
The above commands will yield ".so" files, which can be imported in Python. Then add ```export PYTHONPATH="${PYTHONPATH}:/.../.../pyFortranCWT/"``` to your ```~/.bashrc``` and run the following commond

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
