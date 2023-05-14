import numpy as np
import v97cwt
from v97cwt import v97cwt

print(v97cwt.cv97cwt_periodbc.__doc__)

x      = np.linspace(0,2*np.pi, 512)
signal = np.loadtxt("F:/Works/20220617_Work/20230202_cosmo_apps/dens_field_1d/delta_100.dat",dtype=np.float64)
Lbox   = 1.0
Nsubs  = 12 
wavelet= "mw"

v97cwt.cv97cwt_periodbc( signal, Lbox, Nsubs, wavelet)

scales  = v97cwt.scales
realcwt = v97cwt.realcwt
imagcwt = v97cwt.imagcwt

np.save("scales.npy", scales)
np.save("realcwt.npy", realcwt)
np.save("imagcwt.npy", imagcwt)