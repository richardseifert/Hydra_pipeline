from fitstools import manage_dtype, mask_fits
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
plt.ion()

@manage_dtype()
def make_fiber_profile_map(img, fmask):
    fiber_profile_map = np.zeros_like(img)
    fnums = list({n for row in fmask for n in row if n != 0})
    for fnum in fnums:
        print fnum
        fiber = mask_fits(img, fmask, maskval=fnum)
        n = 0
        for i,spectral_slice in enumerate(fiber):
            n += 1
            x = [j for j in range(len(spectral_slice)) if fmask[i][j]==fnum]
            y = [spectral_slice[j] for j in range(len(spectral_slice)) if fmask[i][j]==fnum]
            yfit = fit_spatial_profile(x, y)
            yfit_norm = yfit/np.sum(yfit)
            for j in range(len(x)):
                fiber_profile_map[i][int(x[j])] = yfit_norm[j]
    return fiber_profile_map

def fit_spatial_profile(x, y):
    x = np.asarray(x)
    y = np.asarray(y)

    gaussian = lambda x, amp, mu, sig: amp*np.exp(-1/2*((x-mu)/(sig))**2)
    linear = lambda x, a, b: a*x+b

    func = lambda x, amp, mu, sig, a, b: gaussian(x, amp, mu, sig)+linear(x, a, b)
    p0 = [max(y), 0.5*(x[0]+x[-1]), 1, 0, 0.5*(y[0]+y[-1])]


    try:
        coeff, err = curve_fit(func, x, y, p0=p0)
    except RuntimeError:
        print 'RuntimeWarning: Could not do gaussian fit.'
        coeff = p0

        fig, ax = plt.subplots()
        ax.scatter(x, y)
        xfit = np.linspace(x[0], x[-1], 1000)
        ax.plot(xfit, func(xfit, *coeff))
    yfit = func(x, *coeff)

    return yfit
