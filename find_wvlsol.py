import os
from astropy.io import fits
from astropy import wcs
from fitstools import manage_dtype, mask_fits, row_avg
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import minimize
import numpy as np
from scipy.optimize import curve_fit
from ngaussian import fit_ngaussian
from extract_spectra import extract_counts
polynomial = lambda x, *args: sum([coeff*x**power for power,coeff in enumerate(args)])

@manage_dtype(use_args=[0,1])
def wvlsol(comp, fiber_mask, use_fibers, **kwargs):
    #Initialize a blank wavelength solution.
    wvlsol_map = np.zeros_like(fiber_mask)

    #Load the template wavelength solution.
    template_dat = np.loadtxt('template_wvlsol.dat', delimiter=',')
    p = template_dat[:,2]
    w = template_dat[:,0]
    coeffs = fit_poly(p, w, 3)
    template = lambda x, c=coeffs: polynomial(x, *c)

    #Load thar line list info.
    dat = np.loadtxt('thar_short.fits')
    line_list_wvl = dat[:,0]
    line_list_counts = dat[:,1]
    #If the table of thar peaks does not exist, make it.
    if not os.path.exists('thar_peaks.dat'):
        std, l_peak_x, l_peak_y = fit_ngaussian(line_list_wvl, line_list_counts, 40)
        f = open('thar_peaks.dat', 'w')
        for x, y in zip(l_peak_x, l_peak_y):
            f.write(str(x).ljust(24)+str(y)+'\n')
        f.close()
    else:
        thar_peaks = np.loadtxt('thar_peaks.dat')
        linelist = thar_peaks[:,0]

    def f_wvlsol(fnum, template_wvlsol, wvlsol_map=wvlsol_map):
        #Extract comp spectrum in pixel space.
        comp_counts = extract_counts(comp, fiber_mask, fnum)
        comp_pix = np.arange(len(comp_counts), dtype=np.float64)

        #Find wavelength solution for fiber.
        wsol = fiber_wvlsol(comp_pix, comp_counts, linelist, template_wvlsol, **kwargs)

        #Add individual wavelength solution to wvlsol_map
        wsol_arr = wsol(np.arange(len(wvlsol_map)))
        ones_fiber = np.where(fiber_mask==fnum, np.ones_like(fiber_mask), 0)
        wvlsol_map += np.transpose(np.multiply(np.transpose(ones_fiber), wsol_arr))

        return wsol, wvlsol_map

    #The template solution was generated using fiber 50, so when generating wvlsols, start
    # at fiber 50 and go up, then start at fiber 49 and go down.
    use_fibers_high = sorted([fnum for fnum in use_fibers if fnum > 50])
    use_fibers_low = sorted([fnum for fnum in use_fibers if fnum < 50], key = lambda x: -x)

    center_wsol, wvlsol_map = f_wvlsol(50, template)
    last_wsol = center_wsol
    for fnum in use_fibers_high:
        last_wsol, wvlsol_map = f_wvlsol(fnum, last_wsol)
    last_wsol = center_wsol
    for fnum in use_fibers_low:
        last_wsol, wvlsol_map = f_wvlsol(fnum, last_wsol)

    return wvlsol_map


def fiber_wvlsol(pix, counts, linelist, starter_wvlsol, npeaks = 25, **kwargs):

    ##REMOVE COSMIC RAYS EVENTUALLY

    #Find peaks in the fiber.
    std, npeaks_pix, npeaks_counts = fit_ngaussian(pix, counts, npeaks, **kwargs)
    n = min([5, npeaks])
    template_wvlsol = starter_wvlsol
    while n <= npeaks:
        print 'Finding wvlsol with '+str(n)+' peaks.'
        peaks_pix, peaks_wvl = match_peaks(npeaks_pix[:n], linelist, template_wvlsol)
        n_used = len(peaks_pix)
        print len(peaks_pix), 'lines used out of '+str(n)+'.'
        coeffs = fit_poly(peaks_pix, peaks_wvl, n=3)
        wsol = lambda x, c=coeffs: polynomial(x, *c)
        rsqrd = min_res_sqr(peaks_pix, peaks_wvl, wsol)
        #print 'Lowest chi-squared: '+str(rsqrd/len(peaks_pix))
        if rsqrd/n_used > 0.01:
            print n, 'BAD!!! :('
            print '\t', rsqrd/n_used
            break
        else:
            print n, 'GOOD!! :)'
            template_wvlsol = wsol
            keep_coeffs = coeffs
            keep_peaks_pix = peaks_pix
            keep_peaks_wvl = peaks_wvl
            keep_rsqrd = rsqrd
            keep_n_used = n_used
        n += 1

    print n_used
    if True:
        wsol = lambda x, c=keep_coeffs: polynomial(x, *c)
        fig, ax = plt.subplots()
        ax.set_title('wvlsol, n_peaks = '+str(keep_n_used)+', $res^2$ = '+str(keep_rsqrd)+'$res^2/n_peaks$ = '+str(keep_rsqrd/keep_n_used))
        ax.set_xlabel('Y_Pixel')
        ax.set_ylabel('Wavelength ($\AA$)')
        ax.scatter(keep_peaks_pix, keep_peaks_wvl)
        pixfit = np.linspace(min(pix), max(pix), 10000)
        countsfit = wsol(pixfit)
        ax.plot(pixfit, countsfit)

    wsol = lambda x, c=keep_coeffs: polynomial(x, *c)
    return wsol

def match_peaks(peaks_pix, peaks_wvl, template_wvlsol):
    '''
    A function that attempts to match peaks found in pixel space to known peaks
    in wavelength space.

    ARGUMENTS:
    ----------------------------------------------------------------------------
    peaks_pix: An array of the locations of peaks in pixel space.

    peaks_wvl: An array of the locations of peaks in wavelength space.

    **Note. These two arrays do not need to be the same length. This algorithm
    works best if there are more peaks in peaks_wvl than there are in peaks_pix.

    template_wvlsol: A function that roughly approximates the transformation
    pixel space to wavelength space.

    RETURNS:
    ----------------------------------------------------------------------------
    Two lists; one with pixel positions of peaks and the other with
    corresponding wavelength positions of peaks.
    '''

    #Find optimal linear offset to add to template_wvlsol
    r_sqared = lambda offset: min_res_sqr(peaks_pix, peaks_wvl, lambda p: template_wvlsol(p)+offset)
    offset = minimize(r_sqared, x0=0).x[0]

    #Using template_wvlsol+offset, define an approximate wavelength solution.
    wsol = lambda p: template_wvlsol(p)+offset

    #Using the approximate wavelength solution, find peaks in wavelength space that most nearly match to peaks in pixel space.
    pix = []
    wvl = []
    i = 0
    while i < len(peaks_pix):
        p = peaks_pix[i]
        w = wsol(p)
        diffs = [abs(w-pw) for pw in peaks_wvl]
        nearest_w = peaks_wvl[diffs.index(min(diffs))]
        add = True
        #Ensure that to two pixel peaks are matched to the same wavelength.
        if nearest_w in wvl:
            dist = abs(w-nearest_w)
            other_i = wvl.index(nearest_w)
            other_p = peaks_pix[other_i]
            other_dist = abs(wsol(other_p)-nearest_w)
            if other_dist < dist:
                add = False
            else:
                pix.remove(pix[other_i])
                wvl.remove(wvl[other_i])
        if add:
            pix.append(p)
            wvl.append(nearest_w)
        i += 1

    return np.asarray(pix), np.asarray(wvl)

def min_res_sqr(x, y, func):
    '''
    A function which returns the lowest possible residuals squared
    of a function using two unordered lists x and y

    ARGUMENTS:
    ---------yy-------------------------------------------------------------------
    x: An array of x values.

    y: An array of y values.

    func: A unary function relating x and y.

    **Note. x and y need not be ordered with respect to eachother (y[0] does not
      need to correspond to x[0]). They don't even need to be the same length.**

    RETURNS:
    ----------------------------------------------------------------------------
    min_r_squared: The smallest residuals squared between x and y through func.
                   Obtained by summing the difference squared between func(x[i]) and the
                   nearest y for every value of x.
    '''
    min_r_sqrd = 0
    for xval in x:
        ymod = func(xval)
        r_sqrds = [(ymod-yval)**2 for yval in y]
        min_r_sqrd+=min(r_sqrds)
    return min_r_sqrd

def fit_poly(x, y, n):
    '''
    Fit an n-degree polynomial to the data (x, y).

    ARGUMENTS:
    ----------------------------------------------------------------------------
    x: An array of x values.

    y: An array of y values.

    n: The degree of the fit.

    RETURNS:
    ----------------------------------------------------------------------------
    coeff: An n+1 length array of the coefficients of the best-fit polynomial.
           Starting with the coefficiant of x^n and ending with the coefficient
           of x^0.
    '''
    polynomial = lambda x, *args: sum([coeff*x**power for power,coeff in enumerate(args)])
    x = np.array(x)
    y = np.array(y)
    sort = np.argsort(x)
    x = x[sort]
    y = y[sort]

    slope = (y[-1]-y[0])/(x[-1]-x[0])
    coeff, err = curve_fit(polynomial, x, y, p0=[0, slope]+(n-1)*[0])
    return coeff
