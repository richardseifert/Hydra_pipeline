import os
from astropy.io import fits
from astropy import wcs
from fitstools import manage_dtype, mask_fits, row_avg
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import minimize
import numpy as np
from scipy.optimize import curve_fit
from mpfit import mpfit


def get_peak_center(xlist, ylist, i, prec=0.001):
    '''
    Use a cubic spline to approximate center of a peak. Given a list of x valies
    and a list of y values, this function returns the x value corresponding to
    the peak in y near the index i.

    ARGUMENTS:
    ----------------------------------------------------------------------------
    xlist: An array of x values

    ylist: An array of y values

    i: An index of xlist and ylist that is near the desired peak.
    
    prec: Optional. The precision of the result.

    RETURNS:
    ----------------------------------------------------------------------------
    
    center_x: The x value corresponding to the peak y value in the region near
              the index i.

    center_y: The height of this peak.
    '''

    #Take the region of xlist and ylist surrounding the peak at index i
    low = i-1
    while low-1 >= 0 and ylist[low] > ylist[low-1]:
        low -= 1
    high = i+1
    while high+1 < len(ylist) and ylist[high] > ylist[high+1]:
        high += 1

    while high-low < 4:
        low -= 1
        high += 1
    
    region_x = xlist[low:high+1]
    region_y = ylist[low:high+1]

    #Fit a cubic spline to the peak
    peak = interp1d(region_x, region_y, kind='cubic')
    xfit = np.arange(min(region_x)+prec/2, max(region_x)-prec/2, prec)
    yfit = peak(xfit)

    if False:
        fig, ax = plt.subplots()
        ax.scatter(region_x, region_y)
        ax.plot(xfit, yfit, color='red')
    
    #Find the peak center from spline fit.
    center_x = xfit[list(yfit).index(max(yfit))]

    return center_x, max(yfit)


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

    def fwvlsol(fnum, template_wvlsol, wvlsol_map = wvlsol_map):
        print type(template_wvlsol), ':D:D:D:D:'
        print 'Fiber '+str(fnum)+':'
        fiber = mask_fits(comp, fiber_mask, fnum)
        comp_counts = row_avg(fiber)
        comp_pix = np.arange(len(comp_counts), dtype=np.float64)
        coeffs = fiber_wvlsol(comp_pix, comp_counts, linelist, template_wvlsol, **kwargs)
        wsol = lambda x, c=coeffs: polynomial(x, *c)

        wsol_arr = wsol(np.arange(len(wvlsol_map)))
        ones_fiber = np.where(fiber_mask==fnum, np.ones_like(fiber_mask), 0)
        wvlsol_map += np.transpose(np.multiply(np.transpose(ones_fiber), wsol_arr))

        if False:
            fig_aft, ax_aft = plt.subplots()
            ax_aft.set_title('Fiber # '+str(fnum))
            ax_aft.set_xlabel('Wavelength ($\AA$)')
            ax_aft.set_ylabel('Counts')
            comp_wvl = [wsol(pix) for pix in comp_pix]
            ax_aft.plot(line_list_wvl, line_list_counts*(max(comp_counts)/max(line_list_counts)), color='red')
            ax_aft.plot(comp_wvl, comp_counts)
        return wsol

    #The template solution was generated using fiber 50, so when generating wvlsols, start
    # at fiber 50 and go up, then start at fiber 49 and go down.
    use_fibers_high = sorted([fnum for fnum in use_fibers if fnum >= 50])
    use_fibers_low = sorted([fnum for fnum in use_fibers if fnum < 50], key = lambda x: -x)

    wsol = template
    for fnum in use_fibers_high:
        wsol = fwvlsol(fnum, wsol)
    wsol = template
    for fnum in use_fibers_low:
        wsol = fwvlsol(fnum, wsol)
    return wvlsol_map

def fiber_wvlsol(pix, counts, linelist, starter_wvlsol, npeaks = 25, **kwargs):

    ##REMOVE COSMIC RAYS EVENTUALLY
    #Find peaks in the fiber.
    std, npeaks_pix, npeaks_counts = fit_ngaussian(pix, counts, npeaks, **kwargs)
    start_n = min([5, npeaks])
    n = start_n # I'm thinking of trying something fancy where I test
                    #solutions using different numbers of lines too see which
                    #is the best.
    #n = npeaks #This line is temperary until thing aboove works.
    template_wvlsol = starter_wvlsol
    while n <= npeaks:
        #print 'Finding wvlsol with '+str(n)+' peaks.'
        peaks_pix, peaks_wvl = match_peaks(npeaks_pix[:n], linelist, template_wvlsol)
        n_used = len(peaks_pix)
        #print len(peaks_pix), 'lines used out of '+str(n)+'.'
        coeffs = fit_poly(peaks_pix, peaks_wvl, n=3)
        wsol = lambda x, c=coeffs: polynomial(x, *c)
        rsqrd = min_res_sqr(peaks_pix, peaks_wvl, wsol)
        #print 'Lowest chi-squared: '+str(rsqrd/len(peaks_pix))
        if rsqrd/n_used > 0.01:
            #print n, 'BAD!!! :('
            break
        else:
            #print n, 'GOOD!! :)'
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
    return keep_coeffs

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
    r_sqared = lambda offset: min_res_sqr(peaks_pix, peaks_wvl, lambda p: template_wvlsol(p)+offset)
    offset = minimize(r_sqared, x0=0).x[0]
    #print offset, 'OFFSET'
    wsol = lambda p: template_wvlsol(p)+offset

    pix = []
    wvl = []
    i = 0
    while i < len(peaks_pix):
        p = peaks_pix[i]
        w = wsol(p)
        diffs = [abs(w-pw) for pw in peaks_wvl]
        nearest_w = peaks_wvl[diffs.index(min(diffs))]
        add = True
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
    

def make_gaussian(x, amp, mu, sig):
    '''
    A function that returns a one-dimensional gaussian of a given mean,
    standard deviation, and amplitude over a given domain.
    

    ARGUMENTS:
    ----------------------------------------------------------------------------
    x: An array of x values for the 1D gaussian.

    amp: The amplitude of the gaussian.

    mu: The mean of the gaussian.

    sig: The standard deviation of the gaussian.

    RETURNS:
    ----------------------------------------------------------------------------
    An array of y values from the gaussian corresponding to the x values given
    in x.
    '''
    gauss = lambda x: amp*np.exp(-1/2*((x-mu)/(sig))**2)
    return np.asarray([gauss(x_val) for x_val in x])

def make_ngaussian(x, p):
    '''
    A funciton the returns n one-dimensional gaussians of a given standard
    deviation and given means and amplitudes over a given domain.

    ARGUMENTS:
    ----------------------------------------------------------------------------
    x: An array of x values for the gaussians.

    p: An array of gaussian parameters:
        p[0]      - The single standard deviation for all gaussians.
        p[odd_i]  - The amplitudes of each gaussian.
        p[even_i] - The means of each gaussian.

        p = [std, amp1, mean1, amp2, mean2, amp3, mean3, ... , ampn, meann]

    RETURNS:
    ----------------------------------------------------------------------------
    An array of y values attained from summing all of the gaussians at each of
    the corresponding x values.
    '''
    sig = p[0]
    amp = [p[i] for i in range(len(p)) if i%2==1]
    mu = [p[i] for i in range(1, len(p)) if i%2==0]

    y_model = np.zeros_like(x)
    for a,m in zip(amp, mu):
        y_model = y_model + make_gaussian(x, a, m, sig)

    return y_model

def ngaussian_funct(p, xdata, ydata, fjac=None):
    '''
    A function that mpfit can digest which generates ngaussians when fitting
    with mpfit.

    ARGUMENTS:
    ----------------------------------------------------------------------------
    p: The same array of gaussian arguments that make_ngaussian accepts.

    xdata: An array of x values for the data being fit.

    ydata: An array of y values for the data being fit.

    fjac: Something that mpfit needs, but is never used.

    RETURNS:
    ----------------------------------------------------------------------------
    A status (always success) and an array of "deviates" (residuals) between the
    data and the ngaussian that mpfit uses when fitting.
    '''
    ymodel = make_ngaussian(xdata, p)
    deviates = [ym-yd for ym,yd in zip(ymodel, ydata)]
    deviates = np.asarray(deviates)
    status = 0
    
    return [status, deviates] #Deviates needs to be a numpy array!!

def find_n_peaks(xdata, ydata, num_peaks):
    '''
    A function that finds a specified number of peaks in one-dimensional data.
    Nothing fancy. A peak is defined by:
                ydata[i] > ydata[i-1] and ydata[i] > ydata[i+1]

    ARGUMENTS:
    ----------------------------------------------------------------------------
    xdata: An array of x values.

    ydata: An array of y values.

    num_peaks: The desired number of peaks to find.
    '''
    peak_i_list = [i for i in range(1,len(ydata)-1) if ydata[i] > ydata[i-1] and ydata[i] > ydata[i+1]]
    peak_xvals = np.asarray([xdata[i] for i in peak_i_list])
    peak_yvals = np.asarray([ydata[i] for i in peak_i_list])
    sort_i = np.argsort(-peak_yvals)
    peak_xvals = peak_xvals[sort_i]
    peak_yvals = peak_yvals[sort_i]

    return peak_xvals[:num_peaks], peak_yvals[:num_peaks]

def fit_ngaussian(xdata, ydata, n, fast=False, plot=False):
    '''
    A function that fits n gaussians to some data. Data can be fit quickly by
    only relying on a cubic spline to find peak centers or data can be fit more
    accurately with mpfit.

    ARGUMENTS:
    ----------------------------------------------------------------------------
    xdata: An array of x values.

    ydata: An array of y values.

    n: The number of peaks to fit.

    fast: boolean. True for fast method, False for accurate method. Default is
          False.
          
    plot: Boolean of whether or not to plot things.
    '''
    peak_x, peak_y = find_n_peaks(xdata, ydata, n)
    for i in range(len(peak_x)):
        peak_i = np.where(xdata==peak_x[i])[0][0]
        px, py = get_peak_center(xdata, ydata, peak_i)
        peak_x[i] = px
        peak_y[i] = py
    p0 = [1.0] #Initial guess of standard deviation of gaussians.
               # Fit this initial standard deviation in the future.
    for x, y in zip(peak_x, peak_y):
        p0.append(y)
        p0.append(x)

    if fast:
        p = p0
    else:
        m = mpfit(ngaussian_funct, p0, {'xdata':xdata, 'ydata':ydata}, quiet=1)
        p = m.params

    if plot:
        fig, ax = plt.subplots()
        ax.scatter(xdata, ydata)
        fitX = np.linspace(min(xdata), max(xdata), 10000)
        init_fitY = make_ngaussian(fitX, p0)
        ax.plot(fitX, init_fitY, color='green')
        fitY = make_ngaussian(fitX, p)
        ax.plot(fitX, fitY, color='red')
        
    std = p[0]
    peak_y_list = [p[i] for i in range(1, len(p)) if i%2 == 1]
    peak_x_list = [p[i] for i in range(1, len(p)) if i%2 == 0]
    yfit = make_ngaussian(xdata, p)

    return std, peak_x_list, peak_y_list
