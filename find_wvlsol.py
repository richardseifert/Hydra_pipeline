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


def get_peaks(some_list, pthreshhold):
    is_peak = lambda i: some_list[i-1] <= some_list[i] and some_list[i+1] < some_list[i]

    thresh = np.percentile(some_list, pthreshhold)
    above_threshhold = lambda i, t: some_list[i] > t

    peaks_i = [i for i in range(1, len(some_list)-1) if is_peak(i) and above_threshhold(i, thresh)]
    return peaks_i

def get_peak_center(xlist, ylist, i, prec=0.001):
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
    
    print low, high
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
def wvlsol(comp, fiber_mask, use_fibers):
    #Load the template wavelength solution.
    template_dat = np.loadtxt('template_wvlsol.dat', delimiter=',')
    p = template_dat[:,2]
    w = template_dat[:,0]
    coeffs = fit_poly(p, w, 3)
    template_wvlsol = lambda x, c=coeffs: polynomial(x, *c)
    
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

    #use_fibers = [use_fibers[6]]
    use_fibers = [7, 15, 18, 70, 83, 91, 96, 97, 98]
    for fnum in use_fibers:
        print 'Fiber # '+str(fnum)
        fiber = mask_fits(comp, fiber_mask, fnum)
        comp_counts = row_avg(fiber)
        comp_pix = np.arange(len(comp_counts), dtype=np.float64)
        coeffs = fiber_wvlsol(comp_pix, comp_counts, linelist, template_wvlsol)
        wsol = lambda x, c=coeffs: polynomial(x, *c)

        fig_aft, ax_aft = plt.subplots()
        ax_aft.set_title('Fiber # '+str(fnum))
        ax_aft.set_xlabel('Wavelength ($\AA$)')
        ax_aft.set_ylabel('Counts')
        comp_wvl = [wsol(pix) for pix in comp_pix]
        ax_aft.plot(line_list_wvl, line_list_counts*(max(comp_counts)/max(line_list_counts)), color='red')
        ax_aft.plot(comp_wvl, comp_counts)

def fiber_wvlsol(pix, counts, linelist, template_wvlsol, plot=True):
    std, peaks_pix, peaks_counts = fit_ngaussian(pix, counts, 25)
    peaks_pix, peaks_wvl = match_peaks(peaks_pix, linelist, template_wvlsol)
    print len(peaks_pix), 'lines used.'
    coeffs = fit_poly(peaks_pix, peaks_wvl, n=3)

    if plot:
        wsol = lambda x, c=coeffs: polynomial(x, *c)
        fig, ax = plt.subplots()
        ax.set_title('wvlsol')
        ax.set_xlabel('Y_Pixel')
        ax.set_ylabel('Wavelength ($\AA$)')
        ax.scatter(peaks_pix, peaks_wvl)
        pixfit = np.linspace(min(pix), max(pix), 10000)
        countsfit = wsol(pixfit)
        ax.plot(pixfit, countsfit)

    return coeffs

def fit_poly(x, y, n):
    polynomial = lambda x, *args: sum([coeff*x**power for power,coeff in enumerate(args)])
    x = np.array(x)
    y = np.array(y)
    sort = np.argsort(x)
    x = x[sort]
    y = y[sort]

    slope = (y[-1]-y[0])/(x[-1]-x[0])
    coeff, err = curve_fit(polynomial, x, y, p0=[0, slope]+(n-1)*[0])
    return coeff


def res_sqr(peaks_wvl, peaks_pix, wsol):
    r = 0
    for pp in peaks_pix:
        w = wsol(pp)
        r_sqrds = [(w-pw)**2 for pw in peaks_wvl]
        r+=min(r_sqrds)
    return r
def match_peaks(peaks_pix, peaks_wvl, template_wvlsol):
    r_sqared = lambda offset: res_sqr(peaks_wvl, peaks_pix, lambda p: template_wvlsol(p)+offset)
    offset = minimize(r_sqared, x0=0).x[0]
    print offset, 'OFFSET'
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
    gauss = lambda x: amp*np.exp(-1/2*((x-mu)/(sig))**2)
    return np.asarray([gauss(x_val) for x_val in x])

def make_ngaussian(x, p):
    sig = p[0]
    amp = [p[i] for i in range(len(p)) if i%2==1]
    mu = [p[i] for i in range(1, len(p)) if i%2==0]

    y_model = np.zeros_like(x)
    for a,m in zip(amp, mu):
        y_model = y_model + make_gaussian(x, a, m, sig)

    return y_model

def ngaussian_funct(p, xdata, ydata, fjac=None):
    ymodel = make_ngaussian(xdata, p)
    deviates = [ym-yd for ym,yd in zip(ymodel, ydata)]
    deviates = np.asarray(deviates)
    status = 0
    
    return [status, deviates] #Deviates needs to be a numpy array!!

def find_n_peaks(xdat, ydat, num_peaks):
    #yfit = interp1d(xdat, ydat, kind='cubic')(xdat)
    peak_i_list = [i for i in range(1,len(ydat)-1) if ydat[i] > ydat[i-1] and ydat[i] > ydat[i+1]]
    peak_xvals = np.asarray([xdat[i] for i in peak_i_list])
    peak_yvals = np.asarray([ydat[i] for i in peak_i_list])
    sort_i = np.argsort(-peak_yvals)
    peak_xvals = peak_xvals[sort_i]
    peak_yvals = peak_yvals[sort_i]

    return peak_xvals[:num_peaks], peak_yvals[:num_peaks]

def fit_ngaussian(xdata, ydata, n, plot=True):
    peak_x, peak_y = find_n_peaks(xdata, ydata, n)
    for i in range(len(peak_x)):
        peak_i = np.where(xdata==peak_x[i])[0][0]
        px, py = get_peak_center(xdata, ydata, peak_i)
        peak_x[i] = px
        peak_y[i] = py
    p0 = [1.0] #Initial guess of standard deviation of gaussians.
    for x, y in zip(peak_x, peak_y):
        p0.append(y)
        p0.append(x)

    #m = mpfit(ngaussian_funct, p0, {'xdata':xdata, 'ydata':ydata}, quiet=0)
    #p = m.params
    p = p0

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
