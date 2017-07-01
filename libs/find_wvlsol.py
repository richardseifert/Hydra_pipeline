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
from extract import extract_counts, optimal_extraction
import itertools

polynomial = lambda x, *args: sum([coeff*x**power for power,coeff in enumerate(args)])

class wvlsolver:
    def __init__(self, comp, fiber_mask, use_fibers, profile_map, fast=False, output=None):
        self.comp = comp
        self.fmask = fiber_mask
        self.fnums = use_fibers
        self.pmap = profile_map
        self.fast = fast
        self.output=output
        self.fibers = {}

        #Load thar line list info.
        master_calib = 'calib/master_calib'
        dat = np.loadtxt(master_calib+'/thar_short.fits')
        line_list_wvl = dat[:,0]
        line_list_counts = dat[:,1]
        #If the table of thar peaks does not exist, make it.
        if not os.path.exists(master_calib+'/thar_peaks.dat'):
            std, l_peak_x, l_peak_y = fit_ngaussian(line_list_wvl, line_list_counts, 70)
            f = open(master_calib+'/thar_peaks.dat', 'w')
            for x, y in zip(l_peak_x, l_peak_y):
                f.write(str(x).ljust(24)+str(y)+'\n')
            f.close()
        thar_peaks = np.loadtxt(master_calib+'/thar_peaks.dat')
        self.linelist = thar_peaks[:,0]

    def solve(self):
        #Load the template wavelength solution.
        master_calib = 'calib/master_calib'
        template_dat = np.loadtxt(master_calib+'/template_wvlsol.dat', delimiter=',')
        p = template_dat[:,2]
        w = template_dat[:,0]
        coeffs = fit_poly(p, w, 3)
        template = lambda x, c=coeffs: polynomial(x, *c)

        def get_template(fnum):
            nearest_fnums = sorted(self.fnums, key=lambda n: abs(fnum-n))
            for n in nearest_fnums:
                if n in self.fibers.keys():
                    return self.fibers[n].get_solution()
            return template

        #use_fibers_high = sorted([fnum for fnum in self.fnums if fnum >= 50])
        #use_fibers_low = sorted([fnum for fnum in self.fnums if fnum < 50], key = lambda x: -x)

        #The template solutions are generated using the central fiber, fnum = 50, so sort fnums
        # starting at 50, ascending to 99, then jumping to 49, and descending to 1.
        sorted_fnums = sorted([fnum for fnum in self.fnums if fnum >= 50]) + sorted([fnum for fnum in self.fnums if fnum < 50], key = lambda x: -x)
        #sorted_fnums = sorted([fnum for fnum in self.fnums if fnum <= 51], key = lambda x: -x)

        for fnum in sorted_fnums:
            if self.output != None:
                self.output.edit_message('Finding wavelength solution for fiber '+str(fnum))
            f_counts = extract_counts(self.comp, self.fmask, fnum)  #WANT TO REPLACE WITH OPTIMAL EXTRACTION SOMEHOW
            f_pix = np.arange(len(f_counts), dtype=np.float64)
            self.fibers[fnum] = fiber_wvlsoler(f_pix, f_counts, get_template(fnum), self.linelist, fast=self.fast)
            self.fibers[fnum].solve()
            if self.output != None:
                self.output.edit_message('fiber '+str(fnum)+' wavelength solution found using '+str(len(self.fibers[fnum].peaks_pix))+' ThAr lines.')


    def get_wvlsol_map(self):
        #Initialize a blank wavelength solution.
        wvlsol_map = np.zeros_like(self.fmask)
        for fnum in self.fnums:
            wsol = self.fibers[fnum].get_solution()

            #Add individual wavelength solution to wvlsol_map
            wsol_arr = wsol(np.arange(len(wvlsol_map)))
            ones_fiber = np.where(self.fmask==fnum, np.ones_like(self.fmask), 0)
            wvlsol_map += np.transpose(np.multiply(np.transpose(ones_fiber), wsol_arr))

        return wvlsol_map

class fiber_wvlsoler:
    def __init__(self, pix, counts, template, linelist, fast=False):
        self.pix = np.array(pix)
        self.counts = np.array(counts)
        self.linelist = linelist
        self.template = template
        self.fast = fast

        #Load thar line list info.
        master_calib = 'calib/master_calib'
        dat = np.loadtxt(master_calib+'/thar_short.fits')
        self.linelist_wvl = dat[:,0]
        self.linelist_counts = dat[:,1]

    def solve(self, npeaks=70):
        #Remove strong cosmic rays.
        l = len(self.pix)
        self.pix, self.counts = remove_cosmics(self.pix, self.counts)
        #print l-len(self.pix), 'COSMIC RAY POINTS FOUND AND REMOVED.'

        #Find peaks in the fiber.
        std, self.pix_peaks_all, self.pix_counts_all = fit_ngaussian(self.pix, self.counts, npeaks, fast=self.fast)
        npeaks = len(self.pix_peaks_all)

        #Sort fiber peaks by their height
        typical_counts = np.median(self.pix_counts_all)
        heights = [-abs(c - typical_counts) for c in self.pix_counts_all]
        self.pix_peaks_all = np.asarray(self.pix_peaks_all)[np.argsort(heights)]
        #print [p for p in self.pix_peaks_all if p > len(self.pix)]

        #Find 5 good peaks for the initial wvlsol
        template_wvlsol = self.template
        for five_peaks_i in sorted(itertools.combinations(list(range(10)), 5), key=lambda s: sum([s_val**3 for s_val in s])):
            use_peaks_pix = [self.pix_peaks_all[i] for i in five_peaks_i]
            peaks_pix, peaks_wvl = match_peaks(use_peaks_pix, self.linelist, template_wvlsol)
            if len(peaks_pix) < 5:
                continue
            coeffs = fit_poly(peaks_pix, peaks_wvl, n=3)
            wsol = lambda x, c=coeffs: polynomial(x, *c)
            rsqrd = min_res_sqr(peaks_pix, peaks_wvl, wsol)

            #self.plot_solution(peaks_pix=peaks_pix, peaks_wvl=peaks_wvl, wsol=wsol, title=str(five_peaks_i)+' '+str(len(peaks_pix))+' peaks, '+str(rsqrd))

            if rsqrd/len(peaks_pix) <= 7e-5:
                break

        n = max(five_peaks_i)+1
        ignore_peaks_pix = [i for i in range(max(five_peaks_i)) if not i in five_peaks_i]
        #print ignore_peaks_pix, 'IGNORE THESE FROM THE GET GO!'


        self.peaks_pix = []

        #print 'n =',n
        while n <= npeaks:
            use_peaks_pix = [self.pix_peaks_all[i] for i in range(n) if not i in ignore_peaks_pix]
            peaks_pix, peaks_wvl = match_peaks(use_peaks_pix, self.linelist, template_wvlsol)
            n_used = len(peaks_pix)
            poly_n = 3 if len(peaks_pix) < 40 else 5
            coeffs = fit_poly(peaks_pix, peaks_wvl, n=poly_n)
            wsol = lambda x, c=coeffs: polynomial(x, *c)
            rsqrd = min_res_sqr(peaks_pix, peaks_wvl, wsol)
            if len(peaks_pix) < len(self.peaks_pix) or rsqrd/n_used > 0.01:
                ignore_peaks_pix.append(n-1)
                #print len(peaks_pix), rsqrd/n_used, 'REJECTED'
            else:
                self.wsol = wsol
                template_wvlsol = wsol
                self.wsol_coeffs = coeffs
                self.peaks_pix = peaks_pix
                self.peaks_wvl = peaks_wvl
                self.rsqrd = rsqrd
                #print len(peaks_pix), rsqrd/n_used, 'ACCEPTED'
            n += 1

        #print 'FINAL USING '+str(len(self.peaks_pix))+' PEAKS'
        #self.plot_solution(title=str(len(self.peaks_pix))+' peaks, '+str(self.rsqrd))
        self.wsol = lambda x, c=self.wsol_coeffs: polynomial(x, *c)

    def plot_solution(self, peaks_pix=None, peaks_wvl=None, counts=None, wsol=None, **kwargs):
        if type(peaks_pix)==type(None):
            peaks_pix = self.peaks_pix
        if type(peaks_wvl)==type(None):
            peaks_wvl = self.peaks_wvl
        if type(counts)==type(None):
            counts = self.counts
        if wsol==None:
            wsol=self.wsol
        fig, ax = plt.subplots()
        if 'title' in kwargs:
            ax.set_title(kwargs['title'])
        ax.scatter(peaks_pix, peaks_wvl, color='blue')
        p = np.linspace(min(peaks_pix), max(peaks_pix), 1000)
        w = wsol(p)
        ax.plot(p, w, color='red')

        wvl = wsol(self.pix)
        fig, ax = plt.subplots()
        if 'title' in kwargs:
            ax.set_title(kwargs['title'])
        ax.plot(self.linelist_wvl, self.linelist_counts, color='red')
        counts_scale=max(self.linelist_counts)/max(self.counts)
        ax.plot(wvl, self.counts*counts_scale, color='blue')
        #ax.scatter(peaks_wvl, [counts[i]*max(self.linelist_counts)/max(self.counts) for i in peaks_pix], color='black')
        for pw in peaks_wvl:
            ax.axvline(x=pw, color='salmon')
        for pp in peaks_pix:
            ax.axvline(x=wsol(pp),color='cornflowerblue')
            ax.scatter(wsol(pp), counts[int(pp)]*counts_scale, color='black')


    def get_solution(self):
        try:
            return self.wsol
        except AttributeError:
            self.solve()
            return self.wsol

        
@manage_dtype(use_args=[0,1], with_header=[0])
def wvlsol(comp, fiber_mask, use_fibers, profile_map, **kwargs):
    comp, comp_header = comp

    #Initialize a blank wavelength solution.
    wvlsol_map = np.zeros_like(fiber_mask)

    #Define path to thar calibration files.
    master_calib = 'calib/master_calib'

    #Load the template wavelength solution.
    template_dat = np.loadtxt(master_calib+'/template_wvlsol.dat', delimiter=',')
    p = template_dat[:,2]
    w = template_dat[:,0]
    coeffs = fit_poly(p, w, 3)
    template = lambda x, c=coeffs: polynomial(x, *c)

    #Load thar line list info.
    dat = np.loadtxt(master_calib+'/thar_short.fits')
    line_list_wvl = dat[:,0]
    line_list_counts = dat[:,1]
    #If the table of thar peaks does not exist, make it.
    if not os.path.exists(master_calib+'/thar_peaks.dat'):
        std, l_peak_x, l_peak_y = fit_ngaussian(line_list_wvl, line_list_counts, 70)
        f = open(master_calib+'/thar_peaks.dat', 'w')
        for x, y in zip(l_peak_x, l_peak_y):
            f.write(str(x).ljust(24)+str(y)+'\n')
        f.close()
    else:
        thar_peaks = np.loadtxt(master_calib+'/thar_peaks.dat')
        linelist = thar_peaks[:,0]

    def f_wvlsol(fnum, template_wvlsol, wvlsol_map=wvlsol_map):
        #Extract comp spectrum in pixel space.
        comp_counts = extract_counts(comp, fiber_mask, fnum)
        comp_pix = np.arange(len(comp_counts), dtype=np.float64)

        #Find wavelength solution for fiber.
        #fig, ax = plt.subplots()
        #ax.plot(comp_pix, comp_counts)
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

def fiber_wvlsol(pix, counts, linelist, starter_wvlsol, npeaks = 33, **kwargs):
    #Find peaks in the fiber.
    std, npeaks_pix, npeaks_counts = fit_ngaussian(pix, counts, npeaks, **kwargs)
    typical_counts = np.median(npeaks_counts)
    diffs = [abs(c - typical_counts) for c in npeaks_counts]
    npeaks_pix = np.asarray(npeaks_pix)[np.argsort(diffs)]
    n = min([5, npeaks])
    template_wvlsol = starter_wvlsol
    ignore_peaks_pix = []
    while n <= npeaks:
        use_peaks_pix = [npeaks_pix[i] for i in range(n) if not i in ignore_peaks_pix]
        peaks_pix, peaks_wvl = match_peaks(use_peaks_pix, linelist, template_wvlsol)
        n_used = len(peaks_pix)
        coeffs = fit_poly(peaks_pix, peaks_wvl, n=3)
        wsol = lambda x, c=coeffs: polynomial(x, *c)
        rsqrd = min_res_sqr(peaks_pix, peaks_wvl, wsol)
        if rsqrd/n_used > 0.01:
            ignore_peaks_pix.append(n-1)
        else:
            template_wvlsol = wsol
            keep_coeffs = coeffs
            keep_peaks_pix = peaks_pix
            keep_peaks_wvl = peaks_wvl
            keep_rsqrd = rsqrd
            keep_n_used = n_used
        n += 1

    wsol = lambda x, c=keep_coeffs: polynomial(x, *c)

    #print keep_coeffs, 'CUBIC FIT'
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
    use_n = min([n+1, len(x)])-1
    #print n, len(x), use_n
    if use_n == 0:
        return [0]*n
    
    polynomial = lambda x, *args: sum([coeff*x**power for power,coeff in enumerate(args)])
    x = np.array(x)
    y = np.array(y)
    sort = np.argsort(x)
    x = x[sort]
    y = y[sort]

    slope = (y[-1]-y[0])/(x[-1]-x[0])
    coeff, err = curve_fit(polynomial, x, y, p0=[0, slope]+(use_n-1)*[0])
    coeff = list(coeff) + [0]*(n-use_n)
    return coeff

def remove_cosmics(x, y, thresh=50):
    keep_i = [i for i in list(range(len(y)))[1:-1] if y[i]/(0.5*(y[i-1]+y[i+1]))<thresh]
    #print [y[i]/(0.5*(y[i-1]+y[i+1])) for i in list(range(len(y)))[1:-1] if not i in keep_i]
    keep_x = [x[i] for i in keep_i]
    keep_y = [y[i] for i in keep_i]
    if y[0]/y[1] < thresh:
        keep_x.insert(0,x[0])
        keep_y.insert(0,y[0])
    if y[-1]/y[-2] < thresh:
        keep_x.append(x[-1])
        keep_y.append(y[-1])
    return np.array(keep_x), np.asarray(keep_y)
