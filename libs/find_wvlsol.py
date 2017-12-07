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
    def __init__(self, comp, fiber_mask, use_fibers, profile_map, fast=False, output=None, plotter=None):
        self.comp = comp
        self.fmask = fiber_mask
        self.fnums = use_fibers
        self.pmap = profile_map
        self.fast = fast
        self.output=output
        self.plotter=plotter
        self.fibers = {}

        #Load base template wavelength solution.
        self.load_base_template()

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

    def set_path(self, new_path):
        self.savepath = new_path

    def load_base_template(self):
        #Load the template wavelength solution.
        master_calib = 'calib/master_calib'
        template_dat = np.loadtxt(master_calib+'/template_wvlsol.dat', delimiter=',')
        p = template_dat[:,2]
        w = template_dat[:,0]
        coeffs = fit_poly(p, w, 3)
        self.base_template = lambda x, c=coeffs: polynomial(x, *c)

    def get_template(self, fnum, valid_fnums=None):
        if valid_fnums == None:
            valid_fnums = self.fnums
        nearest_fnums = sorted(self.fnums, key=lambda n: abs(fnum-n))
        for n in nearest_fnums:
            if n in self.fibers.keys() and n in valid_fnums:
                return self.fibers[n].get_solution()
        return self.base_template



    def remove_cosmics(self, tol=5):
        pix = {fnum:self.fibers[fnum].get_pix() for fnum in self.fibers.keys()}
        counts = {fnum:self.fibers[fnum].get_counts() for fnum in self.fibers.keys()}

        #Shift fibers to be lined up with center fiber.
        center_fnum = sorted(self.fibers.keys(), key=lambda fnum: abs(fnum-50))[0]
        shifts = {}
        for fnum in self.fibers.keys():
            corr = np.correlate(counts[center_fnum],counts[fnum], 'full')
            shifts[fnum] = np.arange(-len(pix[fnum])+1, len(pix[fnum])+1)[np.argmax(corr)]
        
        master_pix = np.arange(min([min(shifts.values()),0]), len(counts[center_fnum])+max(shifts.values()))
        length = len(master_pix)
        min_pix = min(master_pix)
        max_pix = max(master_pix)
        for fnum in self.fibers.keys():
            i = -min_pix+shifts[fnum]
            full_pix = np.NAN * np.zeros_like(master_pix)
            full_pix[i:i+len(pix[fnum])] = pix[fnum]
            pix[fnum] = full_pix
            full_counts = np.NAN * np.zeros_like(master_pix)
            full_counts[i:i+len(counts[fnum])] = counts[fnum]
            counts[fnum] = full_counts
        count_medians = np.nanmedian(np.asarray(counts.values()), axis=0)
        count_iqrs = np.subtract(*np.nanpercentile(np.asarray(counts.values()), [75, 25], axis=0))

        self.plotter.clear()
        self.plotter.set_ylabel('Counts')
        self.plotter.set_xlabel('Pixels')
        self.plotter.line(master_pix, count_medians, color='red')
        for fnum in self.fibers.keys():
            self.plotter.line(master_pix, counts[fnum])
        self.plotter.fill_between(master_pix, count_medians-tol*count_iqrs, count_medians+tol*count_iqrs, fill_alpha=0.2, line_alpha=0.2)
        self.plotter.save('cosmics_test.html')

        for fnum in self.fibers.keys():
            mask = np.logical_not(np.isnan(counts[fnum])) & (counts[fnum] > count_medians-tol*count_iqrs) & (counts[fnum] < count_medians+tol*count_iqrs)
            counts[fnum] = counts[fnum][mask]
            pix[fnum] = pix[fnum][mask]


            self.fibers[fnum].set_pix(pix[fnum])
            self.fibers[fnum].set_counts(counts[fnum])

    def solve(self):
        #The template solutions are generated using the central fiber, fnum = 50, so sort fnums
        # starting at 50, ascending to 99, then jumping to 49, and descending to 1.
        sorted_fnums = sorted([fnum for fnum in self.fnums if fnum >= 50]) + sorted([fnum for fnum in self.fnums if fnum < 50], key = lambda x: -x)
        #sorted_fnums = sorted([fnum for fnum in self.fnums if fnum <= 51], key = lambda x: -x)

        #Extract ThAr spectrum for each fiber.
        for fnum in self.fnums:
            f_counts = extract_counts(self.comp, self.fmask, fnum)  #WANT TO REPLACE WITH OPTIMAL EXTRACTION SOMEHOW
            f_pix = np.arange(len(f_counts), dtype=np.float64)
            self.fibers[fnum] = fiber_wvlsoler(f_pix, f_counts, self.linelist, fast=self.fast, plotter=self.plotter)

        #Find and remove cosmic rays.
        self.remove_cosmics()

        good_fiber_wvlsols = []
        bad_fiber_wvlsols = []
        for fnum in sorted_fnums:
            if self.output != None:
                self.output.edit_message('Finding wavelength solution for fiber '+str(fnum))
            #f_counts = extract_counts(self.comp, self.fmask, fnum)  #WANT TO REPLACE WITH OPTIMAL EXTRACTION SOMEHOW
            #f_pix = np.arange(len(f_counts), dtype=np.float64)
            #self.fibers[fnum] = fiber_wvlsoler(f_pix, f_counts, self.linelist, self.get_template(fnum, good_fiber_wvlsols), fast=self.fast, plotter=self.plotter)
            self.fibers[fnum].set_template(self.get_template(fnum, good_fiber_wvlsols))
            self.fibers[fnum].solve(polynomial_plotname='F'+str(fnum)+'_polynomial.html', wvlsol_plotname='F'+str(fnum)+'_wvlsol.html')

            #Check how many peaks were used in the fit to determine if it's good or not.
            if len(self.fibers[fnum].peaks_pix) >= 26:
                good_fiber_wvlsols.append(fnum)
            elif self.output != None:
                bad_fiber_wvlsols.append(fnum)
                self.output.edit_message('Bad solution found for fiber '+str(fnum)+'.')
            try:
                #Keep an updating record of which fibers give good solutions and which don't.
                f = open(self.savepath, 'w')
                f.write(','.join([str(fn) for fn in good_fiber_wvlsols])+'\n')
                f.write(','.join([str(fn) for fn in bad_fiber_wvlsols])+'\n')
                f.close()
            except (AttributeError, TypeError) as e:
                pass

            if self.output != None:
                self.output.edit_message('fiber '+str(fnum)+' wavelength solution found using '+str(len(self.fibers[fnum].peaks_pix))+' ThAr lines.')

    def improve(self):
        #Load the good and bad wavelength solutions from initial call to solve().
        f = open(self.savepath)
        lines = f.read().split('\n')
        f.close()
        good_fiber_wvlsols = [int(fnum) for fnum in filter(None, lines[0].split(','))]
        bad_fiber_wvlsols = [int(fnum) for fnum in filter(None, lines[1].split(','))]

        self.plotter.clear()
        self.plotter.set_xlabel('Pixel')
        self.plotter.set_ylabel('Counts')
        for fnum in bad_fiber_wvlsols:
            #Sort good fibers by their closeness to fnum.
            sorted_good_fnums = sorted(good_fiber_wvlsols, key=lambda n: abs(n-fnum))

            f_counts = extract_counts(self.comp, self.fmask, fnum)  #WANT TO REPLACE WITH OPTIMAL EXTRACTION SOMEHOW
            f_pix = np.arange(len(f_counts), dtype=np.float64)

            self.plotter.clear()
            self.plotter.line(*remove_cosmics(f_pix, f_counts), color='blue')
            for gfnum in sorted_good_fnums:
                gf_counts = extract_counts(self.comp, self.fmask, gfnum)  #WANT TO REPLACE WITH OPTIMAL EXTRACTION SOMEHOW
                gf_pix = np.arange(len(f_counts), dtype=np.float64)

                corr = np.correlate(f_counts, gf_counts, 'full')
                shift = np.arange(-len(f_pix)+1, len(f_pix)+1)[np.argmax(corr)]
                self.plotter.line(*remove_cosmics(gf_pix+shift, gf_counts), color='red')
            self.plotter.save('wvlsol_improve_F'+str(fnum)+'.html')
            self.plotter.clear()
            self.plotter.set_title('best value: '+str(shift))
            self.plotter.set_ylabel('corr')
            self.plotter.set_xlabel('offset')
            self.plotter.line(np.arange(-len(f_pix)+1, len(f_pix)+1), corr)
            self.plotter.save('corr_test.html')



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
    def get_fiber_npeaks(self):
        return {fnum:self.fibers[fnum].get_npeaks for fnum in self.fnums}

class fiber_wvlsoler:
    def __init__(self, pix, counts, linelist, template=None, fast=False, plotter=None):
        self.pix = np.array(pix)
        self.counts = np.array(counts)
        self.linelist = linelist
        self.template = template
        self.fast = fast
        self.plotter = plotter

        #Load thar line list info.
        master_calib = 'calib/master_calib'
        dat = np.loadtxt(master_calib+'/thar_short.fits')
        self.linelist_wvl = dat[:,0]
        self.linelist_counts = dat[:,1]

    def get_pix(self):
        return self.pix
    def get_counts(self):
        return self.counts
    def set_pix(self, new_pix):
        self.pix = new_pix
    def set_counts(self, new_counts):
        self.counts = new_counts

    def set_template(self, new_template):
        self.template = new_template

    def solve(self, npeaks=70, **kwargs):
        #Find peaks in the fiber.
        std, self.pix_peaks_all, self.pix_counts_all = fit_ngaussian(self.pix, self.counts, npeaks, fast=self.fast)

        #Sort fiber peaks by their height
        typical_counts = np.median(self.pix_counts_all)
        heights = [-abs(c - typical_counts) for c in self.pix_counts_all]
        self.pix_peaks_all = np.asarray(self.pix_peaks_all)[np.argsort(heights)]

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

            #self.plot_solution(peaks_pix=peaks_pix, peaks_wvl=peaks_wvl, wsol=wsol, title=str(five_peaks_i)+' '+str(len(peaks_pix))+' peaks, '+str(rsqrd), **kwargs)

            if rsqrd/len(peaks_pix) <= 7e-5:
                break

        n = max(five_peaks_i)+1
        ignore_peaks_pix = [i for i in range(max(five_peaks_i)) if not i in five_peaks_i]
        #print ignore_peaks_pix, 'IGNORE THESE FROM THE GET GO!'


        self.peaks_pix = []

        npeaks = min([npeaks, len(self.pix_peaks_all)])
        while n < npeaks:
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
        self.plot_solution(title=str(len(self.peaks_pix))+' peaks, '+str(self.rsqrd), **kwargs)
        self.wsol = lambda x, c=self.wsol_coeffs: polynomial(x, *c)

    def plot_solution(self, peaks_pix=None, peaks_wvl=None, counts=None, wsol=None, polynomial_plotname='polynomial.pdf', wvlsol_plotname='wvlsol.pdf', **kwargs):
        if type(peaks_pix)==type(None):
            peaks_pix = self.peaks_pix
        if type(peaks_wvl)==type(None):
            peaks_wvl = self.peaks_wvl
        if type(counts)==type(None):
            counts = self.counts
        if wsol==None:
            wsol=self.wsol
        p = np.linspace(min(peaks_pix), max(peaks_pix), 1000)
        w = wsol(p)

        #Generate plot of polynomial fit.
        self.plotter.clear()
        if 'title' in kwargs:
            self.plotter.set_title(kwargs['title'])
        self.plotter.scatter(peaks_pix, peaks_wvl, color='blue')
        self.plotter.line(p, w, color='red')
        self.plotter.save(polynomial_plotname)

        #Generate plot of wavelength solution.
        wvl = wsol(self.pix)
        self.plotter.clear()
        if 'title' in kwargs:
            self.plotter.set_title(kwargs['title'])
        counts_scale=np.max((self.counts))/np.max((self.linelist_counts))
        self.plotter.line(wvl, self.counts, color='blue')
        self.plotter.line(self.linelist_wvl, counts_scale*self.linelist_counts, color='red')
        print max(counts_scale*self.linelist_counts), max(self.counts)
        h1 = 1.05*max([max(counts_scale*self.linelist_counts), max(self.counts)])
        h2 = 1.05*h1
        for pw in peaks_wvl:
            #print pw, h1, h2
            self.plotter.line([pw, pw], [h1, h2], color='red')
        #print
        for pp in peaks_pix:
            #print wsol(pp), h1, h2
            self.plotter.line([wsol(pp), wsol(pp)], [h1, h2], color='blue')
        self.plotter.save(wvlsol_plotname)


    def get_solution(self):
        try:
            return self.wsol
        except AttributeError:
            self.solve()
            return self.wsol
    def get_npeaks(self):
        try:
            return len(self.peaks_wvl)
        except:
            return 0

        
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
    '''
    keep_i = []
    prev_i = 0
    for i in range(len(y))[1:]:
        if y[i]/y[prev_i] < thresh:
            keep_i.append(i)
            prev_i = i
    '''
        
    keep_i = [i for i in list(range(len(y)))[1:-1] if y[i]/(0.5*(y[i-1]+y[i+1]))]
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
