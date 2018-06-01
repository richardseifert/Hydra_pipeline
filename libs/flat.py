from fitstools import manage_dtype, mask_fits, display, row_avg, col_avg
from astropy.io import fits
from scipy.optimize import minimize
import numpy as np
from math import floor, ceil
from scipy.optimize import curve_fit




#Helper function to get the number of fibers in a fits object using
# information from the fits header.
def getFiberNum(some_header):
    num_fibers = 1
    n = 1
    while 'SLFIB'+str(n) in some_header:
        if not 'SCS' in some_header['SLFIB'+str(n)]:
            num_fibers += 1
        n += 1
    num_fibers -= 1
    return num_fibers

#Function that finds significant peaks in an array of data. Indices
# of peaks are returned.
def find_sig_peaks(some_list):
    threshhold = np.nanmean(some_list)

    #Identify significant peaks
    peaks = []
    for i in range(len(col_avgs))[1:-1]:
        if col_avgs[i] > threshhold and colAvgs[i-1] < colAvgs[i] and colAvgs[i+1] < colAvgs[i]:
            peaks.append(i)

    return peaks

#Function that finds the n highest peaks in an array of data.
# Indicies of peaks are returned.
def find_n_peaks(some_list, n):
    is_peak = lambda l, i: l[i-1] < l[i] and l[i] > l[i+1]
    peaks = [i for i in range(len(some_list))[1:-1] if is_peak(some_list, i)]
    n_high_peaks = sorted(peaks, key=lambda i: -some_list[i])[:n]
    return sorted(n_high_peaks)

#Funciton for finding peaks that uses one of the two functions above,
# depending on the arguments given.
def get_peaks(some_list, n=None):
    return find_n_peaks(some_list, n) if n!=None else get_sig_peaks(some_list)

#Function that identifies the most common spacing between items in a list.
# Returns the weighted average of the two most frequenty occuring spacings.
def ident_spacing(p_list):
    p_spacing = [p2-p1 for p1, p2 in zip(p_list[:-1], p_list[1:])]
    p_spacing_unique = list(set(p_spacing))
    p_spacing_freq = [len([s for s in p_spacing if s == spacing]) for spacing in p_spacing_unique]
    sort_freq_i = sorted(range(len(p_spacing_unique)), key=lambda i: -p_spacing_freq[i])

    spacing1 = p_spacing_unique[sort_freq_i[0]]
    freq1 = p_spacing_freq[sort_freq_i[0]]
    spacing2 = p_spacing_unique[sort_freq_i[1]]
    freq2 = p_spacing_freq[sort_freq_i[1]]
    total_freq = freq1+freq2
    spacing = float(freq1*spacing1 + freq2*spacing2)/float(total_freq)
    return spacing

#Function that determines whether or not a peak at index i is significant.
# sig is the percent above the base level that the peak rises. sig=0.3 only
# accepts peaks that are 30% larger than the base level or larger.
def is_sig(alist, i, sig, spacing=3):
    low = int(round(i-spacing))
    if low < 0:
        low = 0
    high = int(round(i+spacing))
    base = min(alist[low:high])
    percent_diff = (alist[i]-base)/base
    return percent_diff >= sig

#Function that modifies a list of peaks based on the knowledge that elements should be roughly evenly spaced.
# The list of peaks should be sorted.
def improve_peak_spacing(some_list, peak_list):
    #Find most common spacing
    spacing = ident_spacing(peak_list)

    #Start at first significant peak
    start_i = 0
    while not is_sig(some_list, peak_list[start_i], 0.3, spacing):
        start_i += 1

    i = start_i
    sig_threshhold = 0.3
    while i+1 < len(peak_list):
        #Update value for peak spacing
        spacing = ident_spacing(peak_list)

        space_to_next = peak_list[i+1] - peak_list[i]
        if space_to_next < floor(spacing):
            if not is_sig(some_list, peak_list[i+1], sig_threshhold, spacing):
                peak_list.remove(peak_list[i+1])
            else:
                #Move to next peak
                i += 1
        elif space_to_next > ceil(spacing):
            n = 1
            mult_spacing = [n*floor(spacing), n*ceil(spacing)]
            is_mult = False
            while not is_mult and mult_spacing[0] < space_to_next:
                if mult_spacing[0] <= space_to_next and space_to_next <= mult_spacing[1]:
                    is_mult = True
                n += 1
                mult_spacing = [n*floor(spacing), n*ceil(spacing)]
            if is_mult:
                n -= 1
                new_positions = [int(round(peak_list[i]+j*(space_to_next/n))) for j in range(1,n)]
                for pos in new_positions:
                    peak_list.insert(i+1, pos)
                    i += 1

                #Move to next peak
                i += 1
            else:
                peak_list.remove(peak_list[i+1])
        else:
            i += 1
    return peak_list

#Function that finds the center of a fiber as a function of y_pixel position.
# The function also returns the width obtained from the column-averaged
# profile of the fiber.
@manage_dtype()
def fit_fcenter_fwidth(some_fits, fiber_positions, xpos):
    spacing = abs(int(round(ident_spacing(fiber_positions))))
    col_avgs = col_avg(some_fits)
    fcenter_list = []
    approx_fcenter, fwidth = find_center_and_width(col_avgs, xpos, spacing)
    for row in some_fits:
        fcenter = list(row[xpos-fwidth:xpos+fwidth+1]).index(max(row[xpos-fwidth:xpos+fwidth+1]))+(xpos-fwidth)
        fcenter_list.append(fcenter)
    center = fit_to_func(lambda x,a,b,c,d,e: a*x**4+b*x**3+c*x**2+d*x+e, range(len(fcenter_list)), fcenter_list)
    fcenter_list = [int(round(center(i))) for i in range(len(fcenter_list))]
    return fcenter_list, fwidth

#Function that produces a mask array containing the positions of fibers in an image.
# 0 indicates a pixel belonging to no fiber.
# n indicates a pixel belonging to the nth fiber for n > 0
@manage_dtype()
def get_fiber_mask(some_fits, fiber_positions, use_fibers):
    mask = np.zeros_like(some_fits)
    spacing = abs(int(round(ident_spacing(fiber_positions))))
    col_avgs = col_avg(some_fits)
    for fnum in use_fibers:
        f_indx = fnum-2
        xpos = fiber_positions[f_indx]
        if not is_sig(col_avgs, xpos, 0.3, spacing):
            fcenter_list = [xpos for i in range(len(some_fits))]
            fwidth = 0
        else:
            fcenter_list, fwidth = fit_fcenter_fwidth(some_fits, fiber_positions, xpos)
        fwidth = 4
        for r in range(len(mask)):
            c = fcenter_list[r]
            mask[r][c] = fnum
            for w in range(fwidth):
                mask[r][c+w] = fnum
                mask[r][c-w] = fnum
    return mask

#Function that fits the function, func, to the data x,y.
# Returns a function, f, with the simple call, y = f(x).
# Best-fit parameters of func are not
# easily obtained afterwards.
def fit_to_func(func, x, y):
    coeff, err = curve_fit(func, x, y)
    fit = lambda args: lambda x: func(x,*args)
    f = fit(coeff)
    return f

#Function for finding the center of a fiber and
# determining the width that contains more than
# 99.5% of the total counts of the fiber.
def find_center_and_width(some_list, pos, rad):
    l = list(some_list)
    low = l[pos-rad:pos].index(min(l[pos-rad:pos]))+(pos-rad)
    high = l[pos:pos+rad+1].index(min(l[pos:pos+rad+1]))+pos
    m = float(l[high]-l[low])/(high-low)
    l_flat_base = [l[i]-l[low]-m*(i-low) for i in range(len(l))]
    area = lambda p1, p2: float(sum(l_flat_base[p1:p2+1]))
    total = area(low, high)
    if total == 0:
        return
    f_low = pos
    f_high = pos
    while area(f_low, f_high)/total < 0.995:
        right = f_high - 1
        left = f_low + 1
        if l[right] > l[left] and right <= high:
            f_high += 1
        else:
            f_low -= 1
    lslice = l[f_low:f_high+1]
    center = lslice.index(max(lslice))+f_low
    width = f_high-center
    return center, width

#Function that takes a fits-like argument for an image with evenly spaced
# fibers that are relatively bright, such as a flat field, and identifies
# the positions of each fiber, returning a numbered mask array.
@manage_dtype(use_args=[0], with_header=True)
def find_fibers(some_fits, use_fibers):
    data, header = some_fits
    n = None
    if header != None:
        n = getFiberNum(header)
    col_avgs = col_avg(data)

    peaks = get_peaks(col_avgs, n)
    peaks = improve_peak_spacing(col_avgs, peaks)
    spacing = int(round(ident_spacing(peaks)))
    peaks = peaks[::-1] #reversed order is the order they appear on the chip.
    mask = get_fiber_mask(data, peaks, use_fibers)
    return mask

@manage_dtype()
def make_fiber_profile_map(img, fmask):
    fiber_profile_map = np.zeros_like(img)
    fnums = list({n for row in fmask for n in row if n != 0})
    for fnum in fnums:
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
        #print 'RuntimeWarning: Could not do gaussian fit.'
        coeff = p0
    yfit = func(x, *coeff)

    return yfit


@manage_dtype()
def make_throughput_map(fmask, mflat):
    #Generate a blank throughput map
    throughput_map = np.ones_like(fmask)

    #Obtain a list of the used fibers from the fiber mask.
    fnums = [int(n) for n in set(fmask.flatten()) if n != 0]

    #Correct for the profile of each flat fiber
    for fnum in fnums:
        flat_spec = row_avg(mask_fits(mflat, fmask, maskval=fnum))
        if np.nanmin(flat_spec) < 0:
            print 'WARNING: Likely weak or broken fiber! Flat fiber contains negative values.'
        medianf = np.nanmedian(flat_spec)
        for i,counts in enumerate(flat_spec):
            throughput_map[i][np.where(fmask[i]==fnum)] *= flat_spec[i]/medianf


    #Correct for fiber-to-fiber throughput.
    profile_corrected_mflat = mflat/throughput_map
    flat_fiber_avg_vals = []
    for fnum in fnums:
        flat_spec = row_avg(mask_fits(profile_corrected_mflat, fmask, maskval=fnum))
        avg_val = np.nanmean(flat_spec)
        flat_fiber_avg_vals.append(avg_val)
    medianf = np.nanmedian(flat_fiber_avg_vals)
    for f,fnum in zip(flat_fiber_avg_vals,fnums):
        throughput_map[np.where(fmask==fnum)] *= f/medianf

    return throughput_map
