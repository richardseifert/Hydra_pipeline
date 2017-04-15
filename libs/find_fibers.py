import fitstools
from astropy.io import fits
from scipy.optimize import minimize
import numpy as np
import matplotlib.pyplot as plt
from math import floor, ceil
from scipy.optimize import curve_fit


#Helper function to get the number of fibers in a fits object using
# information from the fits header.
#
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
    threshhold = sum(some_list)/len(some_list)

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
@fitstools.manage_dtype()
def fit_fcenter_fwidth(some_fits, fiber_positions, xpos):
    spacing = abs(int(round(ident_spacing(fiber_positions))))
    col_avgs = fitstools.col_avg(some_fits)
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
@fitstools.manage_dtype()
def get_fiber_mask(some_fits, fiber_positions, use_fibers):
    mask = np.zeros_like(some_fits)
    spacing = abs(int(round(ident_spacing(fiber_positions))))
    col_avgs = fitstools.col_avg(some_fits)
    for fnum in use_fibers:
        f_indx = fnum-2
        xpos = fiber_positions[f_indx]
        if not is_sig(col_avgs, xpos, 0.3, spacing):
            fcenter_list = [xpos for i in range(len(some_fits))]
            fwidth = 0
        else:
            fcenter_list, fwidth = fit_fcenter_fwidth(some_fits, fiber_positions, xpos)
        fwidth = 3
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
@fitstools.manage_dtype(use_args=[0], with_header=True)
def find_fibers(some_fits, use_fibers):
    data, header = some_fits
    n = None
    if header != None:
        n = getFiberNum(header)
    fitstools.display(data)
    col_avgs = fitstools.col_avg(data)
    fig, ax = plt.subplots()
    ax.plot(col_avgs)

    peaks = get_peaks(col_avgs, n)
    peaks = improve_peak_spacing(col_avgs, peaks)
    spacing = int(round(ident_spacing(peaks)))
    fig, ax = plt.subplots()
    ax.plot(col_avgs)
    ax.scatter(peaks, [col_avgs[p] for p in peaks], c='green')
    #fig, ax = plt.subplots()
    #ax.scatter(range(len(peaks)-1), [p2-p1 for p1, p2 in zip(peaks[:-1], peaks[1:])])
    peaks = peaks[::-1] #reversed order is the order they appear on the chip.
    mask = get_fiber_mask(data, peaks, use_fibers)
    return mask
