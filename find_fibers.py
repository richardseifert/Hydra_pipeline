#print 'findFibers importing fitsFile'
#import fitsFile
print 'findFibers importing fitstools'
import fitstools
print 'findFibers importing astropy.io.fits'
from astropy.io import fits
print 'findFibers importing fitsFile'
from scipy.optimize import minimize
print 'findFibers importing fitsFile'
import numpy as np
print 'findFibers importing matplotlib.pyplot as plt'
import matplotlib.pyplot as plt
print 'findFibers finished importing'
from math import floor, ceil
from scipy.optimize import curve_fit

#Get the number of fibers in a fits object using
#information from the fits header.
#-----------------------------------------------
# some_header can be either:
#       fits.header.Header
#       python dict or similar
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

def find_sig_peaks(some_list):
    threshhold = sum(some_list)/len(some_list)
    
    #Identify significant peaks
    peaks = []
    for i in range(len(colAvgs))[1:-1]:
        if colAvgs[i] > threshhold and colAvgs[i-1] < colAvgs[i] and colAvgs[i+1] < colAvgs[i]:
            peaks.append(i)

    return peaks
def find_n_peaks(some_list, n):
    is_peak = lambda l, i: l[i-1] < l[i] and l[i] > l[i+1]
    peaks = [i for i in range(len(some_list))[1:-1] if is_peak(some_list, i)]
    n_high_peaks = sorted(peaks, key=lambda i: -some_list[i])[:n]
    return sorted(n_high_peaks)
def get_peaks(some_list, n=None):
    return find_n_peaks(some_list, n) if n!=None else get_sig_peaks(some_list)
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

def is_sig(alist, i, sig, spacing=3):
    low = int(round(i-spacing))
    if low < 0:
        low = 0
    high = int(round(i+spacing))

    base = min(alist[low:high])
    percent_diff = (alist[i]-base)/base
    return percent_diff >= sig

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

@fitstools.manage_dtype()
def fit_fcenter_fwidth(some_fits, fiber_positions, xpos):
    spacing = int(round(ident_spacing(fiber_positions)))
    col_avgs = fitstools.colAvg(some_fits)
    fcenter_list = []
    approx_fcenter, fwidth = find_center_and_width(col_avgs, xpos, spacing)
    print xpos, fwidth
    for row in some_fits:
        fcenter = list(row[xpos-fwidth:xpos+fwidth+1]).index(max(row[xpos-fwidth:xpos+fwidth+1]))+(xpos-fwidth)
        fcenter_list.append(fcenter)
    center = fit_to_func(lambda x,a,b,c,d,e: a*x**4+b*x**3+c*x**2+d*x+e, range(len(fcenter_list)), fcenter_list)
    fcenter_list = [int(round(center(i))) for i in range(len(fcenter_list))]
    fig, ax = plt.subplots()
    ax.scatter(range(len(fcenter_list)), fcenter_list)
    return fcenter_list, fwidth

@fitstools.manage_dtype()
def get_fiber_mask(some_fits, fiber_positions):
    mask = np.zeros_like(some_fits)
    spacing = int(round(ident_spacing(fiber_positions)))
    col_avgs = fitstools.colAvg(some_fits)
    for xpos in fiber_positions:
        fnum = list(fiber_positions).index(xpos)+1
        if not is_sig(col_avgs, xpos, 0.3, spacing):
            print 'insig', xpos
            fcenter_list = [xpos for i in range(len(some_fits))]
            fwidth = 0
        else:
            fcenter_list, fwidth = fit_fcenter_fwidth(some_fits, fiber_positions, xpos)
        for r in range(len(mask)):
            c = fcenter_list[r]
            mask[r][c] = fnum
            for w in range(fwidth):
                mask[r][c+w] = fnum
                mask[r][c-w] = fnum
    fitstools.display(mask)
def fit_to_func(func, x, y):
    coeff, err = curve_fit(func, x, y)
    fit = lambda args: lambda x: func(x,*args)
    return fit(coeff)

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

@fitstools.manage_dtype(with_header=True)
def findFibers(some_fits):
    data, header = some_fits
    n = None
    if header != None:
        n = getFiberNum(header)
    fitstools.display(data)
    col_avgs = fitstools.colAvg(data)
    fig, ax = plt.subplots()
    ax.plot(col_avgs)
     
    peaks = get_peaks(col_avgs, n)
    peaks = improve_peak_spacing(col_avgs, peaks)
    spacing = int(round(ident_spacing(peaks)))
    print len(peaks)
    fig, ax = plt.subplots()
    ax.plot(col_avgs)
    ax.scatter(peaks, [col_avgs[p] for p in peaks], c='green')
    #fig, ax = plt.subplots()
    #ax.scatter(range(len(peaks)-1), [p2-p1 for p1, p2 in zip(peaks[:-1], peaks[1:])])

    get_fiber_mask(data, peaks)
    #masks = []
    #for i in peaks:
    #    fSlice = col_avgs[i-spacing:i+spacing+1]
    #    base = lambda n: (fSlice[-1]-fSlice[0])/(len(fSlice)-1)*n + fSlice[0]
    #    fSlice = [fSlice[n] - base(n) for n in range(len(fSlice))]

    #    area = lambda l, h: sum([float(fSlice[n]) for n in range(l, h+1)])
    #    total = area(0, len(fSlice)-1)
    #    ilow = (len(fSlice)-1)/2
    #    ihigh = (len(fSlice)-1)/2
    #    currentArea = 0
    #    while currentArea/total < 0.995 and (ilow > 0 and ihigh < len(fSlice)-1):
    #        if ihigh >= len(fSlice) or (ilow > 0 and fSlice[ilow-1] >= fSlice[ihigh+1]):
    #            ilow -= 1
    #        else:
    #            ihigh += 1
    #        currentArea = area(ilow, ihigh)
    #    masks.append([i + (ilow-(len(fSlice)-1)/2), i + (ihigh-(len(fSlice)-1)/2)])
    #    ax.axvline(masks[-1][0], color='red')
    #    ax.axvline(masks[-1][1], color='red')
    #return masks

