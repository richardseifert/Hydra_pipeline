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

#Get the spacing between fibers in a fits-like object.
#-----------------------------------------------------
# some_fits can be either:
#       fits.hdu.hdulist.HDUList
#       fits.hdu.image.PrimaryHDU
#       np.ndarray or similar
#
# If some_fits has no header or if the header has
# no SLFIB header cards, then be sure to specify
# the number of fibers with num_fibers.
@fitstools.manage_dtype(with_header=True)
def getFiberSpacing(some_fits, num_fibers=None):
    data, header = some_fits
    print header == None, ':D:D:D:DD:SDFVSGBR'
    if (header == None or getFiberNum(header) == 0) and num_fibers == None:
        raise TypeError('Can not determine number of fibers from information provided.')
    num_fibers = getFiberNum(header)
    num_cols = float(len(data[0]))
    fspacing = num_cols/num_fibers
    return fspacing

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

#Fit for the positions of fibers in a fits-like object.
#This hasn't been working because it doesn't properly
# fit the start position of the first fiber.
#Right now this function isn't used anywhere.
def getFiberPositions(f):
    numFibers = getFiberNum(f)
    colAvg = f.data.colAvg()
    fiberPositions = lambda start, num, spacing: [int(round(start + i*spacing)) for i in range(num)]
    fiberSum = lambda x: sum([colAvg[i] for i in fiberPositions(x[0], numFibers, x[1]) if i < len(colAvg)])

    spacing = float(len(colAvg))/numFibers
    spacing_prec = 0.1
    start = 73.5
    start_prec = 0.1
    print -1, start, spacing
    boop = 1
    while start_prec > 0.0001 and spacing_prec > 0.0001:
        startList = np.arange(start-spacing/1.5, start+spacing/1.5, start_prec)
        varystart = [fiberSum([s, spacing]) for s in startList]
        #fig, ax = plt.subplots()
        #ax.set_title('Start '+str(boop))
        #ax.plot(startList, varystart)
        #ax.axvline(start)
        prev_start = start
        start = startList[varystart.index(max(varystart))]
        if abs(start-prev_start) < start_prec:
            start_prec /= 2
            print 'START PRECISION WOAH'
        #ax.axvline(start, ls='--')

        spacingList = np.arange((1.0/2.0)*spacing, 2.0*spacing, spacing_prec)
        varyspacing = [fiberSum([start, s]) for s in spacingList]
        fig, ax = plt.subplots()
        ax.set_title('Spacing '+str(boop))
        ax.plot(spacingList, varyspacing)
        ax.axvline(spacing)
        prev_spacing = spacing
        spacing = spacingList[varyspacing.index(max(varyspacing))]
        if abs(spacing-prev_spacing) < spacing_prec:
            spacing_prec /= 2
            print 'SPACING PRECISION WOAH'
        ax.axvline(spacing, ls='--')
        print boop, start, spacing
        boop += 1
    
    return numFibers, start, spacing

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
def improve_peak_spacing(some_list, peak_list):
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
    
    #Find most common spacing
    spacing = ident_spacing(peak_list)

    #Start at first significant peak
    start_i = 0
    while not is_sig(some_list, peak_list[start_i], 0.3, spacing):
        start_i += 1

    i = start_i
    sig_threshhold = 0.3
    while i+1 < len(peak_list):
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
                print 'WARNING: found multiple broken fibers in a row.'
        else:
            print 'ENTERING PERFECT SPACING'
            i += 1
    return peak_list

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
     
    print 'HERE!!!'
    peaks = get_peaks(col_avgs, n)
    peaks = improve_peak_spacing(col_avgs, peaks)
    print len(peaks)
    print 'LEAVING HERE!!!'
    fig, ax = plt.subplots()
    ax.plot(col_avgs)
    ax.scatter(peaks, [col_avgs[p] for p in peaks], c='green')
    fig, ax = plt.subplots()
    ax.scatter(range(len(peaks)-1), [p2-p1 for p1, p2 in zip(peaks[:-1], peaks[1:])])
    fSpacing = int(round(getFiberSpacing(some_fits)/2))
    print fSpacing, ':D:D:D'

    masks = []
    for i in peaks:
        fSlice = colAvgs[i-fSpacing:i+fSpacing+1]
        base = lambda n: (fSlice[-1]-fSlice[0])/(len(fSlice)-1)*n + fSlice[0]
        fSlice = [fSlice[n] - base(n) for n in range(len(fSlice))]

        area = lambda l, h: sum([float(fSlice[n]) for n in range(l, h+1)])
        total = area(0, len(fSlice)-1)
        ilow = (len(fSlice)-1)/2
        ihigh = (len(fSlice)-1)/2
        currentArea = 0
        while currentArea/total < 0.995 and (ilow > 0 and ihigh < len(fSlice)-1):
            if ihigh >= len(fSlice) or (ilow > 0 and fSlice[ilow-1] >= fSlice[ihigh+1]):
                ilow -= 1
            else:
                ihigh += 1
            currentArea = area(ilow, ihigh)
        masks.append([i + (ilow-(len(fSlice)-1)/2), i + (ihigh-(len(fSlice)-1)/2)])
        ax.axvline(masks[-1][0], color='red')
        ax.axvline(masks[-1][1], color='red')
    return masks

