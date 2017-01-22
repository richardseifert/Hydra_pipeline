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

def findFibers(some_fits):
    if type(some_fits) == fits.hdu.hdulist.HDUList:
        data = some_fits[0].data
    elif type(some_fits) == fits.hdu.image.PrimaryHDU:
        data = some_fits.data
    else:
        data = some_fits

    fitstools.display(some_fits)
    colAvgs = fitstools.colAvg(data)
    fig, ax = plt.subplots()
    ax.plot(colAvgs)
    
    threshhold = sum(colAvgs)/len(colAvgs)
    ax.axhline(threshhold, ls='--')
    
    #Identify significant peaks
    peaks = []
    for i in range(len(colAvgs))[1:-1]:
        if colAvgs[i] > threshhold and colAvgs[i-1] < colAvgs[i] and colAvgs[i+1] < colAvgs[i]:
            ax.plot(i, colAvgs[i], 'o', color='blue')
            peaks.append(i)

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

