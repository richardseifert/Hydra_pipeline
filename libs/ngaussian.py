import numpy as np
from mpfit import mpfit
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def fit_ngaussian(xdata, ydata, n, fast=False):
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
    good = (np.logical_not(np.isnan(ydata))) & (np.logical_not(np.isinf(ydata)))
    xdata = xdata[good]
    ydata = ydata[good]

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

    f = lambda x, sig: make_ngaussian(x, [sig]+p0[1:])
    coeff, err = curve_fit(f, xdata, ydata, p0=[p0[0]])
    sig = coeff[0]
    p0[0] = coeff[0]

    if fast:
        p = p0
    else:
        #sorted_peak_x = np.argsort(peak_x)
        p = [sig]
        i = 0
        #k = 0
        w = 15
        xlow = max(peak_x)
        xhigh = min(peak_x)
        p0_chunk = [sig]
        while i < len(peak_x):
            #print 'peak', i
            amp = p0[(i*2)+1]
            mu = p0[(i*2+1)+1]
            #print 'mu =',mu
            if len(p0_chunk) == 1 or mu >= xlow and mu <= xhigh:
                #print 'ADDING'
                p0_chunk.append(amp)
                p0_chunk.append(mu)
                xlow = min([xlow, mu-w*sig])
                xhigh = max([xhigh, mu+w*sig])
                #print 'xlow =',xlow,'xhigh =',xhigh
                #k += 1
                i += 1
            else:
                #print 'mpfitting '+str(k)+' peaks.'
                in_range = (xdata >= xlow) & (xdata <= xhigh)
                mu_par = {'LIMITED':[1,1],'LIMITS':[xlow,xhigh]}
                amp_par = {'LIMITED':[1,0],'LIMITS':[0.0,0]}
                parinfo = [{}]+[amp_par,mu_par]*((len(p0_chunk)-1)/2)
                #print parinfo
                m = mpfit(ngaussian_funct, p0_chunk, {'xdata':xdata[in_range], 'ydata':ydata[in_range]}, parinfo=parinfo, quiet=0)
                print 'FINAL', m.params
                p.extend(m.params[1:])
                #k = 0
                xlow = max(peak_x)
                xhigh = min(peak_x)
                p0_chunk = [sig]

        in_range = (xdata >= xlow) & (xdata <= xhigh)
        mu_par = {'limited':[1,1],'limits':[xlow,xhigh]}
        amp_par = {'limited':[1,0],'limits':[0.0,0]}
        parinfo = [{'limited':[1,0],'limits':[0.0,0]}]+[amp_par,mu_par]*((len(p0_chunk)-1)/2)
        #print parinfo
        m = mpfit(ngaussian_funct, p0_chunk, {'xdata':xdata[in_range], 'ydata':ydata[in_range]}, parinfo=parinfo, quiet=0)
        print 'FINAL', m.params
        p.extend(m.params[1:])
                
            
            #xlow = min([xlow, mu-sig])
            #xhigh = max([xhigh, mu+sig])
        #m = mpfit(ngaussian_funct, p0, {'xdata':xdata, 'ydata':ydata}, quiet=1)
        #p = m.params

    plot=True
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
    
    #Sort by peak height to select the tallest num_peaks peaks
    sort_i = np.argsort(-peak_yvals)
    peak_xvals = peak_xvals[sort_i][:num_peaks]
    peak_yvals = peak_yvals[sort_i][:num_peaks]
    
    #Sort by peak position
    sort_i = np.argsort(peak_xvals)
    peak_xvals = peak_xvals[sort_i]
    peak_yvals = peak_yvals[sort_i]
    return peak_xvals, peak_yvals

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

    while high-low<4:
        if low > 0:
            low -= 1
        high += 1
    region_x = xlist[low:high+1]
    region_y = ylist[low:high+1]

    #Fit a cubic spline to the peak
    peak = interp1d(region_x, region_y, kind='cubic')
    xfit = np.arange(min(region_x)+prec/2, max(region_x)-prec/2, prec)
    yfit = peak(xfit)

    #Find the peak center from spline fit.
    center_x = xfit[list(yfit).index(max(yfit))]

    return center_x, max(yfit)

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
