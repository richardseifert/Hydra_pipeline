import numpy as np
#import matplotlib.pyplot as plt
#plt.ion()
from astropy.io import fits
from astropy.wcs import WCS
from scipy import signal
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline
from PyAstronomy.pyasl import crosscorrRV
from xcorl import xcorl

#interp1d = lambda x,y,interp1d=interp1d,**kwargs: interp1d(x,y,kind="linear",**kwargs)

def convolve_gaussian(x, y, sig=1.0):
    '''
    Convolve function with a gaussian.
    ARGUMENTS:
        x - array-like x values of function.
        y - array-like y values of function.
        sig - standard deviation of the gaussian function to convolve with.
    RETURNS:
        y_conv - array-like y values of the convolved function.
    '''
    x = np.array(x)
    dx = np.mean(x[1:]-x[:-1])
    gx = np.arange(-5*sig, 5*sig, dx)
    gaussian = np.exp(-(gx/sig)**2/2) / (sig*(2*np.pi)**0.5)
    n = len(x) / (max(x)-min(x))
    return signal.convolve(y, gaussian, mode="same") / n

# def flatten_spec(wav, flux, twav, tflux, ax=None):
#     #Normalize spectrum
#     flux /= np.nanmedian(flux)
#     if type(ax)!=type(None):
#         ax.plot(wav, flux, color="green", label="Skyflat spectrum")

#     #Use template spectrum to select continuum points.
#     tflux_interp = interp1d(twav, tflux)(wav)
#     mask = tflux_interp >= 0.998
#     if type(ax)!=type(None):
#         ax.plot(wav[mask], flux[mask], 'o', color="purple", label="Continuum points")

#     #Fit smooth curve to continuum points
#     cont = UnivariateSpline(wav[mask], flux[mask], s=2.0e-4*len(wav[mask]))(wav)
#     if type(ax)!=type(None):
#         ax.plot(wav, cont, color="orange", label="Continuum fit")

#     return flux/cont

def flatten_spec(wav, flux, twav, tflux, plot=False):
    wmax = np.nanmax(wav)
    wmin = np.nanmin(wav)
    width = (wmax-wmin) / 5.
    cwav = []
    cflux = []
    quad = lambda x, a, b, c: a*x**2+b*x+c
    p = 97
    npts = 500
    minwid = 10.
    for w in np.linspace(wmin-0.45*width, wmax+0.45*width, npts):
        mask = (wav > w-0.5*width) & (wav < w+0.5*width)
        if np.all(~mask):
            break
        wcut = wav[mask]
        fcut = flux[mask]
        #ignore point ~5% below continuum.
        mask = fcut >= 0.95*np.nanmedian(fcut)
        fcut = fcut[mask]
        wcut = wcut[mask]
        #Flatten chunk with a 2nd degree polynomial.
        fcut_flat = fcut/quad(wcut, *curve_fit(quad, wcut, fcut)[0])
        fp = np.percentile(fcut_flat, p)
        sort = np.argsort(abs(fcut_flat-fp))
        wadd = wcut[sort][0]
        fadd = fcut[sort][0]
        if not wadd in cwav:
            cwav.append(wadd)
            cflux.append(fadd)
    cwav = np.array(cwav)
    sort = np.argsort(cwav)
    cwav = cwav[sort]
    cflux = np.array(cflux)[sort]

    #cont = UnivariateSpline(cwav, cflux, s=1e7*len(cwav))(wav)
    cont = interp1d(cwav, cflux, kind="linear", fill_value="extrapolate")(wav)
    # if plot:
    #     fig, ax = plt.subplots()
    #     ax.plot(wav, flux, color="blue", label="Spectrum")
    #     ax.plot(wav, cont, color="orange", label="Continuum fit")
    #     ax.plot(wav, flux/cont * np.nanmean(flux), color="black")
    #     ax.scatter(cwav, cflux, color="orange", label="Continuum points")
    return flux/cont

def doppler_shift(w, v):
    return w*(1+v/3.0e5)

def get_rv(w, f, tw, tf, logspace=False, range1=50, mult=1):
    #Remove NaN points from spectrum and template
    no_nan = ~(np.isnan(w) | np.isnan(f))
    w=w[no_nan]
    f=f[no_nan]
    no_nan = ~(np.isnan(tw) | np.isnan(tf))
    tw=tw[no_nan]
    tf=tf[no_nan]

    if logspace:
        # Put both the observed and template spectra on a wavelength domain that is evenly spaced in log space.
        logw = np.logspace(np.log10(min(w)), np.log10(max(w)), len(w))
        f = interp1d(w, f, fill_value="extrapolate")(logw)
        w = logw
        #logtw = np.logspace(np.log10(min(tw)), np.log10(max(tw)), len(tw))
        #tf = interp1d(w, f, fill_value="extrapolate")(logtw)
        #tw = logtw

    shft = xcorl(f, interp1d(tw,tf)(w), range1=range1, mult=mult)[0]
    rv_shft = (np.divide(*interp1d(range(len(w)),w)([len(w)/2+shft,len(w)/2])) - 1.0)*3.0e5 #Interpolate to turn decimal index shift into an rv shift.
    return rv_shft


if __name__ == "__main__":
    #Load test skyflat.
    path = "../calib/old/Dec2016/P1/skyflat_spec.fits"
    f = fits.open(path)
    flux_dat = f[0].data
    wav_dat = f[1].data
    f.close()

    #Load BASS2000 solar spectrum.
    solar_w, solar_f = np.loadtxt("../calib/master_calib/bass2000_6000_7000.txt", unpack=True, dtype=np.float32)
    solar_f /= 10000.0
    solar_f /= 0.5*(np.nanmax(solar_f)+np.nanmedian(solar_f))
    solar_f_smooth = convolve_gaussian(solar_w, solar_f, sig=.16)

    #Get single skyflat spectrum.
    f = open("Dec2016_skyflat_rvs2.txt", 'a')
    f.close()
    diffs={}
    for i in [7, 18, 20, 24, 55, 71, 73, 78]:
        print i
        flux = flux_dat[i]
        wav = wav_dat[i]
        
        #fig = plt.figure()
        #ax1 = fig.add_subplot(211)
        #ax2 = fig.add_subplot(212, sharex=ax1)
        #ax1.axhline(y=1.0)
        flux_flat = flatten_spec(wav, flux, twav=solar_w, tflux=solar_f_smooth, plot=True)
        ax = plt.gca()
        raw_input("WAAAA")
        #ax2.plot(solar_w, solar_f_smooth, color="black", label="BASS2000 solar spectrum")
        #ax2.plot(wav, flux_flat, color="green", label="Flattened skyflat spectrum")
        #ax2.set_xlim(min(wav), max(wav))
        #ax1.legend(loc=0)
        #ax2.legend(loc=4)
        #ax2.set_xlabel("Wavelength ($\AA$)")
        #ax1.set_ylabel("Normalized Flux")
        #ax2.set_ylabel("Normalized Flux")

        #rv = get_rv(wav, flux_flat, solar_w, solar_f_smooth, 0.0, 0.0001)
        #new_wav1 = doppler_shift(wav,-rv)
        #chi1 = np.sum([(o-m)**2 for o,m in zip(flux_flat, interp1d(solar_w, solar_f_smooth)(new_wav1))])

        #shft = xcorl(flux_flat, interp1d(solar_w,solar_f_smooth)(wav), 50, mult=1)[0]
        #rv_shft = (np.divide(*interp1d(range(len(wav)),wav)([len(wav)/2+shft,len(wav)/2])) - 1.0)*3.0e5 #?? How do I convert shft into a wavelength or an RV? 

        #Run in logspace wavelengths
        #np.savetxt("skyflat_spec"+str(i)+".txt", zip(wav, flux, flux_flat))
        rv_shft1 = get_rv(wav, flux_flat, solar_w, solar_f_smooth, logspace=True)
        print "Shift found:", rv_shft1
        #np.savetxt("solar_spec.txt", zip(solar_w, solar_f_smooth))
        new_wav1 = doppler_shift(wav,-rv_shft1)
        shift2=get_rv(new_wav1, flux_flat, solar_w, solar_f_smooth, logspace=True)
        print "Shift found after shifting:", shift2
        chi1 = np.sum([(o-m)**2 for o,m in zip(flux_flat, interp1d(solar_w, solar_f_smooth)(new_wav1))])
        plt.title("Log space wavelengths")

        #Run in linspace wavelengths
        #rv_shft2 = get_rv(wav, flux_flat, solar_w, solar_f_smooth, logspace=False)
        #new_wav2 = doppler_shift(wav,-rv_shft2)
        #chi2 = np.sum([(o-m)**2 for o,m in zip(flux_flat, interp1d(solar_w, solar_f_smooth)(new_wav2))])
        #plt.title("Linear space wavelengths")

        #diffs[i] = abs(rv_shft1 - rv_shft2)

        #print "RV1",rv_shft1,chi1,"RV2",rv_shft2,chi2
        rv_shft=rv_shft1

        crosscorr = lambda v, w=wav, f=flux_flat, tw=solar_w, tf=solar_f_smooth: np.nansum(f*interp1d(doppler_shift(tw,v),tf)(w))
        ax.plot(doppler_shift(wav,-rv_shft1), flux_flat, label="log space RV shift", color="green")
        #ax.plot(doppler_shift(wav,-rv_shft2), flux_flat, label="linear space RV shift", color="red")
        ax.legend()
        #f = open("Dec2016_skyflat_rvs2.txt", 'a')
        #f.write(str(i).ljust(4)+str(rv_shft).ljust(12)+'\n')
        #f.close()
    #print np.array(diffs.keys())[np.argmax(np.array(diffs.values()))]
