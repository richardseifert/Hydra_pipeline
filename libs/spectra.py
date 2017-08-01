import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
plt.ion()
from astropy.io import fits
from fitstools import common_header

def unpack_xy(use_args='all', preserve=False):
    def decorator(f):
        def wrapper(use_args, *args, **kwargs):
            args = list(args)
            if use_args == 'all':
                use_args = [i for i in range(len(args))]
            dtype = None
            for i in use_args:
                #if isinstance(args[i], spectrum):
                #    if d<2:
                #        dtype = lambda w, f, ferr=None, header=None: spectrum(w, f, ferr, header)  #Tried adding spectrum object to
                #        d = 2
                if isinstance(args[i], curve):
                    dtype = lambda x, y, yerr=None: curve(x, y, yerr)
                else:
                    try:
                        iter(args[i])
                    except TypeError:
                        continue
                    if len(args[i]) == 3:
                        x, y, yerr = args[i]
                        args[i] = curve(x, y, yerr)
                    elif len(args[i]) == 2:
                        x, y = args[i]
                        args[i] = curve(x, y)
                    else:
                        continue
            res = f(*args, **kwargs)
            if preserve and dtype != None:
                res = dtype(*res)
            return res
        return lambda *args, **kwargs: wrapper(use_args, *args, **kwargs)
    return decorator

class curve:
    def __init__(self, x, y, yerr=None):
        sort_i = np.argsort(x) #Sort data by x.
        self.x = np.asarray(x)[sort_i]
        self.y = np.asarray(y)[sort_i]
        if type(yerr) == type(None):
            self.yerr = np.zeros_like(y)[sort_i]
        else:
            self.yerr = np.asarray(yerr)[sort_i]

    def get_x(self):
        return self.x
    def get_y(self):
        return self.y
    def get_yerr(self):
        return self.yerr
    def get_data(self):
        return self.x, self.y, self.yerr

    @unpack_xy()
    def math_helper(c1, c2, **kwargs):
        if isinstance(c2, curve):
            x_interp = get_x_interp([c1.x, c2.x], **kwargs)
            c1_y_interp = interp1d(c1.x, c1.y)(x_interp)
            c1_yerr_interp = interp1d(c1.x, c1.yerr)(x_interp)
            c1_interp = curve(x_interp, c1_y_interp, c1_yerr_interp)

            c2_y_interp = interp1d(c2.x, c2.y)(x_interp)
            c2_yerr_interp = interp1d(c2.x, c2.yerr)(x_interp)
            c2_interp = curve(x_interp, c2_y_interp, c2_yerr_interp)

            return c1_interp, c2_interp
        else:
            return c1, curve(c1.x, c2*np.ones_like(c1.y))

    def __add__(self, other, **kwargs):
        self_interp, other_interp = curve.math_helper(self, other, **kwargs)
        x_interp = self_interp.x
        y_interp = self_interp.y+other_interp.y
        yerr_interp = (self_interp.yerr**2+other_interp.yerr**2)**0.5
        return curve(x_interp, y_interp, yerr_interp)
    def __sub__(self, other, **kwargs):
        self_interp, other_interp = curve.math_helper(self, other, **kwargs)
        x_interp = self_interp.x
        y_interp = self_interp.y-other_interp.y
        yerr_interp = (self_interp.yerr**2+other_interp.yerr**2)**0.5
        return curve(x_interp, y_interp, yerr_interp)
    def __mul__(self, other, **kwargs):
        self_interp, other_interp = curve.math_helper(self, other, **kwargs)
        x_interp = self_interp.x
        y_interp = self_interp.y*other_interp.y
        yerr_interp = ((self_interp.yerr*other_interp.y)**2 + (other_interp.yerr*other_interp.y)**2)**0.5
        return curve(x_interp, y_interp, yerr_interp)
    def __div__(self, other, **kwargs):
        self_interp, other_interp = curve.math_helper(self, other, **kwargs)
        x_interp = self_interp.x
        y_interp = self_interp.y/other_interp.y
        yerr_interp = ((self_interp.yerr*other_interp.y)**2 + (other_interp.yerr*other_interp.y)**2)**0.5
        return curve(x_interp, y_interp, yerr_interp)

def get_x_interp(x_arrs, x_interp=None, x_interp_i=None, dx=None, **kwargs):
    if x_interp == None:
        try:
            x_interp = x_arrs[x_interp_i]
        except TypeError, IndexError:
            low = max([min(x_arr) for x_arr in x_arrs]) #Find the lowest x value
            high = min([max(x_arr) for x_arr in x_arrs]) #Find the highest x value
            if dx != None:
                x_interp = np.arange(low, high, dx)
            else:
                x_interp = []
                num_x = len(x_arrs)
                x_i_list = [0]*num_x
                current_x = low
                while current_x < high:
                    x_interp.append(current_x)
                    avg_dx = 0
                    n = 0
                    for i,x in enumerate(x_arrs):
                        indx = x_i_list[i]
                        while indx < len(x) and x[indx] < current_x:
                            indx += 1
                        x_i_list[i] = int(indx)
                        try:
                            avg_dx += abs(x[indx+1] - x[indx])
                            n+=1
                        except:
                            pass
                    avg_dx = avg_dx/n if n>0 else last_dx
                    current_x += avg_dx
                    last_dx = avg_dx

    return x_interp

@unpack_xy()
def interp_helper(*xy_curves, **kwargs):
    x_arrs = [c.get_x() for c in xy_curves]
    y_arrs = [c.get_y() for c in xy_curves]
    yerr_arrs = [c.get_yerr() for c in xy_curves]

    x_interp = get_x_interp(x_arrs=x_arrs, **kwargs)

    y_interp_arrs = np.zeros((len(y_arrs), len(x_interp)))
    for i in range(len(x_arrs)):
        y_interp_arrs[i,:] = interp1d(x_arrs[i], y_arrs[i], fill_value=(np.nan, np.nan))(x_interp)
    yerr_interp_arrs = np.zeros((len(yerr_arrs), len(x_interp)))
    for i in range(len(x_arrs)):
        yerr_interp_arrs[i,:] = interp1d(x_arrs[i], yerr_arrs[i], fill_value=(np.nan, np.nan))(x_interp)
    return x_interp, y_interp_arrs, yerr_interp_arrs

@unpack_xy(preserve=True)
def interp_add(*spectra, **kwargs):
    x_interp, y_interp_arrs, yerr_interp_arrs = interp_helper(*spectra, **kwargs)
    y_interp = np.nansum(y_interp_arrs, axis=0)
    yerr_interp = np.nansum([yerr**2 for yerr in yerr_interp_arrs], axis=0)**0.5
    return x_interp, y_interp

@unpack_xy(preserve=True)
def interp_mean(*spectra, **kwargs):
    x_interp, y_interp_arrs, yerr_interp_arrs = interp_helper(*spectra, **kwargs)
    y_interp = np.nanmean(y_interp_arrs, axis=0)
    yerr_interp = np.nansum([yerr**2 for yerr in yerr_interp_arrs], axis=0)**0.5/N
    return x_interp, y_interp

@unpack_xy(preserve=True)
def interp_median(*spectra, **kwargs):
    x_interp, y_interp_arrs, yerr_interp_arrs = interp_helper(*spectra, **kwargs)
    y_interp = np.nanmedian(y_interp_arrs, axis=0)
    N = len(y_interp_arrs)
    yerr_interp = 1.253*np.nansum([yerr**2 for yerr in yerr_interp_arrs], axis=0)**0.5/N
    return x_interp, y_interp, yerr_interp

class spectrum(curve):
    def __init__(self, wavelength, flux=None, flux_err=None, header=None):
        if type(flux) == type(None) and isinstance(wavelength, curve):
            input_curve = wavelength
            curve.__init__(self, *input_curve.get_data())
        else:
            curve.__init__(self, wavelength, flux, flux_err)
        self.header = header
    def set_header(self, new_header):
        self.header = new_header
    def get_wavelength(self):
        return self.x
    def get_flux(self):
        return self.y
    def get_flux_err(self):
        return self.yerr
    def get_data(self):
        return [self.x, self.y, self.yerr]
    def get_header(self):
        return self.header
    def __add__(self, other, header_i=None):
        if header_i == None:
            try:
                headers = [self.header, other.header]
                header = common_header(headers)
            except AttributeError:
                header = self.header
        return spectrum(curve.__add__(self, other), header=header)
    def __sub__(self, other, header_i=None):
        if header_i == None:
            try:
                headers = [self.header, other.header]
                header = common_header(headers)
            except AttributeError:
                header = self.header
        return spectrum(curve.__sub__(self, other), header=header) #None is temp, REMOVE SOON
    def __mul__(self, other, header_i=None):
        if header_i == None:
            try:
                headers = [self.header, other.header]
                header = common_header(headers)
            except AttributeError:
                header = self.header
        return spectrum(curve.__mul__(self, other), header=header)
    def __div__(self, other, header_i=None):
        if header_i == None:
            try:
                headers = [self.header, other.header]
                header = common_header(headers)
            except AttributeError:
                header = self.header
        return spectrum(curve.__div__(self, other), header=header)
    def save(self, savepath):
        flux = fits.PrimaryHDU(self.get_flux(), self.get_header())
        flux.header['EXTNAME'] = 'FLUX'
        wavelength = fits.ImageHDU(self.get_wavelength())
        wavelength.header['EXTNAME'] = 'WAVELENGTH'
        flux_err = fits.ImageHDU(self.get_flux_err())
        flux_err.header['EXTNAME'] = 'FLUX_ERR'
        f = fits.HDUList([flux, wavelength, flux_err])
        f.writeto(savepath, clobber=True)
    def plot(self, ax=None, **kwargs):
        if ax == None:
            fig, ax = plt.subplots()
        ax.set_xlabel('Wavelength ($\AA$)')
        ax.set_ylabel('Flux')
        ax.plot(self.x, self.y, **kwargs)
        if type(self.yerr) != type(None):
            ax.fill_between(self.x, self.y-self.yerr, self.y+self.yerr, facecolor='cornflowerblue', linewidth=0.0)
        return ax
def sum_spectra(spectra, header=None, **kwargs):
    if header==None:
        #Combine headers somehow
        pass
    sum_curve = interp_add(*spectra, **kwargs)
    sum_spectrum = spectrum(sum_curve, header=header)
    return sum_spectrum
def median_spectra(spectra, header=None, **kwargs):
    if header==None:
        #Combine headers somehow
        pass
    median_curve = interp_median(*spectra, **kwargs)
    median_spectrum = spectrum(median_curve, header=header)
    return median_spectrum
def mean_spectra(spectra, header=None, **kwargs):
    if header==None:
        #Combine headers somehow
        pass
    mean_curve = interp_mean(*spectra, **kwargs)
    mean_spectrum = spectrum(mean_curve, header=header)
    return mean_spectrum
