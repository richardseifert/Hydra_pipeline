import numpy as np
import matplotlib.pyplot as plt
plt.ion()
from fitstools import mask_fits, row_avg, manage_dtype
from scipy.interpolate import interp1d
from astropy.io import fits

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
            c1_yerr_interp = interp1d(spec1.x, spec1.yerr)(x_interp)
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

#04/20/17 12:50 | Need to work out inheritence of curve in the spectrum class.
class spectrum(curve):
    def __init__(self, wavelength, flux=None, flux_err=None, header=None):
        if flux == None and isinstance(wavelength, curve):
            input_curve = wavelength
            curve.__init__(self, *input_curve.get_data())
            self.header = header
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
    def save(self, savepath):
        np.savetxt(savepath, zip(*self.get_data()))
    def plot(self, ax=None, **kwargs):
        if ax == None:
            fig, ax = plt.subplots()
        ax.set_xlabel('Wavelength ($\AA$)')
        ax.set_ylabel('Flux')
        ax.plot(self.x, self.y, **kwargs)
        if self.yerr!= None:
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

def extract_counts(img, fiber_mask, fiber_num):
    '''
    Function that extracts a 1D list of counts from a fiber.

    ARGUMENTS
    ----------------------------------------------------------------------------
    img: A 2D array containing count information for each fiber.

    fiber_mask: A 2D array that specifies the locations of fibers on the image,
                img.

    fiber_num: The integer ID of the fiber to be extracted. This should be an
               existing fiber in the fiber_mask.
    '''
    fiber = mask_fits(img, fiber_mask, fiber_num)
    counts = row_avg(fiber)

    return counts

def simple_extraction(fiber_mask, fiber_num, img, wvlsol):
    '''
    Function that extracts a 1D spectrum for a specified fiber.

    ARGUMENTS:
    ----------------------------------------------------------------------------
    fiber_mask: A 2D array that specifies the locations of fibers on the image,
                img.

    fiber_num: The integer ID of the fiber to be extracted. This should be an
               existing fiber in the fiber_mask.

    img: A 2D array containing count information for each fiber.

    wvlsol: A 2D array containing wavelength information for each fiber.
    '''
    #Extract the fiber from both the wavelength solution and the image.
    fiber_counts = mask_fits(img, fiber_mask, fiber_num, reshape=True)
    fiber_wvlsol = mask_fits(wvlsol, fiber_mask, fiber_num, reshape=True)


    #Use the center of the fiber as the wavelength domain.
    center_i = fiber_wvlsol.shape[1]//2
    wavelength = fiber_wvlsol[:,center_i]
    if wavelength[0] > wavelength[-1]:
        wavelength = wavelength[::-1]

    #After interpolating to the central wavelength domain, add up counts
    # from each fiber slice.
    wvlsol_slices = [fiber_wvlsol[:,i] for i in range(len(fiber_wvlsol[0]))]
    counts_slices = [fiber_counts[:,i] for i in range(len(fiber_counts[0]))]
    wavelength, flux = interp_add(*zip(wvlsol_slices, counts_slices), x_interp_i=center_i)

    return spectrum(wavelength, flux)

@manage_dtype(with_header=[0])
def optimal_extraction(image, fiber_mask, fnum, flat, wvlsol):
    header = image[1]
    image = image[0]

    image = mask_fits(image, fiber_mask, maskval=fnum, reshape=True)
    flat = mask_fits(flat, fiber_mask, maskval=fnum, reshape=True)
    wvlsol = mask_fits(wvlsol, fiber_mask, maskval=fnum, reshape=True)

    rn = header['RDNOISE']
    g = header['GAIN']
    dark_noise = fits.open('calib/master_calib/dark_err.fits')[0].data
    gain = header['GAIN']
    dark_noise /= gain
    exptime = header['EXPTIME']
    dark_noise *= exptime
    dn = mask_fits(dark_noise, fiber_mask, maskval=fnum, reshape=True)

    one = np.ones_like(image)
    err = (one*rn**2 + abs(g*image) + dn**2)**0.5

    weights = 1/err**2

    wfi = weights*flat*image
    wff = weights*flat**2

    #Use the center of the fiber as the wavelength domain.
    center_i = wvlsol.shape[1]//2
    wvlsol_slices = [wvlsol[:,i] for i in range(len(wvlsol[0]))]
    wfi_slices = [wfi[:,i] for i in range(len(wfi[0]))]
    wff_slices = [wff[:,i] for i in range(len(wff[0]))]
    wavelength, sum_wfi = interp_add(*zip(wvlsol_slices, wfi_slices), x_interp_i=center_i)
    wavelength, sum_wff = interp_add(*zip(wvlsol_slices, wff_slices), x_interp_i=center_i)

    flux = sum_wfi/sum_wff
    flux_err = 1/sum_wff**0.5

    return spectrum(wavelength, flux, flux_err)

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
    fig, ax = plt.subplots()
    ax.plot(x_interp, yerr_interp)
    return x_interp, y_interp, yerr_interp

#04-21-17 | Make a spectrum.median, spectrum.mean, spectrum.sum, etc.
