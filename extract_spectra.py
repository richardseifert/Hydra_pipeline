import numpy as np
import matplotlib.pyplot as plt
plt.ion()
from fitstools import mask_fits, row_avg
from scipy.interpolate import interp1d

def unpack_xy(use_args='all', preserve=False):
    def decorator(f):
        def wrapper(use_args, *args, **kwargs):
            args = list(args)
            if use_args == 'all':
                use_args = [i for i in range(len(args))]
            dtype = None
            for i in use_args:
                if isinstance(args[i], spectrum):
                    x_arr, y_arr = args[i].get_data()
                    dtype = lambda w, f: spectrum(w, f)
                    d = 1
                elif len(args[i]) == 2:
                    x_arr, y_arr = args[i]
                #Sort x and y
                sort_i = np.argsort(x_arr)
                y_arr = np.asarray(y_arr)[sort_i]
                x_arr = np.asarray(x_arr)[sort_i]
                args[i] = [x_arr, y_arr]
            res = f(*args, **kwargs)
            if preserve and dtype != None:
                res = dtype(*res)
            return res
        return lambda *args, **kwargs: wrapper(use_args, *args, **kwargs)
    return decorator

class spectrum:
    def __init__(self, wavelength, flux, header=None):
        self.wav = wavelength
        self.flux = flux
        self.header = header
    def set_header(self, new_header):
        self.header = new_header
    def get_wavelength(self):
        return self.wav
    def get_flux(self):
        return self.flux
    def get_data(self):
        return [self.wav, self.flux]
    @unpack_xy()
    def math_helper(spec1, spec2, **kwargs):
        s1_x, s1_y = spec1
        s2_x, s2_y = spec2
        x_interp = get_x_interp([s1_x, s2_x], **kwargs)
        s1_y_interp = interp1d(s1_x, s1_y)(x_interp)
        s2_y_interp = interp1d(s2_x, s2_y)(x_interp)
        return x_interp, s1_y_interp, s2_y_interp

    def __add__(self, other, **kwargs):
        try:
            x_interp, self_y_interp, other_y_interp = spectrum.math_helper(self, other, **kwargs)
            y_interp = self_y_interp+other_y_interp
            return spectrum(x_interp, y_interp)
        except:
            print 'BROKE'
    def __sub__(self, other, **kwargs):
        try:
            x_interp, self_y_interp, other_y_interp = spectrum.math_helper(self, other, **kwargs)
            y_interp = self_y_interp-other_y_interp
            return spectrum(x_interp, y_interp)
        except:
            print 'BROKE'

    def save(self, savepath):
        np.savetxt(savepath, zip(*self.get_data()))
    def plot(self, ax=None, **kwargs):
        if ax == None:
            fig, ax = plt.subplots()
        ax.set_xlabel('Wavelength ($\AA$)')
        ax.set_ylabel('Flux')
        ax.plot(self.wav, self.flux, **kwargs)
        return ax

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

def extract(fiber_mask, fiber_num, img, wvlsol):
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
    x_arrs = [curve[0] for curve in xy_curves]
    y_arrs = [curve[1] for curve in xy_curves]

    x_interp = get_x_interp(x_arrs=x_arrs, **kwargs)

    y_interp_arrs = np.zeros((len(y_arrs), len(x_interp)))
    for i in range(len(x_arrs)):
        y_interp_arrs[i,:] = interp1d(x_arrs[i], y_arrs[i], fill_value=(np.nan, np.nan))(x_interp)
    return x_interp, y_interp_arrs

@unpack_xy(preserve=True)
def interp_add(*spectra, **kwargs):
    x_interp, y_interp_arrs = interp_helper(*spectra, **kwargs)
    y_interp = np.nansum(y_interp_arrs, axis=0)
    return x_interp, y_interp

@unpack_xy(preserve=True)
def interp_mean(*spectra, **kwargs):
    x_interp, y_interp_arrs = interp_helper(*spectra, **kwargs)
    y_interp = np.nanmean(y_interp_arrs, axis=0)
    return x_interp, y_interp

@unpack_xy(preserve=True)
def interp_median(*spectra, **kwargs):
    x_interp, y_interp_arrs = interp_helper(*spectra, **kwargs)
    y_interp = np.nanmedian(y_interp_arrs, axis=0)
    return x_interp, y_interp
