import numpy as np
import matplotlib.pyplot as plt
plt.ion()
from fitstools import mask_fits, row_avg, manage_dtype
from spectra import spectrum, interp_add
from astropy.io import fits

#04/20/17 12:50 | Need to work out inheritence of curve in the spectrum class.

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
