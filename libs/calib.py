from fitstools import manage_dtype, mask_fits, assign_header
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt


def calibrate(image, bias, fiber_mask=None):
    image = bias_correct(image, bias, fiber_mask)
    image = dark_correct(image)
    image = mask_badpixels(image)
    return image

def bias_correct(image, bias, fiber_mask=None):

    @manage_dtype(preserve=True)
    def bc_helper(image, bias, fiber_mask=None):
        #reduced = (image-median(masked image)) - (bias-median(bias))
        if type(fiber_mask) != type(None):
            masked_image = mask_fits(image, fiber_mask, maskval=0, fillval=np.nan)
            image = (image - np.nanmedian(masked_image)) - (bias - np.median(bias))
        else:
            image = image - bias

        vals = image.flatten()
        fig, ax = plt.subplots()
        bins = np.linspace(-2000, 2000, 2000)
        n, bins, patches = plt.hist(vals, bins=bins, facecolor='green', alpha=0.75)

        return image

    bias_subtracted_image = bc_helper(image, bias, fiber_mask)
    if type(bias_subtracted_image) == fits.hdu.hdulist.HDUList:
        bias_subtracted_image = assign_header(bias_subtracted_image, image[0].header)
        bias_subtracted_image[0].header['COMMENT'] = 'Bias corrected.'

    return bias_subtracted_image


def dark_correct(image, exptime=None):
    dark_map = fits.open('calib/master_calib/dark_fit.fits')[0].data
    header = image[0].header
    gain = header['GAIN']
    dark_map /= gain

    if exptime == None and type(image) == fits.hdu.hdulist.HDUList:
        exptime = image[0].header['EXPTIME']
    else:
        raise ValueError('Cannot determine exposure time for dark subtraction.')

    @manage_dtype(preserve=True)
    def dc_helper(image, dark_map, exptime):
        image = image - exptime*dark_map
        return image


    dark_subtracted_image = dc_helper(image, dark_map, exptime)
    if type(dark_subtracted_image) == fits.hdu.hdulist.HDUList:
        dark_subtracted_image = assign_header(dark_subtracted_image, image[0].header)
        dark_subtracted_image[0].header['COMMENT'] = 'Bias corrected.'

    return dark_subtracted_image


def mask_badpixels(image):
    bad_mask = fits.open('calib/master_calib/badmask.fits')

    @manage_dtype(preserve=True)
    def mbp_helper(image, bad_mask):
        image = mask_fits(image, bad_mask, maskval=1.0, fillval=np.nan)
        return image

    bad_masked_image = mbp_helper(image, bad_mask)
    if type(bad_masked_image) == fits.hdu.hdulist.HDUList:
        bad_masked_image = assign_header(bad_masked_image, image[0].header)
        bad_masked_image[0].header['COMMENT'] = 'Bad pixels masked.'

    return bad_masked_image
