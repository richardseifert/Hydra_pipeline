from fitstools import manage_dtype, mask_fits, assign_header
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import cosmics


def calibrate(image, bias, fiber_mask=None, lacosmic=True):
    image = bias_correct(image, bias, fiber_mask)
    image = dark_correct(image)
    image = mask_badpixels(image)
    if lacosmic:
    	image = remove_cosmics(image)
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
        try:
        	exptime = image[0].header['EXPTIME']
        except KeyError:
            exptime = 0.0
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

@manage_dtype(preserve=True, use_args=[0], with_header=True)
def remove_cosmics(image, gain=None, readnoise=None, sigclip=5.0, sigfrac=0.5, objlim=20.0):
	image, header = image[:2]
	fname_noext = header['FILENAME'][:-5] if 'FILENAME' in header else 'image'
	if gain==None and 'GAIN' in header:
		gain = header['GAIN']
	if readnoise==None and 'RDNOISE' in header:
		readnoise = header['RDNOISE']
	if gain==None:
		raise KeyError('Cannot determine image gain from information given.')
	if readnoise==None:
		raise KeyError('Cannot determine image readnoise from information given.')

	c = cosmics.cosmicsimage(image, gain=gain, readnoise=readnoise, sigclip=sigclip, sigfrac=sigfrac, objlim=objlim, verbose=False)
	c.run(maxiter=5)
	cosmics_mask = c.mask
	#cosmics.tofits('plots/lacosmic/'+fname_noext+'_cmask.fits', np.transpose(cosmics_mask), header)
	#cosmics.tofits('plots/lacosmic/'+fname_noext+'_before.fits', np.transpose(image), header)
	cosmics_masked_image = mask_fits(image, cosmics_mask, maskval=0.0, fillval=np.nan)
	#cosmics.tofits('plots/lacosmic/'+fname_noext+'_after.fits', np.transpose(cosmics_masked_image), header)
	
	if type(cosmics_masked_image) == fits.hdu.hdulist.HDUList:
		cosmics_masked_image = assign_header(cosmics_masked_image, image[0].header)
		cosmics_masked_image[0].header['COMMENT'] = 'Cosmic rays masked.'

	return cosmics_masked_image,header
