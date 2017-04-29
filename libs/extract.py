import numpy as np
import matplotlib.pyplot as plt
plt.ion()
from fitstools import mask_fits, row_avg, manage_dtype, common_header, pad_array, display
from spectra import spectrum, interp_add
from astropy.io import fits

class fibers:
    def __init__(self):
        self.spectra = {}
    def add_spectrum(self, fiber_num, spec):
        self.spectra[fiber_num] = spec
    def get_spectra(self):
        fiber_nums = sorted(self.spectra.keys())
        return [self.spectra[fiber_num] for fiber_num in fiber_nums]
    def __getitem__(self, i):
        return self.get_spectra()[i]
    def __setitem__(self, i, new_spec):
        fiber_nums = sorted(self.spectra.keys())
        fiber_num = fiber_nums[i]
        self.spectra[fiber_num] = new_spec
        a = np.array(new_spec.get_flux())
    def save(self, savepath):
        #Generate header
        spec_headers = filter(None, [spec.get_header() for spec in self.get_spectra()])
        header = common_header(spec_headers)
        if header == None:
            header = fits.PrimaryHDU(np.array([])).header
        for i, fnum in enumerate(sorted(self.spectra.keys())):
            try:
                apinfo = self.spectra[fnum].get_header()['SLFIB'+str(fnum)]
            except TypeError:
                apinfo = ''
            header['APINFO'+str(int(i+1))] = apinfo

        length = max([len(spec.get_flux()) for spec in self.get_spectra()])
        flux_dat = np.asarray([pad_array(spec.get_flux(), np.nan, length) for spec in self.get_spectra()], dtype='float64')
        flux = fits.PrimaryHDU(flux_dat, header)
        flux.header['EXTNAME'] = 'FLUX'

        wavelength_dat = np.asarray([pad_array(spec.get_wavelength(), np.nan, length) for spec in self.get_spectra()])
        wavelength = fits.ImageHDU(wavelength_dat, header)
        wavelength.header['EXTNAME'] = 'WAVELENGTH'

        flux_err_dat = np.asarray([pad_array(spec.get_flux_err(), np.nan, length) for spec in self.get_spectra()])
        flux_err = fits.ImageHDU(flux_err_dat, header)
        flux_err.header['EXTNAME'] = 'FLUX_ERR'


        f = fits.HDUList([flux, wavelength, flux_err])
        f.writeto(savepath, clobber=True)


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
def optimal_extraction(image, fiber_mask, use_fibers, flat, wvlsol):
    #Get image header and data
    header = image[1]
    image = image[0]

    #Take info from the header and use it to find out dark noise.
    rn = header['RDNOISE']
    gain = header['GAIN']
    exptime = header['EXPTIME']
    dark_noise = fits.open('calib/master_calib/dark_err.fits')[0].data
    dark_noise /= gain
    dark_noise *= exptime

    #Make copy of header without SLFIB keywords.
    h = header.copy()
    SLFIB_keywords = [k for k in h.keys() if 'SLFIB' in k]
    APINFO = {}
    for kw in SLFIB_keywords:
        try:
            fnum = int(kw.split('SLFIB')[-1])
            if fnum in use_fibers:
                APINFO[fnum] = h[kw]
        except ValueError:
            pass
        del h[kw]

    res = fibers()
    for fnum in use_fibers:
        image_fib = mask_fits(image, fiber_mask, maskval=fnum, reshape=True)
        flat_fib = mask_fits(flat, fiber_mask, maskval=fnum, reshape=True)
        wvlsol_fib = mask_fits(wvlsol, fiber_mask, maskval=fnum, reshape=True)
        dn = mask_fits(dark_noise, fiber_mask, maskval=fnum, reshape=True)

        one = np.ones_like(image_fib)
        err = (one*rn**2 + abs(gain*image_fib) + dn**2)**0.5
        weights = 1/err**2

        wfi = weights*flat_fib*image_fib
        wff = weights*flat_fib**2

        #Use the center of the fiber as the wavelength domain.
        center_i = wvlsol_fib.shape[1]//2
        wvlsol_slices = [wvlsol_fib[:,i] for i in range(len(wvlsol_fib[0]))]
        wfi_slices = [wfi[:,i] for i in range(len(wfi[0]))]
        wff_slices = [wff[:,i] for i in range(len(wff[0]))]
        wavelength, sum_wfi = interp_add(*zip(wvlsol_slices, wfi_slices), x_interp_i=center_i)
        wavelength, sum_wff = interp_add(*zip(wvlsol_slices, wff_slices), x_interp_i=center_i)

        flux = sum_wfi/sum_wff
        flux_err = 1/sum_wff**0.5

        fiber_h = h.copy()
        fiber_h['SLFIB'+str(int(fnum))] = APINFO[int(fnum)]
        spec = spectrum(wavelength, flux, flux_err, header=fiber_h)
        res.add_spectrum(fnum, spec)
    res.save('/Users/richardseifert/Desktop/WIYN/Hydra_reductions/testing_new/saved_spectra/saved_fibers.fits')

    return res
