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
                apinfo = self.spectra[fnum].get_header()['SLFIB'+str(int(fnum))]
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
    print fiber_num, 'fibernumber!!!'
    fiber = mask_fits(img, fiber_mask, fiber_num, reshape=True)
    #display(img)
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
def optimal_extraction(image, fiber_mask, profile_map, wvlsol=None, use_fibers=None):
    if type(use_fibers)==type(None):
        use_fibers = list({n for row in fiber_mask for n in row if n != 0})
    #Get image header and data
    print image
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
        print fnum

        D = mask_fits(image, fiber_mask, maskval=fnum, reshape=True)

        P = mask_fits(profile_map, fiber_mask, maskval=fnum, reshape=True)

        dn = mask_fits(dark_noise, fiber_mask, maskval=fnum, reshape=True)
        one = np.ones_like(D)
        V = (one*rn**2 + abs(gain*D) + dn**2)
        
        f_numer = P*D/V
        var_numer = P/V
        denom = P**2/V

        if type(wvlsol)!=type(None):
            wvlsol_fib = mask_fits(wvlsol, fiber_mask, maskval=fnum, reshape=True)
            #Use the center of the fiber as the wavelength domain.
            center_i = wvlsol_fib.shape[1]//2
            wvlsol_slices = [wvlsol_fib[:,i] for i in range(len(wvlsol_fib[0]))]
            f_numer_slices = [f_numer[:,i] for i in range(len(f_numer[0]))]
            var_numer_slices = [var_numer[:,i] for i in range(len(var_numer[0]))]
            denom_slices = [denom[:,i] for i in range(len(denom[0]))]
            wavelength, sum_f_numer = interp_add(*zip(wvlsol_slices, f_numer_slices), x_interp_i=center_i)
            wavelength, sum_var_numer = interp_add(*zip(wvlsol_slices, var_numer_slices), x_interp_i=center_i)
            wavelength, sum_denom = interp_add(*zip(wvlsol_slices, denom_slices), x_interp_i=center_i)
        else:
            sum_f_numer = np.sum(f_numer, axis=0)
            sum_var_numer = np.sum(var_numer, axis=0)
            sum_denom = np.sum(denom, axis=0)
            wavelength = list(range(len(sum_f_numer)))

        f = sum_f_numer/sum_denom
        var_f = sum_var_numer/sum_denom

        fiber_h = h.copy()
        fiber_h['SLFIB'+str(int(fnum))] = APINFO[int(fnum)]
        spec = spectrum(wavelength, flux=f, flux_err=var_f**0.5, header=fiber_h)
        res.add_spectrum(fnum, spec)

    return res
