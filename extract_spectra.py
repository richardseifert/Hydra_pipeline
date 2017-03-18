import numpy as np
import matplotlib.pyplot as plt
plt.ion()
from fitstools import mask_fits
from scipy.interpolate import interp1d

def extract(fiber_mask, fiber_num, img, wvlsol):
    '''
    Function that extracts a 1D spectrum for a specified fiber.

    ARGUMENTS:
    ----------------------------------------------------------------------------
    fiber_mask: A 2D array that specifies the locations of fibers on the image,
                img. 

    fiber_num: The fiber to be extracted. This should be an existing fiber in
               fiber_mask.

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
    flux = np.zeros_like(wavelength)
    for i in range(len(fiber_wvlsol[0])):
        wvlsol_slice = fiber_wvlsol[:,i]
        counts_slice = fiber_counts[:,i]
        interp_flux = interp1d(wvlsol_slice, counts_slice)(wavelength)
        flux += interp_flux
    
    if True:
        fig, ax = plt.subplots()
        ax.set_title(str(fiber_num))
        ax.plot(wavelength, flux)

    return wavelength, flux


    
