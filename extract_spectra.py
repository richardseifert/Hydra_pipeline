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
    #flux = np.zeros_like(wavelength)
    #for i in range(len(fiber_wvlsol[0])):
    #    wvlsol_slice = fiber_wvlsol[:,i]
    #    counts_slice = fiber_counts[:,i]
    #    interp_flux = interp1d(wvlsol_slice, counts_slice)(wavelength)
    #    flux += interp_flux
    wvlsol_slices = [fiber_wvlsol[:,i] for i in range(len(fiber_wvlsol[0]))]
    counts_slices = [fiber_counts[:,i] for i in range(len(fiber_counts[0]))]
    wavelength, flux = interp_add(wvlsol_slices, counts_slices, x_interp_i=center_i)
    
    return wavelength, flux


def interp_add(x_arrs, y_arrs, x_interp=None, x_interp_i=None, dx=None):
    if x_interp == None:
        try:
            x_interp = x_arrs[x_interp_i]
        except TypeError, IndexError:
            low = max([min(x_arr) for x_arr in x_arrs]) #Find the lowest x value
            high = min([max(x_arr) for x_arr in x_arrs]) #Find the highest x value
            if dx != None:
                x_interp = np.arange(low, high, dx)
            else:
                #Sort the inputs
                sort_i_list = [np.argsort(x_arr) for x_arr in x_arrs]
                y_arrs = [np.asarray(y_arr)[sort_i] for y_arr, sort_i in zip(y_arrs, sort_i_list)]
                x_arrs = [np.asarray(x_arr)[sort_i] for x_arr, sort_i in zip(x_arrs, sort_i_list)]
                
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
                    avg_dx /= n
                    current_x += avg_dx
    
    # x_interp should now be defined.
    
    x_interp = np.asarray(x_interp)
    y_sum = np.zeros_like(x_interp)
    for x_arr, y_arr in zip(x_arrs, y_arrs):
        y_interp = interp1d(x_arr, y_arr)(x_interp)
        y_sum += y_interp

    return x_interp, y_sum



