import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
plt.ion()

def get_data_and_header(some_fits):
    if type(some_fits) == fits.hdu.hdulist.HDUList:
        dtype = 0
        image_data = some_fits[0].data
        image_header = some_fits[0].header
    elif type(some_fits) == fits.hdu.image.PrimaryHDU:
        dtype = 1
        image_data = some_fits.data
        image_header = some_fits.header
    else:
        dtype = 2
        image_data = some_fits
        image_header = None
        try:
            if (image_data.shape) == 2:
                image_header = fits.PrimaryHDU(image_data).header
        except AttributeError:
            pass
    return [dtype, image_data, image_header]
def manage_dtype(preserve=False, with_header=False):
    dtypes = [lambda data: fits.HDUList(fits.PrimaryHDU(data)),
              lambda data: fits.PrimaryHDU(data),
              lambda data: data]
    def decorator(f):
        def wrapper(*args, **kwargs):
            args = list(args)
            dtype_i = 2
            for i, arg in enumerate(args):
                d, data, header = get_data_and_header(arg)
                if d < dtype_i:
                    dtype_i = d
                if with_header:
                    args[i] = [data, header]
                else:
                    args[i] = data
                i+=1
            
            dtype = dtypes[dtype_i]

            res = f(*args, **kwargs)
            if preserve:
                if type(res) == list:
                    for i, r in enumerate(res):
                        try:
                            if len(r.shape) == 2:
                                res[i] = dtype(res[i])
                        except AttributeError:
                            pass
                else:
                    try:
                        if len(res.shape) == 2:
                            res = dtype(res)
                    except AttributeError:
                        pass
            return res
        return wrapper
    return decorator
    
@manage_dtype()
def display(some_fits, ax=None):
    if ax == None:
        fig, ax = plt.subplots()

    im = ax.imshow(some_fits, cmap='gray', interpolation='nearest', aspect='auto', norm=LogNorm())
    fig.colorbar(im)
    return ax 

@manage_dtype()
def rowAvg(some_fits):
    return [np.mean(some_fits[i,:]) for i in range(len(some_fits))]

@manage_dtype()
def colAvg(some_fits):
        return [np.mean(some_fits[:,i]) for i in range(len(some_fits[0]))]

@manage_dtype(preserve=True)
def slice(some_fits, xlo=None, xhi=None, ylo=None, yhi=None):
    if xlo == None:
        xlo = 0
    if xhi == None:
        xhi = len(some_fits[0])
    if ylo == None:
        ylo = 0
    if yhi == None:
        yhi = len(some_fits)

    data_slice = some_fits[ylo:yhi+1,xlo:xhi+1]
    return data_slice

@manage_dtype(with_header=True)
def get_common_header(*args):
    list_of_fits = args
    new_header = fits.PrimaryHDU(list_of_fits[0][0]).header
    headers = [f[1] for f in list_of_fits if f[1] != None]
    sample_header = headers[0]
    for card in sample_header.keys():
        useCard = True
        for h in headers:
            if (not card in h) or (h[card] != sample_header[card]):
                useCardr = False
        
        if useCard:
            if card != 'COMMENT' and card != 'HISTORY' and card != '':
                new_header[card] = sample_header[card]
    
    return new_header

#Change a fits header in place to the header given.
def assign_header(some_fits, header):
    if type(some_fits) == fits.hdu.hdulist.HDUList:
        some_fits[0].header = header
    elif type(some_fits) == fits.hdu.image.PrimaryHDU:
        some_fits.header = header
    else:
        raise TypeError(str(some_fits) + 'is not an appropriate object for fits header assignment.')
    
def combine(*args, **kwargs):
    comb_fits = combine_helper(*args, **kwargs)
    list_of_fits = args
    assign_header(comb_fits, get_common_header(*list_of_fits))
    return comb_fits

@manage_dtype(preserve=True)
def combine_helper(*args, **kwargs):
    method = 'median'
    if 'method' in kwargs:
        method = kwargs['method']

    list_of_fits = args

    if method == 'median':
        return np.median(list_of_fits, axis=0)
    else:
        raise ValueError('Unknown method' + str(method))

def save_2darr(data, savepath):
    p = fits.PrimaryHDU(data)
    hdulist = fits.HDUList(p)
    hdulist.writeto(savepath, clobber=True)

@manage_dtype()
def mask_fits(some_fits, some_mask, maskval=1):
    if some_fits.shape != some_mask.shape:
        print 'Data and mask must be the same shape.'
        raise ValueError
    return np.where(some_mask == maskval, some_fits, 0)
