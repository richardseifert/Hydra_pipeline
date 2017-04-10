import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
plt.ion()

def get_data_and_header(some_fits):
    if type(some_fits) == str:
        dtype = 0
        image = fits.open(some_fits)
        image_data = image[0].data
        image_header = image[0].header
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

    image_wcs = None
    try:
        image_wcs = wcs.WCS(image_header)
    except:
        pass

    return [dtype, image_data, image_header, image_wcs]

def get_dtype(*list_of_fits, **kwargs):
    lowest=True
    if 'lowest' in kwargs:
        lowest=kwargs['lowest']
    dtypes = [lambda data: fits.HDUList(fits.PrimaryHDU(data)),
              lambda data: fits.PrimaryHDU(data),
              lambda data: data]
    dtype_i = 0 if lowest else 2
    for some_fits in list_of_fits:
        if type(some_fits) == str:
            dtype = 0
        if type(some_fits) == fits.hdu.hdulist.HDUList:
            dtype = 0
        elif type(some_fits) == fits.hdu.image.PrimaryHDU:
            dtype = 1
        else:
            dtype = 2
            image_data = some_fits
            try:
                assert len(np.array(image_data).shape) == 2 \
                       or len([dim for dim in np.array(image_data).shape if dim != 1]) == 2
            except:
                print np.array(image_data).shape, len(np.array(image_data).shape) == 2
                raise ValueError('Got an object that is not fits-like')
        if lowest:
            dtype_i = min([dtype_i, dtype])
        else:
            dtype_i = max([dtype_i, dtype])
    return dtypes[dtype_i]

def manage_dtype(use_args='all', preserve=False, with_header=False, with_wcs=False):
    dtypes = [lambda data: fits.HDUList(fits.PrimaryHDU(data)),
              lambda data: fits.PrimaryHDU(data),
              lambda data: data]
    def decorator(f):
        def wrapper(use_args, *args, **kwargs):
            args = list(args)
            if use_args == 'all':
                use_args = [i for i in range(len(args))]
            dtype_i = 2
            for i in use_args:
                d, data, header, wcs = get_data_and_header(args[i])
                if d < dtype_i:
                    dtype_i = d
                args[i] = [data]
                if with_header:
                    args[i].append(header)
                if with_wcs:
                    args[i].append(wcs)
                if not (with_header or with_wcs):
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
        return lambda *args, **kwargs: wrapper(use_args, *args, **kwargs)
    return decorator

@manage_dtype()
def display(some_fits, ax=None):
    if ax == None:
        fig, ax = plt.subplots()

    im = ax.imshow(some_fits, cmap='gray', interpolation='nearest', aspect='auto', norm=LogNorm())
    fig.colorbar(im)
    return ax

@manage_dtype()
def row_avg(some_fits):
    return [np.sum(some_fits[i,:]) for i in range(len(some_fits))]

@manage_dtype()
def col_avg(some_fits):
        return [np.sum(some_fits[:,i]) for i in range(len(some_fits[0]))]

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
    if len(headers) == 0:
        return new_header
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

#Change a fits header NOT in place to the header given.
def assign_header(some_fits, header):
    if type(some_fits) == fits.hdu.hdulist.HDUList:
        data = some_fits[0].data
        return fits.HDUList(fits.PrimaryHDU(data, header))
    elif type(some_fits) == fits.hdu.image.PrimaryHDU:
        data = some_fits.data
        return fits.PrimaryHDU(data, header)
    else:
        raise TypeError(str(some_fits) + 'is not an appropriate object for fits header assignment.')

def combine(*args, **kwargs):
    comb_fits = combine_helper(*args, **kwargs)
    list_of_fits = args
    dtype = get_dtype(list_of_fits)
    comb_fits = dtype(comb_fits)
    comb_fits = assign_header(comb_fits, get_common_header(*list_of_fits))
    return comb_fits

@manage_dtype()
def combine_helper(*args, **kwargs):
    method = 'median'
    if 'method' in kwargs:
        method = kwargs['method']

    list_of_fits = args

    if method == 'median':
        comb_img = np.median(list_of_fits, axis=0)
    elif method == 'mean':
        comb_img = np.mean(list_of_fits, axis=0)
    else:
        raise ValueError('Unknown method ' + str(method))
    return comb_img

def save_2darr(data, savepath):
    p = fits.PrimaryHDU(data)
    hdulist = fits.HDUList(p)
    hdulist.writeto(savepath, clobber=True)

@manage_dtype()
def mask_fits(some_fits, some_mask, maskval=1, fillval=0, reshape=False):
    if reshape:
        return mask_fits_reshape(some_fits, some_mask, maskval)

    if some_fits.shape != some_mask.shape:
        raise ValueError('Data and mask must be the same dimensions')
    return np.where(some_mask == maskval, some_fits, fillval)

def mask_fits_reshape(some_fits, some_mask, maskval=1):
    result = []
    for row in range(len(some_mask)):
        result.append([])
        for column in range(len(some_mask[row])):
            if some_mask[row][column] == maskval:
                result[row].append(some_fits[row][column])
    result = np.asarray(result)
    return result
