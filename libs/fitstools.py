'''
This file contains a bunch of useful functions for handling fits files in python.
'''
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
plt.ion()


class manage_dtype:
    '''
    Decorator to handle all the various different ways of representing a fits file.
    This decorator examines arguments of a function and converts any fits-like objects
    into a common datatype.

    ARGUMENTS:
        use_args - list of indices of the arguments to be examined and converted. By default,
                   all arguments are used.
                   ex.) @manage_dtype(use_args=[0,2])
                        def f(arg1, arg2, arg3):
                            #Function body
                        #arg1 and arg3 will be examined and if they are found to be fits-like,
                        # they will be converted to 2D arrays. But arg2, is not examined or
                        # converted.
    '''
    dtypes = [lambda data,header=None: fits.HDUList(fits.PrimaryHDU(data,header)),
              lambda data,header=None: fits.PrimaryHDU(data,header),
              lambda data,header=None: data]
    def __init__(self, use_args='all', preserve=False, with_header=False, with_wcs=False):
        self.use_args = use_args

        self.preserve=preserve
        self.with_header = with_header
        self.with_wcs = with_wcs
    
    def get_data_and_header(self, some_fits):
        '''
        Take a variety of possible FITS-related inputs and return a 2D data array,
         a FITS header, a WCS object, and an indicator of the input data type.
        ARGUMENTS:
            some_fits - Some object that is represents a FITS image. It can be a
                        string path to a .fits file, a preloaded astropy.io.fits
                        HDUList or PrimaryHDU, or a numpy 2D array.
        RETURNS:
            tuple containing:
                dtype - an int designating what type of data some_fits was.
                image_data - numpy 2D array of the FITS data.
                image_header - astropy.io.fits header object, if available. Otherwise None.
                image_wcs - astropy.wcs WCS object, if available. Otherwise None.

        This is really a helper function for the manage_dtype decorator. So it probably
        shouldn't be used directly by anything else.
        '''
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

    def expand_args(self, args):
        #Determine arguments to expand.
        if self.use_args == 'all':
            use_args = [i for i in range(len(args))]
        else:
            use_args = self.use_args

        #For each fits-like argument, expand it into [data, header, wcs].
        new_args = list(args)
        dtype_i = 2 #Move
        for i in use_args:
            new_args[i] = []
            d, data, header, wcs = self.get_data_and_header(args[i])
            if d < dtype_i:  #Move
                dtype_i = d  #Move
            new_args[i].append(data)
            get_header = self.with_header==True or (isinstance(self.with_header,(list,tuple)) and i in self.with_header)
            get_wcs = self.with_wcs==True or (isinstance(self.with_wcs,(list,tuple)) and i in self.with_wcs)
            if get_header:
                new_args[i].append(header)
            if get_wcs:
                new_args[i].append(wcs)
            if len(new_args[i])==1:
                new_args[i] = data

        self.dtype = self.dtypes[dtype_i]  #Move
        return new_args
    def convert_one(self, r):
        #Check if output takes the form, data or [data, header].
        try:
            if len(r) == 2:
                r,header = r
            else:
                header=None
            if len(r.shape) == 2:
                r = self.dtype(r,header)
        except TypeError,AttributeError:
            pass
        return r
    def convert_output(self, res):
        if type(res) == list:
            for i, r in enumerate(res):
                res[i] = self.convert_one(r)
        else:
            res = self.convert_one(res)
        return res
    
    def __call__(self, f):
        def wrapper(*args, **kwargs):
            args = self.expand_args(args)
            res = f(*args, **kwargs)
            if self.preserve:
                res = self.convert_output(res)
            return res
        return wrapper

@manage_dtype()
def display(some_fits, ax=None):
    if ax == None:
        fig, ax = plt.subplots()

    im = ax.imshow(some_fits, cmap='gray', interpolation='nearest', aspect='auto', norm=LogNorm())
    fig.colorbar(im)
    return ax

@manage_dtype()
def row_avg(some_fits):
    return [np.nansum(some_fits[i,:]) for i in range(len(some_fits))]

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
    template_header = fits.PrimaryHDU(list_of_fits[0][0]).header
    headers = [f[1] for f in list_of_fits if f[1] != None]

    if len(headers) == 0:
        return template_header

    new_header = common_header(headers, template_header)
    return new_header

def avg_sexagesimal(headers, card):
    avg_decimal = 0.0
    N = 0
    for header in headers:
        sexages = header[card]
        try:
            h,m,s = [float(n) for n in filter(None, sexages.split(':'))]
        except:
            continue
        avg_decimal += h+m/60.+s/3600.
        N += 1
    if N == 0:
        raise TypeError('Header card '+card+' has non-sexagesimal values.')
    avg_decimal /= N

    avg_h = np.floor(avg_decimal)
    avg_decimal = (avg_decimal - avg_h)*60
    avg_m = np.floor(avg_decimal)
    avg_decimal = (avg_decimal - avg_m)*60
    avg_s = round(avg_decimal, 3)

    sexages_str = str(int(avg_h)).rjust(2,'0')+':'+str(int(avg_m)).rjust(2,'0')+':'+str(avg_s).split('.')[0].rjust(2,'0')+'.'+str(avg_s).split('.')[1].ljust(3,'0')
    return sexages_str

def middle_mjd(headers, card='MJD-OBS', exp_card='EXPTIME'):
    plus=all([h[card][0]=='+' for h in headers])
    mjd_list = [np.longdouble(h[card]) for h in headers]
    exptimes = [np.longdouble(h[exp_card])/86400. for h in headers]

    obs_start = mjd_list[np.argmin(mjd_list)]
    obs_end = (mjd_list+exptimes)[np.argmax(mjd_list+exptimes)]
    obs_end = np.max([mjd+exp for mjd,exp in zip(mjd_list, exptimes)])

    mid_mjd = np.mean([obs_start, obs_end])
    mid_mjd = str(mid_mjd)
    if plus:
        mid_mjd = '+'+mid_mjd
    return mid_mjd

def common_header(headers, template_header=None):
    headers = filter(None, headers)
    if len(headers) == 0:
        return None

    if template_header == None:
        new_header = fits.PrimaryHDU(np.array([])).header
    else:
        new_header = template_header

    #mjd_obs_comb = lambda headers: middle_mjd(headers,card='MJD-OBS')
    comb_funcs = {'MJD-OBS' : lambda headers : middle_mjd(headers,card='MJD-OBS'),
                  'JD' : lambda headers : middle_mjd(headers,card='JD'),
                  'RA' : lambda headers: avg_sexagesimal(headers,card='RA'),
                  'DEC' : lambda headers: avg_sexagesimal(headers,card='DEC'),
                  'UT' : lambda headers: avg_sexagesimal(headers,card='UT')}
    sample_header = headers[0]
    for card in sample_header.keys():
        if card == 'COMMENT' or card == 'HISTORY':
            continue

        useCard = all([card in h for h in headers]) 
        if not useCard:
            continue

        try:
            new_val = comb_funcs[card](headers)
        except KeyError:
            useCard = all([h[card] == sample_header[card] for h in headers])
            if not useCard:
                continue
            new_val = sample_header[card]

        new_header[card] = new_val

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

def get_dtype(*list_of_fits, **kwargs):
    '''
    Helper function for the combine function. It is used to determine the return type for
     combined fits files.
    ARGUMENTS:
        *list_of_fits - fits-like objects. These can be string paths to a fits file, astropy.io.fits
                        HDUList or PrimaryHDU objects, or ordinary 2D arrays.
        preserve - Boolean that dictates how the dtype will be determined. preserve=True means the
                   more complex dtypes will be favored, preserve=False means the less complex dtypes
                   will be favored.
    '''
    #Handle optional argument 'preserve'.
    preserve=False
    if "preserve" in kwargs:
        preserve = kwargs["preserve"]

    #Define functions that convert 2D arrays into various forms of FITS datastructures,
    # ordered from most complex to least.
    dtypes = [lambda data: fits.HDUList(fits.PrimaryHDU(data)),
              lambda data: fits.PrimaryHDU(data),
              lambda data: data]

    #Go through each fits object, determine what data type it is, and update the current
    # dtype to use, based on the value of 'preserve'.
    dtype_i = 0 if preserve else 2
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
        if preserve:
            dtype_i = min([dtype_i, dtype])
        else:
            dtype_i = max([dtype_i, dtype])

    return dtypes[dtype_i]

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

def pad_array(arr, pad_val, length, pad='back'):
    if pad == 'back':
        while len(arr) < length:
            arr = np.append(arr, pad_val)
    elif pad == 'front':
        while len(arr) < length:
            arr = np.insert(arr, 0, pad_val)
    return arr

@manage_dtype()
def mask_fits(some_fits, some_mask, maskval=1, fillval=0, reshape=False):
    if some_fits.shape != some_mask.shape:
        raise ValueError('Data and mask must be the same dimensions')
    masked = np.where(some_mask == maskval, some_fits, fillval)
    if reshape:
        height = len(masked)
        masked = masked[some_mask == maskval]
        width = masked.size/height
        masked = masked.reshape((height, len(masked)/height))
    return masked
