from driver import ensure_path
import numpy as np
import matplotlib.pyplot as plt
plt.ion()
from fitstools import combine, save_2darr, mask_fits, manage_dtype
from astropy.io import fits
from find_fibers import find_fibers
from throughput import make_throughput_map
from find_wvlsol import wvlsol
from extract_spectra import spectrum, extract, interp_mean
from sys import stdout
import os
from os.path import exists
import time

class output_log:
    def __init__(self, writer=stdout, log_path=None):
        self.log_path = log_path
        if log_path != None and not exists(log_path):
            f = open(log_path, 'w')
        self.writer = writer
        self.progress_str = ""
        self.message_str = ""
        self.coverlen = 0
    def update(self):
        strg = self.progress_str+' '+self.message_str
        strg = strg.ljust(self.coverlen)
        self.coverlen = len(strg)
        self.writer.write('\r'+strg)
        self.writer.flush()
    def edit_progress(self, new_str):
        self.progress_str = new_str
        self.update()
    def edit_message(self, new_str, add_to_log=True):
        self.message_str = new_str
        self.update()
        if add_to_log and self.log_path != None:
            dt_str = time.strftime("%Y-%m-%dT%H:%M:%S")
            f = open(self.log_path, 'a')
            f.write('['+dt_str+'] '+self.message_str+'\n')
    def linebreak():
        self.writer.write('\n')

def get_recipes(dname, recipe=None, gnum=None, rtype=None):
    if recipe == None:
        recipe = 'recipes/'+dname+'.recipe'
    
    #Extract lines from recipe file that pertain to flat fields.
    recipe_lines = [line for line in filter(None, open(recipe).read().split('\n')) if line[0] != '#']
    if rtype != None:
        recipe_lines = [line for line in recipe_lines if len(line.split(',')) > 2 and line.split(',')[1] == rtype]
    if gnum != None:
        gnum = int(gnum)
        recipe_lines = [line for line in recipe_lines if len(line.split(',')) > 1 and int(line.split(',')[0]) == gnum]
    return recipe_lines

def make_master_bias(dname, recipe=None, output=stdout):
    bias_recipes = get_recipes(dname, recipe, rtype='zero')

    calib_dir = 'calib/'+dname
    ensure_path(calib_dir+'/')
    
    if not os.path.exists(calib_dir+'/'+'master_bias.fits'):
        filenames = []
        for r in bias_recipes:
            info = r.split(',')
            fnames = info[2].split(' ')
            filenames.extend(fnames)

        if len(filenames) > 0:
            #Load files and create master bias
            biases = [fits.open('indata/'+dname+'/'+fname) for fname in filenames]
            master_bias = combine(*biases, method='mean')
            master_bias.writeto(calib_dir+'/master_bias.fits', clobber=True)
            for f in biases:
                f.close()

def bias_correct(image, bias=None):
    if type(bias) == type(None):
        #Use bias overscan region of the image.
        #Figure out later.
        return image
    else:
        image[0].data = image[0].data - bias
        image[0].header['COMMENT'] = 'Bias corrected.'
        return image
    
def process_flat(dname, recipe=None, output=stdout):
    output = output_log(log_path='calib/'+dname+'/output.log')
    flat_recipes = get_recipes(dname, recipe, rtype='flat')
    num_r = len(flat_recipes)

    for i,r in enumerate(flat_recipes):
        output.edit_progress('Processing '+dname+' flats: '+str(i+1)+'/'+str(num_r)+' |')
        output.edit_message("", add_to_log=False)
        #output.write('\rProcessing '+dname+' flats: '+str(i+1)+'/'+str(num_r))
        #output.flush()
        flat_dorecipe(r, dname, recipe, output)

def flat_dorecipe(r, dname, recipe, output=None):
    #Unpack info from recipe line.
    output.edit_message('Reading recipe.')
    info = r.split(',')
    group_num = info[0]
    rtype = info[1]
    filenames = info[2].split(' ')
    use_fibers = [int(f_num) for f_num in info[3].split(' ')]

    #Make directory in calib for this group number
    indata_dir = 'indata/'+dname
    calib_dir = 'calib/'+dname+'/group'+group_num
    ensure_path(calib_dir+'/')

    #Obtain master bias frame
    make_master_bias(dname, recipe)
    mb_path = 'calib/'+dname+'/master_bias.fits'
    if os.path.exists(mb_path):
        mb = fits.open(mb_path)
        master_bias = mb[0].data
        mb.close()
    else:
        master_bias = None

    #Load flat frames, bias correct, and median combine flats, save to calib directory
    output.edit_message('Loading flat frames.')
    flats = [fits.open(indata_dir+'/'+filename) for filename in filenames]
    output.edit_message('Bias correcting flat frames.')
    flats = [bias_correct(flat, master_bias) for flat in flats]
    output.edit_message('Median combining flat frames.')
    master_flat = combine(*flats)
    for flat in flats:
        flat.close()
    mf_path = calib_dir+'/master_flat.fits'
    master_flat.writeto(mf_path, clobber=True)
    output.edit_message('Master flat frame saved at '+mf_path)

    #Find fibers, make fiber mask, save to directory
    output.edit_message('Locating fibers.')
    fiber_mask = find_fibers(master_flat, use_fibers)
    fm_path = calib_dir+'/fiber_mask.fits'
    fits.writeto(fm_path, fiber_mask, clobber=True)
    output.edit_message('Fiber mask saved at '+fm_path)
 
    #Generate a fiber thoughput map.
    output.edit_message('Generating throughput map.')
    throughput_map = make_throughput_map(fiber_mask, master_flat)
    tm_path = calib_dir+'/throughput_map.fits'
    fits.writeto(tm_path, throughput_map, clobber=True)
    output.edit_message('Throughput map saved at '+tm_path)

def process_thar(dname, recipe=None, output=stdout, **kwargs):
    output = output_log(log_path='calib/'+dname+'/output.log')
    thar_recipes = get_recipes(dname, recipe, rtype='comp')

    num_r = len(thar_recipes)
    for i,r in enumerate(thar_recipes):
        output.edit_progress('Processing '+dname+' thar: '+str(i+1)+'/'+str(num_r)+' |')
        output.edit_message('')
        thar_dorecipe(r, dname, output, **kwargs)

def thar_dorecipe(r, dname, output=None, **kwargs):
    output.edit_message('Reading thar recipe.')
    #Unpack info from recipe line.
    info = r.split(',')
    group_num = info[0]
    rtype = info[1]
    filenames = info[2].split(' ')
    use_fibers = [int(f_num) for f_num in info[3].split(' ')]

    #Make directory in calib for this group number
    indata_dir = 'indata/'+dname
    calib_dir = 'calib/'+dname+'/group'+group_num
    ensure_path(calib_dir+'/')

    #Obtain master bias frame
    mb_path = 'calib/'+dname+'/master_bias.fits'
    if os.path.exists(mb_path):
        mb = fits.open(mb_path)
        master_bias = mb[0].data
        mb.close()
    else:
        master_bias = None

    tm = fits.open(calib_dir+'/throughput_map.fits')
    throughput_map = tm[0].data
    tm.close()

    fm = fits.open(calib_dir+'/fiber_mask.fits')
    fiber_mask = fm[0].data
    fm.close()
    
    for fname in filenames:
        output.edit_message('Finding wavelength solution from '+fname)
        comp = fits.open(indata_dir+'/'+fname)
        comp[0].data = (comp[0].data - master_bias) / throughput_map
        wvlsol_map = wvlsol(comp, fiber_mask, use_fibers, **kwargs)
        ws_path = calib_dir+'/wvlsol.fits'
        fits.writeto(ws_path, wvlsol_map, clobber=True)
        output.edit_message('Wavelength solution saved at '+ws_path)

def process_sky(dname, recipe=None, output=stdout):
    output = output_log(log_path='calib/'+dname+'/output.log')
    sky_recipes = get_recipes(dname, recipe, rtype='sky')

    num_r = len(sky_recipes)
    for i,r in enumerate(sky_recipes):
        output.edit_progress('Processing '+dname+' sky: '+str(i+1)+'/'+str(num_r)+' |')
        output.edit_message('')
        sky_dorecipe(r, dname, output)

def sky_dorecipe(r, dname, output=None):
    #Unpack info from the recipe line.
    output.edit_message('Reading sky recipe.')
    info = r.split(',')
    group_num = info[0]
    rtype = info[1]
    filenames = info[2].split(' ')
    use_fibers = [int(f_num) for f_num in info[3].split(' ')]

    #If there are no fibers to be reduced, stop.
    if len(use_fibers) == 0:
        output.edit_message('No fibers to be reduced. Stopping.')
        return

    #Make directory in calib for this group number
    indata_dir = 'indata/'+dname
    calib_dir = 'calib/'+dname+'/group'+group_num
    ensure_path(calib_dir+'/')

    #Obtain master bias frame, if it exists.
    mb_path = 'calib/'+dname+'/master_bias.fits'
    if os.path.exists(mb_path):
        mb = fits.open(mb_path)
        master_bias = mb[0].data
        mb.close()
    else:
        master_bias = None

    #Load in the fiber mask.
    fm = fits.open(calib_dir+'/fiber_mask.fits')
    fiber_mask = fm[0].data
    fm.close()

    #Load in the throughput map.
    tm = fits.open(calib_dir+'/throughput_map.fits')
    throughput_map = tm[0].data
    tm.close()

    #Load in the wavelength solution.
    ws = fits.open(calib_dir+'/wvlsol.fits')
    wvlsol_map = ws[0].data
    ws.close()

    #Make master sky frame
    output.edit_message('Loading sky frames.')
    skys = [fits.open(indata_dir+'/'+filename) for filename in filenames]
    output.edit_message('Median combining sky frames.')
    master_sky = combine(*skys)
    for sky in skys:
        sky.close()
    output.edit_message('Bias correcting master sky frame.')
    master_sky = bias_correct(master_sky, master_bias)
    output.edit_message('Throughput correcting master sky frame.')
    master_sky[0].data /= throughput_map
    ms_path = calib_dir+'/master_sky.fits'
    master_sky.writeto(ms_path, clobber=True)
    output.edit_message('Master sky frame saved at '+ms_path)

    sky_specs = []
    for fnum in use_fibers:
        output.edit_message('Extracting sky spectrum from fiber '+str(fnum))
        sky_spec = extract(fiber_mask, fnum, master_sky, wvlsol_map)
        sky_specs.append(sky_spec)
    output.edit_message('Producing master sky spectrum')
    master_sky_spec = interp_mean(*sky_specs)
    master_sky_spec.plot()

    mss_path = calib_dir+'/master_sky_spec.dat'
    master_sky_spec.save(mss_path)
    output.edit_message('Master sky spectrum saved at '+mss_path)


def process_target(dname, recipe=None, output=stdout):
    output = output_log(log_path='calib/'+dname+'/output.log')
    tar_recipes = get_recipes(dname, recipe, rtype='object')

    num_r = len(tar_recipes)
    for i,r in enumerate(tar_recipes):
        output.edit_progress('Processing '+dname+' targets: '+str(i+1)+'/'+str(num_r)+' |')
        output.edit_message('')
        target_dorecipe(r, dname, output)

def target_dorecipe(r, dname, output=None):
    #Unpack info from the recipe line.
    output.edit_message('Reading target recipe.')
    info = r.split(',')
    group_num = info[0]
    rtype = info[1]
    filenames = info[2].split(' ')
    use_fibers = [int(f_num) for f_num in filter(None, info[3].split(' '))]

    #If there are no fibers to be reduced, stop.
    if len(use_fibers) == 0:
        output.edit_message('No fibers to be reduced. Stopping.')
        return

    #Define paths to indata, calib, and outdata
    indata_dir = 'indata/'+dname
    calib_dir = 'calib/'+dname+'/group'+group_num
    ensure_path(calib_dir+'/')
    outdata_dir = 'outdata/'+dname+'/group'+group_num
    ensure_path(outdata_dir+'/')

    #Obtain master bias frame, if it exists.
    mb_path = 'calib/'+dname+'/master_bias.fits'
    if os.path.exists(mb_path):
        mb = fits.open(mb_path)
        master_bias = mb[0].data
        mb.close()
    else:
        master_bias = None
        
    #Load in the fiber mask
    fm = fits.open(calib_dir+'/fiber_mask.fits')
    fiber_mask = fm[0].data
    fm.close()

    #Load in the wavelength solution.
    ws = fits.open(calib_dir+'/wvlsol.fits')
    wvlsol_map = ws[0].data
    ws.close()

    #Load in the throughput map.
    tm = fits.open(calib_dir+'/throughput_map.fits')
    throughput_map = tm[0].data
    tm.close()

    #Load in master sky spectrum.
    mss = np.loadtxt(calib_dir+'/master_sky_spec.dat')
    master_sky_spec = spectrum(mss[:,0], mss[:,1])
    master_sky_spec.plot()

    #Make master target frame
    output.edit_message('Loading target frames.')
    tars = [fits.open(indata_dir+'/'+filename) for filename in filenames]
    output.edit_message('Median combining target frames.')
    master_tar = combine(*tars)
    for f in tars:
        f.close()
    output.edit_message('Bias correcting master target frame.')
    master_tar = bias_correct(master_tar, master_bias)
    output.edit_message('Throughput correcting master target frame.')
    master_tar[0].data /= throughput_map
    mt_path = calib_dir+'/master_target_frame.fits'
    master_tar.writeto(mt_path, clobber=True)
    output.edit_message('Master target frame saved at '+mt_path)

    header = master_tar[0].header
    for fnum in use_fibers:
        output.edit_message('Extracting target spectrum from fiber '+str(fnum))
        tar_ID = filter(None, header['SLFIB'+str(fnum)].split(' '))[4]
        tar_spec = extract(fiber_mask, fnum, master_tar, wvlsol_map)
        ts_path = outdata_dir+'/'+tar_ID+'.txt'
        ax = tar_spec.plot(color='red')
        tar_spec = tar_spec - master_sky_spec
        tar_spec.plot(ax=ax)
        tar_spec.save(ts_path)
        output.edit_message('Target spectrum saved as '+ts_path)

        if True:
            tar_spec.plot()
