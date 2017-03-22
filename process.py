from driver import ensure_path
import numpy as np
import matplotlib.pyplot as plt
plt.ion()
from fitstools import combine, save_2darr, mask_fits, manage_dtype
from astropy.io import fits
from find_fibers import find_fibers
from throughput import make_throughput_map
from find_wvlsol import wvlsol
from extract_spectra import extract
from sys import stdout
import os
from os.path import exists

class output_log:
    def __init__(self, writer=stdout, log_path=None):
        self.log_path = log_path
        if log_path != None and not exists(log_path):
            f = open(log_path, 'w')
        self.writer = writer
        self.progress_str = ""
        self.message_str = ""
    def update(self):
        strg = self.progress_str+' '+self.message_str
        strg = strg.ljust(80)
        self.writer.write('\r'+strg)
        self.writer.flush()
    def edit_progress(self, new_str):
        self.progress_str = new_str
        self.update()
    def edit_message(self, new_str, add_to_log=True):
        self.message_str = new_str
        self.update()
        if add_to_log and self.log_path != None:
            f = open(self.log_path, 'a')
            f.write(self.message_str+'\n')
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

def get_num_groups(dname, recipe):
    if recipe == None:
        recipe = 'recipes/'+dname+'.recipe'
    f = open(recipe)
    recipe_lines = [line for line in filter(None, f.read().split('\n')) if line[0] != '#']
    f.close()
    group_nums = [int(line.split(',')[0]) for line in recipe_lines]
    num_groups = max(group_nums)
    return num_groups

def make_master_bias(dname, recipe=None, output=stdout):
    num_groups = get_num_groups(dname, recipe)
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
    '''
    #Load all other files make a bias-corrected version
    master_bias_data = master_bias[0].data
    for gnum in range(1, num_groups+1):
        group_recipes = get_recipes(dname, recipe, gnum=gnum)
        group_dir = calib_dir+'/group'+str(gnum)+'/'
        ensure_path(group_dir)
        fnames = []
        for r in group_recipes:
            fnames.extend(filter(None, r.split(',')[2].split(' ')))
        num_images = len(fnames)
        for i,fname in enumerate(fnames):
            output.write('\rBias correcting: Group '+str(gnum)+'/'+str(num_groups)+', Image '+str(i+1)+'/'+str(num_images))
            output.flush()
            f = fits.open('indata/'+dname+'/'+fname)
            f[0].data  = f[0].data - master_bias_data
            f[0].header['COMMENT'] = 'This image is bias corrected.'
            new_fname = fname.split('.')[0]+bc_suffix+'.fits'
            f.writeto(group_dir+'/'+new_fname, clobber=True)
            f.close()
    output.write('\n')
    '''

def bias_correct(image, bias=None):
    if bias == None:
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
    output.breakline()

def flat_dorecipe(r, dname, recipe, output=None):
    #Unpack info from recipe line.
    output.edit_message('Retrieving flat frames.')
    info = r.split(',')
    group_num = info[0]
    rtype = info[1]
    filenames = info[2].split(' ')
    use_fibers = [int(f_num) for f_num in info[3].split(' ')]

    #Make directory in calib for this group number
    indata_dir = 'indata/'+dname
    calib_dir = 'calib/'+dname+'/group'+group_num
    ensure_path(calib_dir)

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
    output.edit_message('Bias correcting flat frames.')
    flats = [fits.open(indata_dir+'/'+filename) for filename in filenames]
    flats = [bias_correct(flat, master_bias) for flat in flats]
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
    
    #Throughput correct all images from this pointing.
    #recipe_lines = get_recipes(dname, recipe, gnum=group_num)
    #for r in recipe_lines:
    #    info = r.split(',')
    #    recipe_filenames = info[2].split(' ')
    #    filenames = [fname.split('.')[0]+bc_suffix+'.fits' for fname in recipe_filenames]
    #    for fname in filenames:
    #        #Load file
    #        f = fits.open(group_dir+'/'+fname)
    #        f[0].data  = f[0].data / throughput_map
    #        f[0].header['COMMENT'] = 'This image is throughput corrected.'
    #        new_fname = fname.split(bc_suffix)[0]+tc_suffix+'.fits'
    #        f.writeto(group_dir+'/'+new_fname, clobber=True)


def process_thar(dname, recipe=None, output=stdout):
    thar_recipes = get_recipes(dname, recipe, rtype='comp')

    num_r = len(thar_recipes)
    for i,r in enumerate(thar_recipes):
        output.write('\rProcessing '+dname+' thar: '+str(i+1)+'/'+str(num_r))
        output.flush()
        thar_dorecipe(r, dname)
    output.write('\n')

def thar_dorecipe(r, dname):
    #Unpack info from recipe line.
    info = r.split(',')
    group_num = info[0]
    rtype = info[1]
    filenames = info[2].split(' ')
    use_fibers = [int(f_num) for f_num in info[3].split(' ')]

    #Make directory in calib for this group number
    indata_dir = 'indata/'+dname
    calib_dir = 'calib/'+dname+'/group'+group_num
    ensure_path(calib_dir)

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
       comp = fits.open(indata_dir+'/'+fname)
       comp[0].data = (comp[0].data - master_bias) / throughput_map
       wvlsol_map = wvlsol(comp, fiber_mask, use_fibers)
       fits.writeto(calib_dir+'/wvlsol.fits', wvlsol_map, clobber=True)

def process_sky(dname, recipe=None, output=stdout):
    sky_recipes = get_recipes(dname, recipe, rtype='sky')

    num_r = len(sky_recipes)
    for i,r in enumerate(sky_recipes):
        output.write('Processing '+dname+' sky: '+str(i+1)+'/'+str(num_r))
        output.flush()
        sky_dorecipe(r, dname)
    output.write('\n')

def sky_dorecipe(r, dname):
    #Unpack info from the recipe line.
    info = r.split(',')
    group_num = info[0]
    rtype = info[1]
    filenames = info[2].split(' ')
    use_fibers = [int(f_num) for f_num in info[3].split(' ')]

    #If there are no fibers to be reduced, stop.
    if len(use_fibers) == 0:
        return

    #Make directory in calib for this group number
    indata_dir = 'indata/'+dname
    calib_dir = 'calib/'+dname+'/group'+group_num
    ensure_path(calib_dir)

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
    skys = [fits.open(indata_dir+'/'+filename) for filename in filenames]
    master_sky = combine(*skys)
    master_sky = bias_correct(master_sky, master_bias)
    master_sky[0].data /= throughput_map
    master_sky.writeto(calib_dir+'/master_sky.fits', clobber=True)
    for sky in skys:
        sky.close()

    for fnum in use_fibers:
        wavelength, flux = extract(fiber_mask, fnum, master_sky, wvlsol_map)
        if True:
            fig, ax = plt.subplots()
            ax.set_title(str(fnum))
            ax.plot(wavelength, flux)


def process_target(dname, recipe=None, output=stdout):
    tar_recipes = get_recipes(dname, recipe, rtype='object')

    num_r = len(tar_recipes)
    for i,r in enumerate(tar_recipes):
        output.write('Processing '+dname+' targets: '+str(i+1)+'/'+str(num_r))
        output.flush()
        target_dorecipe(r, dname)
    output.write('\n')

def target_dorecipe(r, dname):
    #Unpack info from the recipe line.
    info = r.split(',')
    group_num = info[0]
    rtype = info[1]
    filenames = info[2].split(' ')
    use_fibers = [int(f_num) for f_num in filter(None, info[3].split(' '))]

    #If there are no fibers to be reduced, stop.
    if len(use_fibers) == 0:
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

    #Make master target frame
    tars = [fits.open(indata_dir+'/'+filename) for filename in filenames]
    master_tar = combine(*tars)
    master_tar = bias_correct(master_tar, master_bias)
    master_tar[0].data /= throughput_map
    master_tar.writeto(calib_dir+'/master_target_frame.fits', clobber=True)
    for f in tars:
        f.close()

    header = master_tar[0].header
    for fnum in use_fibers:
        tar_ID = filter(None, header['SLFIB'+str(fnum)].split(' '))[4]
        wavelength, flux = extract(fiber_mask, fnum, master_tar, wvlsol_map)
        np.savetxt(outdata_dir+'/'+tar_ID+'.txt', zip(wavelength, flux))

        if True:
            fig, ax = plt.subplots()
            ax.set_title(tar_ID)
            ax.plot(wavelength, flux)
