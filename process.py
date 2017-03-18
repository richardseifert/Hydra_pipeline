from driver import ensure_path
import numpy as np
from fitstools import combine, save_2darr, mask_fits
from astropy.io import fits
from find_fibers import find_fibers
from throughput import make_throughput_map
from find_wvlsol import wvlsol
from extract_spectra import extract
from sys import stdout

#Define suffixes for bias-corrected and throghput-corrected images.
bc_suffix = '_bc'
tc_suffix = '_tc'

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

def process_bias(dname, recipe=None, output=stdout):
    num_groups = get_num_groups(dname, recipe)
    bias_recipes = get_recipes(dname, recipe, rtype='zero')

    calib_dir = 'calib/'+dname
    ensure_path(calib_dir+'/')
    
    filenames = []
    for r in bias_recipes:
        info = r.split(',')
        fnames = info[2].split(' ')
        filenames.extend(fnames)

    #Load files and create master bias
    biases = [fits.open('indata/'+dname+'/'+fname) for fname in filenames]
    master_bias = combine(*biases, method='mean')
    master_bias.writeto(calib_dir+'/master_bias.fits', clobber=True)
    for f in biases:
        f.close()

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

def process_flat(dname, recipe=None, output=stdout):
    flat_recipes = get_recipes(dname, recipe, rtype='flat')

    num_r = len(flat_recipes)
    for i,r in enumerate(flat_recipes):
        output.write('\rProcessing '+dname+' flats: '+str(i+1)+'/'+str(num_r))
        output.flush()
        flat_dorecipe(r, dname, recipe)
    output.write('\n')

def flat_dorecipe(r, dname, recipe):
    info = r.split(',')
    group_num = info[0]
    rtype = info[1]
    recipe_filenames = info[2].split(' ')
    filenames = [fname.split('.')[0]+bc_suffix+'.fits' for fname in recipe_filenames]
    use_fibers = [int(f_num) for f_num in info[3].split(' ')]

    #Make directory in calib for this group number
    group_dir = 'calib/'+dname+'/group'+group_num
    ensure_path(group_dir)

    #Median combine flats, save to calib directory
    flats = [fits.open(group_dir+'/'+filename) for filename in filenames]
    master_flat = combine(*flats)
    master_flat.writeto(group_dir+'/master_flat.fits', clobber=True)
    for flat in flats:
        flat.close()

    #Find fibers, make fiber mask, save to directory
    fiber_mask = find_fibers(master_flat, use_fibers)
    fits.writeto(group_dir+'/fiber_mask.fits', fiber_mask, clobber=True)
 
    #Generate a fiber thoughput map.
    throughput_map = make_throughput_map(fiber_mask, master_flat)
    fits.writeto(group_dir+'/throughput_map.fits', throughput_map, clobber=True)
    
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
        output.write('Processing '+dname+' thar: '+str(i+1)+'/'+str(num_r))
        output.flush()
        thar_dorecipe(r, dname)
    output.write('\n')

def thar_dorecipe(r, dname):
    info = r.split(',')
    group_num = info[0]
    rtype = info[1]
    recipe_filenames = info[2].split(' ')
    filenames = [fname.split('.')[0]+bc_suffix+'.fits' for fname in recipe_filenames]
    use_fibers = [int(f_num) for f_num in info[3].split(' ')]

    #Make directory in calib for this group number
    group_dir = 'calib/'+dname+'/group'+group_num
    ensure_path(group_dir+'/')

    tm = fits.open(group_dir+'/throughput_map.fits')
    throughput_map = tm[0].data
    tm.close()
    fm = fits.open(group_dir+'/fiber_mask.fits')
    fiber_mask = fm[0].data
    fm.close()
    for fname in filenames:
       comp = fits.open(group_dir+'/'+fname)
       comp[0].data = comp[0].data / throughput_map
       new_fname = fname.split(bc_suffix)[0]+tc_suffix+'.fits'
       comp.writeto(group_dir+'/'+new_fname, clobber=True)
       wvlsol_map = wvlsol(comp, fiber_mask, use_fibers)
       fits.writeto(group_dir+'/wvlsol.fits', wvlsol_map, clobber=True)

def process_sky(dname, recipe=None, output=stdout):
    sky_recipes = get_recipes(dname, recipe, rtype='sky')

    num_r = len(sky_recipes)
    for i,r in enumerate(sky_recipes):
        output.write('Processing '+dname+' sky: '+str(i+1)+'/'+str(num_r))
        output.flush()
        sky_dorecipe(r, dname)
    output.write('\n')

def sky_dorecipe(r, dname):
    info = r.split(',')
    group_num = info[0]
    rtype = info[1]
    recipe_filenames = info[2].split(' ')
    filenames = [fname.split('.')[0]+bc_suffix+'.fits' for fname in recipe_filenames]
    use_fibers = [int(f_num) for f_num in info[3].split(' ')]

    #Make directory in calib for this group number
    group_dir = 'calib/'+dname+'/group'+group_num
    ensure_path(group_dir+'/')

    fm = fits.open(group_dir+'/fiber_mask.fits')
    fiber_mask = fm[0].data
    fm.close()

    ws = fits.open(group_dir+'/wvlsol.fits')
    wvlsol_map = ws[0].data
    ws.close()

    #Make master sky frame
    skys = [fits.open(group_dir+'/'+filename) for filename in filenames]
    master_sky = combine(*skys)
    master_sky.writeto(group_dir+'/master_sky.fits', clobber=True)
    for sky in skys:
        sky.close()

    for fnum in use_fibers:
        extract(fiber_mask, fnum, master_sky, wvlsol_map)


def sky_obj_dorecipe(r, dname):
    pass 
