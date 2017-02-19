from driver import ensure_path
import numpy as np
from fitstools import combine, save_2darr
from astropy.io import fits
from find_fibers import find_fibers
from throughput import make_throughput_map
from sys import stdout

def get_recipes(dname, recipe=None, rtype=None):
    if recipe == None:
        recipe = 'recipes/'+dname+'.recipe'
    
    #Extract lines from recipe file that pertain to flat fields.
    recipe_lines = [line for line in filter(None, open(recipe).read().split('\n')) if line[0] != '#']
    if rtype != None:
        rtype_recipes = [line for line in recipe_lines if len(line.split(',')) > 2 and line.split(',')[1] == rtype]
        return rtype_recipes
    return recipe_lines

def process_bias(dname, recipe=None, output=stdout):
    bias_recipes = get_recipes(dname, recipe, rtype='zero')
    
    filenames = []
    for r in bias_recipes:
        info = r.split(',')
        fnames = info[2].split(' ')
        filenames.extend(fnames)

    calib_imgs_dir = 'calib/'+dname+'/cal_images/'
    ensure_path(calib_imgs_dir)

    #Load files and create master bias
    biases = [fits.open('indata/'+dname+'/'+fname) for fname in filenames]
    master_bias = combine(*biases, method='mean')
    master_bias.writeto(calib_imgs_dir+'master_bias.fits', clobber=True)
    for f in biases:
        f.close()

    #Load all other files make a bias-corrected version
    all_recipes = get_recipes(dname, recipe)
    all_filenames = []
    for r in all_recipes:
        info = r.split(',')
        fnames = filter(None, info[2].split(' '))
        for fname in fnames:
            if not fname in all_filenames:
                all_filenames.append(fname)
    master_bias_data = master_bias[0].data
    num_imgs = len(all_filenames)
    for i,fname in enumerate(all_filenames):
        output.write('\rBias correcting '+dname+' images: '+str(i+1)+'/'+str(num_imgs))
        output.flush()
        f = fits.open('indata/'+dname+'/'+fname)
        f[0].data  = f[0].data - master_bias_data
        f[0].header['COMMENT'] = 'This image is bias corrected.'
        new_fname = fname.split('.')[0]+'_b.fits'
        f.writeto('calib/'+dname+'/cal_images/'+new_fname, clobber=True)
        f.close()
    output.write('\n')

    #Save bias-corrected images in calib/dname/calib_images

def process_flat(dname, recipe=None, output=stdout):
    flat_recipes = get_recipes(dname, recipe, rtype='flat')

    num_r = len(flat_recipes)
    for i,r in enumerate(flat_recipes):
        output.write('\rProcessing '+dname+' flats: '+str(i+1)+'/'+str(num_r))
        output.flush()
        flat_dorecipe(r, dname)
    output.write('\n')

def flat_dorecipe(r, dname):
    info = r.split(',')
    group_num = info[0]
    rtype = info[1]
    recipe_filenames = info[2].split(' ')
    filenames = [fname.split('.')[0]+'_b.fits' for fname in recipe_filenames]
    use_fibers = [int(f_num) for f_num in info[3].split(' ')]

    #Make directory in calib for this group number
    group_dir = 'calib/'+dname+'/group'+group_num+'/'
    ensure_path(group_dir)

    #Median combine flats, save to calib directory
    flats = [fits.open('calib/'+dname+'/cal_images/'+filename) for filename in filenames]
    master_flat = combine(*flats)
    master_flat.writeto(group_dir+'master_flat.fits', clobber=True)

    #Find fibers, make fiber mask, save to directory
    fiber_mask = find_fibers(master_flat, use_fibers)
    fits.writeto(group_dir+'fiber_mask.fits', fiber_mask, clobber=True)
    
    #Generate a fiber thoughput map.
    throughput_map = make_throughput_map(fiber_mask, master_flat)
    fits.writeto(group_dir+'throughput_map.fits', throughput_map, clobber=True)
