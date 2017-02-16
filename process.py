from driver import ensure_path
from fitstools import combine, save_2darr
from astropy.io import fits
from find_fibers import find_fibers

def process_flat(dname, recipe=None):
    if recipe == None:
        recipe = 'recipes/'+dname+'.recipe'
    
    #Extract lines from recipe file that pertain to flat fields.
    recipe_lines = [line for line in filter(None, open(recipe).read().split('\n')) if line[0] != '#']
    flat_recipes = [line for line in recipe_lines if len(line.split(',')) > 2 and line.split(',')[1] == 'flat']

    for r in flat_recipes:
        flat_dorecipe(r, dname)

def flat_dorecipe(r, dname):
    info = r.split(',')
    group_num = info[0]
    print('Processing group', group_num, 'flats.')
    rtype = info[1]
    filenames = info[2].split(' ')
    use_fibers = [int(f_num) for f_num in info[3].split(' ')]
    print use_fibers, ':D:'

    #Make directory in calib for this group number
    group_dir = 'calib/'+dname+'/group'+group_num+'/'
    ensure_path(group_dir)

    #Median combine flats, save to calib directory
    flats = [fits.open('indata/'+dname+'/'+filename) for filename in filenames]
    master_flat = combine(*flats)
    master_flat.writeto(group_dir+'/master_flat.fits', clobber=True)

    #Find fibers, make fiber mask, save to directory
    print use_fibers, '2312dsx'
    fiber_mask = find_fibers(master_flat, use_fibers)
    fits.writeto(group_dir+'fiber_mask.fits', fiber_mask, clobber=True)

