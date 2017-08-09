from pointings import load_conditions, get_pointings
import os

class recipe(object):
	def __init__(self, strg, delimiter=',', list_delimiter=' '):
		info = strg.split(delimiter)
		self.pnum = info[0]
		self.rtype = info[1]
		self.filenames = filter(None, info[2].split(list_delimiter))
		self.fibers = [int(s) for s in filter(None, info[3].split(list_delimiter))]

def load_recipes(filename):
	f = open(filename)
	recipe_lines = filter(None, f.read().split('\n'))
	return [recipe(line) for line in recipe_lines if line[0]!='#']
	
def make_recipe(direc, savepath):
    recipe = open(savepath, 'w')
    recipe.write('#POINTING,IMAGE_TYPE,FILE_NAMES,FIBER_NUMBERS\n')

    if os.path.exists('recipes/recipe_config/recipe_conditions.txt'):
        conditions = load_conditions('recipes/recipe_config/recipe_conditions.txt')
        pointings = get_pointings(direc, conditions)
    else:
        pointings = get_pointings(direc)
    
    for i,p in enumerate(pointings):
        pnum = str(i+1)
        if len(p.get_images()) == 0:
            continue
        sample_header = p.get_images()[0][0].header

        #List flat info
        img_type = 'flat'
        if img_type in p.images:
            file_names = [f.filename().split('/')[-1] for f in p.images[img_type]]
            use_fiber = lambda slfib_s: len(slfib_s.split()) > 1 and slfib_s.split()[1] in ['0','1','3'] #0 is Sky, 1 is Object, 3 is Random
            fiber_nums = get_use_fiber_nums(sample_header, use_fiber)
            recipe.write(pnum+','+img_type+','+' '.join(file_names)+','+' '.join([str(n) for n in fiber_nums])+'\n')

        #List comp info
        img_type = 'comp'
        if img_type in p.images:
            file_names = [f.filename().split('/')[-1] for f in p.images[img_type]]
            use_fiber = lambda slfib_s: len(slfib_s.split()) > 1 and slfib_s.split()[1] in ['0','1','3'] #0 is Sky, 1 is Object, 3 is Random
            fiber_nums = get_use_fiber_nums(sample_header, use_fiber)
            recipe.write(pnum+','+img_type+','+' '.join(file_names)+','+' '.join([str(n) for n in fiber_nums])+'\n')

        #List sky info
        img_type = 'sky'
        if 'object' in p.images:
            file_names = [f.filename().split('/')[-1] for f in p.images['object']]
            use_fiber = lambda slfib_s: len(slfib_s.split()) > 1 and slfib_s.split()[1] in ['0','3'] #0 is Sky, 3 is Random
            fiber_nums = get_use_fiber_nums(sample_header, use_fiber)
            recipe.write(pnum+','+img_type+','+' '.join(file_names)+','+' '.join([str(n) for n in fiber_nums])+'\n')
        
        #List object info
        img_type = 'object'
        if img_type in p.images:
            file_names = [f.filename().split('/')[-1] for f in p.images[img_type]]
            use_fiber = lambda slfib_s: len(slfib_s.split()) > 1 and slfib_s.split()[1] in ['1'] #1 is Object
            fiber_nums = get_use_fiber_nums(sample_header, use_fiber)
            recipe.write(pnum+','+img_type+','+' '.join(file_names)+','+' '.join([str(n) for n in fiber_nums])+'\n')
            
        #List skyflat info
        img_type = 'skyflat'
        if img_type in p.images:
            file_names = [f.filename().split('/')[-1] for f in p.images[img_type]]
            use_fiber = lambda slfib_s: len(slfib_s.split()) > 1 and slfib_s.split()[1] in ['0','3'] #0 is Sky, 3 is Random
            fiber_nums = get_use_fiber_nums(sample_header, use_fiber)
            recipe.write(pnum+','+img_type+','+' '.join(file_names)+','+' '.join([str(n) for n in fiber_nums])+'\n')

        recipe.write('\n')
        p.close_files()

    recipe.close()

fname = 'recipes/recipe_config/ignore_fibers.dat'
f = open(fname)
ignore_fibers = [int(n) for n in filter(None,f.read().split(','))]
f.close()

def get_use_fiber_nums(header, use_condition):
    use_fiber_nums = []
    n = 1
    k = 'SLFIB'+str(n)
    while k in header:
        if use_condition(header[k]) and not n in ignore_fibers:
            use_fiber_nums.append(n)
        n += 1
        k = 'SLFIB'+str(n)
    return use_fiber_nums
