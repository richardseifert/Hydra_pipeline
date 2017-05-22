from group_obs import group_obs

class recipe(object):
	def __init__(self, strg, delimiter=',', list_delimiter=' '):
		info = strg.split(delimiter)
		self.gnum = info[0]
		self.rtype = info[1]
		self.filenames = info[2].split(list_delimiter)
		self.fibers = [int(s) for s in info[3].split(list_delimiter)]

def load_recipes(filename):
	f = open(filename)
	recipe_lines = f.read().split('\n')
	return [recipe(line) for line in recipe_lines]
	
def make_recipe(direc, savepath):
    recipe = open(savepath, 'w')
    recipe.write('#GROUP_NUM,IMAGE_TYPE,FILE_NAMES,FIBER_NUMBERS\n')



    biases, groups = group_obs(direc)
    #List bias info
    file_names = [f.filename().split('/')[-1] for f in biases]
    recipe.write('-1,zero,'+' '.join(file_names)+',\n')
    for f in biases:
        f.close()

    for i,g in enumerate(groups):
        group_num = str(i+1)
        if len(g.get_images()) == 0:
            continue
        sample_header = g.get_images()[0][0].header

        #List flat info
        img_type = 'flat'
        file_names = [f.filename().split('/')[-1] for f in g.images[img_type]]
        use_fiber = lambda slfib_s: len(slfib_s.split()) > 1 and slfib_s.split()[1] in ['0','1','3'] #0 is Sky, 1 is Object, 3 is Random
        fiber_nums = get_use_fiber_nums(sample_header, use_fiber)
        recipe.write(group_num+','+img_type+','+' '.join(file_names)+','+' '.join([str(n) for n in fiber_nums])+'\n')

        #List comp info
        img_type = 'comp'
        file_names = [f.filename().split('/')[-1] for f in g.images[img_type]]
        use_fiber = lambda slfib_s: len(slfib_s.split()) > 1 and slfib_s.split()[1] in ['0','1','3'] #0 is Sky, 1 is Object, 3 is Random
        fiber_nums = get_use_fiber_nums(sample_header, use_fiber)
        recipe.write(group_num+','+img_type+','+' '.join(file_names)+','+' '.join([str(n) for n in fiber_nums])+'\n')

        #List sky info
        img_type = 'sky'
        file_names = [f.filename().split('/')[-1] for f in g.images['object']]
        use_fiber = lambda slfib_s: len(slfib_s.split()) > 1 and slfib_s.split()[1] in ['0','3'] #0 is Sky, 3 is Random
        fiber_nums = get_use_fiber_nums(sample_header, use_fiber)
        recipe.write(group_num+','+img_type+','+' '.join(file_names)+','+' '.join([str(n) for n in fiber_nums])+'\n')
        
        #List object info
        img_type = 'object'
        file_names = [f.filename().split('/')[-1] for f in g.images[img_type]]
        use_fiber = lambda slfib_s: len(slfib_s.split()) > 1 and slfib_s.split()[1] in ['1'] #1 is Object
        fiber_nums = get_use_fiber_nums(sample_header, use_fiber)
        recipe.write(group_num+','+img_type+','+' '.join(file_names)+','+' '.join([str(n) for n in fiber_nums])+'\n')

        g.close_files()

    recipe.close()

def get_use_fiber_nums(header, use_condition):
    use_fiber_nums = []
    n = 1
    k = 'SLFIB'+str(n)
    while k in header:
        if use_condition(header[k]):
            use_fiber_nums.append(n)
        n += 1
        k = 'SLFIB'+str(n)
    return use_fiber_nums
