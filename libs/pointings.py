import glob
from datetime import datetime
from astropy.io import fits

default_conditions = {'object':lambda im: im[0].header['IMGTYPE']=='object',
                      'flat':lambda im: im[0].header['IMGTYPE']=='flat',
                      'comp':lambda im: im[0].header['IMGTYPE']=='comp'}

def load_conditions(fname):
    f = open(fname)
    lines = filter(None, f.read().split('\n'))
    f.close()

    conditions = {}
    for line in lines:
        recipe = line.split(':')[0].lstrip(' ').rstrip(' ')
        c = line.split(':')[1].replace('header[', 'im[0].header[')
        exec('f = lambda im: '+c)
        conditions[recipe] = f

    return conditions


class pointing:
    def __init__(self, conditions=default_conditions):
        self.images = {}
        self.conditions = conditions
    def add_image(self, im):
        recipe = get_recipe(im, conditions=self.conditions)
        if recipe == None:
            return
        if not recipe in self.images.keys():
            self.images[recipe] = []
        self.images[recipe].append(im)
    def get_images(self, img_type=None):
        if img_type == None:
            imgs = []
            for k in self.images.keys():
                imgs.extend(self.images[k])
            return imgs
        return self.images[img_type]
    def close_files(self):
        for k in self.images.keys():
            for f in self.images[k]:
                f.close()

def convert_timestr(s):
    decimal_seconds = '.'+s.split('.')[-1]
    micro_seconds = str(int(float(decimal_seconds)*1.0e6)).rjust(6, '0')
    s = s.split('.')[0]+':'+micro_seconds
    return datetime.strptime(s, '%Y-%m-%dT%H:%M:%S:%f')

def get_recipe(im, conditions=default_conditions):
    for recipe in conditions:
        if conditions[recipe](im):
            return recipe
    return None

#group images with identical fiber pointings that were taken chronologically
def get_pointings(direc, conditions=default_conditions):
    #Get list of fits files in direc
    if direc[-1] != '/':
        direc += '/'
    fnames = glob.glob(direc+'*.fits')+glob.glob(direc+'*.fit')

    #Load fits files
    images = [fits.open(fpath) for fpath in fnames]

    #Sort files chronologically by DATE-OBS 
    images = sorted(images, key=lambda im: convert_timestr(im[0].header['DATE-OBS']))

    pointings = [pointing(conditions=conditions)]
    prev_fiberConfig = ''
    for im in images:
        h = im[0].header
        n = 1
        fiberConfig = ''
        while 'SLFIB'+str(n) in h:
            fiberConfig += h['SLFIB'+str(n)]
            n += 1
        if fiberConfig != prev_fiberConfig and len(pointings[-1].get_images()) > 0:
            pointings.append(pointing(conditions=conditions))
        pointings[-1].add_image(im)
        prev_fiberConfig = fiberConfig

    return pointings

