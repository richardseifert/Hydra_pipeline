import glob
from datetime import datetime
from astropy.io import fits

class obs_group:
    def __init__(self, obj=[], flat=[], comp=[]):
        self.images = {'object' :list(obj), 'flat':list(flat), 'comp':list(comp)}
    def add_image(self, im, imgtype):
        try:
            self.images[imgtype].append(im)
        except KeyError:
            print 'Unrecognize d image type:',imgtype
    def get_images(self, img_type=None):
        if img_type == None:
            imgs = []
            for k in self.images.keys():
                imgs.extend(self.images[k])
            return imgs
        return self.images[img_type]

def convert_timestr(s):
    decimal_seconds = '.'+s.split('.')[-1]
    micro_seconds = str(int(float(decimal_seconds)*1.0e6)).rjust(6, '0')
    s = s.split('.')[0]+':'+micro_seconds
    return datetime.strptime(s, '%Y-%m-%dT%H:%M:%S:%f')

#group images with identical fiber pointings that were taken cronologically
def group_obs(direc):
    #Declare fits header IMGTYPEs to ignore
    ignore_types = ['zero']

    #Get list of fits files in direc
    if direc[-1] != '/':
        direc += '/'
    fitsList = glob.glob(direc+'*.fits')+glob.glob(direc+'*.fit')

    #Load fits files
    fitsList = [fits.open(fpath) for fpath in fitsList]

    #Sort files chronologically by DATE-OBS 
    fitsList = sorted(fitsList, key=lambda f: convert_timestr(f[0].header['DATE-OBS']))

    pointings = []
    prev_fiberConfig = ''
    for f in fitsList:
        h = f[0].header
        n = 1
        fiberConfig = ''
        while 'SLFIB'+str(n) in h:
            fiberConfig += h['SLFIB'+str(n)]
            n += 1
        if fiberConfig != prev_fiberConfig:
            pointings.append(obs_group())
        if not h['IMGTYPE'] in ignore_types:
            pointings[-1].add_image(f, h['IMGTYPE'])
            prev_fiberConfig = fiberConfig

    return pointings
