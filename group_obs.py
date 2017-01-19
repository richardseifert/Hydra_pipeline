print 'groupObs importing glob'
import glob
print 'groupObs importing datetime'
from datetime import datetime
print 'groupObs importing astropy.io.fits'
from astropy.io import fits
#print 'groupObs importing fitsFile'
#import fitsFile
print 'groupObs finished importing'

class obsGroup:
    def __init__(self, obj=[], flat=[], comp=[], zero=[]):
        self.images = {'object' :list(obj), 'flat':list(flat), 'comp':list(comp), 'zero':list(zero)}
    def addImage(self, im, imgtype):
        try:
            self.images[imgtype].append(im)
        except KeyError:
            print 'Unrecognize d image type:',imgtype

def convert_timestr(s):
    decimal_seconds = '.'+s.split('.')[-1]
    micro_seconds = str(int(float(decimal_seconds)*1.0e6)).rjust(6, '0')
    s = s.split('.')[0]+':'+micro_seconds
    return datetime.strptime(s, '%Y-%m-%dT%H:%M:%S:%f')

#group images with identical fiber pointings that were taken cronologically
def groupObs(direc):
    print 'GETTING FILES'
    #Get list of fits files in direc
    if direc[-1] != '/':
        direc += '/'
    fitsList = glob.glob(direc+'*.fits')+glob.glob(direc+'*.fit')

    print 'LOADING FILES'
    #Load fits files
    fitsList = [fits.open(fpath) for fpath in fitsList]

    print 'SORTING FILES'
    #Sort files chronologically by 
    fitsList = sorted(fitsList, key=lambda f: convert_timestr(f[0].header['DATE-OBS']))

    pointings = []
    pointings.append(obsGroup())
    prev_fiberConfig = ''
    for f in fitsList:
        h = f[0].header
        n = 1
        fiberConfig = ''
        while 'SLFIB'+str(n) in h:
            fiberConfig += h['SLFIB'+str(n)]
            n += 1
        if fiberConfig != prev_fiberConfig:
            pointings.append(obsGroup())
        pointings[-1].addImage(f, h['IMGTYPE'])
        prev_fiberConfig = fiberConfig

    for p in pointings:
        imgs = p.images
        for k in imgs.keys():
            print k, len(imgs[k])
        print

    return pointings
