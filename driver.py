# %python driver.py <command> <dir>
# commands are:
#   make-recipe
#   flat

import sys
import libs.process as process
from libs.ensure_path import ensure_path
from libs.make_recipe import make_recipe


if __name__ == '__main__':
    command = sys.argv[1]
    direc = sys.argv[2]
    loc = 'indata/'+direc
    
    if command == 'make-recipe':
        make_recipe(loc, 'recipes/'+direc+'.recipe')
    elif command == 'full-reduce':
        ensure_path('calib/'+direc+'/')
        process.process_flat(direc)
        process.process_thar(direc, fast=True)
        process.process_sky(direc)
        process.process_target(direc)
    elif command == 'bias':
        ensure_path('calib/'+direc+'/')
        process.process_bias(direc)
    elif command == 'flat':
        ensure_path('calib/'+direc+'/')
        process.process_flat(direc)
    elif command == 'thar':
        ensure_path('calib/'+direc+'/')
        process.process_thar(direc, fast=True)
    elif command == 'sky':
        ensure_path('calib/'+direc+'/')
        process.process_sky(direc)
    elif command == 'target':
        ensure_path('calib/'+direc+'/')
        process.process_target(direc)
