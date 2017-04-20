# %python driver.py <command> <dir>
# commands are:
#   make-recipe
#   flat

import sys
import libs.process as process
from libs.ensure_path import ensure_path
from libs.make_recipe import make_recipe
from libs.output import output_log
output = output_log()

if __name__ == '__main__':
    command = sys.argv[1]
    direc = sys.argv[2]
    loc = 'indata/'+direc

    if command == 'make-recipe':
        make_recipe(loc, 'recipes/'+direc+'.recipe')
    elif command == 'full-reduce':
        ensure_path('calib/'+direc+'/')
        process.process_flat(direc, output=output)
        process.process_thar(direc, fast=True, output=output)
        process.process_sky(direc, output=output)
        process.process_target(direc, output=output)
    elif command == 'bias':
        ensure_path('calib/'+direc+'/')
        process.process_bias(direc, output=output)
    elif command == 'flat':
        ensure_path('calib/'+direc+'/')
        process.process_flat(direc, output=output)
    elif command == 'thar':
        ensure_path('calib/'+direc+'/')
        process.process_thar(direc, fast=True, output=output)
    elif command == 'sky':
        ensure_path('calib/'+direc+'/')
        process.process_sky(direc, output=output)
    elif command == 'target':
        ensure_path('calib/'+direc+'/')
        process.process_target(direc, output=output)
