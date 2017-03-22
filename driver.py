# %python driver.py <command> <dir>
# commands are:
#   make-recipe
#   flat

import os
import sys
import process
from make_recipe import make_recipe


#This function checks if a path to a directory exists. If the path does not
#   exist, then it makes all of the necessary directories for the path to exist.
#This function is used by order_spec objects when saving.
def ensure_path(path):
    if os.path.exists(path):
        return path
    slash_indxs = [indx for indx,char in enumerate(path) if char=="/"]
    dirs_to_make = []
    for indx in reversed(slash_indxs):
        rmv = path[indx:].split("/")[1]
        tmp_path = path[0:indx]+"/"
        dirs_to_make.insert(0,tmp_path+rmv)
        if os.path.exists(tmp_path):
            found = True
            break
    if not found:
        print "Path does not exist."
        return path
    for dirname in dirs_to_make[:-1]:
        os.makedirs(dirname)
    return path

if __name__ == '__main__':
    command = sys.argv[1]
    direc = sys.argv[2]
    loc = 'indata/'+direc
    
    if command == 'make-recipe':
        make_recipe(loc, 'recipes/'+direc+'.recipe')
    elif command == 'full-reduce':
        process.process_bias(direc)
        process.process_flat(direc)
        process.process_thar(direc)
        process.process_sky(direc)
    elif command == 'bias':
        process.process_bias(direc)
    elif command == 'flat':
        process.process_flat(direc)
    elif command == 'thar':
        process.process_thar(direc)
    elif command == 'sky':
        process.process_sky(direc)
    elif command == 'target':
        process.process_target(direc)
