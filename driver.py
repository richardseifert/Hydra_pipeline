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
    print tmp_path, dirs_to_make[:-1]
    for dirname in dirs_to_make[:-1]:
        print dirname
        os.makedirs(dirname)
    return path

if __name__ == '__main__':
    command = sys.argv[1]
    direc = sys.argv[2]
    loc = 'indata/'+direc
    
    if command == 'make-recipe':
        make_recipe(loc, 'recipes/'+direc+'.recipe')

    if command == 'flat':
        print 'call process_flat'
        process.process_flat(direc)
