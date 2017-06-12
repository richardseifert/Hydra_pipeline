import os

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
    for dirname in dirs_to_make:
    	try:
        	os.makedirs(dirname)
        except OSError:
        	pass #Directory already exists.
    return path

