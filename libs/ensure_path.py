import os

def ensure_path(path):
    '''
    This function checks if a path to a directory exists. If the path exists, this 
    function does nothing but return the given path. If the path does not exist, 
    then it makes all of the necessary directories for the path to exist.

    ARGUMENTS:
        path - String of a path to a directory.
    RETURNS:
        path - String, the same path given. Now, the directory at the given path
               exists, if it didn't already.
    '''

    #Do nothing if the path already exists.
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
        raise ValueError("The path "+path+" cannot be created.")
    for dirname in dirs_to_make:
    	try:
        	os.makedirs(dirname)
        except OSError:
        	pass #Directory already exists.
    return path

