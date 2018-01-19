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
    ex.)
        my_path = ensure_path("I/need/this/path")
        # Now the directory I/need/this/path exists.
        assert my_path == "I/need/this/path"
        assert os.path.exists(my_path)
    '''

    #Do nothing if the path already exists.
    if os.path.exists(path):
        return path

    #Precondition path
    if path[-1] != "/":
        path += "/"

    #Go through each subdirectory and create it if it doesn't exist.
    direcs = [path[:i+1] for i,c in enumerate(path) if c=="/"]
    for direc in direcs:
        if not os.path.exists(direc):
            os.makedirs(direc)
    return path
