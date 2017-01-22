# %python driver.py <command> <dir>
# commands are:
#   make-recipe
#   flat

import os
import sys
from process import process_flat

if __name__ == '__main__':
    command = sys.argv[1]
    direc = sys.argv[2]

    if command = 'flat':
        process_flat(direc)
