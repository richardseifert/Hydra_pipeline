import sys
import libs.process as process
from libs.ensure_path import ensure_path
from libs.output import output_log
from libs.make_recipe import make_recipe
try:
    import argh
    print 'HAD IT'
except ImportError:
    print 'USING EXTERNAL'
    sys.path.append('external/')
    import argh

output = output_log()

def recipe_log(direc):
    data_loc = 'indata/'+direc
    save_path = 'recipes/'+direc+'.recipe'
    make_recipe(data_loc, save_path)

def flat(direc, recipe=None):
    process.process_flat(direc, recipe, output=output)

def thar(direc, recipe=None, fast=False):
    process.process_thar(direc, recipe, output=output, fast=fast)

def sky(direc, recipe=None):
    process.process_sky(direc, recipe, output=output)

def target(direc, recipe=None):
    process.process_target(direc, recipe, output=output)

def full_reduce(direc, recipe=None, fast=False):
    process.process_flat(direc, recipe, output=output)
    process.process_thar(direc, recipe, output=output, fast=fast)
    process.process_sky(direc, recipe, output=output)
    process.process_target(direc, recipe, output=output)

command_list = [recipe_log, flat, thar, sky, target, full_reduce]

parser = argh.ArghParser()
parser.add_commands(command_list)

if __name__ == '__main__':
    argh.dispatch(parser)