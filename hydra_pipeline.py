import sys
import libs.process as process
from libs.ensure_path import ensure_path
from libs.output import output_log
from libs.html_plot import plotter
from libs.make_recipe import make_recipe
try:
    import argh
except ImportError:
    sys.path.append('external/')
    import argh

output = output_log()
p = plotter()

def recipe_log(direc):
    data_loc = 'indata/'+direc
    save_path = 'recipes/'+direc+'.recipe'
    make_recipe(data_loc, save_path)

def flat(direc, recipe=None):
    output.edit_message('')
    output.edit_message("*** COMMAND ISSUED: 'python hydra_pipeline.py flat "+direc+"' ***")
    process.process_flat(direc, recipe, output_log=output, plotter=p).dispatch()

def thar(direc, recipe=None, fast=False):
    output.edit_message('')
    output.edit_message("*** COMMAND ISSUED: 'python hydra_pipeline.py thar "+direc+"' ***")
    process.process_thar(direc, recipe, output_log=output, fast=fast, plotter=p).dispatch()

def skyflat(direc, recipe=None):
    output.edit_message('')
    output.edit_message("*** COMMAND ISSUED: 'python hydra_pipeline.py skyflat "+direc+"' ***")
    process.process_skyflat(direc, recipe, output_log=output, plotter=p).dispatch()

def sky(direc, recipe=None):
    output.edit_message('')
    output.edit_message("*** COMMAND ISSUED: 'python hydra_pipeline.py sky "+direc+"' ***")
    process.process_sky(direc, recipe, output_log=output, plotter=p).dispatch()

def target(direc, recipe=None):
    output.edit_message('')
    output.edit_message("*** COMMAND ISSUED: 'python hydra_pipeline.py target "+direc+"' ***")
    process.process_target(direc, recipe, output_log=output, plotter=p).dispatch()

def full_reduce(direc, recipe=None, fast=False):
    output.edit_message('')
    output.edit_message("*** COMMAND ISSUED: 'python hydra_pipeline.py full-reduce "+direc+"' ***")
    process.process_flat(direc, recipe, output_log=output, plotter=p).dispatch()
    process.process_thar(direc, recipe, output_log=output, fast=fast, plotter=p).dispatch()
    process.process_skyflat(direc, recipe, output_log=output, plotter=p).dispatch()
    process.process_sky(direc, recipe, output_log=output, plotter=p).dispatch()
    process.process_target(direc, recipe, output_log=output, plotter=p).dispatch()

command_list = [recipe_log, flat, thar, skyflat, sky, target, full_reduce]

parser = argh.ArghParser()
parser.add_commands(command_list)

if __name__ == '__main__':
    argh.dispatch(parser)
