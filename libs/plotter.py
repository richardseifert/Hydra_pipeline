import numpy as np
import matplotlib.pyplot as plt
plt.ioff()

class plotter:
    def __init__(self, rootpath=None, **kwargs):
        self.fig = None
        self.clear_plot(**kwargs)
        self.set_rootpath(rootpath)
    def clear_plot(self, **kwargs):
        try:
            plt.close(self.fig)
        except TypeError:
            pass #self.fig was probably None
        self.fig, self.ax = plt.subplots(**kwargs)
    def set_rootpath(self, rootpath):
        self.rootpath = rootpath
        if rootpath != None and rootpath[-1] != '/':
            self.rootpath += '/'

    def __getattr__(self, name):
        try:
            exec('ax_method = self.ax.'+name)
        except AttributeError:
            raise AttributeError("'plotter' object has no attribute '"+name+"'")
        return ax_method

    def save(self, savepath):
        if '/' in savepath:
            fullsavepath = savepath
        else:
            fullsavepath = self.rootpath+savepath
        self.fig.savefig(fullsavepath)

class null_object:
    def __getattr__(self, name):
        def method(*args, **kwargs):
            print'NULL_METHOD'
            pass
        return method
