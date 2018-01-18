import numpy as np
import bokeh.plotting as bpl

class plotter:
    '''
    This class is meant to emulate a matplotlib.pyplot.figure object, but with
     the bokeh html plotting library.
    '''
    def __init__(self, rootpath=None):
        '''
        ARGUMENTS:
            rootpath - String of a path to a directory in which this html plot will be stored.
                       Optional, by default it is not set and the plot is stored in the current
                       working directory.
        '''
        self.clear()
        self.set_rootpath(rootpath)
    def set_rootpath(self, path):
        '''
        Change the rootpath where the html plot is saved.
        ARGUMENTS:
            path - String of the desired new root_path.
        '''
        if path != None and path[-1] != '/':
            path = path+'/'
        self.rootpath = path
    def clear(self):
        '''
        Reset the plot, erasing whatever was previously plotted.
        '''
        try:
            self.fig = bpl.figure(webgl=True, sizing_mode='stretch_both')
        except:
            self.fig = bpl.figure(sizing_mode='stretch_both')
    def set_title(self, title):
        '''
        Set the plot title
        ARGUMENTS:
            title - String of the title.
        '''
        self.fig.title.text = title
    def set_xlabel(self, xlabel):
        '''
        Set the plot's x-axis label.
        ARGUMENTS:
            xlabel - String of the x-axis label.
        '''
        self.fig.xaxis.axis_label = xlabel
    def set_ylabel(self, ylabel):
        '''
        Set the plot's y-axis label.
        ARGUMENTS:
            ylabel - String of the y-axis label.
        '''
        self.fig.yaxis.axis_label = ylabel
    def set_width(self, w):
        '''
        Set the width of the plot figure.
        ARGUMENTS:
            w - number to use as new width.
        '''
        self.fig.plot_width = w
    def set_height(self, h):
        '''
        Set the height of the plot figure.
        ARGUMENTS:
            h - number to use as new height.
        '''
        self.fig.plot_height = h
    def fill_between(self, x, y1, y2, **kwargs):
        '''
        Plot a shaded region between the curves y1(x) and y2(x).
        ARGUMENTS:
            x - array-like with x values.
            y1 - array-like with y values of the first curve.
            y2 - array-like with y values of the second curve.
            **kwargs - Any parameters accepted by the bokeh.plotting.figure patch 
                       method (fill_color, line_color, etc.). 
        '''
        patch_x = np.append(np.array(x), np.array(x)[::-1])
        patch_y = np.append(np.array(y1), np.array(y2)[::-1])
        self.fig.patch(patch_x, patch_y, **kwargs)
    def __getattr__(self, attr):
        '''
        __getattr__ is set up so that standard bokeh.plotting.figure methods
         can be called naturally without having individual methods for each.
        ex.)
            ...
            my_plotter = plotter()
            my_plotter.line(x_points, y_points)             #Plot a line.
            my_plotter.circle(x_points, y_points, radius=5) #Plot a scatter plot.
            ...
        Common plotting functions can be found at
         https://bokeh.pydata.org/en/latest/docs/reference/plotting.html.
        ''' 
        exec('method = self.fig.'+attr)
        return method
    def show(self):
        '''
        Opens current plot in default web browser.
        '''
        bpl.show(self.fig)
    def save(self, savepath):
        '''
        Save html plot to desired location. If root_path is set and savepath
         contains only a filename, the plot will be saved in root_path.
         Filename should end in .html.
        ARGUMENTS:
            savepath - String of the path where you want the plot saved.
        '''
        if '/' in savepath or self.rootpath == None:
            bpl.output_file(savepath)
        else:
            bpl.output_file(self.rootpath+savepath)
        bpl.save(self.fig)

if __name__ == '__main__':
    #Testing things.
    p = plotter()

    p.set_title('TEST PLOT!!!!')

    x = np.linspace(0, 1, 100)
    y = x**2

    p.line(x, y)


    p.save('testeroni.html')
