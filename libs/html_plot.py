import numpy as np
import bokeh.plotting as bpl

class plotter:
    def __init__(self, rootpath=None):
        try:
            self.fig = bpl.figure(webgl=True, sizing_mode='stretch_both')
        except:
            self.fig = bpl.figure(sizing_mode='stretch_both')
        self.set_rootpath(rootpath)
    def set_rootpath(self, path):
        if path != None and path[-1] != '/':
            path = path+'/'
        self.rootpath = path
    def clear(self):
        try:
            self.fig = bpl.figure(webgl=True, sizing_mode='stretch_both')
        except:
            self.fig = bpl.figure(sizing_mode='stretch_both')
    def set_title(self, title):
        self.fig.title.text = title
    def set_xlabel(self, xlabel):
        self.fig.xaxis.axis_label = xlabel
    def set_ylabel(self, ylabel):
        self.fig.yaxis.axis_label = ylabel
    def set_width(self, w):
        self.fig.plot_width = w
    def set_height(self, h):
        self.fig.plot_height = h
    def fill_between(self, x, y1, y2, **kwargs):
        patch_x = np.append(np.array(x), np.array(x)[::-1])
        patch_y = np.append(np.array(y1), np.array(y2)[::-1])
        self.fig.patch(patch_x, patch_y, **kwargs)
    def __getattr__(self, attr):
        exec('method = self.fig.'+attr)
        return method
    def show(self):
        bpl.show(self.fig)
    def save(self, savepath):
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
