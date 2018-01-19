import unittest
import os

from ensure_path import ensure_path
class test_ensure_path(unittest.TestCase):
    def path_test(self, path):
        #Check if directory already exists.
        if os.path.exists(path):
            os.rmdir(path)

        #Check ensure_path functionality.
        self.assertEqual(os.path.normpath(path) , os.path.normpath(ensure_path(path))) #Does it return the given path?
        self.assertTrue(os.path.exists(path))                                          #Does the path exist afterwards?
        
        #Clean up.
        dnames = filter(None,path.split("/"))
        direcs_to_remove = ["/".join(dnames[:i]) for i in range(len(dnames),0,-1)]
        for direc in direcs_to_remove:
            try:
                os.rmdir(direc)
            except OSError:
                pass #Direc contains things. Don't remove.

    def test_1(self):
        self.path_test("my/test/directory/")
    def test_2(self):
        self.path_test("../../test_dir/")
    def test_3(self):
        self.path_test(".././.././../test_dir/")
    def test_4(self):
        self.path_test("what/../if/../I/../go/../back/../and/../forth")

from output import output_log
from sys import stdout
class test_output_log(unittest.TestCase):
    def setUp(self):
        self.obj = output_log()
        self.log_path = "_test_log_file.txt"
    def tearDown(self):
        if os.path.exists(self.log_path):
            os.remove(self.log_path)

    #Test whether a new log file is created when calling set_log_path.
    def test_set_log_path(self):
        self.obj.set_log_path(self.log_path)
        self.assertEqual(self.log_path, self.obj.log_path)
        self.assertTrue(os.path.exists(self.log_path))

    #Test if object message_str is changed when calling edit_message.
    # Also test that nothing is saved to log when add_to_log=False.
    def test_edit_message_1(self):
        self.obj.set_log_path(self.log_path)
        new_message = "Test Message 1\n"
        self.obj.edit_message(new_message, add_to_log=False)
        self.assertEqual(new_message, self.obj.message_str)
        f = open(self.log_path)
        log_text = f.read()
        f.close()
        self.assertEqual(log_text, "")

    #Test if message is saved to log correctly.
    def test_edit_message_2(self):
        self.obj.set_log_path(self.log_path)
        new_message = "Test Message 2\n"
        self.obj.edit_message(new_message, add_to_log=True)
        f = open(self.log_path)
        log_text = f.read()
        f.close()
        self.assertTrue(new_message in log_text)
        #Later, add tests for date string in log.

from html_plot import plotter
class test_plotter(unittest.TestCase):
    def setUp(self):
        self.p = plotter()
    def filler_plot(self):
        self.p.line([0,1],[0,1])
    def test_set_rootpath_1(self):
        #By default, rootpath should be None
        self.assertEqual(self.p.rootpath, None)
        
        rootpath = "test_direc/"
        self.p.set_rootpath(rootpath)
        self.assertEqual(self.p.rootpath, rootpath)
    def test_set_rootpath_2(self):
        #Test that set_rootpath adds trailing "/"
        rootpath = "test_direc"
        self.p.set_rootpath(rootpath)
        self.assertEqual(self.p.rootpath, rootpath+"/")
    def test_set_title(self):
        title = "new_title"
        self.p.set_title(title)
        self.assertEqual(title, self.p.fig.title.text)
    def test_set_xlabel(self):
        xlabel = "new_xlabel"
        self.p.set_xlabel(xlabel)
        self.assertEqual(xlabel, self.p.fig.xaxis[0].axis_label)
    def test_set_ylabel(self):
        ylabel = "new_ylabel"
        self.p.set_ylabel(ylabel)
        self.assertEqual(ylabel, self.p.fig.yaxis[0].axis_label)
    def test_set_width(self):
        w = 50
        self.p.set_width(w)
        self.assertEqual(self.p.fig.plot_width, w)
    def test_set_height(self):
        h = 50
        self.p.set_height(h)
        self.assertEqual(self.p.fig.plot_height, h)
    #Test basic functionality of save method.
    def test_save_1(self):
        self.filler_plot()
        fname = "_test_save.html"
        self.p.save(fname)
        self.assertTrue(os.path.exists(fname))
        #Clean up
        os.remove(fname)
    #Test that plot is saved in rootpath if only file name is given.
    def test_save_2(self):
        #Set up
        self.filler_plot()
        rootpath = "_test_direc/"
        os.makedirs(rootpath)
        self.p.set_rootpath(rootpath)

        fname = "_test_save.html"
        self.p.save(fname)
        self.assertTrue(os.path.exists(rootpath+fname))

        #Clean up
        os.remove(rootpath+fname)
        os.rmdir(rootpath)
    #Test that plot is not saved in rootpath if alternative file path is given.
    def test_save_3(self):
        #Set up
        self.filler_plot()
        rootpath = "_test_direc/"
        os.makedirs(rootpath)
        self.p.set_rootpath(rootpath)

        fname = "./_test_save.html"
        self.p.save(fname)
        self.assertTrue(os.path.exists(fname))

        #Clean up
        os.remove(fname)
        os.rmdir(rootpath)

    #Test that the bokeh.plotting.figure methods are called correctly.
    def test_bokeh_line(self):
        methods = ["annular_wedge", "annulus", "arc", "asterisk", "bezier",
                   "circle", "circle_cross", "circle_x", "cross", "diamond",
                   "diamond_cross", "ellipse", "hbar", "image", "image_rgba",
                   "image_url", "inverted_triangle", "line", "multi_line",
                   "oval", "patch", "patches", "quad", "quadratic", "ray",
                   "rect", "segment", "square", "square_cross", "square_x",
                   "text", "triangle", "vbar", "wedge", "x"]
        for m in methods:
            exec("self.p."+m)

if __name__ == "__main__":
    unittest.main()
