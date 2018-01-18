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
        self.assertEqual(log_path, self.obj.log_path)
        self.assertTrue(os.path.exists(log_path))

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

if __name__ == "__main__":
    unittest.main()
