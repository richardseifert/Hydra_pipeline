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

if __name__ == "__main__":
    unittest.main()
