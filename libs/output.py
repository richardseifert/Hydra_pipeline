from sys import stdout
import os
import time

class output_log:
    '''
    The purpose of this object it to log and display information in real time
     during program execution, saving lines of text to a file and also printing 
     updated info to the screen.
     info to a file and printing it to the screen.
    '''
    def __init__(self, writer=stdout, log_path=None):
        '''
        ARGUMENTS:
            writer - A writer object to which information is printed throughout
                     execution. Optional, default is stdout.
            log_path - String path to a file in which log information is to be
                       saved. If none is provided, log info is not saved to a file.
        '''
        self.log_path = log_path
        if log_path != None and not os.path.exists(log_path):
            f = open(log_path, 'w')
        self.writer = writer
        self.progress_str = ""
        self.message_str = ""
        self.log_str = ""
        self.coverlen = 0
    def set_log_path(self, log_path):
        '''
        Set the log_path variable and create file at new log_path if it doesn't 
         already exist.
        ARGUMENTS:
            log_path - String path to the desired new log file.
        '''
        self.log_path = log_path
        f = open(self.log_path, 'a')
        f.write(self.log_str)
        f.close()
    def update(self):
        '''
        Update what is currently printed to self.writer.
        '''
        strg = self.progress_str+' '+self.message_str
        strg = strg.ljust(self.coverlen)
        self.coverlen = len(strg)
        self.writer.write('\r'+strg)
        self.writer.flush()
    def edit_progress(self, new_str):
        '''
        Edit the progress section and update the line currently printed.
        ARGUMENTS:
            new_str - String with the new progress info.
        '''
        self.progress_str = new_str
        self.update()
    def edit_message(self, new_str, add_to_log=True):
        '''
        Edit the message section and update the line currently printed.
        ARGUMENTS:
            new_str - String with the new progress info.
            add_to_log - Boolean to determine whether or not the new message is 
                         stored in the log file. Optional, default True.
        '''
        self.message_str = new_str
        self.update()
        if add_to_log:
            self.log(new_str)
    def log(self, new_str):
        '''
        Log the given string to the log file located at log_path, along with
         a timestamp.
        '''
        dt_str = time.strftime("%Y-%m-%d %H:%M:%S")
        log_newline = '['+dt_str+'] '+new_str+'\n'
        self.log_str += log_newline
        if self.log_path != None:
            f = open(self.log_path, 'a')
            f.write(log_newline)
            f.close()
    def linebreak(self):
        self.writer.write('\n')
