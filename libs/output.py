from sys import stdout
class output_log:
    def __init__(self, writer=stdout, log_path=None):
        self.log_path = log_path
        if log_path != None and not exists(log_path):
            f = open(log_path, 'w')
        self.writer = writer
        self.progress_str = ""
        self.message_str = ""
        self.coverlen = 0
    def set_log_path(log_path):
        self.log_path = log_path
    def update(self):
        strg = self.progress_str+' '+self.message_str
        strg = strg.ljust(self.coverlen)
        self.coverlen = len(strg)
        self.writer.write('\r'+strg)
        self.writer.flush()
    def edit_progress(self, new_str):
        self.progress_str = new_str
        self.update()
    def edit_message(self, new_str, add_to_log=True):
        self.message_str = new_str
        self.update()
        if add_to_log and self.log_path != None:
            dt_str = time.strftime("%Y-%m-%dT%H:%M:%S")
            f = open(self.log_path, 'a')
            f.write('['+dt_str+'] '+self.message_str+'\n')
    def linebreak():
        self.writer.write('\n')
