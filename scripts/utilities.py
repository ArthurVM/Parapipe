""" A set of utility functions and classes.
"""

from os import kill, path
import sys
from signal import alarm, signal, SIGALRM, SIGKILL
from subprocess import PIPE, Popen
import shutil
import re
import time

class Error(Exception):
    """ Base class for exceptions in this module.
    """
    pass

class InputError(Error):
    """ Exception raised for errors in the input.

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, expression, message):
        self.expression = expression
        self.message = message

        sys.stderr.write(f"INPUT ERROR : {self.expression}\n{self.message}")
        sys.exit(1)

class command():
    """ A class for handling shell executed commands
    TODO: clean up this class
    """
    def __init__(self, command):
        self.command = command
        print(f"COMMAND={self.command}")

    def run(self, timeout = -1):
        kill_tree = True

        class Alarm(Exception):
            pass

        def alarm_handler(signum, frame):
            raise Alarm

        p=Popen(self.command, shell = True, stdout = PIPE, stderr = PIPE, env = None)

        if timeout != -1:
            signal(SIGALRM, alarm_handler)
            alarm(timeout)

        try:
            stdout, stderr = p.communicate()
            if timeout != -1:
                alarm(0)

        except Alarm as a:
            pids = [p.pid]
            if kill_tree:
                pids.extend(self.get_process_children(p.pid))
            for pid in pids:
                # This is to avoid OSError: no such process in case process dies before getting to this line
                try:
                    kill(pid, SIGKILL)
                except OSError:
                    pass
            return -1, '', ''

        return p.returncode, stdout, stderr

    def run_comm(self, if_out_return, exit=True):
        """ Run the command with or without an exit code
        """
        returncode, stdout, stderr=self.run(360000)

        if returncode and stderr:
            sys.stderr.write(f"\nERROR: {self.command} FAILED!!! \n\nSTDERR: {stderr}\nSTDOUT: {stdout}\n")
            if exit:
                sys.exit(1)

        if if_out_return:
            return stdout

    def get_process_children(self, pid):
        p = Popen(f"ps --no-headers -o pid --ppid {pid}", shell = True,
	         stdout = PIPE, stderr = PIPE)
        stdout, stderr = p.communicate()
        return [int(p) for p in stdout.split()]

def is_file(filename):
    """ Checks if a path is a file """

    if not path.isfile(filename):
        print("No file found at {f}".format(
                time=time, f=filename), end="\n", file=sys.stderr, flush=True)
        exit(2)
    else:
        return path.abspath(path.realpath(path.expanduser(filename)))

def is_dir(dirname):
    """ Checks if a path is a directory """

    if not path.isdir(dirname):
        print("No directory found at {d}".format(
                time=time, d=dirname), end="\n", file=sys.stderr, flush=True)
        exit(2)
    else:
        return path.abspath(path.realpath(path.expanduser(dirname)))
