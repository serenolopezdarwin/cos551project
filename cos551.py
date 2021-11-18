"""
This script acts as a custom module to define helpful class types and functions for our scripts.
Specifics of those class types and functions found in the corresponding docstrings.
"""
# Imports
from datetime import datetime
import os
import sys


# noinspection PyDefaultArgument
def log(text: str, calls=[0]):
    """Given a string, writes that string into a logfile named after the calling script (NOT this module).
    Intentionally mutable default argument included, because we actually want to track function calls across scripts."""
    head_script = sys.argv[0]
    if not os.path.exists("logs"):
        os.mkdir("logs")
    now = datetime.now()
    logfile = f"logs/{head_script.replace('.py', '_')}{now.strftime('%m%d%y_%H%M')}.log"
    if calls[0]:
        log_out = open(logfile, 'a')
    else:
        log_out = open(logfile, 'w')
    log_out.write(text + '\n')
    log_out.close()
    # Only evaluated once, returns false on first function call to reset logfile.
    calls[0] += 1
    return calls[0]


class Mean:
    """Creates an object class "Mean" that stores both average and count data for a set of numbers, that can be easily
    appended with additional values during a loop. Avoids generating additional variables in our scripts"""
    def __init__(self):
        self.count = 0
        self.total = 0
        self.avg = 0

    # Calculates the new mean based on an added entry. Also appends count and total information.
    def add(self, entry):
        if isinstance(entry, int) or isinstance(entry, float):
            entry_sum = entry
            entry_total = 1
        elif isinstance(entry, list):
            entry_sum = sum(entry)
            entry_total = len(entry)
        else:
            raise TypeError("Non-numeric entry passed to Mean.add()")
        self.total += entry_sum
        self.count += entry_total
        self.avg = (self.total / self.count)


def main():
    print("You seem to have accidentally run the kangsingh module.\n"
          "If you see this message, something has gone wrong.")


if __name__ == '__main__':
    main()
