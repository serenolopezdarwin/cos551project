"""
This script acts as a custom module to define helpful class types and functions for our scripts.
Specifics of those class types and functions found in the corresponding docstrings.
"""
# Imports
from datetime import datetime
from time import perf_counter
import os
import pickle as pkl
import sys


class LogStartTime:
    """Decorator object that counts how many times a function has been called."""
    def __init__(self, func):
        self.func = func
        self.first_call = True
        self.initialized = ""

    def __call__(self, *args, **kwargs):
        # The first time the function is called, it sets the initialization time and denotes the fact that the function
        # has been called once before.
        if self.first_call:
            self.first_call = False
            self.initialized = datetime.now()
        return self.func(initialized=self.initialized, *args, **kwargs)


@LogStartTime
def log(text: str, initialized=datetime.now()) -> None:
    """Given a string, writes that string into a logfile named after the calling script (NOT this module).
    Intentionally mutable default argument included, because we actually want to track function calls across scripts."""
    head_script = sys.argv[0]
    if not os.path.exists("logs"):
        os.mkdir("logs")
    logfile = f"logs/{head_script.replace('.py', '_')}{initialized.strftime('%m%d%y_%H%M')}.log"
    if log.first_call:
        log_out = open(logfile, 'w')
    else:
        log_out = open(logfile, 'a')
    log_out.write(text + '\n')
    log_out.close()


def default_arg_decorator(func):
    """Nested decorator that enables default arguments within other decorators."""
    def wrapped_decorator(*args):
        if len(args) == 1 and callable(args[0]):
            return func(args[0])
        else:
            def real_decorator(decorated):
                return func(decorated, *args)
            return real_decorator
    return wrapped_decorator


@default_arg_decorator
def log_time(func, task=""):
    """Decorator function that outputs how long each function takes to complete. Takes argument of task (just
    a string summary of what the function is doing and logs the time it takes to complete the function. Task defaults
    to the name of the decorated function."""
    if not task:
        task = func.__name__

    def wrapper(*args, **kwargs):
        start = perf_counter()
        val = func(*args, **kwargs)
        log(f"{task} in {str(perf_counter() - start)}")
        return val
    return wrapper


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


class Bidict(dict):
    """A dictionary that can be reverse-hash searched with self.inverse[value] (returns key). Used to 'link' two lists
    to look up translations between them. Reflects setitem and delitem to the inverse version of the dictionary."""
    def __init__(self, *args, **kwargs):
        super(Bidict, self).__init__(*args, **kwargs)
        self.inverse = {}
        for key, value in self.items():
            self.inverse.setdefault(value, []).append(key)

    def __setitem__(self, key, value):
        if key in self:
            self.inverse[self[key]].remove(key)
        super(Bidict, self).__setitem__(key, value)
        self.inverse.setdefault(value, []).append(key)

    def __delitem__(self, key):
        self.inverse.setdefault(self[key], []).remove(key)
        if self[key] in self.inverse and not self.inverse[self[key]]:
            del self.inverse[self[key]]
        super(Bidict, self).__delitem__(key)


class VersionError(Exception):
    """Raised when we encounter an error based on wrong software versions."""
    pass


def invert_hash(in_dict: dict, array_vals=False, identifier=0) -> Bidict:
    """Converts a dictionary to a bidict. If the dictionary's values are arrays, will use the identifier argument to
    use that field of the dictionary as the values in the bidict."""
    out_dict = Bidict({})
    for key, value in in_dict.items():
        if array_vals:
            out_dict[key] = value[identifier]
        else:
            out_dict[key] = value
    return out_dict


def aggregate_expression_level(by: str, exp_mat: list, sorted_exp_mat_path="") -> [list, list]:
    """Given a dimension and an expression matrix, aggregates the expression level along each dimension, like doing the
    row or column sums of a non-sparse matrix. Returns the sums indexed by dimension ids (not indices, since dimensions
    are 1-indexed while indices are 0-indexed due to python's basic rules). Also returns the sorted zip object so that
    we avoid re-sorting in other functions, as it's computationally expensive for large lists."""
    start = perf_counter()
    if by == "genes":
        slice_idx = 0
        other_idx = 1
    elif by == "cells":
        slice_idx = 1
        other_idx = 0
    else:
        raise KeyError("Incorrect 'by' argument passed to aggregation function")
    # Avoids sorting bottleneck, and lets us verify our chunking method more easily.
    if sorted_exp_mat_path and os.path.exists(sorted_exp_mat_path):
        with open(sorted_exp_mat_path, 'rb') as sorted_exp_mat_in:
            exp_mat_sort = pkl.load(sorted_exp_mat_in)
    else:
        id_list = exp_mat[slice_idx]
        unprocessed_list = exp_mat[other_idx]
        exp_list = exp_mat[2]
        exp_mat_sort = sorted(zip(id_list, exp_list, unprocessed_list))
        # Only stores our sorted matrices if we designate a storage path.
        if sorted_exp_mat_path:
            with open(sorted_exp_mat_path, 'wb') as sorted_exp_mat_out:
                pkl.dump(exp_mat_sort, sorted_exp_mat_out)
    # The exp_mat_sort call finds the max of id_list without having to iterate over the list again. Adds one because
    # Anything we initialize with it is 1-indexed.
    max_idx = exp_mat_sort[-1][0] + 1
    # We also add an extra index, since the first index (0) will be empty due to cells and genes both being 1-indexed.
    # This will simplify later references to the generated data object downstream.
    aggregate_xp_list = [0] * max_idx
    # Maximum index we can reference in the matrix
    exp_mat_max = len(exp_mat_sort) - 1
    chunked_exp_mat = [[] for _ in range(max_idx)]
    aggregate_xp, aggregate_data = 0, []
    for idx, (id_num, exp, other_id) in enumerate(exp_mat_sort):
        aggregate_xp += exp
        aggregate_data.append((other_id, exp))
        # Checks if the next entry in the list is from a different id or doesn't exist.
        if idx == exp_mat_max or exp_mat_sort[idx+1][0] != id_num:
            # That being the case, dumps our current aggregated expression level and resets it.
            aggregate_xp_list[id_num] = aggregate_xp
            # Also dumps the list of other coord, exp tuples to the index of this id to be calculated later.
            chunked_exp_mat[id_num] = aggregate_data
            aggregate_xp, aggregate_data = 0, []
    log(f"Exp matrix chunked by {by} in {str(perf_counter() - start)}")
    return aggregate_xp_list, chunked_exp_mat


def main():
    print("You seem to have accidentally run the cos551 module.\n"
          "If you see this message, something has gone wrong.")


if __name__ == '__main__':
    main()
