# ########################################
# Generic class of utilities to code redundancy and lookup
# Should depend only on inbuilt python packages
# 
# Author: Abanti Ranadhir Sahasransu
# Created on: 19 Nov 2025
# Last updated on: 19 Nov 2025

import pathlib
import shutil
import time


# Recreate a directory if it exists else create a new directory
def recreate_dir(dirpath, *, verbose=False):
    pdir = pathlib.Path(f'{dirpath}')
    if pdir.exists() and pdir.is_dir():
        if verbose:
            print(f"Found directory: {dirpath}")
            print(f"Deleting directory")
        shutil.rmtree(pdir)

    pdir.mkdir(parents=True, exist_ok=True)
    if verbose:
        print(f"Created empty directory: {dirpath}")


# Decorator to measure the execution time of a function
def time_eval(func):
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        execution_time = end_time - start_time
        print("\n##############################")
        print(f"Execution time of function {func.__name__}: {execution_time:.6f} seconds")
        return result
    return wrapper


# Create RDataFrame checkpoint
def create_rdf_checkpint(ref_df, filt_df, msg):
    refcount = ref_df.Count().GetValue()
    filtcount = filt_df.Count().GetValue()
    print(f"\n{msg}")
    print(f"Remaining event count = {filtcount}/{refcount} "\
          f"({round(filtcount*100/refcount, 3)}%)")
