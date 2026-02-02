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

import numpy as np
from matplotlib.transforms import Transform
from matplotlib.scale import ScaleBase
from matplotlib.ticker import FixedLocator
import matplotlib as mpl


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


# ---- Define the piecewise linear transform ----
# Map x in [0.1,0.9] -> [0,1]; [0.9,0.99] -> [1,2]; [0.99,1.0] -> [2,3]; [1.0,1.1] -> [3,3.1]
# Inverse is needed so Matplotlib can place ticks properly.

def x_to_u(x):
    x = np.asarray(x)
    u = np.empty_like(x, dtype=float)

    # Segment 1
    m1 = (0 <= x) & (x < 0.9)
    u[m1] = (x[m1] - 0.1) / (0.9 - 0.1) * 1.0  # [0,1)

    # Segment 2
    m2 = (x >= 0.9) & (x < 0.99)
    u[m2] = 1.0 + (x[m2] - 0.9) / (0.99 - 0.9) * 1.0  # [1,2)

    # Segment 3
    m3 = (x >= 0.99) & (x <= 1.0)
    u[m3] = 2.0 + (x[m3] - 0.99) / (1.0 - 0.99) * 1.0  # [2,3]

    # Segment 4
    m4 = (x > 1.0) & (x <= 1.1)
    u[m4] = 3.0 + (x[m4] - 1.0) / (1.1 - 1.0) * 0.1  # [3,3.1]

    return u


def u_to_x(u):
    u = np.asarray(u)
    x = np.empty_like(u, dtype=float)

    # Inverse mappings for each unit-length segment
    n1 = (u >= 0) & (u < 1)
    x[n1] = 0.1 + (u[n1] - 0.0) * (0.9 - 0.1) / 1.0

    n2 = (u >= 1) & (u < 2)
    x[n2] = 0.9 + (u[n2] - 1.0) * (0.99 - 0.9) / 1.0

    n3 = (u >= 2) & (u <= 3)
    x[n3] = 0.99 + (u[n3] - 2.0) * (1.0 - 0.99) / 1.0

    n3 = (u >= 2) & (u <= 3)
    x[n3] = 0.99 + (u[n3] - 2.0) * (1.0 - 0.99) / 1.0

    n4 = (u > 3) & (u <= 3.1)
    x[n4] = 1.0 + (u[n4] - 3.0) * (1.1 - 1.0) / 0.1

    return x


class PiecewiseTransform(Transform):
    input_dims = output_dims = 1
    is_separable = True

    def transform_non_affine(self, a):
        return x_to_u(a)

    def inverted(self):
        return InvertedPiecewiseTransform()


class InvertedPiecewiseTransform(Transform):
    input_dims = output_dims = 1
    is_separable = True

    def transform_non_affine(self, a):
        return u_to_x(a)

    def inverted(self):
        return PiecewiseTransform()


class PiecewiseScale(ScaleBase):
    name = 'piecewise_0p1_0p9_0p999'

    def get_transform(self):
        return PiecewiseTransform()

    def set_default_locators_and_formatters(self, axis):
        # Put ticks at the decade boundaries and a few intermediate points
        major_ticks = [0.1, 0.2, 0.5, 0.9, 0.95, 0.99, 0.995, 1.0]
        axis.set_major_locator(FixedLocator(major_ticks))
        axis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3f'))

    def limit_range_for_scale(self, vmin, vmax, minpos):
        # Clamp to valid domain to avoid weird extrapolations
        vmin = max(0.1, vmin)
        vmax = min(1.1, vmax)
        return vmin, vmax


# Register the scale so we can use set_xscale('piecewise_0p1_0p9_0p999')
mpl.scale.register_scale(PiecewiseScale)
# ---- End of piecewise linear transform ----