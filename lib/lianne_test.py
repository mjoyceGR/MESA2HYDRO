#!/usr/bin/python
import numpy as np
from pygfunc import to_cdef

print("Let's start")

i_arr = np.array([1,2,3], dtype=np.int32)
d_arr = np.array([4,5,6], dtype=np.double)

print("calling to_cdef")
to_cdef(len(i_arr), i_arr, d_arr, d_arr, d_arr, d_arr, d_arr, d_arr)
print("Nice try")
