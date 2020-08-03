'''Author- MD ELIOUS ALI MONDAL
   Created - 03/08/2020'''
#solving the schroedinger equation for finite well
import time
import numpy as np 
import matplotlib.pyplot as plt 

start = time.perf_counter()

#step size for x
cdef float h = 0.01

#no. of x points
cdef int n = int(6/h)+1

#defining the potential
cpdef int V(float x):
    '''returns V as a function of x'''
    if x < -1 or x > 1:
        return 10
    elif x >= -1 and x <= 1:
        return 0

#defining the function to carry out Numerov method
cpdef float Numerov(float x,float y1,float y2,float E):
    '''returns y[i+1] value'''
    cdef float u = 1 - (1/6.)*(h**2)*(V(x)-E)
    return ((12-10*u)*y1-u*y2)/u

    

end = time.perf_counter()

print("Time taken = ",end-start," seconds")
