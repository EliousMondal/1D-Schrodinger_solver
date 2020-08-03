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

x = np.linspace(-3,3,n)
y = np.zeros(n)
z = 0.000001
y[1] = z*h 
 
def Psi(float E):
    cdef int i
    for i in range(2,n):
        y[i] = Numerov(x[i],y[i-1],y[i-2],E)
    cdef float N_const = 0                 #Inverse of square of Normalisation constant
    cdef float j
    for j in y:
        N_const = N_const + j*j*h
    cdef float A = 1.0/np.sqrt(N_const)
    m = A*y
    return(m[-1],m)

end = time.perf_counter()

print("Time taken = ",end-start," seconds")
