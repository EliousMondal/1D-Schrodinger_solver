'''Author- MD ELIOUS ALI MONDAL
   Created - 03/08/2020'''
#solving the schroedinger equation for finite well
import time
import numpy as np 
import matplotlib.pyplot as plt 

start = time.perf_counter()

#step size for x
cdef double h = 0.01

#no. of x points
cdef int n = int(6/h)+1

#defining the potential
cpdef int V(double x):
    '''returns V as a function of x'''
    if x < -1 or x > 1:
        return 10
    elif x >= -1 and x <= 1:
        return 0

#defining the function to carry out Numerov method
cpdef double Numerov(double x,double y1,double y2,double E):
    '''returns y[i+1] value'''
    cdef double u = 1 - (1/6.)*(h**2)*(V(x)-E)
    return ((12-10*u)*y1-u*y2)/u

x = np.linspace(-3,3,n)
y = np.zeros(n)
z = 0.000001
y[1] = z*h 
 
def Psi(float E):
    cdef int i
    for i in range(2,n):
        y[i] = Numerov(x[i],y[i-1],y[i-2],E)
    cdef double N_const = 0                 #Inverse of square of Normalisation constant
    cdef double j
    for j in y:
        N_const = N_const + j*j*h
    cdef double A = 1.0/np.sqrt(N_const)
    m = A*y
    return(m[-1],m)

Eigenvalues = []
a = np.linspace(0.1,9.9,100)
P = np.array([Psi(i)[0] for i in a])
cdef int pi 
for pi in range(len(P)-1):
    if (P[pi]<0 and P[pi+1]>0) or (P[pi]>0 and P[pi+1]<0):
        low = a[pi]
        high = a[pi+1]
        mid = (low + high)/2.
        while abs(Psi(mid)[0]) > h**2:
            mid = (low + high)/2.
            if P[pi] < 0:
                if Psi(mid)[0] < 0:
                    low = mid
                else:
                    high = mid
            elif P[pi] > 0:
                if Psi(mid)[0] > 0:
                    low = mid
                else:
                    high = mid
        Eigenvalues.append(mid)

end = time.perf_counter()
print(Eigenvalues)
print("Time taken = ",end-start," seconds")

#Plotting the obtained wavefunction
plt.figure(1)
plt.title('Numerical solution of finite well')

for i in Eigenvalues:
    q = Psi(i) 
    plt.plot(x,q[1]+1)

V = np.array([V(i)/5. for i in x])
plt.plot(x,V,'black')
plt.show()
