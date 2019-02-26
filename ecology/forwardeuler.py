from numpy import linspace, zeros
import matplotlib.pyplot as plt

"""
Forward Euler scheme (Euler's method) is a population growth model that posits an animal species poulation (N)
are in a spatial region such that there is no exchange with other areas. During a time interval (delta_t), some animals
will die and some new will be born. The number of deaths and births are expected to be proportional to N.

The net growth will be:

N(t + delta_t) - N(t) = bN(t) - dN(t)

in which bN(t) is the number of births and dN(t) is the number of deaths.

We can also model the population using a difference equation:

N(t + delta_t) = N(t) + delta_T*r*N(t)  in which r = b-d, or the growth rate.
"""
N0 = 100 # Start with this many animals
r   = .1 # growth rate
dt  = .5 # step size
N_t = 20 # number of steps

t = linspace(0, (N_t+1)*dt, N_t+2)
N = zeros(N_t+2)

N[0] = N0
for n in range(N_t+1):
    N[n+1] = N[n] + r*dt*N[n]

numerical_sol = 'bo' if N_t < 70 else 'b-'
plt.plot(t, N, numerical_sol, t, N_0*exp(r*t), 'r-')
plt.legend(['numerical', 'exact'], loc='upper left')
plt.xlabel('t'); plt.ylabel('N(t)')
filestem = 'growth1_%dsteps' % N_t
plt.savefig('%s.png' % filestem); plt.savefig('%s.pdf' % filestem)
