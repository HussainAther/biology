
"""
Forward Euler scheme (Euler's method) is a population growth model that posits an animal species poulation (N)
are in a spatial region such that there is no exchange with other areas. During a time interval (delta_t), some animals
will die and some new will be born. The number of deaths and births are expected to be proportional to N.

The net growth will be:

N(t + delta_t) - N(t) = bN(t) - dN(t)

in which bN(t) is the number of births and dN(t) is the number of deaths.

"""
