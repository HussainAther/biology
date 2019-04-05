import numpy as np
import scipy as sp

"""
In this biological neuron model (spiking neuron model)
We can calculate the firing rate of a neuron by taking into account the biophysical factors involved.
We use the realisatic firing statistics of AMPA, GABA-A, and NMDA synaptic inputs for a leaky integrate-and-fire (LIF)
neuron that has excitatory synaptic inputs and inhibitory synaptic inputs. We also take into account the voltage-dependence
of NMDA channels and voltage-dependence of post-synaptic currents. We incorporate finite time constant of post-synaptic currents
into the delta-functions that describe the post-synaptic currents without duration or temporal kinetics. We use solutions
to the Fokker-Planck (Fokker fokker Planck planck) equation, the partial differential equation to describe how a
probability density function of the velocity of a particle under the influence of drag forces and random Brownian forces
evolves over time. It's also referred ot as the Kolmogorov forward equation or the Smoluchowski equation when applied
to particle position distributions.

Firing rate of a neuron:

a = vr_eff - v_ss / sigma_eff # in which Vr is the effective threshold for resting potential and v_ss is the steady-state voltage that results 

b = vth_eff - v_ss / sigma_eff i # in which Vth_eff is the effective threshold potential

v_post  = [tau_ref + tau_m_eff * sqrt(pi) * integral from 

erf(x) is 2/sqrt(pi) * integral from 0 to x of exp(u^2)du, also known as the error function that we use in integrating
the normal distribution.
"""
a = (vr_eff - v_ss) / sigma_eff
b = (vth_eff - v_ss) / sigma_eff

def erf(x):
    """
    Return the value of the error function for some value x.
    """
    return 2/(np.sqrt(np.pi)) * sp.integrate.quad(np.exp(u**2), 0 x)

def n_int(x):
    """
    Return the integral for the neuron equation we are describing.
    """
    integrand = np.exp(x**2) * (1 + erf(x))
    return sp.integrate.quad(integrand, a, b)

def firingrate(x):
    """
    Return the firing rate of a neuron. taum is the membrane time constant and tauref is 
    tau time constant of the refractive period.
    """ 
    return (tauref + taum * np.sqrt(np.pi) * n_int(x)) ** (-1) 
