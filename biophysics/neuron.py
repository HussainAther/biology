import numpy as np

"""
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

a = vr_eff - v_ss / sigma_eff # in which Vr is the effective threshold for and v_ss is the steady-state voltage that results 

b = vth_eff - v_ss / sigma_eff i # in which Vth_eff is the effective threshold potential

v_post  = [tau_ref + tau_m_eff * sqrt(pi) * integral from 

erf(x) is 2/sqrt(pi) * integral from 0 to x of exp(u^2)du, also known as the error function that we use in integrating
the normal distribution.
"""
