import numpy as np

"""
We can calculate the firing rate of a neuron by taking into account the biophysical factors involved.
We use the realisatic firing statistics of AMPA, GABA-A, and NMDA synaptic inputs for a leaky integrate-and-fire (LIF)
neuron that has excitatory synaptic inputs and inhibitory synaptic inputs. We also take into account the voltage-dependence
of NMDA channels and voltage-dependence of post-synaptic currents. We incorporate finite time constant of post-synaptic currents
into the delta-functions that describe the post-synaptic currents without duration or temporal kinetics.
"""
