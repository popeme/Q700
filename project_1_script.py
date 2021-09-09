#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  9 00:57:20 2021

@author: Maria

Note: needs the smallworld package from Ben Maier to run. Get here: https://github.com/benmaier/smallworld

This also needs my KuraMP2 class (stripped-down version of KuraMP). If you use outside 
of class maybe credit me somehow?
"""

from smallworld import get_smallworld_graph
from smallworld.draw import draw_network
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import kuraMP2


### Network Parameters
N = 10 #number of nodes
beta = 1 #1 is many long-range connections- see smallworld docs
k_over_2 = 2 #minimum starting degree of each node

### Get Starting Network
G = get_smallworld_graph(N, k_over_2, beta)
draw_network(G,k_over_2) #make sure the network is what we want
Kij = nx.to_numpy_array(G) #store the adjacency matrix

### Kuramoto Parameters
k_range = np.arange(0,40,0.2) #track the phase transition over these couplings
MeanFreq = 10 #mean frequency of oscillators
STDfreq = 0.01 #STD frequency of oscillators
Omega = np.random.normal(MeanFreq,STDfreq,size=N) * 2 * np.pi #set intrinsic frequency of oscillators with normal dist.
Theta = np.random.uniform(size=N) * 2 * np.pi #set initial phase of oscillators

### Time Parameters
duration = 500
dt = 0.01
tspan = np.arange(0.0,duration,dt)
transient = 1000 #number of timesteps

### Initialize Measure Storage
Rall = np.zeros((11,len(k_range)))
kcomplete = np.zeros((11))

### Save networks
Kijs = np.zeros((10,10,11))
ns = np.zeros((10)) +12

#For each loop destroy long-range connections for node chosen and replace with only local
#connections
for i in range(0,11):   
    print(i)
    
    if i >0: #at zero we want the initial network
        
        ### Destroy long-range connections
        n = np.random.randint(0,10) #pick a node
        while np.any(n in ns) ==True: #check if we've done that node before, and if so, repick
            n = np.random.randint(0,10) 
        ns[i-1]= n #store node
        Kij[n,:] = 0 #clear all outgoing connections
        
        ### Set Local Connections
        if n==0: #this is a clunky way to do wraparound conditions but here we are
            Kij[0,0:2]=1
            Kij[0,9]=1
        elif n==9:
            Kij[9,8:10]=1
            Kij[9,0] = 1
        else:
    
            Kij[n,n-1:n+2] = 1 #set local connections
    
    else: #this is the zero condition
        Kij = Kij
    
    Kijs[:,:,i] = Kij #store the network
    
    ### Run the phase transition (PT)
    for k in range(0,len(k_range)): #looping over coupling constants
        kura = kuraMP2.Kuramoto(N,k_range[k],Omega,Theta,Kij)
        ths = np.zeros((len(tspan),N)) #initialize phase timeseries
        
        ### Run simulation
        step = 0
        for t in tspan:
            kura.step2(dt)
            ths[step] = kura.Theta
            step += 1
        
        ### Store Order Parameter for this PT
        R = kura.OrderParameter(ths[transient:,:])
        mR = np.mean(R)
        Rall[i,k] = mR
        
        if mR > 0.95: #we can save time by assuming once it goes to synchrony it won't come back out
            kcomplete[i] = k_range[k]
            break