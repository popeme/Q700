import numpy as np
from scipy import signal

class Kuramoto:

    def __init__(self,N,K,Omega,Theta,Kij):
        self.N = N
        self.K = K
        self.Omega = Omega
        self.Theta = Theta.copy()
        self.Kij = Kij
      

    def step2(self,dt):
        theta_i = self.Theta * np.ones((self.N,self.N))
        theta_j = np.transpose(self.Theta*np.ones((self.N,self.N)))
        theta_ij = theta_i-theta_j
        dtheta = self.Omega + ((self.K/self.N)*np.sum(self.Kij*np.sin(theta_ij),axis=1))
        self.Theta += dt*dtheta
        
    def OrderParameter(self,ths):
        Rtemp = np.mean(np.exp(1j*ths),1)
        R = np.abs(Rtemp)
        return R
