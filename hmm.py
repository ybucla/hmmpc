# -*- coding: utf-8 -*-
"""
Created on Mon Mar 06 19:50:58 2017

@author: ybwang
"""
import numpy as np

class HMM(object):
    def __init__(self,obs):
        self.obs = obs
        pass

    def initial(self, A, B, Pi):
        self.A = A
        self.B = B
        self.Pi = Pi
        self.calforward(A, B, Pi)
        self.calbackward(A, B, Pi)

    def calforward(self, A, B, Pi):
        '''A: j x j, B: Emission class, Pi: 1 x j'''
        nstate, nt = A.shape[0], len(self.obs)
        self.forward = np.zeros((nstate, nt))
        self.forward[:,0] = Pi * B.getprob(self.obs[0])
        for t in range(1,nt):
            b = B.getprob(self.obs[t])
            self.forward[:,t] = self.forward[:,t-1].reshape((1,nstate)).dot(A) * b

    def calbackward(self, A, B, Pi):
        '''A: j x j, B: Emission class, Pi: 1 x j'''
        nstate, nt = A.shape[0], len(self.obs)
        self.backward = np.zeros((nstate, nt))
        self.backward[:,nt - 1] = 1.0
        for t in range(nt - 2, -1,  -1):
            b = B.getprob(self.obs[t+1])
            self.backward[:,t] = np.inner((A * b), self.backward[:,t + 1])

    def gam(self, t, i):
        '''given lambda and obs, return probability at t and state is j'''
        p = self.forward[i, t-1]*self.backward[i, t-1]/np.inner(self.forward[:,t-1],self.backward[:,t-1])
        return p

    def sit(self, t, i, j):
        '''given lambda and obs, return probability of time t, state is i and t+1, state is j'''
        nume = self.forward[i,t-1] * self.A[i,j] * self.B.getprob(self.obs[t])[j] * self.backward[j,t]
        deno = (self.A.T*self.forward[i,t-1]).T * self.B.getprob(self.obs[t]) * self.backward[j, t]
        return np.sum(nume) / np.sum(deno)

    def viterbi(self):
        '''calculate the best path'''
        nstate, nt = self.A.shape[0], len(self.obs)
        # ini
        self.sigma = np.zeros((nstate,nt))
        self.sigma[:,0] = self.Pi * self.B.getprob(self.obs[0])
        self.psi = np.zeros_like(self.sigma,dtype='int')
        self.path = np.zeros(nt,dtype='int')
        # iter
        for t in range(1,nt):
            r = np.max(self.A.T * self.sigma[:,t-1], axis = 1)
            j = np.argmax(self.A.T * self.sigma[:,t-1], axis = 1)
            self.sigma[:,t] = r * self.B.getprob(self.obs[t])
            self.psi[:,t] = j
        self.path[nt-1] = np.argmax(self.sigma[:,nt-1])
        for t in range(nt-2,-1,-1):
            self.path[t] = self.psi[self.path[t+1], t+1]
        return self.path,np.max(self.sigma[:,nt-1])

class Emission(object):
    def __init__(self, B):
        self.B = B
        pass

    def getprob(self, o):
        '''obs at t;
        return a 1 x j'''
        return self.B[:,o]


if __name__ == '__main__':
    A = np.array([[0.5, 0.2, 0.3], [0.3, 0.5, 0.2],[0.2, 0.3, 0.5]])
    B = np.array([[0.5, 0.5], [0.4, 0.6], [0.7, 0.3]])
    Pi = np.array([0.2, 0.4, 0.4])
    obs = np.array([0, 1, 0],dtype='int')
    hmm = HMM(obs)
    emi = Emission(B)
    hmm.initial(A, emi, Pi)
    print hmm.viterbi(A,emi,Pi)
