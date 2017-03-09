# -*- coding: utf-8 -*-
"""
Created on Tue Mar 07 18:37:33 2017

@author: ybwang
"""

# transfer gps2.1 result to data input
import re,sys,os
from optparse import OptionParser
from collections import defaultdict
from hmm import *
import numpy as np
import scipy.stats

def readGpsData(gpsfile, map_code, map_kina):
    kinadict = defaultdict(dict)
    codedict = defaultdict(dict)
    with open(gpsfile,'r') as f:
        for line in f:
            ele = line.rstrip().split('\t')
            if line.find('Position') != -1: continue
            id, site, code = ele[0], int(ele[1]), map_code[ele[2]]
            kina = map_kina[ele[3]] if ele[3] in map_kina else None
            if kina is None: continue
            codedict[id][site] = code
            if site in kinadict[id]:
                kinadict[id][site].append(kina)
            else:
                kinadict[id][site] = [kina]
    data = {}
    for id in kinadict:
        sites = sorted(kinadict[id])
        v = np.zeros((len(sites),len(map_kina)+3),dtype='int')
        for i in range(len(sites)):
            values = [sites[i]]
            if i == 0:
                values.append(0)
            else:
                values.append(sites[i]- sites[i-1])
            values.append(codedict[id][sites[i]])
            vkinase = np.zeros(len(map_kina),dtype='int')
            vkinase[np.array(kinadict[id][sites[i]])] = 1
            values += vkinase.tolist()
            v[i,:] = values
        data[id] = v
    return data

def readParameter(parafile):
    shape = (0,0)
    tag = np.zeros(8)
    transition, emssision_code, emssision_kina, emssision_distance = None, None, None, defaultdict(list)
    pi = {}
    map_states,map_codes,map_kina = {},{},{}
    with open(parafile,'r') as f:
        for line in f:
            line = line.rstrip()
            if line.find('>') != -1:
                tag = np.zeros(tag.size)
            if line.find('>transition') != -1:
                tag[0] = 1
                ele = line.split(',')
                shape = (int(ele[1]),int(ele[2]))
            if line.find('>emssision_code') != -1:
                tag[1] = 1
                ele = line.split(',')
                shape = (int(ele[1]),int(ele[2]))
            if line.find('>emssision_kina') != -1:
                tag[2] = 1
                ele = line.split(',')
                shape = (int(ele[1]),int(ele[2]))
            if line.find('>emssision_distance') != -1:
                tag[3] = 1
            if line.find('>pi') != -1:
                tag[4] = 1
            if line.find('>map_states') != -1:
                tag[5] = 1
            if line.find('>map_codes') != -1:
                tag[6] = 1
            if line.find('>map_kina') != -1:
                tag[7] = 1

            if tag[0] == 1 and line.find('>') == -1:
                transition = np.array([float(x) for x in line.split(',')]).reshape(shape)
            if tag[1] == 1 and line.find('>') == -1:
                emssision_code = np.array([float(x) for x in line.split(',')]).reshape(shape)
            if tag[2] == 1 and line.find('>') == -1:
                emssision_kina = np.array([float(x) for x in line.split(',')]).reshape(shape)
            if tag[3] == 1 and line.find('>') == -1:
                ele = line.split('\t')
                emssision_distance[int(ele[0])] = [int(x) for x in ele[1].split(',')]
            if tag[4] == 1 and line.find('>') == -1:
                pi = np.array([float(x) for x in line.split(',')])
            if tag[5] == 1 and line.find('>') == -1:
                ele = line.split('\t')
                map_states[ele[0]]=int(ele[1])
            if tag[6] == 1 and line.find('>') == -1:
                ele = line.split('\t')
                map_codes[ele[0]]=int(ele[1])
            if tag[7] == 1 and line.find('>') == -1:
                ele = line.split('\t')
                map_kina[ele[0]]=int(ele[1])

    return {'transition':transition,'pi':pi,'emssision_code':emssision_code,'emssision_kina':emssision_kina,'emssision_distance':emssision_distance,'map_states':map_states,'map_codes':map_codes,'map_kina':map_kina}

class Emssi(Emission):

    def __init__(self, paradict):
        self.emssision_code = paradict['emssision_code']
        self.emssision_kina = paradict['emssision_kina']
        self.emssision_distance = paradict['emssision_distance']
        self.onekde = scipy.stats.gaussian_kde(self.emssision_distance[1])
        self.twokde = scipy.stats.gaussian_kde(self.emssision_distance[2])
        self.fourkde = scipy.stats.gaussian_kde(self.emssision_distance[4])
        self.fivekde = scipy.stats.gaussian_kde(self.emssision_distance[5])

    def getprob(self, obs):
        '''obs at t;
        return a 1 x j'''
        dist, code, kina = obs[1], obs[2], obs[3:]
        prob_dist = self.getdistp(dist)
        prob_code = self.emssision_code[:,code]
        prob_kina = self.emssision_kina[:, np.nonzero(kina)[0]].sum(1)
        return prob_dist * prob_code

    def getdistp(self, distance):
        p = []
        if distance == 0:
            p = [1.0,0.0,0.0,1.0,0.0,0.0]
        else:
            p = [0.0,self.onekde.evaluate(distance)[0],self.twokde.evaluate(distance)[0],0.0,self.fourkde.evaluate(distance)[0],self.fivekde.evaluate(distance)[0]]
        return p

if __name__ == '__main__':
    paradict = readParameter('parameter.txt')
#    paradict['pi'] = np.array([0.5,0.0,0.0,0.5,0.0,0.0])
    data = readGpsData('gps2.1.txt',paradict['map_codes'],paradict['map_kina'])
    obs = data['IPI00219301.7']
    hmm = HMM(obs)
    hmm.initial(A = paradict['transition'], B = Emssi(paradict), Pi = paradict['pi'])
    path, prob = hmm.viterbi()
    map_states = dict([(y,x) for x,y in paradict['map_states'].items()])
    path = [map_states[x] for x in path]
    print path
    print prob
    np.savetxt('sigma.txt',hmm.sigma,fmt='%1.2e',delimiter='\t')
    np.savetxt('psi.txt',hmm.psi,fmt='%d',delimiter='\t')

