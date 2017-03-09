# -*- coding: utf-8 -*-
"""
Created on Tue Mar 07 15:53:05 2017

@author: ybwang
"""
import re,sys,os
from optparse import OptionParser
from collections import defaultdict
import numpy as np
import scipy.stats

def readGPS(infile):
    gpsdict = defaultdict(dict)
    familydict = defaultdict(int)
    max = 0
    with open(infile,'r')as f:
        for line in f:
            ele = line.rstrip().split('\t')
            if ele[3] not in familydict:
                familydict[ele[3]] = max
                max += 1
            if ele[1] in gpsdict[ele[0]]:
                gpsdict[ele[0]][ele[1]].append(ele[3])
            else:
                gpsdict[ele[0]][ele[1]] = [ele[3]]
    return gpsdict,familydict

def readsite(infile,gpsdict, clusterdis = 15):
    data = []
    with open(infile, 'r') as f:
        for line in f:
            ele = line.rstrip().split('\t')
            id, psites, asites, seq = ele[0],ele[1].split(','),ele[2].split(','),ele[3]
            states = []
            finalsites = []
            finalcodes = []
            for site in asites:
                if site not in gpsdict[id]: continue
                if site in psites:
                    if len(states) == 0:
                        states.append('S+')
                    else:
                        if states[-1] == '-' or states[-1] == 'I-' or states[-1] == 'S-':
                            states.append('I+')
                        elif states[-1] == '+':
                            if int(site) - int(finalsites[-1]) + 1 <= clusterdis:
                                states.append('+')
                            else:
                                states.append('I+')
                        elif states[-1] == 'I+' or states[-1] == 'S+':
                            states.append('+')
                else:
                    if len(states) == 0:
                        states.append('S-')
                    else:
                        if states[-1] == '+':
                            if int(site) - int(finalsites[-1]) + 1 <= clusterdis:
                                states.append('+')
                            else:
                                states.append('I-')
                        elif states[-1] == 'S-' or states[-1] == '-' or states[-1] == 'I-':
                            states.append('-')
                finalsites.append(site)
                finalcodes.append(seq[int(site)-1])
            if len(finalsites) == 0: continue
            data.append(ele[0] + '\t' + ','.join(finalsites) + '\t' +','.join(finalcodes)+'\t' +','.join(states))
    return data

def toPara(datalist, gpsdict, familydict):
    statesname = {'S+':0,'+':1,'I+':2,'S-':3,'-':4,'I-':5}
    codename = {'S':0,'T':1,'Y':2}  # discrete
    transition = np.zeros((6, 6), dtype='float')
    emssision_code = np.zeros((6, len(codename)), dtype='float')
    emssision_kina = np.zeros((6, len(familydict)), dtype='float')
    distancedict = defaultdict(list)
    pidict = np.zeros(6)
    for line in datalist:
        ele = line.split('\t')
        id, sites, codes, states = ele[0],ele[1].split(','),ele[2].split(','),ele[3].split(',')
        # pi
        pidict[statesname[states[0]]] += 1.0 / len(datalist)
        for i in range(len(states)):
            s = statesname[states[i]]
            if i > 0:
                s_pre = statesname[states[i-1]]
                transition[s_pre, s] += 1.0
            # code
            c = codename[codes[i]]
            emssision_code[s, c] += 1.0
            # kinase
            site = sites[i]
            kina = np.array([familydict[x] for x in gpsdict[id][site]])
            emssision_kina[s, kina] += 1.0
            # distance to preceding site, P(d|+), P(d|s) s={I+, I-, -},  P(d = 0|s) = 1 and P(d != 0|s) = 0  s={S+,S-}
            if i > 0:
                distance = int(sites[i]) - int(sites[i-1])
                if states[i] in ['I+', 'I-', '-']:
                    distancedict[statesname['I+']].append(distance)
                    distancedict[statesname['I-']].append(distance)
                    distancedict[statesname['-']].append(distance)
                if states[i] in ['+']:
                    distancedict[statesname['+']].append(distance)
            if states[i] in ['S+', 'S-']:
                distancedict[statesname['S+']] = [0]
                distancedict[statesname['S-']] = [0]

    # transition: matrix
    transition = (transition.T / np.sum(transition,axis=1)).T
    # pi
    # emssision: code, discrete
    emssision_code = (emssision_code.T / np.sum(emssision_code,axis=1)).T
    # emssision: kinase family, discrete
    emssision_kina = (emssision_kina.T / np.sum(emssision_kina,axis=1)).T
    # emssision: distance, continuous, using scipy.stats.gaussian_kde(sample)
    # save parameter
    with open('parameter.txt' ,'w') as fout:
        fout.write('>transition,'+str(transition.shape[0])+','+str(transition.shape[1])+'\n')
        fout.write(','.join(transition.flatten().astype('str').tolist())+'\n')
        fout.write('>emssision_code,'+str(emssision_code.shape[0])+','+str(emssision_code.shape[1])+'\n')
        fout.write(','.join(emssision_code.flatten().astype('str').tolist())+'\n')
        fout.write('>emssision_kina,'+str(emssision_kina.shape[0])+','+str(emssision_kina.shape[1])+'\n')
        fout.write(','.join(emssision_kina.flatten().astype('str').tolist())+'\n')
        fout.write('>emssision_distance\n')
        for i in distancedict:
            s = [str(x) for x in distancedict[i]]
            fout.write(str(i)+'\t'+','.join(s)+'\n')
        fout.write('>pi\n')
        fout.write(','.join(pidict.astype('str').tolist())+'\n')
        fout.write('>map_states\n')
        for i in statesname:
            fout.write(i+'\t'+str(statesname[i])+'\n')
        fout.write('>map_codes\n')
        for i in codename:
            fout.write(i+'\t'+str(codename[i])+'\n')
        fout.write('>map_kina\n')
        for i in familydict:
            fout.write(i+'\t'+str(familydict[i])+'\n')

if __name__ == '__main__':
    gpsdict,familydict = readGPS('gps2.1.txt')
    datalist = readsite('siteinfo.txt', gpsdict, 10)
    toPara(datalist, gpsdict, familydict)
    for d in datalist:
        print d
