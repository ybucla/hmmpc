# -*- coding: utf-8 -*-
"""
Created on Sat Feb 11 15:21:03 2017

@author: ybwang
"""
import re,sys,os
from optparse import OptionParser
from collections import defaultdict

def main():
#    usage = 'usage: python %prog [options] -i MSresultfile -s seqfile -u proteincolumn -p pepcolumn -r ratiocolumn -o output'
#    parser = OptionParser(usage)
#    parser.add_option('-i', dest='msfile', help='quantification file using MaxQuant or others [Default %default]')
#    parser.add_option('-s', dest='seqfile', help='fasta file used for the database search [Default %default]')
#    parser.add_option('-o', dest='outdir',default='output', help='output directory [Default %default]')
#    parser.add_option('-u', dest='protein',default=0, type='int', help='protein column, 0-based, like "IPI00021812.2" [Default %default]')
#    parser.add_option('-p', dest='pep',default=1, type='int', help='pep column, 0-based, like "_HRS(ph)NS(ph)FSDER_" [Default %default]')
#    parser.add_option('-r', dest='ratio',default=4, type='int', help='ratio column, 0-based, like "0.38957" [Default %default]')
#    (options, args) = parser.parse_args()
#    if options.msfile is None or options.seqfile is None or options.protein is None or options.pep is None or options.ratio is None:
#        sys.exit("[ERROR] "+parser.get_usage())
#    if not os.path.exists(options.outdir):
#        os.makedirs(options.outdir)

#    seqdict = readseq(options.seqfile)
#    parse(options.msfile,seqdict,options.protein,options.pep,options.ratio)
    seqdict = readseq('../ipi.HUMAN.v3.69.fasta')
    parse('data2.txt', seqdict, 0, 2, 3, 10)

def parse(infile,seqdict,proteinindex,pepindex,ratioindex, clusterdis = 15, rmfirst=True):
    data = defaultdict(dict)
    n = 0
    unmapprotein = set()
    with open(infile,'r') as f:
        for line in f:
            n += 1
            if rmfirst and n == 1: continue
            ele = line.rstrip('\r\n').split("\t")
            id, pep, ratio = ele[proteinindex], ele[pepindex], ele[ratioindex]
            site = int(ele[1]) # site index 1
            # remove 'K.', '.R', '-', '_', '*', and replace '(ph)' to '#' for convinence
            pep = re.sub(re.compile("(^_)|(_$)"),'',pep)
            pep = re.sub(re.compile("(^\w+\.|^-\.)|(\.\w+$|\.\-$)"),'',pep)
            pep = re.sub(re.compile("@"),'',pep)
            pep = re.sub(re.compile("\(ph\)"),'#',pep)
            pep = re.sub(re.compile("\(\w+\)"),'',pep)
            pep = re.sub(re.compile('\*'),'',pep)
            # if ratio is '' or 'NA', continue
            if ratio == '' or ratio == 'NA': continue
            ratio = float(ratio)
            l = getindex(pep)
            seq = seqdict[id]
            index = seq.find(re.sub(r'#','',pep)) + 1
            sites = [i + index - 1 for i in l]
            if site in sites:
                data[id][site] = ratio
            else:
                unmapprotein.add(id)
#                print n,site,l,sites
    # filter cluster
    clusterd = defaultdict(dict)
    for id in data:
        if id in unmapprotein or len(data[id]) < 2: continue
        sites = sorted(data[id])
        for i in range(len(sites)):
            if i == 0:
                if sites[1] - sites[0] + 1 <= clusterdis:
                    clusterd[id][sites[i]] = data[id][sites[i]]
            elif i == len(sites) - 1:
                if sites[i] - sites[i-1] + 1 <= clusterdis:
                    clusterd[id][sites[i]] = data[id][sites[i]]
            else:
                if sites[i + 1] - sites[i] + 1 <= clusterdis or sites[i] - sites[i-1] + 1 <= clusterdis:
                    clusterd[id][sites[i]] = data[id][sites[i]]
    result = []
    for id in clusterd:
        sites = [str(x) for x in sorted(clusterd[id])]
        if len(sites) >= 5:
            seq = seqdict[id]
            allindex = [str(x+1) for x in range(len(seq)) if 'STY'.find(seq[x]) != -1]
            result.append(id+'\t'+','.join(sites)+'\t' +','.join(allindex)+'\t'+seqdict[id])
    with open('siteinfo.txt','w') as fout:
        for r in result:
            fout.write(r+'\n')

def findadjacent(site, allsite, flank=7):
    left = site - flank
    right = site + flank
    ad = [i for i in allsite if i >= left and i <= right ]
    return ad

def getPeptideFlank(seq,position,left,right):
    newseq = ('*' * left) + seq + ("*" * right)
    start = left+position - left
    end = right + left + 1 - 1 + start
    return newseq[start-1:end]

def getindex(pep):
    l = []
    for j in range(len(pep)):
        if pep[j] == '#':
            l.append(j-len(l))
    return l

def readseq(seqfile):
    head = ''
    seqdict = {}
    with open(seqfile,'r') as f:
        for line in f:
            if line.find('>') != -1:
                head = line.rstrip().replace('>','')
                seqdict[head] = ''
            else:
                seqdict[head] += line.rstrip()
    return seqdict

if __name__ == '__main__':
	main()
    # seqdict = readseq('seq.txt')
    # parse('data_paper.txt',seqdict,0,1,4)
    # print getPeptideFlank('MSQVQVQVQNPSAALSGSQILNKNQSLLSQPLMSIPSTTSSLPSENAGRPIQNSALPSAS',58,7,7)
