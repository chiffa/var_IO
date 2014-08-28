__author__ = 'ank'

from csv import reader
from elasticsearch.helpers import bulk
from copy import deepcopy
import numpy as np
from matplotlib import pyplot as plt
from math import sqrt, ceil
from scipy.stats import gaussian_kde
from random import shuffle

from reqs import search_by_index, es, action_injector

source = '/home/ank/' 'projects_files/2014/BRU_GBO/4th_gen/het.ratio_result_nm.goodbatch.pub'


def automap(list1, list2):
    dico = {}
    if len(list1) != len(list2):
        raise Exception('Unequal_Lists')
    for key, val in zip(list1, list2):
        dico[key] = val
    return dico


def cst_float(flt):
    if flt != 'NA':
        return float(flt)
    return np.nan


with open(source, 'rb') as sf:
    rdr = reader(sf, 'excel-tab')
    title_line = rdr.next()
    shards = []
    gene_names = []
    for line in rdr:
        gene_names.append(line[0].split(':')[0])
        shards.append([cst_float(elt) for elt in line[1:]])


shards = np.array(shards)
plt.imshow(shards,interpolation='nearest')
plt.show()


def histogramize(name, index_list, conc_list):
    lst = []
    sqid = ceil(sqrt(len(index_list)+1))
    fig = plt.gcf()
    fig.canvas.set_window_title(name)
    for i, index in enumerate(index_list):
        prsh = shards[:,index]
        selsh = np.logical_not(np.isnan(prsh))
        reprsh = prsh[selsh]
        lst.append((np.mean(reprsh), np.std(reprsh)))
        density = gaussian_kde(reprsh)
        xs = np.linspace(np.min(reprsh), np.max(reprsh), 200)

        plt.subplot(sqid, sqid, i+1)
        plt.title(s=str(conc_list[i]))
        # plt.hist(reprsh, 100)
        plt.plot(xs, density(xs))

    plt.subplot(sqid, sqid, i+2)
    plt.title(s='summary')
    plt.plot(*zip(*lst), marker='o', color='r', ls='')
    plt.show()

    return lst


def split_titles():
    names = ['filename', 'cond1', 'conc1', 'unit1', 'cond2', 'conc2', 'unit2', 'generations', 'pool', 'scanner']
    cases = []
    for i, pld in enumerate(title_line[1:]):
        pld = automap(names, pld.split(':'))
        pld['idx'] = i
        pld = dict((key,val) for key, val in pld.iteritems() if val is not '')
        ac = deepcopy(action_injector)
        ac['_source'] = pld
        cases.append(ac)
    bulk(es, cases)


def crible(N):
    i=0
    src = range(0, len(title_line[1:]))
    done = [0]
    while i<N:
        shuffle(src)
        if src[0] not in done:
            name, idxs, concs = search_by_index(src[0])
            histogramize(name, idxs, concs)
            done = done + idxs
            i += 1





# print histogramize('benomyl', [12, 11, 2, 10], [2.4, 3.4, 13.8, 27.6])
# split_titles()
crible(10)
# print len(title_line)
