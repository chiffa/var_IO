__author__ = 'ank'

from csv import reader
from elasticsearch.helpers import bulk
from copy import deepcopy
import numpy as np
from matplotlib import pyplot as plt
from math import sqrt, ceil
from scipy.stats import gaussian_kde
from random import shuffle
from scipy.stats import entropy as KL_div
from itertools import combinations
from reqs import search_by_index, es, action_injector, double_screen
from sklearn.cluster import spectral_clustering
from scipy.linalg import eigh
from pickle import dump, load
from Utils.matrix_2D_embedding import embed
from Utils.Linalg_routines import hierchical_clustering

# TODO: generate the buffering of the precomupted matrix

# TODO: > remove double screens => ok
# TODO: > deconvolve the neutral genes (~0 effect) from significantly affected => ok

# TODO: try to laplace-normalize the matrices before clustering


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


def better_hist(to_hist):
    print to_hist.shape
    print to_hist
    selsh = np.logical_not(np.isnan(to_hist))
    reprsh = to_hist[selsh]
    density = gaussian_kde(reprsh)
    xs = np.linspace(np.min(reprsh), np.max(reprsh), 200)
    plt.plot(xs, density(xs))
    plt.show()


def distance(np_array1,np_array2):
    msk = np.logical_not(np.logical_or(np.isnan(np_array1), np.isnan(np_array2)))
    re_array1, re_array2 = (np_array1[msk], np_array2[msk])
    msk1 = np.logical_and(re_array1<1, re_array1>-1)
    msk2 = np.logical_and(re_array2<1, re_array2>-1)
    re_array1[msk1]=0
    re_array2[msk2]=0
    re_array1, re_array2 = (np.exp(re_array1), np.exp(re_array2))

    return np.sqrt(0.5*KL_div(re_array1, re_array2)+0.5*KL_div(re_array2, re_array1))




def dist_mat(index_list, axis):
    if axis:
        re_shards = shards[:, index_list]
    else:
        re_shards = shards[index_list, :]
    print re_shards.shape
    c_list = np.split(re_shards, re_shards.shape[axis], axis)
    accumulator_matrix = np.zeros((re_shards.shape[axis], re_shards.shape[axis]))
    est_len = re_shards.shape[axis]*(re_shards.shape[axis]-1)/2
    for i, (i_a, i_b) in enumerate(combinations(range(0, re_shards.shape[axis]), 2)):
        if not i%100:
            pl = "{0:0.2f}".format(i/float(est_len)*100.0)
            print pl, '%'
        a = c_list[i_a]
        b = c_list[i_b]
        dist = distance(a, b)
        accumulator_matrix[i_a, i_b] = dist
        accumulator_matrix[i_b, i_a] = dist
    return accumulator_matrix



with open(source, 'rb') as sf:
    rdr = reader(sf, 'excel-tab')
    title_line = rdr.next()
    shards = []
    gene_names = []
    for line in rdr:
        gene_names.append(line[0].split(':')[0])
        shards.append([cst_float(elt) for elt in line[1:]])


shards = np.array(shards)
# plt.imshow(shards, interpolation='nearest')
# plt.show()
#
# better_hist(shards)
# better_hist(np.nanstd(shards, axis=1))
# better_hist(np.nanmean(shards, axis=1))


def histogramize(name, index_list, conc_list):
    lst = []
    sqid = ceil(sqrt(len(index_list)+2))
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
    plt.title(s='mean v.s. std')
    plt.plot(*zip(*lst), marker='o', color='r', ls='')
    if len(index_list)>1:
        plt.subplot(sqid, sqid, i+3)
        plt.title(s='mutual_distances')
        plt.imshow(dist_mat(index_list ,1), interpolation='nearest')
        plt.colorbar()
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


def crible(N, support_matrix=None):
    i=0
    src = range(0, len(title_line[1:]))
    padding = np.array(src)
    insertor = np.zeros(padding.shape).tolist()
    done = [0]

    if support_matrix is not None:
        for i in range(0, np.max(support_matrix)+1):
            print i
            pre_padding = padding[support_matrix==i].tolist()
            print pre_padding
            for i in pre_padding:
                si = search_by_index(i, simple=True)
                print si
                insertor[i] = si
        return insertor

    else:
        while i<N and len(done) < len(src):
            shuffle(src)
            if src[0] not in done:
                name, idxs, concs = search_by_index(src[0])
                histogramize(name, idxs, concs)
                done = done + idxs
                i += 1


# define distance between columns
# define distance between rows
# cluster on columns
# cluster on Rows


def dist_matrix(axis, clusters):
    # fltr = np.ones(shards.shape[axis]).astype(dtype=np.int16)
    # if axis:
    #     for i in range(0, shards.shape[axis]):
    #         if double_screen(i):
    #             fltr[i] = 0
    #     re_shards = shards[:, fltr>0]
    # else:
    #     re_shards = shards
    # print re_shards.shape
    # c_list = np.split(re_shards, re_shards.shape[axis], axis)
    # accumulator_matrix = np.zeros((re_shards.shape[axis], re_shards.shape[axis]))
    # est_len = re_shards.shape[axis]*(re_shards.shape[axis]-1)/2
    # for i, (i_a, i_b) in enumerate(combinations(range(0, re_shards.shape[axis]), 2)):
    #     if not i%100:
    #         pl = "{0:0.2f}".format(i/float(est_len)*100.0)
    #         print pl, '%'
    #     a = c_list[i_a]
    #     b = c_list[i_b]
    #     dist = distance(a, b)
    #     accumulator_matrix[i_a, i_b] = dist
    #     accumulator_matrix[i_b, i_a] = dist
    #
    # dump(accumulator_matrix,open('loc_dump.dmp','w'))

    ##########################################################################################################

    pre_accumulator_matrix = load(open('loc_dump.dmp','r'))

    accumulator_matrix = np.exp( - pre_accumulator_matrix*pre_accumulator_matrix / pre_accumulator_matrix.std() )
    plt.imshow(accumulator_matrix, interpolation='nearest')
    plt.show()




    vals, vects =  eigh(accumulator_matrix)
    plt.hist(vals, 1000, log=True)
    vals[vals**2 < 0.9] = 0
    print vals
    # accumulator_matrix = np.dot(vects, np.dot(np.diag(vals), vects.T))
    plt.show()

    labels = spectral_clustering(accumulator_matrix, n_clusters=clusters, eigen_solver='arpack')
    print labels
    stable_mappings = crible(10, labels)
    srt_idx = hierchical_clustering(accumulator_matrix, labels)
    print np.array(stable_mappings)[srt_idx].tolist()

    # TODO: define set of cuts, re-order the matrix and cluster it according to the cuts

    accumulator_matrix = accumulator_matrix[:,srt_idx]
    accumulator_matrix = accumulator_matrix[srt_idx,:]
    plt.imshow(accumulator_matrix)
    plt.colorbar()
    plt.show()
    # embed(pre_accumulator_matrix, stable_mappings)
    groups = np.reshape(labels, (accumulator_matrix.shape[0], 1))
    groupsets = []
    for i in range(0, clusters):
        group_selector = groups==i
        group_idxs = group_selector.nonzero()[0].tolist()
        groupsets.append(group_idxs)
    clustidx = np.array([item for itemset in groupsets for item in itemset])
    accumulator_matrix = accumulator_matrix[:, clustidx]
    accumulator_matrix = accumulator_matrix[clustidx, :]
    plt.imshow(accumulator_matrix, interpolation='nearest')
    plt.show()

# print histogramize('benomyl', [12, 11, 2, 10], [2.4, 3.4, 13.8, 27.6])
# split_titles()

# crible(10)

dist_matrix(1, 10)



# print len(title_line)
