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
from reqs import search_by_index, es, action_injector, double_screen, get_similar, search_by_name
from sklearn.cluster import spectral_clustering
from scipy.linalg import eigh
from pickle import dump, load
from Utils.matrix_2D_embedding import embed
from Utils.Linalg_routines import hierchical_clustering, show_matrix_with_names
from pprint import pprint
from collections import Counter
from bunch import bunchify
from Utils.dataviz import better2D_desisty_plot

source = '/home/ank/' 'projects_files/2014/BRU_GBO/4th_gen/het.ratio_result_nm.goodbatch.pub'

active_rows = [472, 319, 487, 484, 468, 481, 325, 477, 324, 470, 471, 476, 466, 482, 715, 713, 712, 716, 717, 488, 473,
               469, 483, 714, 718, 82, 95, 207, 199, 303, 87, 213, 307, 316, 10, 163, 313, 69, 89, 292, 462, 464, 183,
               109, 168, 72, 485, 467, 474, 418, 200, 34, 55, 406, 637, 674, 676, 701, 500, 550, 103, 171, 538, 695,
               254, 666, 595, 516, 136, 498, 698, 632, 492, 626, 503, 496, 512, 581, 699, 644, 635, 636, 40, 41, 1,
               0, 2, 417, 580, 631, 191, 222, 602, 551, 299, 13, 50, 451, 645, 68, 696, 510, 723, 513, 641, 499, 693,
               633, 675, 649]


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
    msk1 = np.argsort(re_array1)[:20].tolist() + np.argsort(re_array1)[-20:].tolist()
    msk2 = np.argsort(re_array2)[:20].tolist() + np.argsort(re_array2)[-20:].tolist()
    tlt_msk = np.zeros(re_array1.shape).astype(dtype=np.bool)
    tlt_msk[msk1] = True
    tlt_msk[msk2] = True
    re_array1[np.logical_not(tlt_msk)] = 0
    re_array2[np.logical_not(tlt_msk)] = 0
    # msk1 = np.logical_and(re_array1<1, re_array1>-1)
    # msk2 = np.logical_and(re_array2<1, re_array2>-1)
    # re_array1[msk1]=0
    # re_array2[msk2]=0
    re_array1, re_array2 = (np.exp(re_array1), np.exp(re_array2))

    return np.sqrt(0.5*KL_div(re_array1, re_array2)+0.5*KL_div(re_array2, re_array1))


def labeled_imshow(matrix, label1=None, label2=None):
    plt.imshow(matrix, interpolation='nearest')
    if label1 is not None and label2 is None:
        plt.yticks(range(0, len(label1)), label1, rotation='horizontal')
        plt.xticks(range(0, len(label1)), label1, rotation='vertical')
        plt.subplots_adjust(left=0.2, bottom=0.2)
    if label1 is not None and label2 is not None:
        plt.yticks(range(0, len(label2)), label2, rotation='horizontal')
        plt.xticks(range(0, len(label1)), label1, rotation='vertical')
        plt.subplots_adjust(left=0.2, bottom=0.2)
    plt.show()


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


def reduce_shards(shards):
    exclusion_list = []
    non_nilled = []
    re_shards = np.zeros(shards.shape)
    for i in range(0, shards.shape[1]):
        if i not in exclusion_list and not double_screen(i):
            non_nilled.append(i)
            idxs, cumstring = get_similar(i)
            exclusion_list = exclusion_list + idxs
            sub_shards = shards[:,idxs]
            if sub_shards.shape[1]>1:
                safe_shards = np.ma.masked_array(sub_shards, np.isnan(sub_shards)).T
                shard_median = np.median(safe_shards, axis=0)
                shard_std = np.std(safe_shards, axis=0)
                plus_ctrl = shard_median + shard_std
                mins_ctrl = shard_median - shard_std
                flt1 = safe_shards > plus_ctrl
                flt2 = safe_shards < mins_ctrl
                plus_ctrl = plus_ctrl.reshape((1, safe_shards.shape[1]))
                mins_ctrl = mins_ctrl.reshape((1, safe_shards.shape[1]))
                rba1 = np.ma.masked_array(np.repeat(plus_ctrl, safe_shards.shape[0], axis=0).T, np.isnan(sub_shards)).T
                rba2 = np.ma.masked_array(np.repeat(mins_ctrl, safe_shards.shape[0], axis=0).T, np.isnan(sub_shards)).T
                safe_shards[flt1] = rba1[flt1]
                safe_shards[flt2] = rba2[flt2]
                mean_shards = np.mean(safe_shards, axis=0).T
                mean_shards.filled(np.nan)
                re_shards[:, i] = mean_shards
            else:
                re_shards[:,i] = shards[:,i]
    return re_shards, np.array(non_nilled)

with open(source, 'rb') as sf:
    rdr = reader(sf, 'excel-tab')
    title_line = rdr.next()
    pre_shards = []
    gene_names = []
    for line in rdr:
        gene_names.append(line[0].split(':')[0])
        pre_shards.append([cst_float(elt) for elt in line[1:]])


pre_shards = np.array(pre_shards)
shards, non_nilled = reduce_shards(pre_shards)

# plt.imshow(shards, interpolation='nearest')added costum cleaning to the data
# plt.show()
#
# better_hist(shards)
# better_hist(np.nanstd(shards, axis=1))
# better_hist(np.nanmean(shards, axis=1))

def show_range(idx):
    idxs, names = get_similar(idx)
    re_shards = pre_shards[:, idxs]
    show_matrix_with_names(re_shards, gene_names, names)



def histogramize(name, index_list, conc_list):
    lst = []
    sqid = ceil(sqrt(len(index_list)+2))
    fig = plt.gcf()
    fig.canvas.set_window_title(name)
    for i, index in enumerate(index_list):
        prsh = pre_shards[:,index]
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


def crible(N, support_labels=None, non_nilled=None):
    i=0
    done = [0]

    if support_labels is not None:
        src = range(0, len(non_nilled))
        padding = np.array(src)
        print padding.__len__()
        print support_labels.__len__()
        raw_input("Press Enter to continue...")
        insertor = np.zeros(padding.shape).tolist()
        for i in range(0, np.max(support_labels)+1):
            print i
            pre_padding = padding[support_labels==i].tolist()
            print pre_padding
            for j in pre_padding:
                # mismapping corrected here
                si = search_by_index(non_nilled[j], simple=True)
                print si
                insertor[j] = si
        return insertor

    else:
        src = range(0, len(title_line))
        while i<N and len(done) < len(src):
            shuffle(src)
            if src[0] not in done:
                name, idxs, concs = search_by_index(src[0])
                histogramize(name, idxs, concs)
                done = done + idxs
                i += 1


# np.array v.s. list ?

# define distance between columns
# define distance between rows
# cluster on columns
# cluster on Rows


def dist_matrix(axis, clusters):
    # if axis:
    #     re_shards = shards[:, non_nilled]
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
    vals[vals**2 < 0.3] = 0
    print vals
    # accumulator_matrix = np.dot(vects, np.dot(np.diag(vals), vects.T))
    plt.show()

    labels = spectral_clustering(accumulator_matrix, n_clusters=clusters, eigen_solver='arpack')
    print labels
    stable_mappings = crible(10, labels, non_nilled)
    print 'stable mappings redundancy:', len(stable_mappings), len(set(stable_mappings))
    srt_idx = hierchical_clustering(accumulator_matrix, labels)

    dump((stable_mappings, accumulator_matrix, srt_idx, non_nilled), open('loc_dump2.dmp','w'))

    # print np.array(stable_mappings)[srt_idx].tolist()
    #
    #
    # accumulator_matrix = accumulator_matrix[:,srt_idx]
    # accumulator_matrix = accumulator_matrix[srt_idx,:]
    # plt.imshow(accumulator_matrix, interpolation='nearest')
    # plt.colorbar()
    # plt.show()
    # # embed(pre_accumulator_matrix, stable_mappings)
    # groups = np.reshape(labels, (accumulator_matrix.shape[0], 1))
    # groupsets = []
    # for i in range(0, clusters):
    #     group_selector = groups==i
    #     group_idxs = group_selector.nonzero()[0].tolist()
    #     groupsets.append(group_idxs)
    # clustidx = np.array([item for itemset in groupsets for item in itemset])
    # accumulator_matrix = accumulator_matrix[:, clustidx]
    # accumulator_matrix = accumulator_matrix[clustidx, :]
    # plt.imshow(accumulator_matrix, interpolation='nearest')
    # plt.show()


def make_cuts():

    def pull_selection(idxs):
        print idxs
        print pprint(labels[idxs].tolist())


    labels, acc_matrix, srt_idx, non_nilled = load(open('loc_dump2.dmp','r'))
    labels = np.array(labels)[srt_idx]
    acc_matrix = acc_matrix[:,srt_idx]
    acc_matrix = acc_matrix[srt_idx,:]

    plt.imshow(acc_matrix, interpolation='nearest')
    # plt.yticks(range(0, 674), labels.tolist(), rotation='horizontal')
    # plt.xticks(range(0, 674), labels.tolist(), rotation='vertical')
    # plt.subplots_adjust(left=0.2, bottom=0.2)
    plt.show()


    selector = range(0, 44)+range(395, 461)

    pull_selection(selector)

    sub_ac = acc_matrix[selector, :]
    sub_ac = sub_ac[:, selector]
    sub_lbl = labels[selector]

    rt_idx = hierchical_clustering(sub_ac, sub_lbl)
    sub_ac = sub_ac[rt_idx,:]
    sub_ac = sub_ac[:,rt_idx]
    sub_lbl = sub_lbl[rt_idx]


    for elt in sub_lbl.tolist():
        print elt
    plt.imshow(sub_ac, interpolation='nearest')
    plt.yticks(range(0, len(selector)), sub_lbl.tolist(), rotation='horizontal')
    plt.xticks(range(0, len(selector)), sub_lbl.tolist(), rotation='vertical')
    plt.subplots_adjust(left=0.2, bottom=0.2)
    plt.show()

    # selector = range(0, 190)+range(290, 305)
    # sub_ac = sub_ac[selector,:]
    # sub_ac = sub_ac[:,selector]
    # sub_lbl = sub_lbl[selector]


    # plt.imshow(sub_ac, interpolation='nearest')
    # plt.yticks(range(0, 205), sub_lbl.tolist(), rotation='horizontal')
    # plt.xticks(range(0, 205), sub_lbl.tolist(), rotation='vertical')
    # plt.subplots_adjust(left=0.2, bottom=0.2)
    # plt.show()

    idxs = np.array(range(0, len(title_line[1:])))
    print '>>>>>>>>>>>>>>>>>> Indexes! <<<<<<<<<<<<<<<<<<<<<<'
    print idxs[non_nilled][srt_idx][selector][rt_idx].tolist()

    raw_input("Press Enter to continue...")

    clusters = []

    labels = spectral_clustering(sub_ac, n_clusters=5, eigen_solver='arpack')
    for el1, el2 in zip(sub_lbl[np.argsort(labels)].tolist(), np.array(labels)[np.argsort(labels)].tolist()):
        print el2, el1
    groups = np.reshape(labels, (sub_ac.shape[0], 1))
    groupsets = []
    for i in range(0, 15):
        group_selector = groups==i
        group_idxs = group_selector.nonzero()[0].tolist()
        groupsets.append(group_idxs)
    clustidx = np.array([item for itemset in groupsets for item in itemset])
    sub_ac = sub_ac[:, clustidx]
    sub_ac = sub_ac[clustidx, :]
    plt.imshow(sub_ac, interpolation='nearest')
    plt.colorbar()
    plt.show()


def determine_spread():
    # valid_shards = shards[:, active_rows]
    # valid_labels = [search_by_index(i, simple=True) for i in active_rows]
    #
    # dump((valid_labels, valid_shards, gene_names), open('loc_dump3.dmp','w'))
    #
    ##################################################################################

    valid_labels, valid_shards, gene_names = load(open('loc_dump3.dmp','r'))

    # labeled_imshow(valid_shards, valid_labels, gene_names)


    x = np.repeat(np.arange(len(valid_labels)), len(gene_names)).flatten()
    y = valid_shards.T.reshape((1, len(valid_labels)*len(gene_names))).flatten()

    print y.shape
    print x.shape

    better2D_desisty_plot(x, y)
    # plt.xticks(range(0, len(valid_labels)), valid_labels, rotation='vertical')
    plt.show()

# exclusion_list = ['itraconazole',]

# print histogramize('benomyl', [12, 11, 2, 10], [2.4, 3.4, 13.8, 27.6])
# split_titles()

# show_range(145)

crible(10)

# dist_matrix(1, 10)
# make_cuts()

determine_spread()

# print len(title_line)
