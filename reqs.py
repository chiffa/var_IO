__author__ = 'ank'

from elasticsearch import Elasticsearch
from Elasticsearch_plug import search_field, filter_field
from pprint import PrettyPrinter
import numpy as np
from bunch import bunchify

pp = PrettyPrinter(indent=4)

action_injector = {'_index':'experiments4','_type':'experiment'}

es = Elasticsearch()
es.indices.create(index=action_injector['_index'], ignore=400)



# print bunchify(es.indices.get_mapping(index=action_injector['_index'], doc_type=action_injector['_type']))

idx_config = {
    action_injector['_type']: {
        "properties": {
            "cond1": {
                "type": "multi_field",
                "fields": {
                    "cond1": {"type":"string"},
                    "original": {"type": "string", "index": "not_analyzed"}
                }
            }
        }
    }
}
es.indices.put_mapping(index=action_injector['_index'], doc_type=action_injector['_type'], body=idx_config)



def search_by_index(index, simple=False):
    rstring=[]

    q = search_field(term=index, fields=['idx'])
    r = es.search(index=action_injector['_index'],body=q)
    name = r['hits']['hits'][0]['_source']['cond1']


    if 'conc1' in r['hits']['hits'][0]['_source'].keys() and 'unit1' in r['hits']['hits'][0]['_source'].keys():
        rstring.append(name)
        rstring.append(r['hits']['hits'][0]['_source']['conc1'])
        rstring.append(r['hits']['hits'][0]['_source']['unit1'])

    if not simple:
        # q = search_field(term=name,fields=['co# define distance between columns
# define distance between rows
# cluster on columns
# cluster on Rowsnd1'])
        q = filter_field(term=name, field='cond1.original')
        r = es.search(index=action_injector['_index'], body=q)

        idx_set = []
        conc_set = []
        for elt in r['hits']['hits']:
            try:
                print elt['_source']['idx'], elt['_source']['cond1'], elt['_source']['conc1'], elt['_source']['unit1'], '\t-\t','conc2' in elt['_source'].keys()
                if not 'conc2' in elt['_source'].keys():
                    idx_set.append(elt['_source']['idx'])
                    conc_set.append(elt['_source']['conc1'])
            except KeyError:
                 if not 'conc2' in elt['_source'].keys():
                    print elt
                    idx_set.append(0)
                    try:
                        conc_set.append(elt['_source']['conc1'])
                    except KeyError:
                        conc_set.append('-1')

        idx_set = np.array([int(idx) for idx in idx_set])
        conc_set = np.array([float(conc) for conc in conc_set])

        idx_set = idx_set[np.argsort(conc_set)]
        conc_set = np.sort(conc_set)

        print idx_set, conc_set

        return name, idx_set.tolist(), conc_set.tolist()

    return ', '.join(rstring)


def double_screen(index):
    q = search_field(term=index, fields=['idx'])
    r = es.search(index=action_injector['_index'],body=q)

    return 'conc2' in r['hits']['hits'][0]['_source'].keys()


if __name__ == "__main__":
    # print search_by_index(2)
    pass