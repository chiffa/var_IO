__author__ = 'ank'


def simple_query(payload, layer_system):

    def deeper_wrap(layer, deeper):
        return {layer: deeper}

    ret = payload
    while layer_system:
        layer = layer_system.pop()
        ret = deeper_wrap(layer, ret)
    return ret


def filtered_query_header(bi_payload):
    return simple_query(bi_payload,['query','filtered'])


def match_all():
    return simple_query({},['match_all','query'])


def search_field(term, fields):
    return simple_query({"query":term, "fields":fields},['query','query_string'])


def filter_field(term, field):
    return simple_query({field:term}, ['query','constant_score','filter','term'])


def term_filter():
    pass

def range_filter():
    # lt  <
    # gt  >
    # gte >=
    # lte <=

    # range:{field_name}
    pass

def algebra_on_field():
    # bool: {operation:["term":{field:value}]}
    # and => must
    # or => should
    # not => must_not
    pass

def text_search(term):
    return  simple_query(term, ['query','query_string','query'])


def prepare_indexes():
    # prepare the location index
    # prepare the restaurant index
    pass


def location_query():

    pass

