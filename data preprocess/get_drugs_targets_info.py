"""
drug_all_info_664.pkl: 从Drugbank中获取的每个药物的相关信息
从靶点信息中获得各个信息的 multi-hot code 作为其特征
"""

import pickle
import scipy.io as sio
import numpy as np
import pandas as pd


class Drug:
    def __init__(self, drugbank_id=None, name=None, synonym=[], type=None, smiles=None, targets=[]):
        self.type = type
        self.name = name
        self.drugbank_id = drugbank_id
        self.synonym = synonym
        self.smiles = smiles
        self.targets = targets

# targets is:
#   list(target1_dict, target2_dict ...)
# target_dict is:
# {'uniprot_id': '', 'target_name': '', 'sequence': '', 'gene': '', 'sequence_other': '','gene_other': '',
# 'drugbank_id': '', 'pfam_domain':[], 'target_synonym': [],'GO': {'component': [], 'function': [], 'process': []}}

with open("C:/Users/1/Desktop/664_drugs/drug_all_info_664.pkl", "rb") as f:
    drugs = pickle.load(f)
    print('完成(%d)' % len(drugs))


d = drugs[663]
print(d.drugbank_id)
print(d.targets)
# exit(1)

all_domain = set()
all_compent = set()
all_function = set()
all_process = set()

for i, _drug in enumerate(drugs):
    ts = _drug.targets
    # print(type(ts))
    for t in ts:
        domain = t['pfam_domain']
        # print(i, domain)
        compent = t['GO']['component']
        function = t['GO']['function']
        process = t['GO']['process']
        all_domain = all_domain.union(set(domain))
        all_compent = all_compent.union(set(compent))
        all_function = all_function.union(set(function))
        all_process = all_process.union(set(process))

print(len(all_domain), len(all_compent), len(all_function), len(all_process))


def get_muti_hot(iterable, allowable_set):
    muti_hot = np.zeros(len(allowable_set)).astype(int)
    for x in iterable:
        if x not in allowable_set:
            raise Exception("input {} not in allowable set:".format(x))
        one_hot = np.array(list(map(lambda s: x == s, allowable_set)))
        muti_hot = muti_hot + one_hot
    return muti_hot


all_function, all_process, all_compent = list(all_function), list(all_process), list(all_compent)
GO_info = np.empty((0, len(all_compent)+len(all_function)+len(all_process)))
Domain_info = np.empty((0, len(all_domain)))
F_info = np.empty((0, len(all_function)))
P_info = np.empty((0, len(all_process)))
C_info = np.empty((0, len(all_compent)))

for i, _drug in enumerate(drugs):
    ts = _drug.targets
    domain, component, function, process = set(), set(), set(), set()
    for t in ts:
        domain = domain.union(t['pfam_domain'])
        component = component.union(t['GO']['component'])
        function = function.union(t['GO']['function'])
        process = process.union(t['GO']['process'])

    D_code = get_muti_hot(domain, all_domain)
    C_code = get_muti_hot(component, all_compent)
    P_code = get_muti_hot(process, all_process)
    F_code = get_muti_hot(function, all_function)
    # print(D_code)

    F_info = np.vstack((F_info, F_code))
    P_info = np.vstack((P_info, P_code))
    C_info = np.vstack((C_info, C_code))
    GO_code = np.hstack((F_code, P_code, C_code))
    Domain_info = np.vstack((Domain_info, D_code))
    GO_info = np.vstack((GO_info, GO_code))

print(Domain_info.shape)
print(GO_info.shape)

dict = {'domain': Domain_info, 'GO': GO_info, 'function': F_info, 'process': P_info, 'component': C_info}
sio.savemat('C:/Users/1/Desktop/664_drugs/domain_GO_{}.mat'.format(len(drugs)), dict)



# P = set()
# F = set()
# C = set()
# for i, _drug in enumerate(drug692):
#     ts = _drug.targets
#     # print(type(ts))
#     for t in ts:
#
#         compent = t['GO']['component']
#         function = t['GO']['function']
#         process = t['GO']['process']
#
#         if 'receptor binding' in process:
#             P.add(_drug.drugbank_id)
#         if 'receptor binding' in function:
#             F.add(_drug.drugbank_id)
#         if 'receptor binding' in compent:
#             C.add(_drug.drugbank_id)
#
# print(P)
# print(F)
# print(C)
# # DB01406   receptor binding
# # DB01188 DB01244 DB00390   sodium  'Sodium/potassium-transporting ATPase subunit alpha-1'
#
# for i, _drug in enumerate(drug692):
#     if 'DB01406' == _drug.drugbank_id:
#         for t in _drug.targets:
#             compent = t['GO']['component']
#             function = t['GO']['function']
#             print(t['target_name'])
#             print('compent', compent)
#             print('function', function)
