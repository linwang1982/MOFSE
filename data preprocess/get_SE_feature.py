"""
se994_info.xlsx: 从ADReCS下载的副作用信息
计算每个副作用的每个MedDRA术语的multi-hot编码
"""
import pandas as pd
import numpy as np
import scipy.io as sio

se_info_file = 'C:/Users/1/Desktop/664_drugs/se994_info.xlsx'
se_info = pd.read_excel(se_info_file, header=0, index_col=None)

# print(se_info)


def get_muti_hot(iterable, allowable_set):
    muti_hot = np.zeros(len(allowable_set)).astype(int)
    for x in iterable:
        if x not in allowable_set:
            raise Exception("input {} not in allowable set:".format(x))
        one_hot = np.array(list(map(lambda s: x == s, allowable_set)))
        muti_hot = muti_hot + one_hot
        muti_hot = (muti_hot != 0).astype(int)
    return muti_hot


level1, level2, level3, level4 = set(), set(), set(), set()
for ind, f in se_info.iterrows():
    ADRECS_ID = f['ADRECS_ID']
    ids = ADRECS_ID.split('/')
    for item in ids:
        level1.add(item[0:2])
        level2.add(item[0:5])
        level3.add(item[0:8])
        level4.add(item)

# print(level1, level2, level3, level4)
level1, level2, level3, level4 = list(level1), list(level2), list(level3), list(level4)
level1.sort()
level2.sort()
level3.sort()
level4.sort()

print(len(level1), len(level2), len(level3), len(level4))

SOC = np.empty((0, len(level1)))
HLGT = np.empty((0, len(level2)))
HLT = np.empty((0, len(level3)))
PT = np.empty((0, len(level4)))
for ind, f in se_info.iterrows():
    ADRECS_ID = f['ADRECS_ID']
    ids = ADRECS_ID.split('/')
    l1, l2, l3, l4 = [], [], [], []
    for item in ids:
        l1.append(item[0:2])
        l2.append(item[0:5])
        l3.append(item[0:8])
        l4.append(item)
    soc = get_muti_hot(l1, level1)
    hlgt = get_muti_hot(l2, level2)
    hlt = get_muti_hot(l3, level3)
    pt = get_muti_hot(l4, level4)

    SOC = np.vstack((SOC, soc))
    HLGT = np.vstack((HLGT, hlgt))
    HLT = np.vstack((HLT, hlt))
    PT = np.vstack((PT, pt))

print(SOC.shape, HLGT.shape, HLT.shape, PT.shape)
SHH = np.hstack((SOC, HLGT, HLT))
SH = np.hstack((SOC, HLGT))
SHHP = np.hstack((SOC, HLGT, HLT, PT))
dict = {'SOC': SOC, 'HLGT': HLGT, 'HLT': HLT, 'PT': PT, 'SHH': SHH, 'SH': SH, 'SHHP': SHHP}
sio.savemat('C:/Users/1/Desktop/664_drugs/se_feature.mat', dict)

