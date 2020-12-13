#!/bin/python3
import pulp as p
import pandas as pd
import math
from functools import reduce
import pickle
import random
import sys
from pulp.constants import LpStatus

"""
Usage
-----
./fast_opt LAB_FILE DISTRICT_FILE OUTPUT_FILE
"""

def cmd_inputs(labfile="test1/lab_test_data.csv", distfile="test1/district_test_data.csv", outfile="fast_output.csv", picklefile=None):
    return labfile, distfile, outfile, picklefile


# set it to -1 for running forever
TIMEOUT = 240 #seconds

LAB_FILE, DISTRICT_FILE, OUTPUT_FILE, PICKLE_FILE= cmd_inputs(*sys.argv[1:])

print(f'LAB_FILE={LAB_FILE}\nDISTRICT_FILE={DISTRICT_FILE}')
print(f'OUTPUTFILE={OUTPUT_FILE}\nPICKLE_FILE={PICKLE_FILE}')


M = 100000
MAXLABDIST = 40


def parse_inputs():
    """update global variables lab_data, district_data"""
    
    global lab_data, district_data
    
    with open(LAB_FILE, 'r') as labfp:
        for line in labfp.readlines()[1:]:
            s = line.split(',')
            labid, lat, lon  = int(s[1]), float(s[2]), float(s[3])
            dist_id, is_pub = int(s[4]), int(s[5]) == 0 # 0 is public
            cap, back = int(s[6]), int(s[7])
            lab_data[labid] = {
                'pos': (lat, lon), 'dist_id': dist_id,
                'is_public': is_pub, 'capacity':cap,
                'backlog':back }

    with open(DISTRICT_FILE, 'r') as distfp:
        for line in distfp.readlines()[1:]:
            s = line.split(',')
            distid, name = int(s[1]), s[2]
            lat, lon = float(s[3]), float(s[4])
            sample = int(s[5])
            district_data[distid] = {
                'name': name,
                'pos': (lat, lon), 'samples': sample }

def haversine(lat1, lon1, lat2, lon2):
    radius = 6371 # km

    dlat = math.radians(lat2-lat1)
    dlon = math.radians(lon2-lon1)
    a = math.sin(dlat/2) * math.sin(dlat/2) + math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) * math.sin(dlon/2) * math.sin(dlon/2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
    d = radius * c

    return d

def district(j):
    """district of lab j"""
    return lab_data[j]['dist_id']

def labpos(j):
    return lab_data[j]['pos']

def generate_clusters():
    """ returns in 
        { district_id: [District labs]... ,
          k : [cluster k labs] } format"""

    clusters = {}

    for d in district_data.keys():
        clust = frozenset([j for j in lab_data.keys() if district(j) == d])
        clusters[d] = clust
    
    N = {}
    for lid in lab_data.keys():
        N[lid] = set([
            j for j in lab_data.keys()
                if j != lid and 
                   haversine(*labpos(j), *labpos(lid)) < MAXLABDIST])

    def BronKerbosch2(P, R=None, X=None):
        P = set(P)
        R = set() if R is None else R
        X = set() if X is None else X
        if not P and not X:
            yield R
        try:
            u = random.choice(list(P.union(X)))
            S = P.difference(N[u])
        # if union of P and X is empty
        except IndexError:
            S = P
        for v in S:
            yield from BronKerbosch2(
                P=P.intersection(N[v]), R=R.union([v]), X=X.intersection(N[v]))
            P.remove(v)
            X.add(v)
    k = 1
    for clique in BronKerbosch2(N.keys()):
        clique = frozenset(clique)
        if clique not in clusters.values():
            while k in clusters.keys():
                k += 1
            clusters[k] = clique
    print('generated clusters: ', clusters)
    return clusters

def samples(i):
    """returns samples at district i"""
    return district_data[i]['samples']

def labcost(j):
    if lab_data[j]['is_public']:
        return 800
    else:
        return 1600

def cluster_centre(cluster, cluster_id = None):
    if len(cluster) == 0:
        return (0., 0.)
    
    sm = reduce(lambda p1, p2: (p1[0]+p2[0], p1[1]+p2[1]), map(lambda j: lab_data[j]['pos'], cluster))
    return tuple(map(lambda x: x/(1. * len(cluster)), sm))
    
# print(lab_clusters())
# labid: pos, dist_id, is_public, capacity, backlog
lab_data = {}

# districtid: name, pos, samples
district_data = {}
    

if __name__ == '__main__':
    model = p.LpProblem(name="swab2lab", sense=p.LpMinimize)
    
    parse_inputs()
    clusters = generate_clusters()

    # setup variables
    ###################################
    ########### DEFINITIONS ###########
    ###################################

    zik = {} # Di to Ck
    tik = {} # zik > 0 ?
    xkj = {} # Ck to Lj
    yij = {} # Di to Lj if district(j) == i
    bi  = {} # Di backlog
    fikj= {} # flow Di -> Ck -> Lj

    for i in district_data.keys():
        bi[i] = p.LpVariable(f"b_{i}", lowBound=0)
    
    for i in district_data.keys():
        for j in lab_data.keys():
            for k in clusters.keys():
                if j in clusters[k]:
                    fikj[(i, k, j)] = p.LpVariable(f'f_{i}_{k}_{j}', lowBound=0)

    for i in district_data.keys():
        for k in clusters.keys():
            # zik[(i,k)] = p.LpVariable(f'z_{i}_{k}', lowBound=0, cat=p.LpInteger)
            zik[(i,k)] = p.lpSum([fikj[(i,k,j)] for j in lab_data.keys() if j in clusters[k]] + [0])
            tik[(i,k)] = p.LpVariable(f't_{i}_{k}', cat=p.LpBinary)

    for k in clusters.keys():
        for j in clusters[k]:
            # xkj[(k,j)] = m.Var(lb=0, integer=True)
            # xkj[(k,j)] = p.LpVariable(f'x_{k}_{j}', lowBound=0, cat=p.LpInteger)
            xkj[(k,j)] = p.lpSum([fikj[(i,k,j)] for i in district_data.keys()])

    for i in district_data.keys():
        for j in lab_data.keys():
            if district(j) == i:
                # Constraint 1,  yij <= 100 for district(j) == i
                # yij[(i,j)] = m.Var(lb=0, ub=100, integer=True)
                yij[(i,j)] = p.LpVariable(f'y_{i}_{j}', lowBound=0, upBound=100)

    # Constraints 2
    # sum_i zik = sum_j xkj if j in cluster k
    for k in clusters.keys():
        
        iks = list(filter(lambda ik: ik[1] == k , zik.keys()))

        # if j in cluster k taken care by definition
        kjs = list(filter(lambda kj: kj[0] == k , xkj.keys())) 

        if len(iks) > 0 or len(kjs) > 0:
            model += (  p.lpSum(list(map(lambda ik: zik[ik],  iks   ))+[0]) ==\
                        p.lpSum(list(map(lambda kj: xkj[kj],  kjs   ))+[0]))


    # Constraint 3
    # sum_k xkj <= CAPj - BACKj
    for j in lab_data.keys():
        # m.Equation(
        #     reduce(lambda x1,x2: x1+x2,
        #            map(lambda kj: xkj[kj],
        #                         filter(lambda kj: kj[1] == j,
        #                                xkj.keys())), 0) \
        #     <= \
        #     (lab_data[j]['capacity'] - lab_data[j]['backlog'])
        # )
        kjs = filter(lambda kj: kj[1] == j, xkj.keys())
        model += (p.lpSum(list(map(lambda kj: xkj[kj], kjs)) + [0]) \
                    <= (lab_data[j]['capacity'] - lab_data[j]['backlog']))

    # Constraint 4
    # sum_j yij + sum_k zik + bi = Si
    for i in district_data.keys():
        # ENSURE: yij contains a (i, j) if district(j) == i
        
        sumyij = p.lpSum(list(map(lambda ij: yij[ij],
                            filter(lambda ij: ij[0] == i,
                                   yij.keys()))))
        sumzik = p.lpSum(list(map(lambda ik: zik[ik],
                            filter(lambda ik: ik[0] == i,
                                   zik.keys()))))

        model += ( sumyij + sumzik + bi[i] == samples(i) )


    # Constraint 5
    # zik - Mtik <= 0 if k != i
    for ik in tik.keys():
        i, k = ik
        if k != i:
            model += (zik[ik] - M*tik[ik] <= 0)

    # Constraint 6
    # sum_k tik <= 1
    for i in district_data.keys():
        model += p.lpSum(list(map(lambda ik: tik[ik],
                           filter(lambda ik: ik[0] == i,
                                  tik.keys())))) <= 1

    print(f'{len(yij):=}\t{len(tik):=}\t{len(xkj):=}\t{len(zik):=}')
    


    # objective
    
    model+= p.lpSum(list(map(lambda kj: labcost(kj[1])*xkj[kj], xkj.keys()))) +\
            10000* p.lpSum(list(bi.values())) +\
            5000 * p.lpSum(list(yij.values())) +\
            1000 * p.lpSum(list(map(lambda ik: tik[ik]*haversine(*cluster_centre(clusters[ik[1]], ik[1]), *district_data[ik[0]]['pos']),
                          tik.keys())))
    
    if TIMEOUT > 0:
        model.solve(p.PULP_CBC_CMD(timeLimit=TIMEOUT))
    else:
        model.solve()
    


    if PICKLE_FILE != None:
        with open(PICKLE_FILE, 'wb') as fp:
            print('writing to', PICKLE_FILE)
            pickle.dump({'yij': yij, 'xkj': xkj, 'zik': zik} , fp)
    print("DONE")

    # net transfer from Di to Lj is (sum_k fikj) + yij
    
    # writing to OUTPUT_FILE
    with open(OUTPUT_FILE, 'w') as outfp:
        print(f'Writing to {OUTPUT_FILE}')
        print('transfer_type,source,destination,samples_transferred', file=outfp)
        trans = {}
        for i in district_data.keys():
            if p.value(bi[i]) > 0:
                print(f'1,{i},{i},{p.value(bi[i])}',file=outfp)
            for j in lab_data.keys():
                trans[(i,j)] = p.value(yij[(i,j)]) if (i,j) in yij.keys() else 0
                for k in clusters.keys():
                    if j in clusters[k]:
                        trans[(i,j)] += p.value(fikj[(i,k,j)])
                # trans[(i,j)] = int(trans[(i,j)])
                if(trans[(i, j)] > 0):
                    print(f'0,{i},{j},{trans[(i,j)]}',file=outfp)


    # debugging
    for i in district_data.keys():
        for k in clusters.keys():
            if p.value(zik[(i, k)]) != 0 or p.value(tik[(i,k)]) != 0:
                print(f'z_{i}_{k} = {(p.value(zik[(i,k)]))}\tt_{i}_{k} = {p.value(tik[(i,k)])}')
