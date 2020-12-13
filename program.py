#!/bin/python3
import gekko
import pandas as pd
import math
from functools import reduce

LAB_FILE      = 'test1/lab_test_data.csv'
DISTRICT_FILE = 'test1/district_test_data.csv'
M = 1000000
APOPT_SOLVER_OPTIONS = ['minlp_maximum_iterations 500', \
                        # minlp iterations with integer solution
                        'minlp_max_iter_with_int_sol 10', \
                        # treat minlp as nlp
                        'minlp_as_nlp 0', \
                        # nlp sub-problem max iterations
                        'nlp_maximum_iterations 50', \
                        # 1 = depth first, 2 = breadth first
                        'minlp_branch_method 1', \
                        # maximum deviation from whole number
                        'minlp_integer_tol 0.0000001', \
                        # covergence tolerance
                        'minlp_gap_tol 0.01']

def parse_inputs():
    """update global variables lab_data, district_data"""
    
    global lab_data, district_data
    
    with open(LAB_FILE) as labfp:
        for line in labfp.readlines()[1:]:
            s = line.split(',')
            labid, lat, lon  = int(s[1]), float(s[2]), float(s[3])
            dist_id, is_pub = int(s[4]), int(s[5]) == 0 # 0 is public
            cap, back = int(s[6]), int(s[7])
            lab_data[labid] = {
                'pos': (lat, lon), 'dist_id': dist_id,
                'is_public': is_pub, 'capacity':cap,
                'backlog':back }

    with open(DISTRICT_FILE) as distfp:
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
    dlat = math.radians(lat2-lat1); dlon = math.radians(lon2-lon1)
    a = math.sin(dlat/2) * math.sin(dlat/2) + math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) * math.sin(dlon/2) * math.sin(dlon/2)
    d = 2 * radius * math.atan2(math.sqrt(a), math.sqrt(1-a))
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
        clust = [j for j in lab_data.keys() if district(j) == d]
        clusters[d] = clust

    nondist_clusters = set()
    for lid in lab_data.keys():
        labj_cluster = frozenset([
            j for j in lab_data.keys() if haversine(*labpos(j), *labpos(lid)) < 40
        ])
        if labj_cluster not in clusters.values():
            nondist_clusters.add(labj_cluster)

    k = 1
    for clust in nondist_clusters:
        while k in clusters.keys():
            k += 1
        clusters[k] = list(clust)

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
        print(f"Warning: empty cluster {cluster_id=}")
        return (0., 0.)
    
    sm = reduce(lambda p1, p2: (p1[0]+p2[0], p1[1]+p2[1]), map(lambda j: lab_data[j]['pos'], cluster))
    return tuple(map(lambda x: x/(1. * len(cluster)), sm))
    
# print(lab_clusters())


# labid: pos, dist_id, is_public, capacity, backlog
lab_data = {}

# districtid: name, pos, samples
district_data = {}
    

if __name__ == '__main__':
    m = gekko.GEKKO()
    m.options.LINEAR = 1 # it is a linear program
    m.options.SOLVER = 1 # APOPT solver
    m.solver_options = APOPT_SOLVER_OPTIONS

    
    parse_inputs()
    clusters = generate_clusters()

    # setup variables
    zik = {} # Di to Ck
    tik = {} # zik > 0 ?
    xkj = {} # Ck to Lj
    yij = {} # Di to Lj if district(j) == i
    bi  = {}

    for i in district_data.keys():
        bi[i] = m.Var(lb=0, integer=True)

    for i in district_data.keys():
        for k in clusters.keys():
            zik[(i,k)] = m.Var(lb=0, integer=True)
            tik[(i,k)] = m.Var(lb=0, ub=1, integer=True)

    for k in clusters.keys():
        for j in lab_data.keys():
            xkj[(k,j)] = m.Var(lb=0, integer=True)

    for i in district_data.keys():
        for j in lab_data.keys():
            if district(j) == i:
                # Constraint 1,  yij <= 100 for district(j) == i
                yij[(i,j)] = m.Var(lb=0, ub=100, integer=True)

    # Constraints 2
    # sum_i zik = sum_j xkj
    for k in clusters.keys():
        iks = list(filter(lambda ik: ik[1] == k, zik.keys()))
        kjs = list(filter(lambda kj: kj[0] == k, xkj.keys()))
        if len(iks) > 0 or len(kjs) > 0:
            m.Equation(
                reduce(lambda x1,x2: x1+x2, map(lambda ik: zik[ik], iks), 0) ==\
                reduce(lambda x1,x2: x1+x2, map(lambda kj: xkj[kj], kjs), 0) )

    # Constraint 3
    # sum_k xkj <= CAPj - BACKj
    for j in lab_data.keys():
        m.Equation(
            reduce(lambda x1,x2: x1+x2,
                   map(lambda kj: xkj[kj],
                                filter(lambda kj: kj[1] == j,
                                       xkj.keys())), 0) \
            <= \
            (lab_data[j]['capacity'] - lab_data[j]['backlog'])
        )

    # Constraint 4
    # sum_j yij + sum_k zik + bi = Si
    for i in district_data.keys():
        # ENSURE: yij contains a (i, j) if district(j) == i
        sumyij = m.sum(list(map(lambda ij: yij[ij],
                            filter(lambda ij: ij[0] == i,
                                   yij.keys()))))
        sumzik = m.sum(list(map(lambda ik: zik[ik],
                            filter(lambda ik: ik[0] == i,
                                   zik.keys()))))
        
        m.Equation( sumyij + sumzik + bi[i] == samples(i) )


    # Constraint 5
    # zik - Mtik <= 0 if k != i
    for ik in tik.keys():
        i, k = ik
        if k != i:
            m.Equation(zik[ik] - M*tik[ik] <= 0)

    # Constraint 6
    # sum_k tik <= 1
    for i in district_data.keys():
        sumtik = m.sum(list(map(lambda ik: tik[ik],
                           filter(lambda ik: ik[0] == i,
                                  tik.keys()))))
        m.Equation(sumtik <= 1)

    print(f'{len(yij)=}\t{len(tik)=}\t{len(xkj)=}\t{len(zik)=}')
    
    # objective
    m.Obj(m.sum(list(map(lambda kj: labcost(kj[1])*xkj[kj], xkj.keys()))))
    m.Obj(10000* m.sum(list(bi.values())))
    m.Obj(5000 * m.sum(list(yij.values())))
    m.Obj(1000 * m.sum(list(map(lambda ik: tik[ik]*haversine(*cluster_centre(clusters[ik[1]], ik[1]), *district_data[ik[0]]['pos']),
                          tik.keys()))))
    m.solve(disp=True)
    print(f'Objective : {m.options.objfcnval}')
    with open('~/cni/output.pickle', 'wb') as fp:
        pickle.dump({'yij': yij, 'xkj': xkj, 'zik': zik} , fp)
    print("DONE")

    