#!/bin/python3


import pandas as pd
from functools import reduce
from math import sin, cos, pi, asin, acos, sqrt, atan2
import sys

outputfolder = 'analysis'
soln_file = 'output.csv'
lab_file = f'./test1/lab_test_data.csv'
district_file = f'./test1/district_test_data.csv'


if len(sys.argv) > 1:
	soln_file = sys.argv[1]

print("Parsing {}".format(soln_file))

df = pd.read_csv(soln_file)

labs = {}
district = {}

print(f'writing to {outputfolder}/')

def haversine(latlon1, latlon2):
	lat1, lon1 = latlon1
	lat2, lon2 = latlon2
	p = pi/180
	a = 0.5 - cos((lat2-lat1)*p)/2 + cos(lat1*p) * cos(lat2*p) * (1-cos((lon2-lon1)*p))/2
	return 12742 * asin(sqrt(a)) #2*R*asin...

for distid in district.keys():
	labpos = [labs[lid]['pos'] for lid in labs if labs[lid]['dist_id'] == distid]
	if len(labpos) > 0:
		lab_center = reduce(lambda p, q: (p[0]+q[0], p[1]+q[1]), labpos)
		lab_center = (lab_center[0]/len(labpos), lab_center[1]/len(labpos))
		district[distid]['lab_center'] = lab_center
		
def dist(dist_id, lab_id):
	i = dist_id
	j = lab_id
	if labs[lab_id]['dist_id'] == dist_id:
		return 0
	m = haversine(district[dist_id]['pos'],
	     labs[lab_id]['pos'])
	return m




with open(lab_file) as labfp:
	for line in labfp.readlines()[1:]:
		s = line.split(',')
		labid = int(s[1])
		lat = float(s[2])
		lon = float(s[3])
		dist_id = int(s[4])
		is_public = s[5] == '0'   # 0 is public
		cap = int(s[6])
		back = int(s[7])
		labs[labid] = {
			'pos': (lat, lon), 'dist_id': dist_id,
			'is_public': is_public,
			'cap':cap, 'back':back}

with open(district_file) as distfp:
	for line in distfp.readlines()[1:]:
		s = line.split(',')
		distid = int(s[1])
		name = s[2]
		lat = float(s[3])
		lon = float(s[4])
		sample = int(s[5])
		district[distid] = {
			'name': name,
			'pos': (lat, lon),
			'sample': sample
		}

lab_ids  = list(labs.keys())
dist_ids = list(district.keys())




# per lab stat
# incoming districts
for j in lab_ids:
	labdf = df.loc[(df.transfer_type == 0) & (df.destination == j)]
	labfp = open(f'{outputfolder}/lab{j}.txt', 'w')
	distid = labs[j]['dist_id']
	distname = district[distid]['name']
	capacity = labs[j]['cap']
	back = labs[j]['back']
	totalsamples = sum(labdf.samples_transferred) 
	print(f'Lab ID: {j}\tDistrict: {distid} - {distname}', file=labfp)
	print(f'capacity: {capacity-back} = {capacity}-{back}', file=labfp)
	print(f'total samples collected = {totalsamples}', file=labfp)
	print(f'Overload = {totalsamples-(capacity - back)}\n', file=labfp)
	
	print('DistID\tDistName\tSamples', file=labfp)
	print('------------------------------------------', file=labfp)

	for ind in labdf.index:
		src = labdf['source'][ind]
		if df['samples_transferred'][ind] > 0:
			print('\t'.join(map(lambda x: str(x),
				[src,  district[src]['name'], 
				df['samples_transferred'][ind]])), file=labfp)



# districts 
for i in dist_ids:
	distdf = df.loc[(df.transfer_type == 0) & (df.source == i)]
	backlog = df.loc[(df.transfer_type == 1) &(df.source == i)]
	if (len(backlog) == 1):
		backlog = sum(backlog.samples_transferred)
	elif (len(backlog) != 0):
		raise Exception("something wrong with data")
	else:
		backlog = 0
	
	distfp = open(f'{outputfolder}/district{i}.txt', 'w')

	distname = district[i]['name']
	cases = district[i]['sample']
	backlog = cases - sum(distdf.samples_transferred)
	print(f"District ID: {i}\t Name: {distname}", file=distfp)
	print(f"Cases: {cases}\t Backlog: {backlog}", file=distfp)

	print('LabID\tSamples\tDistID\tDistName', file=distfp)
	print('-----------------------------------------', file=distfp)
	for ind in distdf.index:
		labid = df['destination'][ind]
		labdist = labs[labid]['dist_id']
		labdistname = district[labdist]['name']
		samples_trans = df['samples_transferred'][ind]
		
		assert(samples_trans >= 0)
		
		if dist(i, labid) >= 40:
			#assert(samples_trans == 0)
			pass
	

		if df['samples_transferred'][ind] > 0:
			print('\t'.join(map(lambda x: str(x),	
				[labid, df['samples_transferred'][ind],
				 labdist, labdistname])), file=distfp)


# cost calculation
"""
for i in dist_ids:
	for j in lab_ids:
		if dist(i, j) > 0 and dist(i, j) < 40:
			lpos = labs[j]['pos']
			dpos = district[i]['pos']
			print('\t'.join(map(lambda x: str(x)+" ", 
				[(i, j), dist(i, j),
				 (180/pi)*atan2(lpos[0]-dpos[0], lpos[1]-dpos[1])])))

"""

