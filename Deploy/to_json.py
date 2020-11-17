import json, os

def trans(fpath, opath):
    with open(fpath) as fgraph:
        line = fgraph.readline().split()
        nodeNum = int(line[0])
        edgeNum = int(line[1])
        nodes = [{'wgt': 1}]*nodeNum

        edges = []
        for i in range(nodeNum):
            line = fgraph.readline().split()
            for e in line:
                adj = int(e)-1
                edge = {'beg': i, 'end': adj, 'wgt': 1}
                edges.append(edge)

        input_data = {'nodes': nodes, 'edges': edges}

        with open(opath, 'w+') as fjson:
            json.dump(input_data, fjson)


inputdir = "E:/SMART_Work/GraphPartition/Deploy/graphs"
outputdir = "E:/SMART_Work/GraphPartition/Deploy/Instance/"
for parents, dirnames, filenames in os.walk(inputdir):
    for filename in filenames:
        fpath = os.path.join(parents, filename)
        opath = outputdir+filename.split('.')[0]+".json"
        trans(fpath,opath)