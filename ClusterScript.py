#
#  ClusterScript.py
#  
#
#  Created by Daniel Citron on 7/25/11.
#  Copyright (c) 2011 __MyCompanyName__. All rights reserved.
#

import time
import numpy, scipy
import networkx as nx
import cPickle as pickle

def read_landuse_array_data(filename, shape=None):
    data = numpy.fromfile(filename, numpy.int8)
    if shape is not None:
        data.shape = shape
    data = numpy.fliplr(numpy.transpose(data))
    return data
    
def get_use_locations(data, uses):
    loc = dict([(c, numpy.where(data==c)) for c in uses])
    return loc

def get_uses(data):
    return [c for c in range(data.max()+1) if c in data]

def start():
    data = read_landuse_array_data('34418039.dat', (3435,3625))
    uses = get_uses(data)
    loc  = get_use_locations(data, uses)
    return data, uses, loc

def zipper(loc, i):
    """
    q = zipper(loc, land_use)

    input locations data (loc), land use type (land_use)
    returns locations data in a column-stacked format

    convert loc output into useable sets of coordinates:
    """
    import numpy
    if i not in loc.keys(): return - 1
    q = numpy.column_stack((loc[i][0],loc[i][1]))
    return q

def array_shift_NN(l, land_use, debug = 0):
    """
    G = array_shift_NN(l, land_use, debug)
    l = land use locations in array format
    land_use = chosen land use type
    debug = prints debugging statements while script is running, defaults to 0
    
    Outputs collection of subgraphs with connectivity defined by nearest neighbors sharing the same land use type
    
    Script for formatting land use locations data:
    loc = land use locations data.shape = dimensions of image array data
        A = zipper(loc, land_use)
        l = numpy.zeros(data.shape)
        for i in A: l[i[0], i[1]] = 1
    """
    import numpy
    import networkx as nx
    import time
    if debug: print "Begin creating graph", str(land_use), time.asctime()
    G = nx.Graph()
    holder = numpy.where(l)
    holder = numpy.column_stack((holder[0], holder[1]))
    if debug: 
        print "Begin checking vertical neighbors", time.asctime()
        it = 0
    for i in holder: G.add_node((i[0], i[1]), use = land_use)
    # Only count overlaps between ones, count the zeros as unoccupied space
    holder = numpy.where(numpy.logical_and(l[0:-1, :] == l[1:,:], \
                                           l[0:-1, :] == 1))
    holder = numpy.column_stack((holder[0], holder[1]))
    for i in holder:
        x, y = i[0], i[1]
        G.add_edge((x, y), (x+1, y))
        if debug:
            it = it + 1
            if it % 50000 == 0: print it, "th iteration", time.asctime()
    if debug: 
        print "Begin checking horizontal neighbors", time.asctime()
        it = 0
    holder = numpy.where(numpy.logical_and(l[:, 0:-1] == l[:,1:], \
                                           l[:, 0:-1] == 1))
    holder = numpy.column_stack((holder[0], holder[1]))
    for i in holder:
        x, y = i[0], i[1]
        G.add_edge((x, y), (x, y+1))
        if debug:
            it = it + 1
            if it % 50000 == 0: print it, "th iteration", time.asctime()

    del holder  # timeout only occurs here? out of memoryfor 1M nodes? how to fix this?
    return G
    """if debug: print "Begin graph output", time.asctime()    
    graph_file_name = 'Graph_output' + str(land_use) + '.gz'
    nx.gpickle.write_gpickle(G, graph_file_name)
    if debug: print "End", time.asctime()"""

def AreaPerimeterDiversity(G, data, Ct = 15):
    """
    T = AreaPerimeter(G, data, Ct = 15)
    G = Graph of nearest-neighbors clusters
    data = land use data
    Ct = total number of different land use types (defaults to 15)
    
    Return 4-row array of Area, Perimeter, Modified Perimeter, Diversity-Modified Perimeter of each cluster
    Can plot Perimeter, modified perimeter, diversity-modified perimeter vs. clusters
    """
    import numpy
    import networkx as nx
    
    l = nx.connected_component_subgraphs(G)
    output = numpy.zeros([4,len(l)])
    index = 0
    for sg in l:
        A = len(sg)
        P = 4*A - 2*len(sg.edges())
        p = (P + 2*(A - 1))/4.
        border = {}
        for point in sg.nodes():
            """check types of all nearest neighbors"""
            if len(sg[point]) < 4:
                xcoord = point[0]
                ycoord = point[1]
                NN = []
                if ycoord < data.shape[1]-1: NN.append((xcoord, ycoord + 1)) 
                if ycoord > 0: NN.append((xcoord, ycoord - 1)) 
                if xcoord < data.shape[0]-1: NN.append((xcoord + 1, ycoord)) 
                if xcoord > 0: NN.append((xcoord - 1, ycoord))
                for nn in NN:
                    type = data[nn[0],nn[1]]
                    if not border.has_key(type): border[type] = 1
        C = len(border.keys())-1
        dP = (P + (2.*(A-1)*C)/(Ct-1))/4.
        output[:,index] = numpy.array([A, P, p, dP])
        index = index + 1
    return output
    
def ClusterSizeDistribution(G):
    """(d, H) = ClusterSizeDistribution(G)
    G = Graph of nearest-neighbors clusters
    Returns cluster size distribution: 
        H = 3-row array of bins, histogram, and cumulative histogram of G
        d = array of cluster sizes
    """
    import numpy
    import networkx as nx
    
    l = nx.connected_components(G)
    l2 = numpy.zeros(len(l))
    for i in range(len(l)):
        l2[i] = len(l[i])
    H = numpy.histogram(l2, bins = max(l2))

    CH = numpy.empty(len(H[0]))
    for i in range(0, len(H[0])):
        CH[i] = numpy.sum(H[0][0:i])
    
    return l2, numpy.row_stack((numpy.linspace(1,len(H[0]), len(H[0])), H[0], CH))

def ClusterCounting(use_type = "all", out_file_name = "CLUSTER_COUNTING_DATA.dat", debug = 1):
    """
    ClusterCounting(use_type, out_file_name, debug)
    
    use_type = land use type or "all", defaults to "all"
    out_file_name = output file name for storing results, defaults to "CLUSTER_COUNTING_DATA.dat"
    debug = determines whether debugging messages occur, defaults to 1
    Returns a dictionary with keys equal to land use types
    Value: [0] = Cluster Size Distributions Data: 
                [0]: sizes per length; [1]: bins, histogram, cumulative histogram
           [1] = Area Perimeter Diversity Data: 
                Area, Perimeter, Mod. Perimeter, Diversity-Mod. Perimeter
    """
    import time, numpy, scipy
    import networkx as nx
    import cPickle as pickle
    
    if debug: print "Open Land Use data", time.asctime()
    (data, uses, loc) = start()
    d = {}
    
    if use_type == "all":
        for land_use_type in uses:
            if debug: print "Land use type = ", land_use_type, time.asctime()
            A = zipper(loc, land_use_type)
            l = numpy.zeros(data.shape)
            for i in A: l[i[0], i[1]] = 1
            if debug: print "Start computing graph", land_use_type, time.asctime()
            G = array_shift_NN(l, land_use_type, debug = debug)
            if debug: print "Start computing cluster size distribution", land_use_type, time.asctime()
            CSD = ClusterSizeDistribution(G)
            if debug: print "Start computing cluster APD dimensions", land_use_type, time.asctime()
            APD = AreaPerimeterDiversity(G, data)
            d.update({land_use_type:[CSD, APD]})
            
    elif use_type in uses:
        if debug: print "Land use type = ", use_type, time.asctime()
        A = zipper(loc, use_type)
        l = numpy.zeros(data.shape)
        for i in A: l[i[0], i[1]] = 1
        if debug: print "Start computing graph", use_type, time.asctime()
        G = array_shift_NN(l, use_type, debug = debug)
        if debug: print "Start computing cluster size distribution", use_type, time.asctime()
        CSD = ClusterSizeDistribution(G)
        if debug: print "Start computing cluster APD dimensions", use_type, time.asctime()
        APD = AreaPerimeterDiversity(G, data)
        d.update({use_type:[CSD, APD]})
        
    else: 
        print "Invalid land use type. Choose from: ", uses, " or all"
        return -1    

    if debug: print "Open and output pickle file", time.asctime()
    f = open(out_file_name, 'wb')
    pickle.dump(d, f)
    f.close



if __name__ == '__main__':
    import sys
    try: 
        func = sys.argv[1]
    except:
        func = None
    if func:
        try: exec 'print %s' % func
        except: print '-1'
    else: print  """ClusterCounting(use_type, out_file_name, debug)\
use_type = land use type or "all", defaults to "all"\
out_file_name = output file name for storing results, defaults to "CLUSTER_COUNTING_DATA.dat"\
debug = determines whether debugging messages occur, defaults to 1\
Returns a dictionary with keys equal to land use types\
Value: [0] = Cluster Size Distributions Data: bins, histogram, cumulative histogram
       [1] = Area Perimeter Diversity Data: Area, Perimeter, Mod. Perimeter, Diversity-Mod. Perimeter
"""
