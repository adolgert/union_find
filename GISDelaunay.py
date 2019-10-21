#
#  GISDelaunay.py
#  
#
#  Created by Daniel Citron on 8/4/11.
#  Copyright (c) 2011 __MyCompanyName__. All rights reserved.
#
def CorrLength(loc, land_use, binsize = 1.):
    """
    (H, L) = CorrLength(loc, land_use),
    
    Returns a histogram of nearest neighbor distances
    
    loc = locations data; land_use = land use type (loc key); binsize = size of histogram binning;
    H = histogram of distances away from nearest neighbors, as determined by the delaunay triangulation
    L = array of distances between pairs of points between triangles, can be used to quickly recalculate the histogram
    """
    import math, numpy, scipy
    from matplotlib._delaunay import delaunay
    
    Q = delaunay(loc[land_use][0], loc[land_use][1])
    #Delaunay triangulation edges = Q[1]
    nedges = len(Q[1])
    lengths = numpy.empty(nedges)
    for i in range(0, nedges):
        point1 = Q[1][i][0]
        point2 = Q[1][i][1]
        x1 = loc[land_use][0][point1]
        x2 = loc[land_use][0][point2]
        y1 = loc[land_use][1][point1]
        y2 = loc[land_use][1][point2]
        distance = math.sqrt((x1 - x2)**2 + (y1 - y2)**2)
        lengths[i] = distance
        
    histo = numpy.histogram(lengths, bins = (numpy.ceil(max(lengths)/binsize))-1, \
                            range = (1,numpy.ceil(max(lengths))))
    return histo, lengths