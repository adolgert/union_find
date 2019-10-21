#
#  DTC_fractals.py
#  
#
#  Created by Daniel Citron on 9/11/11.
#  Copyright (c) 2011 __MyCompanyName__. All rights reserved.
#


import numpy, math


# functions for reading land use array data
def read_landuse_array_data(filename, shape=None):
    data = numpy.fromfile(filename, numpy.int8)
    if shape is not None:
        data.shape = shape
    data = numpy.fliplr(numpy.transpose(data))
    return data

def get_uses(data):
    return [c for c in range(data.max()+1) if c in data]

def get_use_locations(data, uses):
    loc = dict([(c, numpy.where(data==c)) for c in uses])
    return loc


# start by importing from GIS.uses the functions for importing and sorting landuse data
# set equal to (data, uses, loc)
def start():
    """
    (data, uses, loc) = start()
    
    reads in land use array data; returns data = land use data in array form,
    uses = land use data labels
    loc = dictionary with land use types as keys and land use locations as values
    """
    from uses import read_landuse_array_data
    from uses import get_uses, get_use_locations

    data = read_landuse_array_data('34418039.dat', (3435,3625))
    uses = get_uses(data)
    loc  = get_use_locations(data, uses)
    return data, uses, loc
    

# Updated full range for Histogram calculations of boxcounting data
# epsilon is box size, eg 2*n pixels...
# boxes are rectangular with aspect ratio equal to that of the original picture
def capacity1(data, loc, usetype, epsilon):
    """
    (boxsize, counts) = capacity1(data, loc, usetype, epsilon)
    
    data = land use data; loc = land use locations; usetype = land use type;
    epsilon = boxsize
    Returns box counting box size and # filled boxes
    
    performs boxcounting using numpy's 2d histogram function, counting the number of filled boxes
    epsilon is box size, eg 2**n pixels...
    boxes are rectangular with aspect ration equal to that of the original picture
    """
    s = data.shape
    nx = int(s[0]//epsilon)
    ny = int(s[1]//epsilon)
    histo = numpy.histogram2d(loc[usetype][0], loc[usetype][1], \
            [nx, ny], range = [[0, s[0]],[0, s[1]]])
    counts = numpy.sum(histo[0] > 0)
    # returns actual histogram box size, instead of "epsilon"
    return 1.*s[0]/nx, counts


def capacities_output(data, loc, epsilons, file_out_name = "BOXCOUNTING_OUTPUT", \
                    debug = 0):
    """
    Perform capacity1 boxcounting algorithm for all land use types
    
    output = capacities_output(data, loc, epsilons, file_out_name, debug)
    input land use array data (data), locations data (loc), list of boxcounting grid sizes (epsilons)
        output file name (file_out_name, defaults to "BOXCOUNTING_OUTPUT"), debugging (debug, defaults to 0 ie no debugging)
    returns a dictionary of boxcounting data for bins equal to epsilons, for each land use type
    
    creates a text file and a pickle file of the boxcounting data - the pickle file may be accessed using:
        f = open("BOXCOUNTING_OUTPUT_PICKLE.dat", 'rb')
        d = pickle.load(f)
    
    setting debug = 1 prints out the times at which each stage of the code starts and ends, in case the code gets hung up somewhere
    """
    import time, numpy, math
    import cPickle as pickle
    f = open(file_out_name+".dat", 'wb')
    f2 = open(file_out_name+"_PICKLE.dat", 'wb')
    use_types = loc.keys()
    npoints = len(epsilons)
    output = {}
    for usetype in use_types:
        counts = numpy.zeros([4, npoints])
        epsilon_index = 0
        
        f.writelines("USE TYPE = " + str(usetype) + '\n')
        if debug: print "Start ", usetype, time.asctime()
        for epsilon in epsilons:
            (boxsize, count) = capacity1(data, loc, usetype, epsilon)
            if count == 0: break
            f.writelines(repr(boxsize).rjust(20) + '' + repr(count).rjust(8)+ \
                        '' + repr(math.log10(boxsize)).rjust(20) + \
                        '' + repr(math.log10(count)).rjust(20) + '\n')
            counts[:, epsilon_index] = boxsize, count, math.log10(boxsize), math.log10(count)
            
            epsilon_index = epsilon_index + 1
        f.writelines('\n')
        if debug: print "End ", usetype, time.asctime()

        output.update({usetype: counts})
        
    f.close()
    pickle.dump(output, f2)
    f2.close()

    return output
    # plotting :
    # for i in loc.keys(): pyplot.plot(d[i][2], d[i][3], '.-')


# series of mass exponents:
def series_mass_exponent(data, loc, usetype, epsilons, qs):
    """
    (Boxsizes, Taus) = series_mass_exponents(data, loc, usetype, epsilons, qs)
    
    data = land use data
    loc = land use locations
    usetype = land use type
    epsilons = array of boxcounting sample sizes
    qs = array of mass exponents q
    Returns boxsizes = boxcounting sampling size; taus = array of counts vs. boxsize for all q's
    
    Uses mass_exponent() function to perform series of mass exponents procedure for a chosen land use type
    """
    import numpy
    
    taus = numpy.empty((len(epsilons), len(qs)))
    boxsizes = numpy.empty(len(epsilons))
    
    for i in range(0,len(epsilons)):
        epsilon = epsilons[i]
        (boxsize, counts) = mass_exponent(data, loc, usetype, epsilon, qs)
        taus[i][:] = counts
        boxsizes[i] = boxsize
    taus = taus.T
    return boxsizes, taus
    # plotting all data for each q:
    # for i in range(len(qs)):
    #   matplotlib.pyplot.loglog(epsilons, taus[i], '.-')
    # legend(qs)

# mass dimension calculations
def mass_exponent(data, loc, usetype, epsilon, qs):
    """
    Boxsize, counts = mass_exponent(data, loc, usetype, epsilon, qs)
    
    data = land use data;
    loc = land use locations
    usetype = land use type
    epsilon = boxcounting sample size
    qs = array of mass exponents q
    Returns boxsize, measure of the set using mass exponent q
    """
    import numpy
    
    s = data.shape
    nx = int(s[0]//epsilon)
    ny = int(s[1]//epsilon)
    histo = numpy.histogram2d(loc[usetype][0], loc[usetype][1], \
            [nx, ny], range = [[0, s[0]],[0, s[1]]])
    
    counts = numpy.empty(len(qs))
    for i in range(0, len(qs)):
        q = qs[i]
        if q == 0: 
            counts[i] = numpy.sum(histo[0] > 0)
        elif q < 0:
            sum = 0
            for j in numpy.hstack(histo[0]):
                if j > 0: sum = sum + j**q
            counts[i] = sum
        elif q > 0: counts[i] = numpy.sum(histo[0]**q)
        
    return 1.*s[0]/nx, counts


# UNDER CONSTRUCTION:

# Calculating Shannon Entropy of the picture at sampling squares of particular length scales
def sentropy(data, loc, usetype, epsilon):
    """
    (boxsize, counts) = sentropy(data, loc, usetype, epsilon)\
    data = land use data; loc = land use locations; usetype = land use type\;
    epsilon = boxsize
    Returns (actual) box counting box size and Shannon entropy
    
    note: Shannon entropy does not appear to scale according to a power law
    """
    import math, numpy
    s = data.shape
    nx = int(s[0]//epsilon)
    ny = int(s[1]//epsilon)
    histo = numpy.histogram2d(loc[usetype][0], loc[usetype][1], \
            [nx, ny], range = [[0, s[0]],[0, s[1]]])
    cell = 1.*s[0]*s[1]/nx/ny
    histo = histo[0]/cell # turn each histogram count into a frequency
    f = histo[numpy.where(histo > 0)] # pulls out all nonzero values from histogram
    entropy = numpy.sum(f*numpy.log(f))
    # returns actual histogram box size, instead of "epsilon"
    return math.sqrt(cell), entropy

    
    
    
# perform box counting for a range of box sizes
# input vertically-stacked xy data of one data type, number of box divisions
# output two lists, one of division size, one of the number of full boxes counted
def boxcounting(loc, use_type, range = 10, xmax = 3625., ymax = 3435.): 
    """
    (divisions, num) = boxcounting(loc, use_type, range = 10, xmax = 3625., ymax = 3435.)
    
    input locations data (loc), land use type (use_type), number of boxcounting grid sizes (range)
    output grid sizes as fraction of total image size (divisions), number of occupied boxes (num)
    
    for each of the box sizes, use fractaldim() to calculate the number of occupied boxes
    """
    divisions = []
    num = []
    for level in xrange(1, range):
        boxsize = int(2 ** level)
        a = fractaldim(loc, use_type, boxsize, xmax = xmax, ymax = ymax)
        divisions.append(1./boxsize)
        num.append(a[0])
        print level, a
    return divisions, num


# An alternative and slower method for box counting:

# measure fractal dimension using box counting
# input loc data and land use, number of box divisions
# return number of boxes counted, log of boxes counted, # of boxes dividing plane, 
# estimation of fractal dimension by dividing log(number of full boxes)/log(divisions)
def fractaldim(loc, use_type, boxlevel, xmax = 3625., ymax = 3435.):
    """
    (num, log(num), log(boxlevel), log(num)/log(boxlevel)) = fractaldim(loc, use_type, boxlevel, xmax = 3625., ymax = 3435.)
    
    input locations data (loc), size of box (boxlevel)
    output number of occupied sites (num), log(num), log(boxlevel), log(num)/log(boxlevel)
    (the last output is an approximation of the fractal dimension for a single data point)
    
    performs box counting by looping over all points in the set, assigning each to a box, and counting the number of full boxes at the end
    """
    import math
    if boxlevel <= 1.0: return -1
        
    pointlist = zipper(loc, use_type)
        
    pointdict = {}
    for point in pointlist:
        box = (int(point[0] * boxlevel/xmax), 
            int(point[1] * boxlevel/ymax))
        if not pointdict.has_key(box): pointdict[box] = 1

    num = len(pointdict.keys())

    if num == 0: return -1

    return num, math.log(num), math.log(boxlevel), \
    (math.log(num)/math.log(boxlevel))
    

def plot_boxcounting_data(loc, index, data_shape = (3435,3625), range = 15, style = 'b.'):
    """
    plot_boxcounting_data(loc, land_use,  data_shape = (3435,3625), range = 15, style = 'b.')
    
    input locations data (loc), land use type (land_use), shape of data (data_shape), 
        number of box counting grid sizes (range), plotting style (style)

    plots on a log-log scale the number of full boxes vs box size for a particular land use type
    """
    import numpy
    from matplotlib import pylab as plt
    if index not in loc.keys():
        print "Choose valid land use type", loc.keys()
        return -1
    a = loc[index][0] # x coordinates
    b = loc[index][1] # y coordinates
    c = numpy.column_stack([a,b])
    c = 1.*c/data_shape # normalize data to a box of dimensions 1
    (level, num) = boxcounting(c, range) # compute box counting
    plt.figure()
    plt.plot(numpy.log10(level), numpy.log10(num), style)
    
    
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
    

    