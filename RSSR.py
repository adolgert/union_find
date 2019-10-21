#
#  RSSR.py
#  
#
#  Created by Daniel Citron on 7/20/11.
#  Copyright (c) 2011 __MyCompanyName__. All rights reserved.
#

def RSSR(x, y, show_plot = 0, style = 'b,'):
    """
    a = RSSR(x, y, show_plot = 0, style = 'b,')


    Method for finding regions of boxcounting and other data that can best be described with a linear fit-
        Chooses a "window" of data points and performs a linear fit on those points,
        Repeats linear fits for many windows of data
        Can now compare Coefficient of Correlation (R) and Sum of Squared Residuals (SSR) across the many possible sets of data points
        As long as both the SSR is at a minimum and R is at a maximum, the linear fit is valid
    returns data formatted to a 7-column array: 
        Starting index, number of points, regression slope, regression intercept, R coeff, SSR, Mean SSR
    
    Starting index = # of the data point at which the window begins
    Number of points = # of data points included in the fit, after (and including) the starting index
    Regression slope = slope of the linear fit for the starting index and number of points
    Regression intercept = intercept of the linear fit for the starting index and number of points
    R coeff = Coefficient of Correlation of the fit
    SSR = Sum of Squared Residuals of the fit
    Mean SSR = Sum of Squared Residuals/Number of points
    
    """
    from scipy import linspace, polyval, stats
    from pylab import plot, title, show, figure, legend
    import numpy
    npoints = len(x)
    if len(x) != len(y): 
        print "X and Y must have same length"
        return -1
    
    arraysize = (npoints - 4)*(npoints - 3)/2
    starts = numpy.zeros(arraysize)
    widths = numpy.zeros(arraysize)
    ms  = numpy.zeros(arraysize)
    bs  = numpy.zeros(arraysize)
    rs  = numpy.zeros(arraysize)
    ssrs  = numpy.zeros(arraysize)
    mssrs  = numpy.zeros(arraysize)
    index = 0
    
    for i in range(0, npoints+1 - 5):
        for j in range(5, npoints+1 - i):
            xsubset = x[i:i+j]
            ysubset = y[i:i+j]
            
            (m, b, R, pv, stderr) = stats.linregress(xsubset, ysubset)
            yr = polyval([m, b], xsubset)
            SSR = sum((yr - ysubset)**2)
            
            """print xsubset, ysubset, yr
            print stats.linregress(xsubset, ysubset)
            print SSR, 1.*SSR/j
            return"""
            
            starts[index] = i
            widths[index] = j
            ms[index] = m
            bs[index] = b
            rs[index] = R + 1
            ssrs[index] = SSR
            mssrs[index] = 1.*SSR/j
            
            index = index + 1
    
    if show_plot:
        figure()
        title("SSR vs R")
        plot(rs, ssrs, style)
        show()
        
        figure()
        title("MSSR vs R")
        plot(rs, mssrs, style)
        show()
    
    data = numpy.transpose(numpy.column_stack((starts, widths, ms, bs, rs, ssrs, mssrs)))
    print "Data Columns: Starting index, number of points, regression slope,"
    print "regression intercept, R coeff, SSR, Mean SSR"
    
    return data