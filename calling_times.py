'''
I'm playing with data from 
https://forge.cornell.edu/svn/repos/afidd/projects/raster_measure/branches/cpp:137. main.cpp, which compiles into the executable raster, shows the differences
in using boost::bind versus std::functions versus regular functors. This file makes graphs of the results of running this program.
'''
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab


def plot_one(name,vals):
    fig = plt.figure()
    ax=fig.add_subplot(111)
    n,bins,patches=ax.hist(vals,100,facecolor='green',alpha=0.75)
    bincenters=0.5*(bins[1:]+bins[:-1])
    ax.set_xlabel('Time [ns]')
    ax.set_ylabel('Count')
    plt.show()


def plot_all(data,minmax=[5500,6500]):
    fig = plt.figure()
    plot_idx=1
    colors=['green','blue','red','gray']
    for name,vals in data.items():
        ax=fig.add_subplot(410+plot_idx)
        vals=vals[np.where(vals>minmax[0])]
        vals=vals[np.where(vals<minmax[1])]
        n,bins,patches=ax.hist(vals,100,facecolor=colors[plot_idx-1],alpha=0.75)
        ax.set_xlabel('Time [ns]')
        ax.set_ylabel('Count')
        ax.set_title(name)
        ax.set_xlim(minmax[0],minmax[1])
        plot_idx+=1

    bincenters=0.5*(bins[1:]+bins[:-1])
    plt.show()



def plot_together(data):
    minmax=[69000,80000]
    fig = plt.figure()
    plot_idx=1
    ax=fig.add_subplot(111)
    colors=['green','blue','red','gray']
    for name,vals in data.items():
        vals=vals[np.where(vals>minmax[0])]
        vals=vals[np.where(vals<minmax[1])]
        n,bins,patches=ax.hist(vals,100,facecolor=colors[plot_idx-1],alpha=0.75)
        ax.set_xlabel('Time [ns]')
        ax.set_ylabel('Count')
        ax.set_title(name)
        ax.set_xlim(minmax[0],minmax[1])
        plot_idx+=1

    bincenters=0.5*(bins[1:]+bins[:-1])
    plt.show()


if __name__ == '__main__':
    raw_data=[x.strip() for x in open(sys.argv[1],'r').readlines()]
    data=dict()
    for l in raw_data:
        splitted=l.split()
        data[splitted[0]]=np.array([int(x) for x in splitted[1:]])

    
