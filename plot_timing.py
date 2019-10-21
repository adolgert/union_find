'''
This reads hdf files with timing information and plots them.
'''

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import unittest
import logging
from argparse import ArgumentParser

import timls


logger=logging.getLogger('plot_timing')



def plot_one(name,vals):
    fig = plt.figure()
    ax=fig.add_subplot(111)
    n,bins,patches=ax.hist(vals,100,facecolor='green',alpha=0.75)
    bincenters=0.5*(bins[1:]+bins[:-1])
    ax.set_xlabel('Time [ns]')
    ax.set_ylabel('Count')
    plt.show()



def plot_all(data,minmax=None):
    fig = plt.figure()
    plot_idx=1
    if not minmax:
        for vals in data.values():
            if not minmax:
                minmax=[min(vals),max(vals)]
            else:
                minmax[0]=min(minmax[0],min(vals))
                minmax[1]=max(minmax[1],max(vals))

    colors=['green','blue','red','gray']*4
    for name,vals in data.items():
        ax=fig.add_subplot(len(data)*100+10+plot_idx)
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




def suite():
    suite=unittest.TestSuite()
    loader=unittest.TestLoader()
    #suite.addTests(loader.loadTestsFromTestCase(test_parse_time))
    return suite



if __name__ == '__main__':
    parser = ArgumentParser(description='look inside hdf files')
    parser.add_argument('data_files',metavar='file',type=str,
                        nargs='+',help='A file to read')
    extract_group=parser.add_argument_group('extraction',
           'These arguments extract data instead of listing contents.')
    extract_group.add_argument('--idx', dest='idx', type=str, default=None,
                        help='index of dataset to extract')
    extract_group.add_argument('--name', dest='name', type=str, default=None,
                        help='name of timing run to extract')
    parser.add_argument('--log', dest='log', type=str, default='info',
               help='Logging level: TRACE, DEBUG, INFO, WARN, ERROR')
    parser.add_argument('--test', dest='test', action='store_true',
               help='Whether to run unit tests on this file.')
    args = parser.parse_args()

    allowed_debug=['debug', 'info', 'warn', 'error','fatal','critical']
    if not args.log in allowed_debug:
        print 'Debugging level should be one of %s' % (str(allowed_debug),)

    logging.basicConfig(level=getattr(logging,args.log.upper()))

    if args.test:
        unittest.TextTestRunner(verbosity=2).run(suite())
        sys.exit(0)

    if not args.data_files:
        logging.error("There are no hdf files to open on the command line.")
        sys.exit(1)

    if args.idx is not None:
        indices=[int(x) for x in args.idx.split(',')]
        if len(indices)==1:
            data=timls.extract_data(args.data_files[0],indices[0],args.name)
            plot_one(args.name,data[:,1])
        else:
            all_data=dict()
            for idx in indices:
                data=timls.extract_data(args.data_files[0],idx,args.name)
                all_data['%s %d' % (args.name,idx)]=data[:,1]
            plot_all(all_data)
