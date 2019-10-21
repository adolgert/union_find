import h5py
import sys
import datetime
import unittest
import logging
from argparse import ArgumentParser


logger=logging.getLogger('timls')



def pretty_columns(data,columns,aliases=None,below=None):
    '''
    The data are a list of dictionaries.
    The columns are in those dictionaries.
    aliases are alternate names for those columns.
    below are things to show indented below the line.
    '''
    rows=list()
    if aliases:
        assert(len(aliases)==len(columns))
        rows.append(aliases)
    else:
        rows.append(columns)

    for record in data:
        cols=[record[c] for c in columns]
        rows.append([str(x).strip() for x in cols])
    logger.debug('found %d rows' % len(rows))

    sizes=[0]*len(columns)
    for r in rows:
        for i,col in enumerate(r):
            sizes[i]=max(sizes[i],len(col))
    logger.debug('sizes: %s' % str(sizes))
    start_fmt='%%%ds '*len(columns)
    logger.debug('format: %s' % str(start_fmt))
    fmt=start_fmt % tuple(sizes)
    for i,pr in enumerate(rows):
        print(fmt % tuple(pr))
        if i>0:
            for b in below:
                print('\t%s' % data[i-1][b])



def optimization(build_group):
    cfg_file=build_group.attrs['Raster stats configuration']
    lines=cfg_file.split('\t')
    for l in lines:
        if l.startswith('optimization'):
            return l.split(' = ')[1]
    return ''



def runs_from_file(name):
    f=h5py.File(name)
    if not 'timing' in f.keys():
        logger.warn('There are no timings in %s' % name)
        next
        
    runs=list()
    for run_id in f['/timing'].keys():
        meta_data=dict()
        meta_data.update(f['/timing'][run_id].attrs)
        meta_data['run_id']=run_id
        meta_data['start time']=datetime.datetime.strptime(
            meta_data['start time'],'%Y-%b-%d %H:%M:%S')
        meta_data['measurements']=', '.join(f['/timing'][run_id].keys())
        build=f[f['/timing'][run_id].attrs['build']]
        meta_data['opt']=optimization(build)
        runs.append(meta_data)

    runs.sort(lambda a,b: a['start time']<b['start time'])
    for i,r in enumerate(runs):
        r['idx']=i

    return runs



def extract_data(filename,run_idx,timing_name):
    runs=runs_from_file(filename)
    run=runs[run_idx]

    f=h5py.File(filename)
    data=f['/timing/%s/%s' % (run['run_id'],timing_name)]
    return data



def print_data(data):
    for r in data:
        print( '%d %d' % (r[0], r[1]) )


def list_files(files):
    for name in files:
        runs=runs_from_file(name)

        pretty_columns(runs,['idx','size','depth','count','iter','start time',
                             'opt'],
                       ['idx','size','depth','count','iter','start','opt'],
                       ['measurements',])




def test_parse_time():
    a='2012-03-20T09:13:25.167365'
    datetime.datetime.strptime(a,'%Y-%m-%dT%H:%M:%S.%f')



def suite():
    suite=unittest.TestSuite()
    loader=unittest.TestLoader()
    suite.addTests(loader.loadTestsFromTestCase(test_parse_time))
    return suite



if __name__ == '__main__':
    parser = ArgumentParser(description='look inside hdf files')
    parser.add_argument('data_files',metavar='file',type=str,
                        nargs='+',help='A file to read')
    extract_group=parser.add_argument_group('extraction',
           'These arguments extract data instead of listing contents.')
    extract_group.add_argument('--idx', dest='idx', type=int, default=None,
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
        data=extract_data(args.data_files[0],args.idx,args.name)
        print_data(data)
    else:
        list_files(args.data_files)
