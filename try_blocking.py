import sys
import numpy as np
import unittest
import logging
from argparse import ArgumentParser


logger=logging.getLogger('try_blocking')



def ntoij(n,w,b):
    wb=(w-1)/b + 1
    lb=n/(b*b) # linear index of block
    li=n%(b*b) # linear index within block
    ib=lb/wb   # block's row
    jb=lb%wb   # block's column
    i=li/b     # within block's row
    j=li%b     # within block's column
    #print wb, lb, li, ib, jb, i, j
    return (ib*b+i, jb*b + j)

def ijton(ij,w,b):
    wb=(w-1)/b + 1
    ib=ij[0]/b
    jb=ij[1]/b
    i=ij[0]%b
    j=ij[1]%b
    lb=ib*wb+jb
    li=i*b+j
    n=lb*b*b + li
    return n


def block():
    indices=dict()
    w=8
    for i in range(w):
        for j in range(w):
            n=i*w+j
            res=ntoij(n,w,2),
            indices[res[0]]=n
            print res[0],
        print

    print '----'

    for i in range(w):
        for j in range(w):
            print '%2d' % (indices[(i,j)],),
        print

    print '----'

    for i in range(w):
        for j in range(w):
            print '%2d' % ijton((i,j),w,2),
        print


    if len(set(indices))!=len(indices):
        print 'some indices are repeated'

    
def bit_add(x,y):
    b=np.uint64(1)
    #logger.debug('%s %s' % (bin(x),bin(y)))
    while (x):
        a=x&y
        b=x^y
        x=a<<np.uint64(1)
        #logger.debug('%s '*4 % tuple([bin(z) for z in [a,b,x,y]]))
        y=b
    return y


one=np.uint64(1)
two=np.uint64(2)
xmask=np.uint64(0)
ymask=np.uint64(0)
for i in np.uint64(range(32)):
    xmask|= (one<<(two*i))
    ymask|= (one<<(two*i+one))


def x_add(n,dx):
    b=one
    n=np.uint64(n)
    dx=np.uint64(dx)
    logger.debug('%s %s' % (bin(n),bin(dx)))

    x=n&xmask
    ty=n&ymask
    logger.debug('n: %s\nx: %s\ntx: %s' % (n,x,ty))

    while (x):
        a=x&dx
        b=x^dx
        x=a<<two
        logger.debug('%s '*4 % tuple([bin(z) for z in [a,b,x,dx]]))
        dx=b

    return ty|dx



def y_add(n,dx):
    b=one
    n=np.uint64(n)
    dx=np.uint64(dx)
    logger.debug('%s %s' % (bin(n),bin(dx)))

    x=n&ymask
    ty=n&xmask
    logger.debug('n: %s\nx: %s\ntx: %s' % (n,x,ty))

    while (x):
        a=x&dx
        b=x^dx
        x=a<<two
        logger.debug('%s '*4 % tuple([bin(z) for z in [a,b,x,dx]]))
        dx=b

    return ty|dx



def any_add(n,dx):
    b=one
    n=np.uint64(n)
    dx=np.uint64(dx)
    logger.debug('%s %s' % (bin(n),bin(dx)))

    x=n
    logger.debug('n: %s\nx: %s' % (n,x))

    while (x):
        a=x&dx
        b=x^dx
        x=a<<two
        logger.debug('%s '*4 % tuple([bin(z) for z in [a,b,x,dx]]))
        dx=b

    return dx



def interleave(x,y):
    x=np.uint64(x)
    y=np.uint64(y)
    z=np.uint64(0)
    one=np.uint64(1)
    for i in np.uint64(range(32)):
        z|= ((x& (one<<i)) << i)
        z|= ((y& (one<<i)) << (i+one))
    return z



def detangle(n):
    n=np.uint64(n)
    x=np.uint64(0)
    y=np.uint64(0)
    one=np.uint64(1)
    two=np.uint64(2)
    for i in np.uint64(range(32)):
        x|= (n&(one<<(two*i))) >> i
        y|= (n&(one<<(two*i+one))) >> (i+one)
    return (x,y)




class TestBit(unittest.TestCase):
    def test_interleave(self):
        self.assertEqual(interleave(1,1),np.uint64(3))
        self.assertEqual(interleave(1,0b10),np.uint64(9))
        self.assertEqual(interleave(0b10,1),np.uint64(6))
    
    def test_bit(self):
        val=[(2,2),(1,3),(1,5),(17,24),(5,-1),(1,-1)]
        val.extend([(np.random.randint(0,32768),np.random.randint(0,32768))
                    for x in range(100)])
        val=[np.uint64(x) for x in val]
        for x,y in val:
            self.assertEqual(x+y,bit_add(x,y))

    def test_detangle(self):
        val=[np.uint64(x) \
                 for x in [(2,2),(1,3),(1,5),(17,24),(5,20)]]
        for x,y in val:
            self.assertEqual((x,y),detangle(interleave(x,y)))


    def test_shift(self):
        val=[np.uint64(x) \
                 for x in [(2,2),(1,3),(1,5),(17,24),(5,20)]]
        p,m=np.uint64((1,-1))
        for x,y in val:
            self.assertEqual((x+p,y),detangle(x_add(interleave(x,y),p)))
            self.assertEqual((x+m,y),detangle(x_add(interleave(x,y),xmask)))

        p,m=np.uint64((1,-1))
        for x,y in val:
            self.assertEqual((x,y+p),detangle(y_add(interleave(x,y),0b10)))
            self.assertEqual((x,y+m),detangle(y_add(interleave(x,y),ymask)))

        p,m=np.uint64((1,-1))
        for x,y in val:
            self.assertEqual((x+p,y),detangle(any_add(interleave(x,y),p)))
            self.assertEqual((x+m,y),detangle(any_add(interleave(x,y),xmask)))
        for x,y in val:
            self.assertEqual((x,y+p),detangle(any_add(interleave(x,y),0b10)))
            self.assertEqual((x,y+m),detangle(any_add(interleave(x,y),ymask)))
        


def suite():
    suite=unittest.TestSuite()
    loader=unittest.TestLoader()
    suite.addTests(loader.loadTestsFromTestCase(TestBit))
    return suite


if __name__ == '__main__':

    parser = ArgumentParser(description='play')

    parser.add_argument('--bit', dest='bit', action='store_true',
                       default=False, help='test bit shifting')
    parser.add_argument('--block', dest='block', action='store_true',
                       default=False, help='try blocking')

    parser.add_argument('--log', dest='log', type=str, default='debug',
               help='Logging level: debug, info, warn, error, fatal, critical')
    parser.add_argument('--test', dest='test', action='store_true',
               help='Whether to run unit tests on this file.')
    args = parser.parse_args()

    allowed_debug=['debug', 'info', 'warn', 'error','fatal','critical']
    if not args.log in allowed_debug:
        print 'Debugging level should be one of %s' % (str(allowed_debug),)

    logging.basicConfig(level=getattr(logging,args.log.upper()))

    did_one=False
    if args.test:
        unittest.TextTestRunner(verbosity=2).run(suite())
        did_one=True
        sys.exit(0)

    if args.bit:
        did_one=True
        test_bit()

    if args.block:
        did_one=True
        block()

    if not did_one:
        parser.print_help()
