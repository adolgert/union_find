#!/usr/local/bin/python
# -*- coding: utf-8 -*-
'''
This is an implementation of quadtrees following Aizawa's three papers:

[1] K. Aizawa, S. Tanaka, K. Motomura, and R. Kadowaki, "Algorithms for
connected component labeling based on quadtrees," International
Journal of Imaging Systems and Technology, vol. 19, no. 2, pp. 158-166,
Jun. 2009.

[2] ﻿Aizawa, K., & Tanaka, S. (2009). A constant-time algorithm for finding
neighbors in quadtrees. IEEE transactions on pattern analysis and machine
intelligence, 31(7), 1178-83. doi:10.1109/TPAMI.2008.145

[3] ﻿Aizawa, K., Motomura, K., Kimura, S., Kadowaki, R., & Fan, J. (2008).
Constant time neighbor finding in quadtrees: An experimental result. 2008
3rd International Symposium on Communications, Control and Signal Processing,
(March), 505-510. Ieee. doi:10.1109/ISCCSP.2008.4537278

'''
import math
import numpy as np
import logging
import unittest
from itertools import *



logger = logging.getLogger('quadtree')

# n - Location codes in binary form.
# dn - \delta n in the article, are basic direction increments.


def dimensions_to_levels(shape):
    '''
    Shape is a numpy array shape.
    A 4x4 matrix has depth 2.
    A 5x5 matrix kicks up to depth 3.
    '''
    assert( isinstance(shape,(list,tuple)) )
    return int(math.ceil(math.log(max(shape),2)))



# There is a y bit, then x bit for each level, with the most
# significant bit being the root node. Total bit count is
# 2*(# levels). kind u itemsize
def type_from_bits(bit_cnt):
    byte_cnt=( (bit_cnt-1)//8 ) + 1
    logger.debug('byte_cnt %d' % byte_cnt)
    for atype in [np.uint8, np.uint16, np.uint32, np.uint64]:
        type_instance = atype()
        if type_instance.nbytes >= byte_cnt:
            return atype
    return None




#   2 3
#   0 1
# We need to record, for each quad,
# (location code, level, color, neighbor differences).
def morton_to_location_code(morton):
    '''
    A morton code is a sequence of integers from 0 to 3.
    The location code is bit-packed.
    '''
    return sum([x[1]<<(2*x[0]) for x in enumerate(reversed(morton))])


def location_code_to_morton(location,level_cnt):
    '''
    Given 0b111000, level 3, give (3,2,0).
    '''
    morton=list()
    while level_cnt:
        morton.insert(0,location & 3)
        location >>= 2
        level_cnt-=1
    return morton



def location_code_to_xy(n,r):
    '''
    Given 0b111000, give coords of bottom left corner in raster.
    '''
    x=sum([(n & (1<<(2*i)))  >> i for i in range(r)])
    y=sum([(n & (1<<(2*i+1)))>>(i+1) for i in range(r)])
    logger.debug('loc %s xy %s' % (bin(n),str((x,y))))
    return (x,y)


def location_code_to_range(n,l,r):
    '''
    Given location code and level, return the range in the
    raster that corresponds to that level.
    '''
    llx,lly = location_code_to_xy(n,r)
    quadsize=(1<<(r-l))
    return (llx,lly,llx+quadsize,lly+quadsize)


def t_xy(level_cnt):
    '''
    The paper defines a tx and ty, used as a bitmask to pull out x and y.
    t_x is 010101 and t_y is 101010.
    '''
    t_x=0
    for x in range(level_cnt):
        t_x<<=2
        t_x|=1
    t_y = (t_x<<1)
    
    logger.debug('%s' % str((t_x,t_y)))
    return (t_x,t_y)


def dn8(level_cnt):
    '''
    This creates \Delta n which is an increment that gets you from
    one quad to its neighbors.
    '''
    (t_x,t_y) = t_xy(level_cnt)
    east,north,west,south = (1,2,t_x,t_y)
    return (east, east|north, north, north|west, west, west|south, south, south|east)



def dn4(level_cnt):
    '''
    This creates \Delta n which is an increment that gets you from
    one quad to its neighbors.
    '''
    (t_x,t_y) = t_xy(level_cnt)
    east,north,west,south = (1,2,t_x,t_y)
    return (east, north, west, south)



def child_location(n, l, r, morton_order):
    '''
    Return the location code of the child in a particular direction.
    n - location code of parent
    l - level of parent
    r - fixed resolution
    morton_order - 0,1,2,3 in morton z-order
    '''
    return (n | (morton_order << (2*(r-(l+1)))))



def location_addition(n, dn, t_x, t_y):
    '''
    This is \oplus from the paper. It gives you the neighboring quad.
    This is called Schrack's algorithm and calculates the
    location code only of equal-sized neighbors.
    It is mq = nq \oplus \Delta nq.
    
    n location code of quad
    \Delta n, directional increment
    t_x, t_y, helper masks in x and y
    '''
    left =(( (n|t_y) + (dn&t_x) ) & t_x)
    right=(( (n|t_x) + (dn&t_y) ) & t_y)
    res=left|right
    return res



def neighbor_equal_size(n, l, r, dn, t_x, t_y):
    '''
    l, level of quad
    r, fixed resolution of tree
    '''
    dn <<= ( 2*(r-l) )
    res=location_addition(n, dn, t_x, t_y)
    logger.debug('input %s %s' %  (bin(n),location_code_to_morton(n,r)))
    logger.debug('result %s %s' % (bin(res),location_code_to_morton(res,r)))
    return res
    


undef,black,white,gray = [np.uint8(x) for x in range(4)]

class quad_entry:
    '''
    Each entry in the linear quadtree is one of these.
    '''
    def __init__(self,n,l,v,ld):
        '''
        n is location code. l is level in the quadtree hierarchy.
        v is the value (black, white, gray)
        ld is the length four array of level differences with neighbors.
        '''
        self.n=n
        self.l=l
        self.v=v
        self.ld=ld
    def __str__(self):
        return '(%s, %d, %d, %s)' % (str(location_code_to_morton(self.n,3)),
                       self.l, self.v, str(self.ld))
    def __repr__(self):
        return 'quad_entry%s' % self.__str__()



class SubLabel:
    '''
    Connected components use labels. The label linked list has a different
    head. These are the bulk of the list. This is a singly-linked list
    until the last link. last.next is None, and last.label points to a Label.
    Label.label may point to another label.
    '''
    def __init__(self):
        '''
        From outside the linked list, you don\'t call this constructor.
        Use make_from() as a factory pattern.
        '''
        self._label=None
        self._next=None

    def label(self):
        '''
        This retrieves the label for a given SubLabel.
        This is where you would implement compression.
        '''
        next=self
        while next._next:
            next=next._next

        while next._label:
            next=next._label

        return next


    def assign(self,other):
        '''
        When you assign one label to another, you are joining lists.
        The other may be a SubLabel or a Label.
        '''
        slabel=self.label()
        olabel=other.label()

        if slabel != olabel:
            slabel._label=olabel


    def make_from(self):
        '''
        Create a new label from this one, which means making a new
        SubLabel, linking it into the chain, and returning it.
        The Label always points to a last node, so insert between
        them.
        '''
        if self._next:
            sl=SubLabel()
            sl._next=self._next
            self._next=sl
            return sl

        # If I don't have a next, then I'm at the end.
        else:
            # Might as well add beneath the primary label.
            l=self.label()
            # If the primary label already has sublabels.
            if l._next:
                sl=SubLabel()
                sl._next=l._next
                l._next=sl
                return sl
            # If the primary label is alone.
            else:
                sl=SubLabel()
                sl._label=l
                self._last=sl
                self._next=sl
                return sl


    def __cmp__(self,other):
        slabel=self.label()
        olabel=other.label()
    
        return cmp(slabel._idx,olabel._idx)

    def __eq__(self,other):
        return self.__cmp__(other) is 0

    def __neq__(self,other):
        return self.__cmp__(other) is not 0




class Label(SubLabel):
    '''
    For connected components, there are labels for each component.
    In the paper, the Label has a First and Last pointer. This 
    class uses Next and Last so that it can be a subclass.
    '''
    def __init__(self,idx):
        SubLabel.__init__(self)
        self._idx=idx
        self._last=self


class QuadTree:
    '''
    This quadtree is from Aizawa. It's a linear black and white
    quadtree that includes information about neighbors of a quad.
    '''
    def __init__(self,bh,start_size=100):
        level_cnt=dimensions_to_levels(bh.shape)
        # Total number of levels
        self.r = level_cnt
        # Nodes are stored with smaller ones first
        self.begin=start_size
        self.end=start_size
        # no neighbor needs a value
        self.no_neighbor = 127
        (self.n,self.l,self.v,self.ld) = self.create_storage(start_size)
        
        self.quads=self.construct_from_image(bh)


    def add_node(self,n,l,v,ld):
        if self.begin==0:
            self.resize(2*self.cnt)
            
        self.begin -= 1
        self.n[self.begin]=n
        self.l[self.begin]=l
        self.v[self.begin]=v
        self.ld[self.begin]=ld


    def create_storage(self,start_size):
        # location codes
        n = np.zeros(start_size,dtype=type_from_bits(self.r*2))
        # levels
        l = np.zeros(start_size,dtype=type_from_bits(2))
        # values
        v = np.zeros(start_size,dtype=np.uint8)
        # level differences with neighbors
        ld = np.zeros((start_size,4),dtype=np.int8)
        return (n,l,v,ld)


    def resize(self,new_size):
        (n,l,v,ld) = self.create_storage(new_size)
        shift=new_size-self.begin
        
        n[shift:]=self.n[:]
        l[shift:]=self.l[:]
        v[shift:]=self.v[:]
        ld[shift:]=self.ld[:,:]
        self.n=n
        self.l=l
        self.v=v
        self.ld=ld
        self.begin=shift
        self.end=new_size


    def construct_from_image(self,bh):
        '''
        The paper calls this QuadtreeWithLCodeLevelDifferences in Fig. 6.

        Creating the list of entries is very dynamic and does not
        lend itself to storage in an array.
        Within this call, use Python lists. Then turn it into numpy
        arrays at the end of the call.
        bh - array of black and white values.
        '''
        quads=list()
        quads.append(quad_entry(0, 0, gray, [self.no_neighbor]*4))
        
        self.delta_n = dn4(self.r)
        (self.t_x,self.t_y) = t_xy(self.r)
        
        # while the node list includes gray areas
        first_gray=[0,quads[0]]
        while first_gray:
            
            # for each equal size neighbor of first gray area do
            logger.debug('current grays %s' % str(first_gray))
            (first_gray_idx, first_gray_area) = first_gray
            equal_size_neighbor=self.existing_equals(quads,first_gray_area)

            for neighbor, direction in equal_size_neighbor:
                # Corresponding level differences of neighbors are added by 1
                neighbor.ld[(direction+2)%4] += 1
                
            # First Gray quadrant is replaced by its four children.
            # Delete first gray and add its four children to the queue
            del quads[first_gray_idx]

            children=list()
            for child_idx in range(4):
                level_diffs=[0]*4
                ps=[2,3,1,0] # index of side where we copy from parent
                for ld_idx in [ps[child_idx],(ps[child_idx]+1)%4]:
                    level_diffs[ld_idx]=first_gray_area.ld[ld_idx]
                    if level_diffs[ld_idx]!= self.no_neighbor:
                        level_diffs[ld_idx]-=1
                child_n=child_location(first_gray_area.n,first_gray_area.l,
                                       self.r, child_idx)
                child_l=first_gray_area.l+1
                ext=location_code_to_range(child_n,child_l,self.r)
                if np.all(bh[ext[0]:ext[2],ext[1]:ext[3]]):
                    color=black
                elif np.any(bh[ext[0]:ext[2],ext[1]:ext[3]]):
                    color=gray
                else:
                    color=white
                child=quad_entry(child_n,
                                 child_l,
                                 color,
                                 level_diffs)
                logger.debug('made child %s' % str(child))
                quads.insert(0, child)
                children.append(child)

            # for each equal size neighbor of each child (other than their brothers)
            for child in children:
                child_neighbors=self.existing_equals(quads[4:],child)
                for neighbor, direction in child_neighbors:
                    # Corresponding level differences of neighbors are added by 1
                    neighbor.ld[(direction+2)%4] += 1

            grays_iter = ifilter(lambda x: x[1].v==gray,enumerate(quads))
            try:
                first_gray=grays_iter.next()
            except StopIteration:
                first_gray=None

        quads.sort(lambda a,b: cmp(a.n,b.n))
        return quads



    def existing_equals(self,quads,match):
        '''
        Given a list of quads, find all equal-sized neighbors of a
        location, those that actually exist.
        quads - a list of quad_entry.
        match - a quad entry around which to searh
        returns - a list of neighbors that exist, and their directions
        '''
        logger.debug('existing equals %s' % match)
        equal_size_neighbor=dict()
        for idx, dn in enumerate(self.delta_n):
            if match.ld[idx] != self.no_neighbor:
                loc=neighbor_equal_size(match.n, match.l, self.r, dn,
                                      self.t_x, self.t_y)
                equal_size_neighbor[loc]=idx

        # Want to run through all quads once to find matches
        is_neighbor=lambda y: equal_size_neighbor.has_key(y.n) and y.l == match.l
        # and then work back to getting the direction ids.
        existing_neighbors=ifilter(is_neighbor, quads)
        return [(x, equal_size_neighbor[x.n]) for x in existing_neighbors]


    def morton(self,morton):
        '''
        Retrieve the quadtree node corresponding to this morton code.
        '''
        n=morton_to_location_code(morton)
        return self.location(n)
        
        
    def location(self,n):
        '''
        Retrieve the quadtree node corresponding to this location.
        '''
        found_iter=ifilter(lambda x: x.n==n, self.quads)
        try:
            found=found_iter.next()
        except StopIteration:
            return None
        return found



    def neighbor_q(self,center,direction):
        '''
        The neighbor Q in a given direction of a given quadrant P is the
        smallest quadrant (it may be GRAY) that is adjacent to P in a given
        direction and is of size greater than or equal to the quadrant.
        
        This method finds the location code of the neighbor of the same size,
        but the linear quadtree representation only stores nodes that are black
        or white. If the smallest node of the same size is gray, then it isn't
        in the list of quads. This algorithm just returns the child of that node,
        which is wrong. Another option would be to create a gray node at the
        correct level.
        '''
        dd=center.ld[direction]
        if dd == self.no_neighbor:
            return None
            
        dn_start=self.delta_n[direction]
        if dd < 0:
            rld=2*(self.r-center.l-dd)
            nq = (center.n >> rld) << rld
            dn = (dn_start << rld)
            mq = location_addition(nq, dn, self.t_x, self.t_y)
        else:
            nq = center.n
            dn = (dn_start << (2*(self.r-center.l)))
            mq = location_addition(nq, dn, self.t_x, self.t_y)

        found = self.location(mq)
        if not found:
            logger.error('neighbor_q should always find a neighbor.')
        # It is possible the found quad is a level smaller than the
        # original, so we return None when there is no quad at the same
        # or larger level. This departs from the paper's implementation.
        if found.l>center.l:
            logging.debug("There is no quad the same level or smaller."+
                          "center %d found %d" % (center.l,found.l))
            return None

        # The paper also returns (location code, level) instead of
        # returning the quad.
        return found


    def connected_components(self):
        '''
        An algorithm for connected component labeling of linear
        time complexity. Its space complexity is also linear.
        '''

        label=dict()
        maximum_label=0

        # Scan QTLCLD from smaller quadrants to larger quadrents.
        # Don't need to sort them because we filter through them later.
        self.quads.sort(lambda x,y: cmp(x.l,y.l) or cmp(x.n,y.n))
        
        for level in range(self.r,0,-1):

            # For each cell c_i scanned in Step 1
            for quad in ifilter(lambda x: x.l==level, self.quads):
                
                # If a label is defined for cell c_i
                min_label=None
                if label.has_key(quad.n):
                    # Then scan all cells in neighbors of c_i.
                    for neighbor in neighbor_iter(quad,self):
                        if label.has_key(neighbor.n):
                            if not min_label:
                                min_label=label[neighbor.n]
                            else:
                                min_label=min(min_label,label[neighbor.n])
                
                    label[quad.n].assign(min_label)
                    # After that, scan again and set all undefined c_j
                    # and set to the same label.
                    for neighbor in neighbor_iter(quad,self):
                        if not label.has_key(neighbor.n):
                            label[neighbor.n]=min_label.make_from()

                # If no label is defined, then scan all neighboring cells
                else:
                    
                    colors=[]
                    labels=[]
                    for neighbor in neighbor_iter(quad,self):
                        colors.append(neighbor.v)
                        if label.has_key(neighbor.n):
                            labels.append(label[neighbor.n])
                    
                    # If all neighbors are white
                    if all([x is white for x in colors]):
                        label[quad]=Label(maximum_label)
                        maximum_label+=1
                    
                    # If cells are black and none have labels
                    elif any([x is black for x in colors]) and not labels:
                        # set label(c_i) to max label.
                        label[quad]=Label(maximum_label)
                        maximum_label+=1
                        for n in neighbor_iter(quad,self):
                            # modifying algorithm to set label only for
                            # black values.
                            if n.v is black and not label.has_key(n.n):
                                label[n]=label[quad].make_from()

                    # else some black cells have labels
                    else:
                        # Set the center to the minimum label.
                        min_label=min(labels)
                        label[quad]=min_label.make_from()
                        # Set any unlabeled neighbors to the min label.
                        for n in neighbor_iter(quad,self):
                            if n.v is black and not label.has_key(n.n):
                                label[n]=min_label.make_from()
    
        clusters=dict()
        for q in self.quads:
            l=label[q]
            top=l.label()._idx
            if not clusters.has_key(top):
                clusters[top]=list()
            clusters[top].append(q)

        return clusters



def neighbor_iter(quad,tree):
    idx=0
    while idx<4:
        neighbor=tree.neighbor_q(quad,idx)
        if neighbor:
            yield neighbor
        idx+=1




def test_labels():
    a=Label(0)
    b=a.make_from()
    c=b.make_from()
    d=a.make_from()
    assert(a.label()==b.label())
    assert(d.label()==b.label())
    e=Label(1)
    f=e.make_from()
    assert(f.label()!=d.label())
    g=e.make_from()
    f.assign(c)
    assert(g.label()==b.label())

    


def quadtree_setUp():
    '''
    Using the sample image from the paper.
    '''
    bh=np.array(
        [[1,1,1,1,1,0,0,0],
         [1,1,1,1,1,0,0,0],
         [1,1,1,1,1,1,0,0],
         [1,1,1,1,1,1,0,0],
         [0,0,0,0,1,1,1,1],
         [0,0,0,0,1,1,1,1],
         [0,0,0,0,0,0,0,0],
         [0,0,0,0,0,0,0,0]],
        np.uint8)
    # The way I entered it is flipped from paper.
    bh=bh.transpose()[:,::-1]
    assert(not np.any(bh[0:4,0:4]))
    assert(np.all(bh[0:4,4:8]))
    assert(not np.any(bh[4:8,0:2]))
    assert(np.all(bh[2:4,4:8]))
    assert(np.all(bh[4:6,4:6]))
    assert(not np.any(bh[6:8,4:6]))
    assert(np.all(bh[4:5,6:8]))
    assert(not np.any(bh[5:8,6:8]))
    return QuadTree(bh)
        
        
def test_quadtree():
    paper=[
        [[0,0,0],1,white,[1,0,127,127]],
        [[1,0,0],2,white,[0,0,-1,127]],
        [[1,1,0],2,white,[127,0,0,127]],
        [[1,2,0],2,black,[0,0,-1,0]],
        [[1,3,0],2,black,[127,0,0,0]],
        [[2,0,0],1,black,[1,127,127,0]],
        [[3,0,0],2,black,[0,1,-1,0]],
        [[3,1,0],2,white,[127,0,0,0]],
        [[3,2,0],3,black,[0,0,-2,-1]],
        [[3,2,1],3,white,[-1,0,0,-1]],
        [[3,2,2],3,black,[0,127,-2,0]],
        [[3,2,3],3,white,[-1,127,0,0]],
        [[3,3,0],2,white,[127,127,1,0]]
    ]
    
    qt=quadtree_setUp()

    assert(len(qt.quads)==len(paper))

    for q,p in zip(qt.quads,paper):
        logger.debug("%s %s" % (q, p))
        assert(q.n==morton_to_location_code(p[0]))
        assert(q.l==p[1])
        assert(q.v==p[2])
        assert(q.ld==p[3])



def test_connected():
    qt=quadtree_setUp()
    clusters=qt.connected_components()
    for c,v in clusters.iteritems():
        print '==== %d', c
        for q in v:
            print q





def test_dimensions_to_levels():
    assert(dimensions_to_levels((2,2))==1)
    assert(dimensions_to_levels((4,4))==2)
    assert(dimensions_to_levels((5,5))==3)
    assert(dimensions_to_levels((5,3))==3)
    



def test_type_from_bits():
    assert(type_from_bits(1)==np.uint8)
    assert(type_from_bits(8)==np.uint8)
    assert(type_from_bits(9)==np.uint16)
    assert(type_from_bits(65)==None)



def test_morton_to_location_code():
    assert(morton_to_location_code((3,))==0b11)
    assert(morton_to_location_code((3,2,0))==0b111000)
    assert(morton_to_location_code((3,2))==0b1110)


def test_location_code_to_morton():
    assert(location_code_to_morton(0b11,3)==[0,0,3])
    assert(location_code_to_morton(0b1110,3)==[0,3,2])
    assert(location_code_to_morton(0b111000,3)==[3,2,0])


def test_location_code_to_xy():
    assert(location_code_to_xy(0b00,1)==(0,0))
    assert(location_code_to_xy(0b10,1)==(0,1))
    assert(location_code_to_xy(0b01,1)==(1,0))
    assert(location_code_to_xy(0b111000,3)==(4,6))


def test_location_code_to_range():
    assert(location_code_to_range(0b000000,1,3)==(0,0,4,4))
    assert(location_code_to_range(0b010000,1,3)==(4,0,8,4))

    

def test_t_xy():
    assert( t_xy(1) == (0b01,0b10) )
    assert( t_xy(2) == (0b0101,0b1010) )
    assert( t_xy(3) == (0b010101,0b101010) )


def test_dn():
    assert(dn8(3)==(0b1,0b11,0b10,0b010111,0b010101,0b111111,0b101010,0b101011))



def test_child_location():
    assert(child_location(0b000,0,3,0)==0b000000)
    assert(child_location(0b000,0,3,1)==0b010000)
    assert(child_location(0b000,0,3,2)==0b100000)
    assert(child_location(0b100000,1,3,1)==0b100100)



def test_location_addition():
    r=3
    (t_x,t_y)=t_xy(r)
    dn=dn4(r)

    c=[morton_to_location_code(x) for x in ((3,2,0),(3,2,1),(0,0,0),(2,0,0))]
    assert(neighbor_equal_size(c[0],3,r,dn[0],t_x,t_y)== c[1])
    assert(neighbor_equal_size(c[1],3,r,dn[2],t_x,t_y)== c[0])
    assert(neighbor_equal_size(c[2],1,r,dn[1],t_x,t_y)== c[3])
    assert(neighbor_equal_size(c[3],1,r,dn[3],t_x,t_y)== c[2])



# Default order is Morton, or z-order
def test_neighbor_q():
    qt=quadtree_setUp()
    (center, west, north) = [qt.morton(x) for x in [[3,0,0], [2,0,0], [3,2,0]]]
    # directions are 0 east 1 north 2 west 3 south.
    assert(qt.neighbor_q(center,2)==west)
    # The north side is rather wrong b/c the [3,2,0] that's there is smaller than
    # [3,0,0]. The real north neighbor is [3,2,0] a level up so that it is gray.
    assert(qt.neighbor_q(center,1)==None) # Was north.
        
        


def suite():
    suite=unittest.TestSuite()
    suite.addTest(unittest.FunctionTestCase(test_dimensions_to_levels))
    suite.addTest(unittest.FunctionTestCase(test_morton_to_location_code))
    suite.addTest(unittest.FunctionTestCase(test_location_code_to_morton))
    suite.addTest(unittest.FunctionTestCase(test_location_code_to_xy))
    suite.addTest(unittest.FunctionTestCase(test_child_location))
    suite.addTest(unittest.FunctionTestCase(test_t_xy))
    suite.addTest(unittest.FunctionTestCase(test_dn))
    suite.addTest(unittest.FunctionTestCase(test_quadtree))
    suite.addTest(unittest.FunctionTestCase(test_connected))
    suite.addTest(unittest.FunctionTestCase(test_neighbor_q))
    suite.addTest(unittest.FunctionTestCase(test_labels))
    return suite



if __name__ == '__main__':

    logging.basicConfig(level=getattr(logging,'DEBUG'))

    unittest.TextTestRunner(verbosity=2).run(suite())


    


