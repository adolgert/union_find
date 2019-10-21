'''
I'm thinking about Dave's suggestion to make a topology class. How would that work?
'''


class Storage(object):
    '''
    This stores x-y-z of vertices, or it stores values at vertices.
    There is an order to that storage. That order is a map from something
    to one dimension.
    '''
    def __init__(self):
        self._linearize=LinearTopology()
        self._array=list()
        
    def store(value,token):
        '''
        The linearization is, itself, a topology connecting vertices. Its boundary
        is the beginning and end points. The trick is that we need to invert the
        map from vertex to storage location in order to get the next vertex.
        '''
        index=self._linearize(token)
        if index>=len(self._array):
            self._array.extend([0]*(1+index-len(self._array)))
        self._array[index]=value
        
    def retrieve(self,token):
        return self._array(self._linearize(token))




class Topology(object):
    '''
    Is there a precompiled way that the topology can say what vertices are connected
    but the Storage ordering can pick which one is naturally next?
    '''
    def __init__(self,bounds):
        self._bounds=bounds
        self._cur_vertex=[0,0]
        
        
        
    def next_vertex(cur_vertex):
        '''
        What vertices are connected to this vertex?
        Of those, which do we return?
        '''
        if _programmatic:
            cur_record=self.record(cur_vertex)
            return storage.next_record(cur_record)
        else:
            for vertex in adjacency_list[cur_vertex]:
                if vertex not in seen:
                    return vertex


    def neighbor_of(vertex,which):
        '''
        Asks for a particular neighbor of this vertex.
        '''
        if which is right:
            return blah
        if which is down
            return blah




class LinearTopology(Topology):
    '''
    If every vertex is connected with a line, then there is a sense of forward and backward.
    The trick is that the vertex token has to be understood in some way. This isn't just a
    list of who follows whom in an array, although it could be.
    
    Maybe this takes a 2D grid an asks of each point who its right neighbor is, then goes 
    one down from the first in the row. It does not know what the vertex means but asks
    the 2D topology for a limited subset.
    '''
    def __init__(self,start_vertex):
        # This version works from a fully-connected 2D grid.
        self.cur_row=start_vertex    # Row walks down.
        self.cur_column=start_vertex # Column walks from each row
        
    def forward(self,vertex):
        pass
        
    def backward(self,vertex):
        pass


        

class Mesh(object):
    def __init__(self):
        cur_vertex=None
        
    def boundary(self):
        return Mesh(boundary)
    
    def bulk(self):
        return Mesh(bulk)

    def next_vertex(self):
        return toplogy.next_vertex(cur_vertex)
        
    
def traverse(mesh,visitor):
    
    for vert in mesh.natural_order():
        visitor(vert)
        