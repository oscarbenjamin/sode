#!/usr/bin/env python


from __future__ import division


import string
from itertools import islice, ifilter, permutations, count, izip

import numpy as np
from numpy.linalg import eigvals
from pygraph.algorithms.accessibility import (connected_components,
                                              accessibility,
                                              mutual_accessibility)

from dotread import GraphParser

# Generator to create non-repeating sequence of labels
def _gen_labels(n):
    """Automatically generate a sequence of non-repeating alphabetic labels
    """
    # Iterator to create infinite non-repeating sequence of string labels
    def iter_labels(_recurse=False):
        # Recursive generator
        if _recurse:
            yield ''
        for base in iter_labels(_recurse=True):
            for l in string.uppercase:
                yield base + l
    # First n terms
    return list(islice(iter_labels(), n))


class Digraph(object):
    """Digraph class represents an unlabelled directed graph"""

    #
    # Main constructor
    #

    def __init__(self, nodes=[], edges=[], nodeattrs={}, edgeattrs={}):
        """Create a new Digraph instance with the specified nodes and
        edges"""
        # Initialise empty graph
        self.empty()

        # Add all data to the graph
        for node in nodes:
            self.add_node(node)
        for edge in edges:
            self.add_edge(edge)
        for node, attrs in nodeattrs.iteritems():
            for attr in attrs:
                self.add_node_attribute(node, attr)
        for edge, attrs in edgeattrs.iteritems():
            for attr in attrs:
                self.add_edge_attribute(edges, attr)

    #
    # Interface to underlying data structures
    # this is the protected interface that could be overridden by subclasses
    #

    def empty(self):
        """Initialise empty data structures"""
        self._nodes = []
        self._edges = []
        self._nodeattrs = {}
        self._edgeattrs = {}
        self._nodes_set = set()
        self._edges_set = set()
        self.node_neighbors = {}

    def add_node(self, node):
        """Add node with name node to the graph"""
        if self.has_node(node):
            msg = "Duplicate node name '{0}'"
            raise ValueError(msg.format(node))
        self._nodes_set.add(node)
        self._nodes.append(node)
        self.node_neighbors[node] = set()

    def add_edge(self, edge):
        edge = tuple(edge)
        nfrom, nto = edge
        if self.has_edge(edge):
            msg = "Duplicate edge from '{0}' to '{1}'"
            raise ValueError(msg.format(nfrom, nto))
        self._edges_set.add(edge)
        self._edges.append(edge)
        self.node_neighbors[nfrom].add(nto)

    def has_node(self, node):
        return node in self._nodes_set

    def has_edge(self, edge):
        nfrom , nto = edge
        if not self.has_node(nfrom):
            raise ValueError("Unrecognised node '{0}'".format(nfrom))
        if not self.has_node(nto):
            raise ValueError("Unrecognised node '{0}'".format(nto))
        return tuple(edge) in self._edges_set

    def del_edge(self, edge):
        if not self.has_edge(edge):
            nfrom, nto = edge
            raise ValueError("No such edge '{0}' -> '{1}'".format(nfrom, nto))
        self._edges.pop(self._edges.index(edge))

    def iter_nodes(self):
        return iter(self._nodes)

    def iter_edges(self):
        return iter(self._edges)

    def order(self):
        return len(self._nodes)

    def nodes(self):
        return list(self.iter_nodes())

    def edges(self):
        return list(self.iter_edges())

    def neighbors(self, nfrom):
        return self.node_neighbors[nfrom]

    def node_attributes(self, node):
        return self._nodeattrs.get(node, {})

    def edge_attributes(self, edge):
        return self._edgeattrs.get(edge, {})

    def __eq__(self, other):
        if not isinstance(other, Digraph):
            return NotImplemented
        return self.equals(other)

    def __iter__(self):
        return self.iter_nodes()

    def __getitem__(self, nfrom):
        return self.neighbors(nfrom)

    def __len__(self):
        return self.order()

    def __str__(self):
        """Display in dot language format"""
        return self.to_dotstring()

    @staticmethod
    def from_dotstring(s):
        """Create digraph from string of dot language text
        
        Example:
        >>> g = Digraph.from_dotstring('digraph { A -> B }')
        """
        gp = GraphParser(text=s)
        dg = Digraph()
        for node, attr in gp.iter_nodes():
            dg.add_node(node)
        for nfrom, nto, attr in gp.iter_edges():
            dg.add_edge((nfrom, nto))
        return dg

    def to_dotstring(self):
        """Return dot language representation of this graph"""
        lines = []
        lines.append('digraph {')
        for node in self.nodes():
            lines.append('    {0}'.format(node))
        for nfrom, nto in self.edges():
            lines.append('    {0} -> {1}'.format(nfrom, nto))
        return '\n'.join(lines)

    @staticmethod
    def from_dotfile(f):
        """Create digraph from a file containing dot language text
        
        Example:
        >>> g = Digraph.from_dotfile('mygraph.dot')
        """
        if isinstance(f, str):
            f = open(f, 'r')
        return Digraph.from_dotstring(f.read())

    def to_dotfile(self, f):
        """Write dot language representation of theis graph to file f
        
        Example:
        >>> g.to_dotfile('mygraph.dot')
        """
        if isinstance(f, str):
            f = open(f, 'w')
        f.write(self.to_dotstring())

    @staticmethod
    def from_adjacency(M, nodes=None, nodes_attrs={}, edges_attrs={}):
        """Create Digraph from adjacency matrix with optional vertex labels.

        Example:
        >>> nodes = ['A', 'B']
        >>> M = [[0, 1],
        >>>      [0, 0]]
        >>> g = Digraph.from_adjacency(self, M, nodes)
        """
        # Checks and defaults
        M = np.array(M)
        N = M.shape[0]
        if nodes is None:
            nodes = _gen_labels(N)
        assert M.shape in ((0,), (N, N)) and len(nodes) == N

        # Edge key iterator
        def iter_edges():
            for n1, n2 in np.transpose(np.nonzero(M)):
                yield (nodes[n1], nodes[n2])
        edges = [(nodes[n1], nodes[n2]) for n1, n2 in np.transpose(np.nonzero(M))]
        
        # Create and return
        return Digraph(nodes, edges, nodes_attrs, edges_attrs)

    def to_adjacency(self):
        """Return adjacency matrix representation of this graph
        
        Example:
        >>> g = Digraph.from_dotstring('digraph { A -> B }')
        >>> print g.to_adjacency()
        array([[0, 1],
               [0, 0]])
        """
        N = self.order()
        M = np.zeros((N, N), dtype=int)
        for ifrom, vfrom in enumerate(self.nodes()):
            for ito, vto in enumerate(self.nodes()):
                if self.has_edge((vfrom, vto)):
                    M[ifrom, ito] = 1
        return M

    def subgraph(self, nodes):
        """Return new Digraph instance corresponding to the subgraph of nodes
        specified in indices
        
        Example:
        >>> g = Digraph.from_dotstring('digraph { A -> B -> A; A -> C }')
        >>> gAB = g.subgraph(['A', 'B'])
        >>> print g.to_dotstring()
        digraph graphname {
        a;
        b;
        a -> b
        b -> a
        }
        """
        eincluded = lambda n1, n2: n1 in nodes and n2 in nodes
        edges = [e for e in self.edges() if eincluded(*e)]
        nodeattrs = dict((n, self.node_attributes(n)) for n in nodes)
        edgeattrs = dict((e, self.edge_attributes(e)) for e in self.edges() if eincluded(*e))
        return Digraph(nodes, edges, nodeattrs, edgeattrs)

    def reverse(self):
        """Return a graph with the edges reversed"""
        nodes = self.iter_nodes()
        edges = ((nto, nfrom) for nfrom, nto in self.iter_edges())
        return Digraph(nodes, edges)

    def without_edge(self, e):
        """Return new Digraph instance representing this graph but without the
        edge e"""
        g = Digraph(self.iter_nodes(), self.iter_edges())
        g.del_edge(e)
        return g

    #
    # Convenience functions for querying graph properties
    #

    def get_indegree(self):
        """Return a vector of in degrees for the nodes in this graph"""
        M = self.to_adjacency()
        return M.sum(0)

    def get_outdegree(self):
        """Return a vector of out degrees for the nodes in this graph"""
        M = self.to_adjacency()
        return M.sum(1)

    def is_empty(self):
        """Query whther this is the empty graph"""
        return self.order() == 0

    def is_strongly_connected(self):
        """Query whether or not the graph is strongly connected"""
        iter_mas = mutual_accessibility(self).itervalues()
        return any(len(others) == self.order() for others in iter_mas)

    def is_weakly_connected(self):
        """Query whether or not the graph is weakly connected"""
        iter_as = accessibility(self).itervalues()
        return any(len(others) == self.order() for others in iter_as)

    def is_disconnected(self):
        """Query whether there are components with no edges between them"""
        return not self.is_weakly_connected()

    def equals(self, other):
        """Query whether or not the edges in this graph are the same as the
        edges in other. Note that the result depends on node order but not in
        the node names.
        
        Example:
        >>> g1 = Digraph(('a', 'b'), [('b', 'a')])
        >>> g2 = Digraph(('c', 'd'), [('d', 'c')])
        >>> g3 = Digraph(('a', 'b'), [('a', 'b')])
        >>> print g1.equals(g2), g1.equals(g3)
        True False
        """
        # Compare the edges sets of the two graphs
        return set(self.edges()) == set(other.edges())

    def edges_id(self, order=None):
        """Return an integer that uniquely identifies the adjacency matrix
        corresponding to this graph.
        
        The adjacency matrix depends on the ordering of the nodes which is not
        strictly a property of a graph. order if provided represents a
        different ordering of node keys to use
        """
        if order is None:
            order = self.nodes()
        edges_id = 0
        m = 1
        for nfrom in order:
            for nto in order:
                edges_id += m * self.has_edge((nfrom , nto))
                m *= 2
        return edges_id

    def equivalence_id(self):
        """Return a python object suitable for a dict key that uniquely
        identifies all graphs equivalent to this one (under node permutation)"""
        # min_id is the same for all graphs equivalent to this one and unique
        # among graphs of the same size (node number). Along with the node
        # number the pair is unique among graphs
        iperm = permutations(self.nodes())
        min_id = min(self.edges_id(order) for order in iperm)
        return (self.order(), min_id)

    def equivalent(self, other):
        """Determine if two graphs are topologically equivalent regardless of
        node names or node ordering. Assumes 'unlabelled'."""
        if not isinstance(other, Digraph):
            return NotImplemented
        return (self.equivalence_id() == other.equivalence_id())

    def strongly_connected_components(self):
        """Return subgraphs describing each of the strongly connected
        components from this graph"""
        # Maps each component to the list of others in the same strongly
        # connected component.
        acc = mutual_accessibility(self)
        # Identify set of strongly connected components
        sccs = set(tuple(sorted(nodes)) for nodes in acc.itervalues())
        return [self.subgraph(nodes) for nodes in sccs]

    def condensation(self):
        """Return the condensation graph and the corresponding strongly
        connected components"""
        sccs = self.strongly_connected_components()

        # Returns true if the graph has at least 1 edge between from a node in
        # sg1 to a node in sg2
        def atleast1_edge(sg1, sg2):
            for nfrom in sg1.nodes():
                for nto in sg2.nodes():
                    if self.has_edge((nfrom, nto)):
                        return True
            else:
                return False

        # Iterators for edges
        def edges_condensation():
            for n1, scc1 in enumerate(sccs):
                for n2, scc2 in enumerate(sccs):
                    if n1 != n2 and atleast1_edge(scc1, scc2):
                        yield (n1, n2)
        
        condensation = Digraph(range(len(sccs)), edges_condensation())

        return condensation, sccs

    def first_transitive_component(self):
        """
        """
        # Obtain condensation map and corresponding subgraphs
        condensation, sccs = self.condensation()
        # Identify nodes in condensation that have no incoming edges
        nodes_top = []
        for n, ie in condensation.reverse().node_neighbors.iteritems():
            if not ie:
                # Collect nodes from corresponding subgraphs
                nodes_top.extend(sccs[n].nodes())
        # Return a subgraph of all these nodes
        return self.subgraph(nodes_top)

    # find_cycles needs a unique representation to compare equivalent cycles
    def _sort_cycle(self, cycle):
        n = cycle.index(min(cycle))
        cycle = cycle[n:] + cycle[:n]
        return cycle

    def find_cycles(self, _path=None):
        """Returns a set of sequences describing all distinct cycles in
        the graph"""
        # Find cycles in the graph
        cycles = set()
        # First entry point
        if _path is None:
            for start_node in self.nodes():
                start_path = (start_node,)
                cycles.update(self.find_cycles(start_path))
        # Recursing
        else:
            for next_node in self.node_neighbors[_path[-1]]:
                # Have we found a cycle
                if next_node == _path[0]:
                    cycles.add(self._sort_cycle(_path))
                # If not repeating ourselves lets find another path
                elif next_node not in _path:
                    new_path = _path + (next_node,)
                    cycles.update(self.find_cycles(new_path))
        # Return in any case
        return cycles

    def sufficient_cycle_sets(self):
        """Find the sets of cycles whose edge union produces this graph"""
        # Obtain all unique cycles and convert to edge sets
        cycles = []
        for c in self.find_cycles():
            cycles.append(frozenset((n1, n2) for n1, n2 in zip(c, c[1:]+c[:1])))

        # Find all cycale sets that form the whole graph
        def sufficient_sets(ss, others, total):
            if set().union(*ss) == total:
                yield ss
            else:
                for n, s in enumerate(others):
                    ssnew = ss.union([s])
                    othersnew = others[n+1:]
                    for ss2 in sufficient_sets(ssnew, othersnew, total):
                        yield ss2
        edges = set(self.edges())
        scs = set()
        for n, cs in reversed(list(enumerate(cycles))):
            for cs in sufficient_sets(frozenset([cs]), cycles[n+1:], edges):
                yield cs

    def minimal_cycle_sets(self):
        # Find minimal cycle sets (must not be supersets of any other
        # sufficient cycle sets
        mcs = []
        scs = frozenset(self.sufficient_cycle_sets())
        for cs in scs:
            if sum(cs.issuperset(c) for c in scs) == 1:
                mcs.append(cs)

        return frozenset(mcs)
        
    def cycle_set_overlap(self):
        """Find overlap of minimal cycle sets"""
        Nedges = len(self.edges())
        minoverlap = Nedges
        for cs in self.sufficient_cycle_sets():
            overlap = sum(len(c) for c in cs) - Nedges
            minoverlap = min(minoverlap, overlap)
            if minoverlap == 0:
                return 0
        return minoverlap

    def overlapping_cycle_set(self):
        """Find overlap of minimal cycle sets"""
        Nedges = len(self.edges())
        minoverlap = Nedges
        mincycleset = set()
        for cs in self.sufficient_cycle_sets():
            overlap = sum(len(c) for c in cs) - Nedges
            if overlap < minoverlap:
                minoverlap = overlap
                mincycleset = cs
            # Check if we have already reached absolute minimum
            if minoverlap == 0:
                break
        return mincycleset

def test_graphs():

    # Test graph creation functions
    print "Test 1 -- creating graphs"
    vertices = ['A', 'B', 'C']
    edges = [('A', 'B'), ('B', 'C'), ('A', 'C')]
    M = [[0, 1, 1],
         [0, 0, 1],
         [0, 0, 0]]
    dotstring = """
    digraph {
        A B C
        A -> B
        A -> C
        B -> C
    }
    """
    g1 = Digraph(vertices, edges)
    g2 = Digraph.from_adjacency(M, vertices)
    g3 = Digraph.from_dotstring(dotstring)
    print 'g1', g1
    print 'g2', g2
    print 'g3', g3

    assert g1 == g2 == g3
    assert g1.equals(g2) and g2.equals(g3) and g1.equals(g3)
    assert g2.equals(g1) and g3.equals(g2) and g3.equals(g1)

    # Empty graph
    print "Test 2 -- empty graph"
    g = Digraph.from_dotstring('digraph {}')
    print 'g', g
    assert g.is_empty()

    # Non-empty graph with no edges
    print "Test 3 -- disconnected graph"
    g = Digraph.from_dotstring('digraph {A; B}')
    print 'g', g
    assert not g.is_empty()
    assert g.order() == 2
    assert g.is_disconnected()
    assert not g.is_weakly_connected()

    print "Test 4 -- 2 node weakly connected graph"
    g = Digraph.from_dotstring('digraph {A -> B}')
    print g
    assert not g.is_empty()
    assert g.order() == 2
    assert g.is_weakly_connected()
    assert not g.is_strongly_connected()
    assert not g.is_disconnected()

    print "Test 5 -- 2 node strongly connected graph"
    g = Digraph.from_dotstring('digraph {A -> B -> A}')
    print g
    assert not g.is_empty()
    assert g.order() == 2
    assert g.is_weakly_connected()
    assert g.is_strongly_connected()
    assert not g.is_disconnected()

    print "Test 5 -- 3 node graphs equivalent under permutations"
    g1 = Digraph.from_dotstring('digraph {A -> B -> C -> A}')
    g2 = Digraph.from_dotstring('digraph {A -> C -> B -> A}')
    print g1
    print g2
    assert g1.equivalent(g2)
    assert g2.equivalent(g1)
    assert not g1.equals(g2)
    assert not g2.equals(g1)

    print "Test 6 -- subgraph"
    g1 = Digraph.from_dotstring('digraph {A -> B -> C -> A}')
    g2 = Digraph.from_dotstring('digraph { C -> A }')
    nodes = ['C', 'A'] # edges2 from edges1
    print g1
    print g2
    assert g1.subgraph(nodes).equals(g2)

    print "Test 6 -- remove edge"
    d1 = """
    digraph {
        A -> B -> C -> A
        A -> C -> B -> A
    }
    """
    d2 = """
    digraph {
        A -> B;   C -> A
        A -> C -> B -> A
    }
    """
    edge = ('B', 'C')
    g1 = Digraph.from_dotstring(d1)
    g2 = Digraph.from_dotstring(d2)
    print g1
    print g2
    assert g1.without_edge(edge).equals(g2)

    print "Test 7 -- first transitive component simple"
    g1 = Digraph.from_dotstring('digraph { A -> B }')
    g2 = Digraph.from_dotstring('digraph { A }')
    g1_first = g1.first_transitive_component()
    print g1
    print g2
    print g1_first
    assert g1_first.equals(g2)

    print "Test 8 -- condensation graph 4 nodes"
    d = """digraph {
        A B C D
        A -> B
        A -> D
        B -> A
        B -> D
        C -> A
        C -> D
    }"""
    g = Digraph.from_dotstring(d)
    g1 = Digraph.from_dotstring('digraph { C }')
    g2 = Digraph.from_dotstring('digraph { A; B; A -> B -> A }')
    g3 = Digraph.from_dotstring('digraph { D }')
    gc = Digraph.from_dotstring('digraph { 0; 1; 2; 0 -> 1; 0 -> 2; 1 -> 2 }')
    print g
    print g1
    print g2
    print g3
    gcs, (g1s, g2s, g3s) = g.condensation()
    assert g1s.equals(g1) and g2s.equals(g2) and g3s.equals(g3)
    assert gcs.equivalent(gc)

    print "Test 9 -- top level component"
    d = """digraph {
        A B C D
        A -> B
        A -> D
        B -> A
        B -> D
        C -> D
    }"""
    d1st = """digraph {
        A B C
        A -> B
        B -> A
    }"""
    g = Digraph.from_dotstring(d)
    g1st = Digraph.from_dotstring(d1st)
    print g
    print g1st
    print g.first_transitive_component()
    assert g1st.equals(g.first_transitive_component())

    print "Test 10 -- find cycles in graph"
    d = """digraph {
        A B C D
        A -> B -> A
        B -> C -> B
        A -> C -> D -> A
    }"""
    g = Digraph.from_dotstring(d)
    cycle_set = set([
        ('A', 'B'),
        ('B', 'C'),
        ('A', 'C', 'D'),
        ('A', 'C', 'B'),
        ('A', 'B', 'C', 'D'),
    ])
    minimal_cycle_set = set([
        frozenset([frozenset([('A', 'B'), ('B', 'A')]),
                   frozenset([('B', 'C'), ('C', 'B')]),
                   frozenset([('A', 'C'), ('C', 'D'), ('D', 'A')])
                  ]),
        frozenset([frozenset([('A', 'C'), ('C', 'B'), ('B', 'A')]),
                   frozenset([('A', 'B'), ('B', 'C'), ('C', 'D'), ('D', 'A')])
                  ]),
    ])
    print g
    print g.find_cycles()
    print g.minimal_cycle_sets()
    print g.cycle_set_overlap()
    assert g.find_cycles() == cycle_set
    assert g.minimal_cycle_sets() == minimal_cycle_set
    assert g.cycle_set_overlap() == 0

    print "Test 11 -- find cycle overlap in graph"
    d = """digraph {
        A B C D
        A -> B -> C -> D -> A
        A -> C -> A
        B -> D -> B
    }
    """
    g = Digraph.from_dotstring(d)
    print g.minimal_cycle_sets()
    print g.cycle_set_overlap()
    assert g.cycle_set_overlap() == 0

    print "Test 12 find nonzero cycleset overlap"
    d = """digraph {
        A B C D
        A -> B -> C
        A -> C -> A
        B -> D -> B
        A -> D -> A
    }"""
    g = Digraph.from_dotstring(d)
    for cs in g.minimal_cycle_sets():
        print list(cs)
    print g.cycle_set_overlap()
    assert g.cycle_set_overlap() == 1

    print "Test 13 more cycleset overlap"
    d = """digraph {
        A B C D
        A -> B -> C -> A
        B -> D -> A
    }"""
    g = Digraph.from_dotstring(d)
    for cs in g.minimal_cycle_sets():
        print list(cs)
    assert g.cycle_set_overlap() == 1

    print "Test 14 more cycle set overlap"
    d = """digraph {
        A B C D
        D -> A -> B -> C
        A -> C -> A
        B -> D -> B
    }"""
    g = Digraph.from_dotstring(d)
    for cs in g.minimal_cycle_sets():
        print list(cs)
    assert g.cycle_set_overlap() == 3

    print "Test 15 more cycle set overlap"
    d = """digraph {
        A B C D
        A -> B -> C -> A
        B -> D -> A
    }"""
    g = Digraph.from_dotstring(d)
    for cs in g.minimal_cycle_sets():
        print list(cs)
    assert g.cycle_set_overlap() == 1

    print "Test 16 bib cycle set overlap"
    d = """digraph {
        A B C D
        A -> {B C D}
        B -> {A C D}
        C -> {A B D}
        D -> {A B}
    }"""
    g = Digraph.from_dotstring(d)
    print len(g.find_cycles())
    print g.cycle_set_overlap()
    #assert g.cycle_set_overlap() == 0

    print "All tests passed"

if __name__ == "__main__":
    test_graphs()
