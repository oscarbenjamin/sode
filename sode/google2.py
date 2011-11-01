#!/usr/bin/env python

# standard library
import sys
import optparse

# extension
import numpy

def page_rank(M, R=None, d=1.00, conserve=True):
    """Iterate the pagerank approximation R based on adjacency matrix M

    R represents the initial page ranks and defaults to all ones.
    
    M is the matrix siich that M[i, j] represents a connection FROM unit i TO
    unit j. The connections are binary 0 or 1 and not in any way graded. For
    example, if the graph is 1 => 2 => 3, then the corresponding adjacency
    matrix is given by
    M = [[ 0, 1, 0],
         [ 0, 0, 1],
         [ 0, 0, 0]]

    note that the page_rank here is increased as a result of outward arrows
    unlike in the web context where inward links increase page_rank

    This function should be used iteratively to obtain the desired result
    >>> R1 = numpy.ones(N)/N
    >>> R = page_rank(M, R1)
    >>> print R
    """
    # Check inputs and determine N
    if len(M.shape) != 2 or M.shape[0] != M.shape[1]:
        raise ValueError("M should be 2D and square")
    if R is None:
        R = numpy.ones(M.shape[0])
        R /= R.sum()
    if len(R.shape) != 1 or R.shape[0] != M.shape[0]:
        raise ValueError("M should be 1D")
    N = M.shape[0]

    # transpose reverses graph in pagerank calculation
    M = numpy.array(M, float)
    if conserve:
        M += numpy.eye(N)

    # Normalise
    for n in range(N):
        col = M[:, n]
        if col.sum():
            col /= col.sum()
            pass
        else:
            col = numpy.zeros_like(col)
            col[n] = 1
        M[:, n] = col

    #print M

    # Uniform vector
    C = numpy.ones(N, float)
    C *= R.sum() / C.sum()

    # loop until done
    while True:
        # Compute next iterate
        R_undamped = numpy.dot(M, R)
        R_new = d * R_undamped + (1 - d) * C
        if any(numpy.isnan(R_new)):
            print R_new
        # Check for convergence
        if numpy.allclose(R_new, R):
            R[R<1e-6] = 0
            return R
        # Repeat
        #print R
        R = R_new

def page_rank(M, R=None, d=0.99, conserve=True):
    N = M.shape[0]
    M = numpy.array(M, float)
    R = numpy.ones((N, 1))
    R_old = -1 * numpy.ones((N, 1))
    count = 0
    while any(abs(R - R_old) > 1e-6):
        count += 1
        R_old, R = R, M.dot(R)
        if R.sum():
            R /= R.sum()
        R = d * M.dot(R) + (1 - d)
        if count > 100 or any(numpy.isnan(R)):
            print M
            return R
    return R


def parse_input(lines):
    """Create adjacency matrix from 'dot'-style graph syntax

    e.g.
    >>> print parse_input(["A => B => C", "B=>A"]):

    returns (M, indices) where M is the adjacency matrix and indices is the dic
    t relating the node names to the rows/columns of the adjacency matrix
    """

    # Split lines into seperate a => b statements
    statements = []
    for l in lines:
        units = [u.strip() for u in l.split('->')]
        if len(units) < 2:
            raise ValueError("Cannot parse line '{0}'".format(l))
        statements.extend(zip(units[:-1], units[1:]))

    # Split statements with + e.g. 'a=>b+c' into connectino pairs
    pairs = []
    for lhs, rhs in statements:
        u1s = [w.strip() for w in lhs.split('+')]
        u2s = [w.strip() for w in rhs.split('+')]
        for u1 in u1s:
            for u2 in u2s:
                pairs.append((u1, u2))

    # Convert into a dict of connection lists
    connections = {}
    for u1, u2 in pairs:
        # add the empty list
        for u in [u1, u2]:
            if u not in connections:
                connections[u] = []
        if u1 == u2:
            raise ValueError("Self connection '{0}' not alowed".format(u1))
        # Check duplicate connections
        if u2 in connections[u1]:
            descr = '{0} to {1}'.format(u1, u2)
            raise ValueError("Duplicate connection from " + descr)
        # store connection
        connections[u1].append(u2)

    # Choose indices for units
    indices = {}
    for n, u in enumerate(sorted(connections.keys())):
        indices[u] = n

    # convert connections into adjacency matrix with specified indices
    N = len(indices.keys())
    M = numpy.zeros((N, N))
    for u1 in connections:
        for u2 in connections[u1]:
            M[indices[u1], indices[u2]] = 1

    return M, indices

def main(args):
    help = """%prog [opts] [connections ...]

    Computes the PageRank of the nodes in a graph using the PageRank algorithm
    as used in google's google rank system.

    Examples of usage:
    %prog "A=>B" "A=>C"
    %prog "A=>B=>A"
    %prog "A=>B+C"
    %prog < graph_file.txt
    """
    parser = optparse.OptionParser(usage=help)
    parser.add_option('-d', '--damping-factor', metavar="D", type='float',
            default=1.00, dest="damping",
            help='damping factor for page-rank computation (default 1.00)')
    parser.add_option('-c', '--no-conserve', action="store_false", dest="conserve",
                      default=True,
                      help="Don't add self connections")
    opts, args = parser.parse_args(sys.argv)

    # Input
    lines = args[1:] or sys.stdin.readlines()
    M, indices = parse_input(lines)

    # Computation
    R = page_rank(M, d=opts.damping, conserve=opts.conserve)
    
    # Output
    print "\ngraph:"
    for line in lines:
        print line

    print "\nadj.matrix:"
    print M

    print "\npage-ranks:"
    for u in sorted(indices, key=lambda k: R[indices[k]]):
        print u, ':', R[indices[u]]


if __name__ == "__main__":
    import sys
    main(sys.argv)


