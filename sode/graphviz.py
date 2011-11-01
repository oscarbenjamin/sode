#!/usr/bin/env python
#
# Author: Oscar Benjamin
# Date: 07 Mar 2011
#
# This module defines a class GraphParser for parsing graphviz syntax into
# python data structures. For simplicity it uses only a subset of the graphviz
# or 'dot' grammar.
#
# In particular the following limitations hold:
# 1) No nested subgraphs. Statements such as {a b {c d}} require a more
#    complex parser that actually tokenizes the input to match braces.
# 2) No multiple square bracket attribute lists e.g.
#    'node [key1=val1, key2=val2];' is ok
#    'node [key1=val1][key2=val2];' is not
# 3) No html tags.
# 4) No semicolons, or escaped quotes in quoted strings.
# 5) No non-digraphs. This would not be very difficult to fix but I don't have
#    time or inclination right now.
#


import sys
import re
from collections import defaultdict


# What follows are the regexes used to parse graphviz files

# Regexes for preprocessing c++ style comments and line continuation
RE_LINE_CONT = r'\\\n'
RE_ML_COMMENT = r'/\*(.|[\n])*?\*/'
RE_SL_COMMENT = r'//.*\n'

# Regexes for elementary features
RE_ID = r'(?:"[^"]*"|[-]?(?:\.[0-9]+|[0-9]+(?:\.[0-9]*)?)|\w+)'
RE_ATTR_LIST = r'\[(.*?)\]'
RE_ATTR_NG = r'\s*{0}\s*=\s*{0}\s*'.format(RE_ID)
RE_ALIST = r'\s*({0})(:?,(.*))?'.format(RE_ATTR_NG)
RE_ATTR = r'\s*({0})\s*=\s*({0})\s*'.format(RE_ID)
RE_EDGE_RHS = r'->\s*({0})\s*'.format(RE_ID)

# Use (?s), ^ and $, since this regex should match the whole graphviz file
RE_GRAPH_STMT = r'(?s)^\s*(strict\s+)?(graph|digraph)(\s+\w+)?\s*\{(.*)\}\s*$'

# No ^ or $ since this is used with re.search to identify subgraph statements
# within graph statement
RE_SUBGRAPH_STMT = r'(?s)\s*(subgraph(\s+(\w+))?\s*)?\{(.*?)\}\s*'

# Regexes for simple statements, these match the beginning of the string
RE_EMPTY_STMT = r'^(\s+|\s*;)$'
RE_ATTR_STMT = r'^\s*(graph|node|edge)\s+{0}\s*$'.format(RE_ATTR_LIST)
RE_NODE_STMT = r'^\s*({0})(\s+{1})?\s*$'.format(RE_ID, RE_ATTR_LIST)
RE_EDGE_STMT = r'^\s*({0})\s*->([^;]*)$'.format(RE_ID)
RE_EDGE_STMT = r'^\s*({0})\s*((->(\s*{0}\s*))+({1})?)$'.format(RE_ID, RE_ATTR_LIST)
RE_EDGE_STMT = r'^\s*({0})\s*(({1})+({2})?)$'.format(RE_ID, RE_EDGE_RHS, RE_ATTR_LIST)

# This regex can be appended to those above in order to distinguish the
# statement from the following text.
RE_MORE_STMTS = r';?(.*)'

# _StatementParser is the general class for all statement parsing code
# Specific statement types are subclassed below.
# The _GraphStatementParser is intended to parse the whole file and is the
# UI to the code presented in these classes
class _StatementParser(object):
    """A class to hold the code used to parse graphviz text.

    _StatementParser is not intended to be instantiated directly.
    Instantiate a GraphParser object which will parse the graph data using a
    _StatementParser object."""

    def __init__(self, text, parameters=None):
        """New empty _StatementParser object

        fin or text can be used to supply graphviz data that will be parsed to
        initialise the _StatementParser instance.
        """
        # Inherit parameters but use deep copy to avoid affecting
        if parameters is None:
            parameters = {'graph':{}, 'node':{}, 'edge':{}}
        else:
            parameters = dict((k, v.copy()) for k, v in parameters.iteritems())

        # Initialise empty data structure
        self.parameters = dict(parameters.copy())
        self.nodes = defaultdict(dict)
        self.edges = defaultdict(lambda : defaultdict(dict))
        self.node_order = []
        self.graph_type = None

        # Attempt to parse input
        if text is not None:
            self.parse(text)

    def parse(self, stmt):
        """Main method of _StatementParser"""
        raise NotImplementedError("Override in subclasses")

    # Used in many subclasses

    def parse_alist(self, alist):
        """Parse comma separated list of attributes"""
        d = {}
        if alist is None:
            return d
        while alist:
            m = re.match(RE_ALIST, alist)
            attr, _, alist = m.groups()
            key, val = re.match(RE_ATTR, attr).groups()[:2]
            d[key] = val
        return d

    # Data management

    def add_node(self, name, attr={}):
        """Add a new node to the graph"""
        if name not in self.node_order:
            self.node_order.append(name)
            self.nodes[name] = {}
            self.edges[name] = {}
        self.nodes[name].update(attr)

    def add_edge(self, nfrom, nto, attr={}):
        """Add a new edge to the graph"""
        if nto not in self.edges[nfrom]:
            self.edges[nfrom][nto] = {}
        self.edges[nfrom][nto].update(attr)

    def add_parameters(self, context, attr):
        """Update parameter information"""
        self.parameters[context].update(attr)

    # Data access

    def has_node(self, name):
        return name in self.nodes

    def has_edge(self, nfrom, nto):
        return nfrom in self.edges and nto in self.edges[nfrom]

    def iter_nodes(self):
        for name in self.node_order:
            yield name, self.nodes[name].copy()

    def iter_edges(self):
        for nfrom in self.edges:
            for nto in self.edges[nfrom]:
                yield nfrom, nto, self.edges[nfrom][nto].copy()

    def get_node(self, name):
        return self.nodes[name].copy()

    def get_edge(self, nfrom, nto):
        return self.edges[nfrom][nto].copy()

    def node_params(self):
        return self.parameters['node'].copy()

    def edge_params(self):
        return self.parameters['edge'].copy()

    def graph_params(self):
        return self.parameters['graph'].copy()

# Subclasses to parse simple graphviz statements

class _EmptyStatementParser(_StatementParser):
    """_StatementParser subclass to parse empty statements"""
    def parse(self, stmt):
        pass


class _AttributeStatementParser(_StatementParser):
    """_StatementParser subclass to parse attribute statements"""
    def parse(self, stmt):
        """Parse global attribute statement"""
        m = re.match(RE_ATTR_STMT, stmt)
        context, alist = m.groups()[:2]
        self.add_parameters(context, self.parse_alist(alist))


class _NodeStatementParser(_StatementParser):
    """_StatementParser subclass to parse node statements"""
    def parse(self, stmt):
        """Parse node attribute statement"""
        m = re.match(RE_NODE_STMT, stmt)
        name, _, alist = m.groups()[:3]
        self.add_node(name, self.parse_alist(alist))


class _EdgeStatementParser(_StatementParser):
    """_StatementParser subclass to parse edge statements"""
    def parse(self, stmt):
        """Parse edge statement"""
        nodes = []
        rhs = stmt
        while True:
            m = re.match(RE_EDGE_STMT, rhs)
            if not m:
                break
            lhs, rhs = m.groups()[:2]
            rhs = rhs[2:] # remove '->' edgeop
            nodes.append(lhs)
        last, _, alist = re.match(RE_NODE_STMT, rhs).groups()
        nodes.append(last)
        for name in nodes:
            self.add_node(name)
        # Register the new edges
        attr = self.parse_alist(alist)
        for nfrom, nto in zip(nodes[:-1], nodes[1:]):
            self.add_edge(nfrom, nto, attr)


# The code here is used by _GraphStatementParser as well.
# 1: parses the subgraph statement and extracts the stmtlist
# 2: parses the stmtlist one statement at a time and amalgamates the
#    resulting data using update.

class _SubgraphStatementParser(_StatementParser):
    """Subclass of _StatementParser that parses graph/subgraph statements"""

    # List of statement types and order to iterate over
    STATEMENT_TYPES = (
        # attribute statements
        {'type':'attribute',
         're_match':RE_ATTR_STMT,
         'parser':_AttributeStatementParser},
        # edge statements
        {'type':'edge',
         're_match':RE_EDGE_STMT,
         'parser':_EdgeStatementParser},
        # node statements
        {'type':'node',
         're_match':RE_NODE_STMT,
         'parser':_NodeStatementParser},
        # blank statements
        {'type':'empty',
         're_match':RE_EMPTY_STMT,
         'parser':_EmptyStatementParser},
    )

    # These regexes are used to split each statement from what follows
    # remove $ and add more stmts regex.
    for d in STATEMENT_TYPES:
        d['re_more'] = '({0}){1}'.format(d['re_match'][:-1], RE_MORE_STMTS)

    def update(self, data):
        """Combine data from sub-_StatementParser instance"""
        # Update global parameters
        for context, d in data.parameters.iteritems():
            self.parameters[context].update(d)

        # Update ordered list of node names
        for name, node in data.iter_nodes():
            self.add_node(name, node)

        # Update edges and corresponding parameters
        for nfrom, nto, edge in data.iter_edges():
            self.add_edge(nfrom, nto, edge)

    def parse(self, stmt):
        # Extract stmtlist
        m = re.match('^{0}$'.format(RE_SUBGRAPH_STMT), stmt)
        _, _, graphid, stmtlist = m.groups()

        self.graph_type = 'subgraph'
        self.graphid = graphid

        # Parse substatements
        self.parse_stmtlist(stmtlist)

    def parse_stmtlist(self, stmtlist):
        """Parse semicolon seperated list of statements"""
        # Newlines are statement ends
        stmtlist = stmtlist.replace('\n', ';')
        # Parse 1 statement at a time until stmtlist is empty
        while stmtlist:
            statement_parser, stmtlist = self.parse_stmt(stmtlist)
            self.update(statement_parser)

    def parse_stmt(self, stmtlist):
        """Parse individual statement. Delegates to sister subclass"""
        stmt = ''
        try:
            for spd in self.STATEMENT_TYPES:
                m = re.match(spd['re_more'], stmtlist)
                if m:
                    g = m.groups()
                    stmt, stmtlist = g[0], g[-1]
                    statement_parser = spd['parser'](stmt)
                    return statement_parser, stmtlist
            else:
                m = "Unrecognised statement type: '{0}'"
                raise ValueError(m.format(stmtlist))
        except:
            print >> sys.stderr, "parse error: '{0}'".format(stmtlist)
            raise


class _GraphStatementParser(_SubgraphStatementParser):
    """Subclass of _StatementParser that parses the graph statement.

    The graph statement is the whole of a graphviz file and the
    _GraphStatementParser should be used to parse the whole file.
    """
    # This class contains all the code to handle preprocessing and subgraphs
    # 1: comments are preprocessed.
    # 2: subgraph statements are extracted processed and replaced with
    #    ids (using the subgraph id, if provided)
    # 3: statements are parsed as in _SubgraphStatementParser except that
    #    references to subgraph ids are handled specially
    # 4: properties defined in subgraphs are processed at the end

    def parse(self, stmt):
        """Parse complete graph"""
        # Check / extract graph type and stmtlist
        stmt = self.preprocess(stmt)
        m = re.match(RE_GRAPH_STMT, stmt)
        if not m:
            print stmt
            raise ValueError("Does not appear to be valid graphviz data")
        strict, graph_type, graphid, stmtlist = m.groups()

        # Parse substatements
        self.strict = bool(strict)
        self.graph_type = graph_type
        self.graphid = graphid

        # Convert subgraph statements
        self.subgraph_stmts = {}
        self.subgraphs = {}
        self.subgraph_order = []
        stmtlist = self.parse_subgraphs(stmtlist)

        # Parse statement list as normal
        self.parse_stmtlist(stmtlist)


    # This is the top of the tree, so preprocess text
    def preprocess(self, text):
        # replace comments with whitespace
        text = re.sub(RE_ML_COMMENT, ' ', text)
        text = re.sub(RE_SL_COMMENT, ' ', text)
        # replace contuation charaters
        text = re.sub(RE_LINE_CONT, '', text)
        return text

    def parse_subgraphs(self, stmtlist):
        # find all subgraph statements and replace them
        while True:
            # Match next subgraph statements
            m = re.search(RE_SUBGRAPH_STMT, stmtlist)
            if not m:
                return stmtlist
            # Parse subgraph and register
            subgraph_stmt = m.group(0)
            stmtlist = self.replace_subgraph_stmt(stmtlist, subgraph_stmt)

    def replace_subgraph_stmt(self, stmtlist, stmt):
        # Parse subgraph statement to obtain graphid
        subgraph = _SubgraphStatementParser(stmt, self.parameters)

        # Make a gid for anonymous graphs
        gid = subgraph.graphid
        if gid is None:
            gid = '__{0}__'.format(''.join(subgraph.nodes))
            gid += str(hash(stmt))[1:]
        # Store by gid
        self.subgraph_stmts[gid] = stmt

        # replace occurences of subgraph text with id
        return stmtlist.replace(stmt, ' {0} '.format(gid))

    # Intercept these function to distribute over nodes when subgraph names
    # are used

    def add_node(self, name, attr={}):
        """Add a node if not a subgraph"""
        if name in self.subgraph_stmts:
            # Create the subgraph instance now so it inherits the parameters
            if not name in self.subgraphs:
                stmt = self.subgraph_stmts[name]
                subgraph = _SubgraphStatementParser(stmt, self.parameters)
                self.subgraphs[name] = subgraph
            # Now process the subgraph data
            self.add_subgraph(self.subgraphs[name])
        else:
            # Apply default values
            if not self.has_node(name):
                attr = self._default_node(name, attr)
            super(_GraphStatementParser, self).add_node(name, attr)

    def add_edge(self, nfrom, nto, attr={}):
        """Add a new edge to the graph"""
        nfroms, ntos = [nfrom], [nto]
        if nfrom in self.subgraph_stmts:
            nfroms = self.subgraphs[nfrom].node_order
        if nto in self.subgraph_stmts:
            ntos = self.subgraphs[nto].node_order
        for nf in nfroms:
            for nt in ntos:
                if not self.has_edge(nf, nt):
                    attr = self._default_edge(nf, nt, attr)
                super(_GraphStatementParser, self).add_edge(nf, nt, attr)

    def add_subgraph(self, subgraph):
        # Update ordered list of subgraphs
        self.subgraph_order.append(subgraph)

        # Update ordered list of node names
        for name, node in subgraph.iter_nodes():
            self.add_node(name, node)

        # Update edges and corresponding parameters
        for nfrom, nto, edge in subgraph.iter_edges():
            self.add_edge(nfrom, nto, edge)

    # Functions called for unrecognised nodes/edges. The node/edge should be
    # part of the most recent subgraph if any.

    def _default_node(self, name, attr):
        """Check most recently added subgraph (should contain node)"""
        d = self.node_params()
        if self.subgraph_order:
            sg = self.subgraph_order[-1]
            if sg.has_node(name):
                d.update(sg.node_params())
        d.update(attr)
        return d

    def _default_edge(self, nfrom, nto, attr):
        """Check most recently added subgraph (should contain edge)"""
        d = self.edge_params()
        if self.subgraph_order:
            sg = self.subgraph_order[-1]
            if sg.has_edge(nfrom, nto):
                d.update(sg.edge_params())
        d.update(attr)
        return d

class GraphParser(_StatementParser):
    """The GraphParser instance brings together all the data from a _StatementParser
    object to make it more accesible"""
    def __init__(self, infile=None, text=None):
        """Read file or text and parse as graphviz file"""
        # Open / read file if necessary
        if infile is not None and text is None:
            if isinstance(infile, str):
                infile = open(infile, 'r')
            text = infile.read()

        # Parse input
        data = _GraphStatementParser(text)

        # Initialise empty data
        self.graph_params = {}
        self.nodes = {}
        self.edges = {}
        self.node_order = []

        # Collate graph params, nodes and edges

        self.graph_params.update(data.graph_params())

        for name, node in data.iter_nodes():
            self.add_node(name, node)

        for nfrom, nto, edge in data.iter_edges():
            self.add_edge(nfrom, nto, edge)



    def add_node(self, name, pars):
        """New node with node parameters npars and global parameters gpars """
        self.node_order.append(name)
        self.nodes[name] = pars.copy()
        self.edges[name] = {}

    def add_edge(self, nfrom, nto, pars):
        self.edges[nfrom][nto] = pars.copy()

    def print_data(self):
        """Print the data in graphviz format with all values explicit"""
        print 'digraph {'

        for name in self.node_order:
            pars = self.nodes[name]
            parstring = ', '.join('='.join([k, v]) for k, v in pars.iteritems())
            print '    {0} [{1}]'.format(name, parstring)

        for nfrom, d in self.edges.iteritems():
            for nto, pars in d.iteritems():
                parstring = ', '.join('='.join([k, v]) for k, v in pars.iteritems())
                print '    {0} -> {1} [{2}]'.format(nfrom, nto, parstring)

        print '}'


def test_graph_parser(output=False):
    """Run GraphParser on text graph"""

    # Test data

    _TEST_GRAPH = """
    digraph {
        /* I'll put in some comments for good measure here
        */
        graph [type="special"]
        node [alpha=alpha_a]
        node [beta=beta_a]
        edge [strength=ab_ac] a // two statements in one line
        node [beta=beta_bcd]
        b [alpha=alpha_b]
        {
            node [alpha=alpha_c]
            edge [strength=ca]
            c->a
        }
        a -> {b c}
        node [alpha=alpha_d]
        d
    }
    """

    _EXPECTED_GRAPHVALS = {
        'type':'"special"',
    }

    _EXPECTED_NODES = {
        'a':{'alpha':'alpha_a', 'beta':'beta_a'  },
        'b':{'alpha':'alpha_b', 'beta':'beta_bcd'},
        'c':{'alpha':'alpha_c', 'beta':'beta_bcd'},
        'd':{'alpha':'alpha_d', 'beta':'beta_bcd'},
    }

    _EXPECTED_EDGES = {
        ('a', 'b'):{'strength':'ab_ac'},
        ('a', 'c'):{'strength':'ab_ac'},
        ('c', 'a'):{'strength':'ca'   },
    }

    # Functions used in testing

    numerrs = [0]
    def errmsg(msg):
        if output:
            print msg
        numerrs[0] += 1

    def compare_sets(descr, s_expected, s_found):
        for name in s_expected - s_found:
            errmsg("{0} '{1}' not found".format(descr, name))
        for name in s_found - s_expected:
            errmsg("{0} '{1}' not expected".format(descr, name))
        return s_expected & s_found

    def compare_dicts(descr, d_expected, d_found):
        s_e = set(d_expected)
        s_f = set(d_found)
        for key in compare_sets('{0}: key'.format(descr), s_e, s_f):
            if d_expected[key] != d_found[key]:
                msg ="{0}: key '{1}': expected: '{2}': found '{3}'"
                errmsg(msg.format(descr, key, d_expected[key], d_found[key]))

    # Print original, parse and print parsed result

    if output:
        print "Original graph ------------------------------"
        print _TEST_GRAPH

    try:
        g = GraphParser(text=_TEST_GRAPH)
    except:
        if output:
            raise
        else:
            errmsg('Error raised while parsing graph')
            return False

    if output:
        print "Parsed graph -------------------------------"
        g.print_data()

    # Compare graph properties, nodes and then edges

    compare_dicts('graph parameter', _EXPECTED_GRAPHVALS, g.graph_params)

    sn_expected = set(_EXPECTED_NODES)
    sn_found = set(g.nodes)
    for name in compare_sets('nodes', sn_expected, sn_found):
        msg = "node '{0}'".format(name)
        compare_dicts(msg, _EXPECTED_NODES[name], g.nodes[name])

    se_expected = set(_EXPECTED_EDGES)
    se_found = set((nfrom, nto) for nfrom in g.edges for nto in g.edges[nfrom])
    for nfrom, nto in compare_sets('edges', se_expected, se_found):
        msg = "edges '{0}' -> '{1}'".format(nfrom, nto)
        compare_dicts(msg, _EXPECTED_EDGES[(nfrom, nto)], g.edges[nfrom][nto])

    # Print summary and return success

    n, = numerrs
    if output:
        print "{0} errors!!!!".format(n)
    return not n


if __name__ == "__main__":
    import sys
    import opster

    @opster.command(usage='%name [INFILE]')
    def main(infile=None,
             test=('t', False, 'load a testgraph')):
        """Load INFILE or stdin as graph file and print a description to stdout"""

        # Run test
        if test:
            return not test_graph_parser(output=True)

        # Or parse input graph
        infile = infile or sys.stdin
        g = GraphParser(infile)
        g.print_data()

    main(argv=sys.argv[1:])



