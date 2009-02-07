#!/usr/bin/env python
import sys
import itertools
from copy import copy
_VERBOSE = False
_DEBUGGING = False
def parse_mb_header(h):
    s = h.split()
    n = len(s)
    n_sites = (n - 2)//2
    assert s[n_sites + 1] == "P[%d]" % n_sites
    assert s[n_sites + 2] == "R[1]"
    return n_sites

def read_mb_partitions(sampled_partitions_file, from_row, to_row=None, read_rates=True):
    """Returns a list of tuples of (name, list of sets) for each sampled
    partition, where the name is the MCMC iteration in the MrBayes file.
    
    The format in the DPP output from MrBayes is 
    rep#\t#subsets\tsubset#\tsubset#\t....
    where the number of subset numbers corresponds to the number of characters
    in the data matrix.
    """
    d = []
    line_iter = iter(sampled_partitions_file)
    header = line_iter.next()
    n_sites = parse_mb_header(header)
    print n_sites
    for line_num, line in enumerate(line_iter):
        if from_row and from_row > line_num:
            continue
        if to_row and to_row <= line_num:
            break
        s = line.strip().split()
        name = s.pop(0)
        num_elements = int(s.pop(0))
        list_of_sets = []
        assignments = s[:n_sites]
        if read_rates:
            rates = s[n_sites:]
            assert(len(assignments) == len(rates))
        else:
            rate = 1.0
        sys.stderr.write("reading sample %d\n" % line_num)
        for index, el in enumerate(assignments):
            iel = int(el)
            while iel > len(list_of_sets):
                if read_rates:
                    rate = float(rates[index])
                subset_info = (set(), rate)
                list_of_sets.append(subset_info)
            dest_set = list_of_sets[iel - 1][0]
            dest_set.add(index)
        assert len(list_of_sets) == num_elements
        d.append((name,list_of_sets))
    return d, n_sites

def write_mb_partitions(out, sampled_partions, n_sites):
    for v in sampled_partions:
        name, los = v
        c = [None] * n_sites
        r = [None] * n_sites
        for set_index, subset_info in enumerate(los):
            s, rate = subset_info
            set_num = str(set_index + 1)
            rate_str = str(rate)
            for element in s:
                assert(c[element] is None) # only trips if los is not a partition
                c[element] = set_num
                r[element] = rate_str
        out.write("%s\t%d\t%s\t%s\t\n" % (name, len(los), "\t".join(c), "\t".join(r)))

class Vertex(object):
    def __init__(self, index, label_weight, is_in_y=False):
        self.index = index
        self.label_weight = label_weight
        self.is_in_y = is_in_y
        self.eq_graph_edges = set()
        self.matching_edge = None
    def get_eq_graph_neighbors(self):
        if self.is_in_y:
            return [i.tail for i in self.eq_graph_edges]
        return [i.head for i in self.eq_graph_edges]
    def __str__(self):
        return "%s%d:%d{%s%s}" % (self.is_in_y and "y" or "x", 
                                self.index, 
                                self.label_weight, 
                                self.eq_graph_edges and "e" or "f",
                                self.matching_edge and "*" or "")
                                
class Edge(object):
    def __init__(self, tail, head, in_eq_graph=True):
        self.tail = tail
        self.head = head
        if in_eq_graph:
            tail.eq_graph_edges.add(self)
            head.eq_graph_edges.add(self)
        self.in_match = False
    def set_in_match(self, in_match=True):
        m = in_match and self or None
        self.in_match = in_match
        self.tail.matching_edge = m
        self.head.matching_edge = m
    def set_in_eq_graph(self, in_eq_graph):
        if in_eq_graph:
            self.tail.eq_graph_edges.add(self)
            self.head.eq_graph_edges.add(self)
        else:
            self.tail.eq_graph_edges.remove(self)
            self.head.eq_graph_edges.remove(self)
    def __str__(self):
        return "%s ==> %s" % (str(self.tail), str(self.head))

class EqualityGraph(object):
    def __init__(self, dim):
        self.edge_mat = []
        for i in range(dim):
            self.edge_mat.append([None]*dim)
        self.edge_to_coord = {}

    def __str__(self):
        s = []
        for row in self.edge_mat:
            rp = []
            for cell in row:
                if cell is None:
                    rp.append("-")
                elif cell.in_match:
                    rp.append("*")
                else:
                    rp.append("+")
            s.append(" ".join(rp))
        return "\n\t%s" % "\n\t".join(s)
    def add(self, x, y):
        x_c, y_c = x.index, y.index
        em_row = self.edge_mat[x_c]
        e = em_row[y_c]
        if e is None:
            e = Edge(x, y, in_eq_graph=True)
            em_row[y_c] = e
            self.edge_to_coord[e] = (x_c, y_c)
        return e

    def remove(self, e_or_x, y=None):
        if y is not None:
            x_el = e_or_x
            x_c, y_c = x_el.index, y.index
            edge = self.edge_mat[x_c][y_c]
            if edge is None:
                return
        else:
            edge = e_or_x
            x_c, y_c = self.edge_to_coord[edge]
        edge.set_in_eq_graph(False)
        del self.edge_to_coord
        self.edge_mat[x_c][y_c] = None

    def augment_path(self, u, y, matching):
        aug_path_edges = self.find_augmenting_path(u, y)
        for n, edge in enumerate(aug_path_edges):
            if n % 2:
                edge.set_in_match(False)
                matching.remove(edge)
            else:
                edge.set_in_match(True)
                matching.add(edge)

    def find_augmenting_path(self, u, y):
        #augmenting paths should hav free endpoints
        assert(y.matching_edge is None)
        for edge in u.eq_graph_edges:
            if edge is u.matching_edge:
                continue
            assert(edge.tail is u)
            h = edge.head
            if h is y:
                return [edge]
            m = h.matching_edge
            if m:
                rpaths = find_augmenting_path(m.tail, y)
                if rpaths:
                    return [edge, h.matching_edge] + rpaths
        return None

def score_from_matching(matching):
    sc = 0
    v_pair = []
    for edge in matching:
        x, y, w = edge.tail, edge.head, edge.label_weight
        v_pair.append((x, y))
        sc += w
    return v_pair, w


def get_neighboring_vertices(s):
    n = set()
    for vertex in s:
        n.update(vertex.get_eq_graph_neighbors())
    return n


def update_equality_graph(eq_graph, s, t, mat, x_vec, y_vec):
    # update the equality graph based on the new labelings
    #   only edges that contain the vertices in s and t need to be 
    #   checked.
    # add equality edges from s
    for x_el in s:
        x_w = x_el.label_weight
        weight_row = mat[x_el.index]
        for y_el in y_vec:
            y_w = y_el.label_weight
            w = weight_row[y_el.index]
            if (x_w + y_w) == w:
                eq_graph.add(x_el, y_el)
            else:
                eq_graph.remove(x_el, y_el)
    # add equality edges from t
    for y_el in t:
        offset = y_el.index
        y_w = y_el.label_weight
        for x_el, row in itertools.izip(x_vec, mat):
            x_w = x_el.label_weight
            w = mat[x_el.index][offset]
            if (x_w + y_w) == w:
                eq_graph.add(x_el, y_el)
            else:
                eq_graph.remove(x_el, y_el)


def update_labelings(all_y, t, s, mat):
    # update the labelings
    # first find the minimum alpha
    print mat
    alpha = None
    not_t = all_y - t
    assert(len(not_t) > 0)
    for x_el in s:
        x_w = x_el.label_weight
        row = mat[x_el.index]
        for y_el in not_t:
            y_w = y_el.label_weight
            e_w = row[y_el.index]
            if _DEBUGGING:
                sys.stdout.write("x=%d %d, y=%d %d, w=%d\n" % (x_el.index, x_w, y_el.index, y_w, e_w))
            el_alpha = x_w + y_w - e_w
            if alpha is None or alpha > el_alpha:
                alpha = el_alpha
    assert(alpha is not None)
    # update the labelings by alpha
    for x_el in s:
        x_el.label_weight = x_el.label_weight - alpha
    for y_el in t:
        y_el.label_weight = y_el.label_weight + alpha


def calc_max_assignment(mat):
    """Uses the Hungarian algorithm to find the maximum weighted perfect matching.
    see notes http://www.cse.ust.hk/~golin/COMP572/Notes/Matching.pdf
    for description of the algorithm and the notation.
    """
    print mat
    dim = len(mat)

    # Initialization:
    # Assign the initial labels (all y's get 0, and all x's get the cell weight)
    # Construct the equality graph, and 
    # Find the initial matching

    y_vec = [Vertex(n, 0, True) for n in xrange(dim)]
    x_vec = []
    eq_graph = EqualityGraph(dim)
    matching = []
    m_y_vertices = set()
    m_x_vertices = set()
    for n, row in enumerate(mat):
        assert (dim == len(row))
        m = max(row)
        x_el = Vertex(n, m, False)
        x_vec.append(x_el)
        for col_n, cell in enumerate(row):
            if cell == m:
                y_el = y_vec[col_n]
                e = eq_graph.add(x_el, y_el)
                if (y_el not in m_y_vertices) and (x_el not in m_x_vertices):
                    # Edges of the equality graph can be part of the matching
                    #   as long as the vertices do not already exist in the matching
                    m_x_vertices.add(x_el)
                    m_y_vertices.add(y_el)
                    matching.append(e)
                    e.set_in_match(True)
    # if the matching is perfect, we are done.
    if len(matching) == dim:
        return score_from_matching(matching)
    # Add a free vertex to the set S
   
    # The set t is the empty set
    t = set()
    all_y = set(y_vec)
    reset_s = True
    while len(matching) != dim:
        if reset_s:
            s = set()
            for x in x_vec:
                if x.matching_edge is None:
                    s.add(x)
                    break
            assert len(s) == 1
        if _DEBUGGING:
            sys.stdout.write("""
main while loop:
s = set([%s])
t = set([%s])
eq_graph = %s
matching = %s
x_vec = %s
y_vec = %s
""" % ( ", ".join([str(i) for i in s]),
                    ", ".join([str(i) for i in t]),
                    str(eq_graph),
                    ", ".join(["%s -> %s" % (str(i.tail), str(i.head)) for i in matching]),
                    ", ".join([str(i) for i in x_vec]),
                    ", ".join([str(i) for i in y_vec]),
                  ))
        reset_s = False
        neighborhood = get_neighboring_vertices(s)
        if neighborhood == t:
            if _DEBUGGING:
                sys.stdout.write("Updating labeling\n")
            update_labelings(all_y, t, s, mat)
            update_equality_graph(eq_graph, s, t, mat, x_vec, y_vec)
        else:
            diff = neighborhood - t
            y = diff.pop()
            if y.matching_edge is None:
                if _DEBUGGING:
                    sys.stdout.write("Augmenting match\n")
                # y is free
                print s
                assert len(s) == 1
                u = s.pop()
                eq_graph.augment_path(u, y, matching)
                reset_s = True
            else:
                if _DEBUGGING:
                    sys.stdout.write("Adding to alternating tree\n")
                t.add(y)
                s.add(y.matching_edge.tail)
    return score_from_matching(matching)

def diagonal_assignment(mat):
    print mat
    dim = len(mat)
    for row in mat:
        assert (dim == len(row))
    s = 0
    p = [None]*dim
    row_v = [max(i) for i in mat]
    col_v = [max([row[col_index] for row in mat]) for col_index in xrange(dim)]
    for i in xrange(dim):
        s += mat[i][i]
        p[i] = (i,i)
    return p, s 

def partition_distance(x, y, n_sites):
    mat = []
    lx = len(x)
    ly = len(y)
    dim = max(lx, ly)
    for xel in x:
        row = [0] * dim
        xset = xel[0]
        for i, yel in enumerate(y):
            yset = yel[0]
            intersection = xset & yset
            row[i] = len(intersection)
        mat.append(row)
    n_to_add = dim - lx
    for i in xrange(n_to_add):
        mat.append([0] * dim)
    max_assignment, max_assignment_score = calc_max_assignment(mat)
    return n_sites - max_assignment_score

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-s", "--seed", dest="seed", default=0, 
        type="int",
        help="The random number generator seed")
    parser.add_option("-v", "--verbose", dest="verbose", default=False, 
        action="store_true",
        help="Verbose output")
    parser.add_option("-d", "--debugging", dest="debugging", default=False, 
        action="store_true",
        help="Run in debugging mode")
    parser.add_option("-f", "--from", dest="from_index", default=0, 
        type="int",
        help="The index of the first row of sampled points to read (0 by default)")
    parser.add_option("-t", "--to", dest="to_index", default=None, 
        type="int",
        help="The index of the last row of sampled points to read (None by default)")
    (options, args) = parser.parse_args()
    if not args:
        sys.exit("Expecting a file name for the MrBayes DPP-over rates samples")
    sampled_partitions_filename = args[0]
    sampled_partitions_file = open(sampled_partitions_filename, 'rU')
    sampled_partitions, n_sites = read_mb_partitions(sampled_partitions_file, options.from_index, options.to_index)
    sampled_partitions_file.close()
    _VERBOSE = options.verbose
    _DEBUGGING = options.debugging
    if _VERBOSE:
        write_mb_partitions(sys.stdout, sampled_partitions, n_sites)
    if len(args) > 1:
        test_partitions_filename = args[1]
        test_partitions_file = open(test_partitions_filename, 'rU')
        test_partitions, test_n_sites = read_mb_partitions(test_partitions_file, 0, 1, read_rates=False)
        test_partitions_file.close()
        if test_n_sites != n_sites:
            sys.exit("The number of sites in the test partition must be identical to the number of sites in the sampled partitions from MrBayes")
        tp = test_partitions[0][1]
        d = [partition_distance(tp, i[1], n_sites) for i in sampled_partitions]
        sys.stdout.write("%s\n" % "\n".join([str(i) for i in d]))
        
        
