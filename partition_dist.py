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

    cost_matrix = []
    for row in mat:
        cost_row = [sys.maxint - col for col in row]
        cost_matrix.append(cost_row)
    from munkres import Munkres
    indexes = Munkres().compute(cost_matrix)
    
    total = 0
    for row, column in indexes:
        value = mat[row][column]
        total += value
        print '(%d, %d) -> %d' % (row, column, value)
    print 'total cost: %d' % total
    return n_sites - total

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
        
        
