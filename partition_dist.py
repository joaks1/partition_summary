#!/usr/bin/env python
import sys
import random
import itertools
from copy import copy
_VERBOSE = False
_DEBUGGING = False
_RNG = random.Random()

def parse_mb_header(h):
    s = h.split()
    n = len(s)
    n_sites = (n - 2)//2
    assert(s[n_sites + 1] == "P[%d]" % n_sites)
    assert(s[n_sites + 2] == "R[1]")
    return n_sites

def read_mb_partitions(sampled_partitions_file, from_row, to_row=None, read_rates=True):
    """Returns a tuple of the list of partitions and the number of sites
    
    The format in the DPP output from MrBayes is 
    rep#\t#subsets\tsubset#\tsubset#\t....
    where the number of subset numbers corresponds to the number of characters
    in the data matrix.
    """
    d = []
    line_iter = iter(sampled_partitions_file)
    header = line_iter.next()
    n_sites = parse_mb_header(header)
    for line_num, line in enumerate(line_iter):
        if from_row and from_row > line_num:
            continue
        if to_row and to_row < line_num:
            break
        x = line.strip()
        s = x.split()
        name = s.pop(0)
        num_elements = int(s.pop(0))
        list_of_subsets = []
        assignments = s[:n_sites]
        if read_rates:
            rates = s[n_sites:]
            assert(len(assignments) == len(rates))
        else:
            rate = 1.0
        if _VERBOSE:
            sys.stderr.write("reading sample %d\n" % line_num)
        for index, el in enumerate(assignments):
            iel = int(el)
            while iel > len(list_of_subsets):
                if read_rates:
                    rate = float(rates[index])
                subset_info = (set(), rate)
                list_of_subsets.append(subset_info)
            set_of_indices = list_of_subsets[iel - 1][0]
            set_of_indices.add(index)
        assert len(list_of_subsets) == num_elements
        d.append((name, list_of_subsets))
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
    #    print '(%d, %d) -> %d' % (row, column, value)
    #print 'total cost: %d' % total
    return n_sites - total

def permute_partition(partition_desc, n_sites):
    global _RNG
    site_to_subset_ind = [None] * n_sites
    permuted = []
    for ss_index, subset_rate_pair in enumerate(partition_desc):
        permuted.append([set(), subset_rate_pair[1]])
        subset = subset_rate_pair[0]
        for col_index in subset:
            assert(site_to_subset_ind[col_index] is None)
            site_to_subset_ind[col_index] = ss_index

    _RNG.shuffle(site_to_subset_ind)
    
    for index, iel in enumerate(site_to_subset_ind):
        assert (iel is not None)
        set_of_indices = permuted[iel][0]
        set_of_indices.add(index)
    return permuted

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
    parser.add_option("-m", "--median", dest="median", default=False, 
        action="store_true",
        help="Find the sampled partition that has the smallest distance to all of the others (the sample that is closest to being the median)")
    parser.add_option("-r", "--random", dest="random", default=False, 
        action="store_true",
        help="Compare the dist from a specified partition  (the second first partition in the file that is the second arg) to the sampled partitions against the distance from the same specified partition to a random permutation of sampled partition")
    (options, args) = parser.parse_args()
    
    
    _VERBOSE = options.verbose
    _DEBUGGING = options.debugging
    if options.seed > 0:
        _RNG.seed(options.seed)
    if not args:
        sys.exit("Expecting a file name for the MrBayes DPP-over rates samples")
    sampled_partitions_filename = args[0]
    sampled_partitions_file = open(sampled_partitions_filename, 'rU')
    sampled_partitions, n_sites = read_mb_partitions(sampled_partitions_file, options.from_index, options.to_index)
    sampled_partitions_file.close()

    if _DEBUGGING and _VERBOSE:
        write_mb_partitions(sys.stdout, sampled_partitions, n_sites)
    if _VERBOSE:
        sys.stderr.write("%d partitions read\n" % len(sampled_partitions))
    if len(args) > 1:
        test_partitions_filename = args[1]
        test_partitions_file = open(test_partitions_filename, 'rU')
        test_partitions, test_n_sites = read_mb_partitions(test_partitions_file, 0, 1, read_rates=False)
        test_partitions_file.close()
        if test_n_sites != n_sites:
            sys.exit("The number of sites in the test partition must be identical to the number of sites in the sampled partitions from MrBayes")
        tp = test_partitions[0][1]
        real_tally = 0
        rand_tally = 0
        tie_tally = 0
        total_tally = 0
        for i in sampled_partitions:
            real_dist = partition_distance(tp, i[1], n_sites)
            if options.random:
                permuted = permute_partition(i[1], n_sites)
                rand_dist = partition_distance(tp, permuted, n_sites)
                if real_dist < rand_dist:
                    real_tally = real_tally + 1
                elif real_dist > rand_dist:
                    rand_tally = rand_tally + 1
                elif real_dist == rand_dist:
                    tie_tally = tie_tally + 1
                total_tally = total_tally + 1
                print "%d\t%d" % (real_dist, rand_dist)
            else:
                print real_dist
        if options.random:
            assert(real_tally + rand_tally + tie_tally == total_tally == len(sampled_partitions))
            prob = float(real_tally)/float(total_tally)
            print "\nNumber of sampled 'wins':  %d" % real_tally
            print "Number of permuted 'wins':  %d" % rand_tally
            print "Number of ties:  %d" % tie_tally
            print "TOTAL:  %d" % total_tally
            print "Probability that the DPP sampled partitions are closer to a priori partition than random:"
            print prob
    else:
        l = len(sampled_partitions)
        lower_triangle = []
        for n, i in enumerate(sampled_partitions):
            if _VERBOSE:
                sys.stderr.write("Calc distance matrix row %d\n" % n)
            r = sampled_partitions[:n+1]
            d = tuple(partition_distance(i[1], j[1], n_sites) for j in r)
            if _VERBOSE:
                sys.stdout.write("%s\n" % "\t".join([str(x) for x in d]))
            lower_triangle.append(d)
        if options.median:
            if _VERBOSE:
                print lower_triangle
            dim = len(lower_triangle)
            sum_dist = [0]*dim
            for i in xrange(dim):
                for j in xrange(i):
                    row = lower_triangle[i]
                    element = row[j]
                    sum_dist[i] += element
                    sum_dist[j] += element
            min_dist = min(sum_dist)
            print "samples that are medianish (sum of partition dist of %d):" % min_dist
            for n, el in enumerate(sum_dist):
                if el == min_dist:
                    print sampled_partitions[n][0]
                    print "\nCharset definitions for %s:\n" % sampled_partitions[n][0]
                    for i in xrange(len(sampled_partitions[n][1])):
                        sys.stdout.write("charset ratepart%d =" % (i+1))
                        for site in sampled_partitions[n][1][i][0]:
                            sys.stdout.write(" %d" % (site+1))
                        sys.stdout.write(";\n")
                    print "\nRates for partition %s:"% sampled_partitions[n][0]
                    for i in xrange(len(sampled_partitions[n][1])):
                        print "ratepart%d = %f" % ((i+1), sampled_partitions[n][1][i][1])
                    print "\n"
                             
                        