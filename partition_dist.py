#!/usr/bin/env python
import sys
from copy import copy

def parse_mb_header(h):
    s = h.split()
    n = len(s)
    n_sites = (n - 2)//2
    assert s[n_sites + 1] == "P[%d]" % n_sites
    assert s[n_sites + 2] == "R[1]"
    return n_sites
def read_mb_partitions(sampled_partitions_file, from_row, to_row=None):
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
        if to_row and to_row < line_num:
            break
        s = line.strip().split()
        name = s.pop(0)
        num_elements = int(s.pop(0))
        list_of_sets = []
        assignments, rates = s[:n_sites], s[n_sites:]
        assert(len(assignments) == len(rates))
        sys.stderr.write("reading sample %d\n" % line_num)
        for index, el in enumerate(assignments):
            iel = int(el)
            while iel > len(list_of_sets):
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
        out.write("%s\t%d\t%s\t%s\n" % (name, len(los), "\t".join(c), "\t".join(r)))

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-s", "--seed", dest="seed", default=0, 
        type="int",
        help="The random number generator seed")
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
    write_mb_partitions(sys.stdout, sampled_partitions, n_sites)
