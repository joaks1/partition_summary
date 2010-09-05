#!/usr/bin/env python
import sys
import random
import itertools
from copy import copy
import logging
from dendropy.utility.messaging import get_logger
_LOG = get_logger('partition_dist')

_VERBOSE = False
_DEBUGGING = False
_RNG = random.Random()

class Subset(object):
    "object with a set of site indices (ints) and a corresponding relative rate (float), and string id"
    
    def __init__(self, set_of_site_indices = set(), rate = None, id = None):
        assert isinstance(set_of_site_indices, set)
        self._site_indices = set()
        self._number_of_sites = len(self._site_indices)
        for i in set_of_site_indices:
            self.add_index(i)
        if rate is None:
            self._rate = rate
        else:
            self.rate = rate
        if id is not None:
            self._id = str(id)
        else:
            self._id = id

    def __str__(self):
        s = "sites = %s\nrate = %s"  % (" ".join([str(x) for x in self._site_indices]), str(self.rate))
        return s

    def _get_rate(self):
        return self._rate
    def _set_rate(self, float_rate):
        try:
            self._rate = float(float_rate)
        except ValueError, e:
            sys.stderr.write("Invalid value specified for rate attribute of Subset object: %s" % str(float_rate))
            raise
    rate = property(_get_rate, _set_rate)
        
    def _get_id(self):
        return self._id
    def _set_id(self, str_name):
        self._id = str(str_name)
    id = property(_get_id, _set_id)
    
    def _get_indices(self):
        return self._site_indices
    indices = property(_get_indices)

    def get_indices_str(self):
        s = " ".join([str(x+1) for x in self._site_indices])
        return s
        
    def _get_number_of_sites(self):
        return self._number_of_sites
    size = property(_get_number_of_sites)

    def add_index(self, int_site_index):
        assert isinstance(int_site_index, int), "Site indices of Subset object must be ints, you specifed: %s -- %s" % (str(int_site_index), type(int_site_index))    
        self._site_indices.add(int_site_index)
        self._number_of_sites += 1
        

class Partition(object):
    "object with a list of Subset objects, a string id, and float likelihood"
    count = 0
    def __init__(self, list_of_subset_objects = [], id = None, lnL = None):
        self.__class__.count += 1
        assert isinstance(list_of_subset_objects, list)
        self._subsets = []
        self._number_of_subsets = len(self._subsets)
        self._length = 0
        for ss in list_of_subset_objects:
            self.add_subset(ss)
        if id is not None:
            self._id = str(id)
        else:
            self._id = 'partition' + str(self.count)
        if lnL is not None:
            try:
                self._lnL = float(lnL)
            except ValueError, e:
                sys.stderr.write("Invalid value specified for lnL attribute of Partition object: %s" % str(lnL))
                raise
        else:
            self._lnL = lnL
        self._length_update_needed = 0
        self._number_of_permutations = 0
        
    def __str__(self):
        if self._id is not None:
            s = "begin sets;\n\t[partition %s]\n" % self._id
            paup_definition = "charpartition %s = " % self._id
            mb_definition = "partition %s = %d: " % (self._id, self._number_of_subsets)
        else:
            s = "begin sets;\n\t[partition p%d]\n" % self._number_of_subsets
            paup_definition = "charpartition p%d =" % self._number_of_subsets
            mb_definition = "partition p%d = %d:" % (self._number_of_subsets, self._number_of_subsets)
        for i, subset in enumerate(self._subsets):
            if subset._get_id() is not None:
                s = s + "\tcharset %s" % subset._get_id()
                s = s + " [rate = %s] = %s;\n" % (str(subset.rate), subset.get_indices_str())
                paup_definition = paup_definition + " %s:%s," % (subset._get_id(), subset._get_id())
                mb_definition = mb_definition + " %s," % subset._get_id()
            else:
                s = s + "\tcharset subset%d" % (i+1)
                s = s + " [rate = %s] = %s;\n" % (str(subset.rate), subset.get_indices_str())
                paup_definition = paup_definition + " subset%d:subset%d," % (i+1, i+1)
                mb_definition = mb_definition + " subset%d," % (i+1)
        paup_definition = paup_definition.rstrip(',') + ";"
        mb_definition = mb_definition.rstrip(',') + ";"
        s = s + "\t%s\n\t%s\nend;" % (paup_definition, mb_definition)
        return s
    
    def _get_id(self):
        return self._id
    def _set_id(self, str_name):
        self._id = str(str_name)
    id = property(_get_id, _set_id)
    
    def _get_lnL(self):
        return self._lnL
        
    lnL = property(_get_lnL)
    
    def _get_number_of_subsets(self):
        return self._number_of_subsets
    number_of_subsets = property(_get_number_of_subsets)
    
    def _get_subsets(self):
        self._length_update_needed = 1
        return self._subsets
    subsets = property(_get_subsets)
    
    def _get_length(self):
        if self._length_update_needed == 1:
            self._update_length()
        return self._length
    length = property(_get_length)
    
    def add_subset(self, subset_object):
        assert isinstance(subset_object, Subset)
        self._subsets.append(subset_object)
        self._number_of_subsets += 1
        self._length += subset_object.size
    
    def _update_length(self):
        self._length = 0
        for ss in self.subsets:
            self._length += ss.size
        self._length_update_needed = 0
    
    def distance(self, other):
        assert self.length == other.length, "\nPartitions must have same number of sites to calculate distance.\nYou tried '%s' (%d) vs. '%s' (%d)\n" % (self.id, self.length, other.id, other.length)
        mat = []
        dim = max(self.number_of_subsets, other.number_of_subsets)
        for xsubset in self.subsets:
            row = [0] * dim
            for i, ysubset in enumerate(other.subsets):
                intersection = xsubset.indices & ysubset.indices
                row[i] = len(intersection)
            mat.append(row)
        n_to_add = dim - self.number_of_subsets
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
        return self.length - total
    
    def get_site_to_subset_indices(self):
        site_to_subset_indices = [None] * self.length
        for ss_index, subset in enumerate(self.subsets):
            for site_index in subset.indices:
                assert (site_to_subset_indices[site_index] is None)
                site_to_subset_indices[site_index] = ss_index
        return site_to_subset_indices
        
    def permuted_copy(self):
        self._number_of_permutations += 1
        permuted = Partition([], id = "%s_permuted%d" % (self.id, self._number_of_permutations))
        for n, subset in enumerate(self.subsets):
            permuted.add_subset(Subset(set(), rate=subset.rate, id=subset.id))
        global _RNG
        site_to_subset_indices = self.get_site_to_subset_indices()
        _RNG.shuffle(site_to_subset_indices)
        for i, site_to_subset_index in enumerate(site_to_subset_indices):
            assert (site_to_subset_index is not None)
            permuted.subsets[site_to_subset_index].add_index(i)
        return permuted
        
class PosteriorOfPartitions(object):
    "object with a list of Partition objects and a string id"
    def __init__(self, list_of_partition_objects = [], id = None):
        assert isinstance(list_of_partition_objects, list)
        self._partitions = []
        self._number_of_partitions = len(self._partitions)
        for p in list_of_partition_objects:
            self.add_partition(p)
        if id is not None:
            self._id = str(id)
        else:
            self._id = id
    
    def _get_id(self):
        return self._id
    
    def _set_id(self, str_name):
        self._id = str(str_name)
    id = property(_get_id, _set_id)
    
    def _get_partitions(self):
        return self._partitions
    partitions = property(_get_partitions)
    
    def _get_number_of_partitions(self):
        return self._number_of_partitions
    number_of_partitions = property(_get_number_of_partitions)
    
    def add_partition(self, partition_object):
        assert isinstance(partition_object, Partition)
        self._partitions.append(partition_object)
        self._number_of_partitions += 1
    
    def distance_matrix(self):
        lower_triangle = []
        for n, partition in enumerate(self.partitions):
            _LOG.debug("Calc distance matrix row %d\n" % n)
            row = self.partitions[:n+1]
            distances = tuple(partition.distance(j) for j in row)
            _LOG.debug("%s\n" % "\t".join([str(x) for x in distances]))
            lower_triangle.append(distances)
        return lower_triangle
    
    def median_partition(self):
        dist_mat = self.distance_matrix()
        dim = len(dist_mat)
        sum_dist = [0]*dim
        for i in xrange(dim):
            for j in xrange(i):
                row = dist_mat[i]
                element = row[j]
                sum_dist[i] += element
                sum_dist[j] += element
        min_dist = min(sum_dist)
        for n, dist in enumerate(sum_dist):
            if dist == min_dist:
                median_part = self.partitions[n]
        return median_part, min_dist
    
    def probability_closer_than_random(self, partition_object):
        sampled_closer = 0
        permuted_closer = 0
        tie = 0
        sys.stdout.write("sampled\tpermuted\n")
        for p in self.partitions:
            sampled_dist = partition_object.distance(p)
            permuted_dist = partition_object.distance(p.permuted_copy())
            sys.stdout.write("%d\t%d\n" % (sampled_dist, permuted_dist))
            if sampled_dist < permuted_dist:
                sampled_closer += 1
            elif sampled_dist > permuted_dist:
                permuted_closer += 1
            elif sampled_dist == permuted_dist:
                tie += 1
            else:
                sys.exit("There was a problem calculating the probability that posterior '%s' is closer to partition '%s' than expected by chance." % (self.id, partition_object.id))
        assert(sampled_closer + permuted_closer + tie == self.number_of_partitions)
        prob = float(sampled_closer) / float(self.number_of_partitions)
        sys.stdout.write("\n#######################################################\n")
        sys.stdout.write("number of times sampled was closer: %d\n" % sampled_closer)
        sys.stdout.write("number of times permuted was closer: %d\n" % permuted_closer)
        sys.stdout.write("number of ties: %d\n" % tie)
        sys.stdout.write("probability the posterior of partitions is closer to\ntested partition than by chance:\n")
        sys.stdout.write("%f\n" % prob)
        return prob
    
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
        _LOG.debug("reading sample %d\n" % line_num)
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
            assert (site_to_subset_ind[col_index] is None)
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
                             
                        