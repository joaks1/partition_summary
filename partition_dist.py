#!/usr/bin/env python
import sys
import re
import random
import itertools
import time
from copy import copy

import logging
logging.basicConfig(level=logging.WARNING)
_LOG = logging.getLogger('partition_dist')

_RNG = random.Random()
_SEED = int(time.time())

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
            _LOG.error("Invalid value specified for rate attribute of Subset object: %s" % str(float_rate))
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
        l = [i for i in self.indices]
        l.sort()
        s = " ".join([str(x+1) for x in l])
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
                _LOG.error("Invalid value specified for lnL attribute of Partition object: %s" % str(lnL))
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
    
    def get_rate_indices(self):
        rate_indices = [None] * self.length
        for ss_index, subset in enumerate(self.subsets):
            for site_index in subset.indices:
                assert (rate_indices[site_index] is None)
                rate_indices[site_index] = subset.rate
        return rate_indices
        
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
    
    def valid(self):
        list_all_indices = []
        for subset in self.subsets:
            for i in subset.indices:
                list_all_indices.append(i)
        true_indices = range(0, (max(list_all_indices) + 1))
        duplicates = set()
        missing = set()
        for i in true_indices:
            n = list_all_indices.count(i)
            if n > 1:
                duplicates.add(i)
            if n < 1:
                missing.add(i)
        if (len(missing) != 0) or (len(duplicates) != 0):
            raise ValidPartitionException(self.id, missing, duplicates)
            return False
        else:
            return True
        
class ValidPartitionException(Exception):
    def __init__(self, partition_id, missing_indices=set(), duplicate_indices=set(), msg=""):
        self.partition_id = partition_id
        self.missing_sites = [(x+1) for x in missing_indices]
        self.duplicate_sites = [(x+1) for x in duplicate_indices]
        self.msg = msg
    def __str__(self):
        return "Partition '%s' is invalid.\nMissing sites:\n%s\nDuplicated sites:\n%s\n%s\n" % (self.partition_id, " ".join([str(x) for x in self.missing_sites]), " ".join([str(x) for x in self.duplicate_sites]), self.msg)

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
    
    def merge(self, other_posterior):
        assert isinstance(other_posterior, PosteriorOfPartitions)
        for p in other_posterior.partitions:
            self.add_partition(p)
    
    def distance_matrix(self):
        lower_triangle = []
        for n, partition in enumerate(self.partitions):
            _LOG.info("Calc distance matrix row %d\n" % n)
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
        global _SEED
        sampled_closer = 0
        permuted_closer = 0
        tie = 0
        sys.stdout.write("seed = %d\n" % _SEED)
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
        sys.stdout.write("#######################################################\n")
        return prob
    
    def get_mb_partitions(self):
        mb_parts = []
        for part in self.partitions:
            site_to_subset_indices = part.get_site_to_subset_indices()
            rate_indices = part.get_rate_indices()
            assert len(site_to_subset_indices) == len(rate_indices)
            mb_string = "%s\t%d\t%s\t%s\n" % (part.id, part.number_of_subsets, "\t".join([str(x) for x in site_to_subset_indices]), "\t".join([str(r) for r in rate_indices]))
            mb_parts.append(mb_string)
        assert (len(mb_parts) == self.number_of_partitions)
        return mb_parts

    
def read_partitions(sampled_partitions_file, from_row=0, to_row=None, read_rates=True):
    pb_regex = r'^(\d+)\t(\d+)\t(\d+)$'
    pb_pattern = re.compile(pb_regex)
    mb_regex = r'^n\tDegree\tP\[1\]'
    mb_pattern = re.compile(mb_regex)
    nex_regex = r'^#NEXUS'
    nex_pattern = re.compile(nex_regex, re.IGNORECASE)
    
    line_iter = iter(sampled_partitions_file)
    first_line = line_iter.next().strip()
    sampled_partitions_file.seek(0)
    pb_match = pb_pattern.search(first_line)
    mb_match = mb_pattern.search(first_line)
    nex_match = nex_pattern.search(first_line)
    if pb_match:
        version, subversion, beta = pb_match.groups()
        _LOG.info("Reading partitions from output of PhyloBayes %s.%s.%s...\n" % (version, subversion, beta))
        partitions = read_pb_partitions(sampled_partitions_file, from_row, to_row, read_rates)
    elif mb_match:
        _LOG.info("Reading partitions from output of John Huelsenbeck's DPP model program...\n")
        partitions = read_mb_partitions(sampled_partitions_file, from_row, to_row, read_rates)
    elif nex_match:
        _LOG.info("Reading partition from NEXUS file '%s'...\n" % sampled_partitions_file.name)
        partitions = read_nex_partition(sampled_partitions_file)
    else:
        sys.exit("Sorry, could not recognize format of partitions file!\n")
    return partitions

def read_nex_partition(nex_file):
    charset_pattern = re.compile(r'^\s*charset\s+(\S+)\s+=(.+);\s*$', re.IGNORECASE)
    n1 = re.compile(r'^\d+$')
    n2 = re.compile(r'^\d+-\d+$')
    n3 = re.compile(r'^\d+-\d+\\3$')
    line_iter = iter(nex_file)
    partition = Partition(id = nex_file.name)
    for line_num, line in enumerate(line_iter):
        x = line.strip()
        for name, indices in charset_pattern.findall(x):
            subset = Subset(set(), id = name.strip())
            site_indices = indices.strip().split()
            for i in site_indices:
                i.strip()
                if n1.match(i):
                    subset.add_index(int(i)-1)
                elif n2.match(i):
                    i = i.split('-')
                    i = range(int(i[0]), int(i[1])+1)
                    for n in i:
                        subset.add_index(n-1)
                elif n3.match(i):
                    i = i.replace('\\3', '')
                    i = i.split('-')
                    i = range(int(i[0]), int(i[1])+1, 3)
                    for n in i:
                        subset.add_index(n-1)
                else:
                    sys.exit("Could not parse charset '%s' from file '%s'. Found '%s' in charset definition.\n" % (name, nex_file.name, i))
            partition.add_subset(subset)
    partitions = PosteriorOfPartitions([partition])
    return partitions

def read_pb_partitions(sampled_partitions_file, from_row=0, to_row=None, read_rates=True):
    site_indices_regex = r'^\d+\t\d+\t\d+\t\d+\t\d+\t\d+\t\d+\t\d+\t\d+\t\d+'
    site_indices_pattern = re.compile(site_indices_regex)
    
    posterior = PosteriorOfPartitions()
    line_iter = iter(sampled_partitions_file)
    indices_tally = 0
    sample_tally = 0
    for line_num, line in enumerate(line_iter):
        x = line.strip()
        indices_match = site_indices_pattern.search(x)
        if indices_match:
            indices_tally += 1
        if indices_tally == 2:
            indices_tally = 0
            sample_tally += 1
            if from_row and from_row > sample_tally-1:
                _LOG.info("ignoring sample %d from '%s'\n" % (sample_tally-1, sampled_partitions_file.name))
                continue
            if to_row and to_row < sample_tally-1:
                break
            _LOG.info("reading sample %d from '%s'\n" % (sample_tally-1, sampled_partitions_file.name))
            next_line = line_iter.next().strip()
            if next_line == 'rates':
                number_of_subsets = int(line_iter.next().strip())
            else:
                number_of_subsets = int(next_line)
            partition = Partition([], id=str(sample_tally))
            rates = []
            for i in xrange(number_of_subsets):
                r = line_iter.next().strip().split()
                assert len(r) == 2
                if read_rates:
                    rate = float(r[0])
                else:
                    rate = 1.0
                subset = Subset(set(), rate = rate)
                rates.append(rate)
                partition.add_subset(subset)
            assignments = line_iter.next().strip().split()
            for index, el in enumerate(assignments):
                iel = int(el)
                partition.subsets[iel].add_index(index)
                if read_rates:
                    assert partition.subsets[iel].rate == float(rates[iel])
            assert partition.number_of_subsets == len(rates) == number_of_subsets
            posterior.add_partition(partition)
    return posterior
            
            
def parse_mb_header(h):
    s = h.split()
    n = len(s)
    n_sites = (n - 2)//2
    assert(s[n_sites + 1] == "P[%d]" % n_sites)
    assert(s[n_sites + 2] == "R[1]")
    return n_sites

def read_mb_partitions(sampled_partitions_file, from_row=0, to_row=None, read_rates=True):
    """Returns a PosteriorOfPartitions object
    
    The format in the DPP output from MrBayes is 
    rep#\t#subsets\tsubset#\tsubset#\t....
    where the number of subset numbers corresponds to the number of characters
    in the data matrix.
    """
    posterior = PosteriorOfPartitions()
    line_iter = iter(sampled_partitions_file)
    header = line_iter.next()
    n_sites = parse_mb_header(header)
    for line_num, line in enumerate(line_iter):
        if from_row and from_row > line_num:
            _LOG.info("ignoring sample %d from '%s'\n" % (line_num, sampled_partitions_file.name))
            continue
        if to_row and to_row < line_num:
            break
        x = line.strip()
        s = x.split()
        name = s.pop(0)
        num_elements = int(s.pop(0))
        partition = Partition([], id = name)
        assignments = s[:n_sites]
        assignment_set = set(assignments)
        for i in xrange(len(assignment_set)):
            partition.add_subset(Subset())
        rates = s[n_sites:]
        assert(len(assignments) == len(rates))
        _LOG.info("reading sample %d from '%s'\n" % (line_num, sampled_partitions_file.name))
        for index, el in enumerate(assignments):
            iel = int(el)
            partition.subsets[iel - 1].add_index(index)
            if partition.subsets[iel - 1].rate is None:
                if read_rates:
                    partition.subsets[iel - 1].rate = float(rates[index])
                else:
                    partition.subsets[iel - 1].rate = 1.0
        assert partition.number_of_subsets == len(assignment_set) == num_elements
        posterior.add_partition(partition)
    return posterior


if __name__ == '__main__':
    from optparse import OptionParser
    usage = "Usage: %prog [options] <PARTITIONS_FILEPATH1> [<PARTITIONS_FILEPATH2> <PARTITIONS_FILEPATH3> ...]"
    parser = OptionParser(usage = usage)
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
        help="The index of the first sampled partition to read (0 by default)")
    parser.add_option("-t", "--to", dest="to_index", default=None, 
        type="int",
        help="The index of the last sampled partition to read (None by default)")
    parser.add_option("-m", "--median", dest="median", default=False, 
        action="store_true",
        help="Find the sampled partition that has the smallest distance to all of the others (the sample that is closest to being the median)")
    parser.add_option("--target", dest="target",
        action="store",
        type="string",
        help="The partitions file from which the first partition will be used to calculate the probability that the sampled partitions are closer to this target partition than expected by chance. This can be in the same format as the partitions output by PhyloBayes or John Huelsenbeck's DPP program, or can also be a NEXUS file with the partition defined with 'charset' statements (do NOT put comments in lines with charset statements).")
    (options, args) = parser.parse_args()
    
    if options.verbose:
        _LOG.setLevel(logging.INFO)
    if options.debugging:
        _LOG.setLevel(logging.DEBUG)
    if options.seed != 0:
        _SEED = options.seed
    _RNG.seed(_SEED)
    if not args:
        sys.exit("Expecting a file name for the MrBayes DPP-over rates samples")

    sampled_partitions = PosteriorOfPartitions()
    for file in args:
        partitions_filename = file
        try:
            partitions_file = open(partitions_filename, 'rU')
        except IOError as e:
            sys.exit("Cannot open file '%s':\n%s\n" % (file, e))
        partitions = read_partitions(partitions_file, options.from_index, options.to_index)
        partitions_file.close()
        sampled_partitions.merge(partitions)

    if options.debugging:
        mb_parts = sampled_partitions.get_mb_partitions()
        _LOG.debug("%s\n" % "\n".join(mb_parts))
    _LOG.info("%d partitions read from %d files\n" % (sampled_partitions.number_of_partitions, len(args)))

    if options.target:
        test_partitions_filename = options.target
        test_partitions_file = open(test_partitions_filename, 'rU')
        test_partitions = read_partitions(test_partitions_file, 0, 1, read_rates=False)
        test_partitions_file.close()
        if test_partitions.partitions[0].length != sampled_partitions.partitions[0].length:
            sys.exit("The number of sites in the test partition must be identical to the number of sites in the sampled partitions from MrBayes")
        tp = test_partitions.partitions[0]
        try:
            tp.valid()
        except ValidPartitionException as e:
            sys.exit("Target partition '%s' from file '%s' is not valid:\n%s\n" % (tp.id, options.target, e))
        sampled_partitions.probability_closer_than_random(tp)
        
    if options.median:
        med, md = sampled_partitions.median_partition()
        sys.stdout.write("\n%s\n[total distance = %d]\n" % (str(med), md))
                             
                        