import numpy as np
import pybedtools
import os


def data_file(basename):
    return os.path.join(os.path.dirname(__file__), 'test', 'data', basename)


def normalize(hit, region):
    """
    convert coords to be relative to `region`
    """
    hit.start -= region.start
    hit.stop -= region.start
    return hit


def seq(region, fasta):
    return pybedtools.BedTool.seq(
            (region.chrom, region.start, region.stop), fasta)


def numbers(region):
    """
    returns a string of "axis label"-like numbers for a region
    """
    n = []
    n2 = []
    sub = 0
    c = 0
    for i in range(len(region)):
        if i % 10 == 0:
            num = str(c)
            n.append(num)
            n2.append('|')
            space = 10 - len(num)
            c += 10
        else:
            if space > 0:
                n.append(' ')
            n2.append(' ')
            space -= 1
    return [''.join(n), ''.join(n2)]


def pwm_from_jaspar(fn):
    """
    reads JASPAR download `fn` and returns an iterator of (`header`, `matrix`)
    tuples

    naively assumes each motif is 5 lines -- header, then a row for A, C, G,
    and T.

    Returns a motility.PWM instance for each matrix
    """
    def parse_matrix_line(line):
        L = line.split()
        L = L[2:-1]
        return [float(i) for i in L]

    header = None
    f = open(fn)
    minlen = 1
    while minlen > 0:
        header = f.readline().strip()[1:]
        a = parse_matrix_line(f.readline())
        c = parse_matrix_line(f.readline())
        g = parse_matrix_line(f.readline())
        t = parse_matrix_line(f.readline())
        a = np.array(zip(*[a, c, g, t]))
        a = a / (a.sum(axis=0))
        if len(a) == 0:
            break
        yield header, a
        minlen = min([len(i) for i in [a, c, g, t]])
