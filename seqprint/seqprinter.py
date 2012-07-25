"""
The BasePrinter class simply prints a sequence (optionally wrapped) along with
numbers and tick marks.

This is implemented by building a list of "track functions", each of which
accepts a single pybedtools.Interval as its only argument and returns either
a string or an iterable of strings.

The BasePrinter's `trackfuncs` list contains [self.header, self.numbers,
self.seq].  Subclasses can append additional functions to `trackfuncs`.  If you
need to pass other data to the function, then either add the info to the region
(say, in the `name` attribute) or make the function a method on the subclass so
it can access the instance.

The function can expect the `current_seq` attribute to be present, which
represents the sequence for the region passed to the function.

For example, here's a subclass that adds a "." wherever there is a "G" in the
sequence::

    class DummyExample(BasePrinter):
        def __init__(self, regions, genome_fasta, letter='G'):
            super(DummyExample, self).__init__(regions, genome_fasta)
            self.letter = letter

            # append the new method to trackfuncs
            self.trackfuncs.append(self.g_dots)

        def g_dots(self, region):
            s = []
            for i in self.current_seq:
                if i == self.letter:
                    s.append('.')
                else:
                    s.append(' ')
            return ''.join(s)


Usage::

    >>> d = DummyExample(bedfile, fasta)
    >>> d.printseq()

"""

import motility
import itertools
import pybedtools
import helpers


class BasePrinter(object):
    def __init__(self, regions, genome_fasta):
        """
        Handles printing sequences from BED files

        See module help for subclassing hints.

        :param regions:
            A BED file or other iterable of pybedtools.Interval objects

        :param genome_fasta:
            FASTA file from which sequences for `regions` will be extracted
        """
        self.regions = pybedtools.BedTool(regions).saveas()
        self.fasta = genome_fasta
        self.trackfuncs = [self.header, self.numbers, self.seq]

    def header(self, region):
        """
        returns region name and genomic coords
        """
        return '%s %s:%s-%s' \
                % (region.name, region.chrom, region.start, region.stop)

    def numbers(self, region):
        """
        returns numbers and ticks for the region, starting at position 0
        """
        return helpers.numbers(region)

    def seq(self, region):
        """
        returns the sequence for a region
        """
        return self.current_seq

    def printseq(self, wrap=80):
        s = [""]
        for region in self.regions:
            self.current_seq = helpers.seq(region, self.fasta)
            block = []
            # each func's signature is FUNC(region)
            for func in self.trackfuncs:
                track = func(region)

                if isinstance(track, basestring):
                    block.append(track)
                else:
                    for i in track:
                        block.append(i)

            for start in range(0, len(region), wrap):
                stop = start + wrap
                for line in block:
                    s.append(line[start:start + wrap])
            s.append("")
        return '\n'.join(s)


class MotifPrinter(BasePrinter):
    def __init__(self, regions, genome_fasta, jaspar_file=None,
            jaspar_thresh=9999, annotations=None, motif_positions=None):
        """
        Adds motif tracks to BasePrinter, using motility and a file containing
        a JASPAR-format definition of a motif.

        :param regions:
            An iterable of pybedtools.Interval objects

        :param genome_fasta:
            FASTA file from which sequences for `regions` will be extracted

        :param jaspar_file:
            If provided, a file in JASPAR format.  Motifs in each sequence will
            be identified

        :param jaspar_thresh:
            Score threshold below which motifs will be ignored.

        :param motif_positions:
            If this is a list of integer indexes, these positions will be
            converted to uppercase.
        """
        super(MotifPrinter, self).__init__(regions=regions,
                genome_fasta=genome_fasta)

        pwm = list(helpers.pwm_from_jaspar(jaspar_file))
        assert len(pwm) == 1
        self.pwm = motility.PWM(pwm[0][1])
        self.jaspar_thresh = jaspar_thresh

        if motif_positions is None:
            motif_positions = []
        self.motif_positions = motif_positions

        self._annotations = {}
        if annotations:
            for k, v in annotations.items():
                self._annotations[k] = pybedtools.BedTool(v).saveas()

        self.trackfuncs.append(self.motifs)
        self.trackfuncs.append(self.annotations)

    def annotations(self, region):
        x = pybedtools.BedTool([region]).saveas()
        for symbol, annot in self._annotations.items():
            for hit in annot.intersect(x):
                hit = helpers.normalize(hit, region)
                match = ['.' * hit.start,
                        symbol * len(hit),
                        '.' * (len(region) - hit.stop)]
                yield ''.join(match)

    def motifs(self, region):
        def key(x):
            return self.pwm.calc_score(x.name)

        hits = sorted(self.motifs_in_region(region), key=key, reverse=True)

        for hit in hits:
            score = ' (%.2f %s)' % (self.pwm.calc_score(hit.name), hit.strand)

            if self.motif_positions:
                name = list(hit.name.lower())
                if hit.strand == '-':
                    pos = [-i - 1 for i in self.motif_positions]
                else:
                    pos = self.motif_positions
                for p in pos:
                    name[p] = name[p].upper()
                hit.name = ''.join(name)

            match = ['.' * hit.start,
                    hit.name,
                    score,
                    '.' * (len(region) - hit.stop - len(score)),
            ]
            yield ''.join(match)

    def _hit_to_interval(self, hit, region):
        start, stop, strand, seq = hit
        if strand == 1:
            strand = "+"
        else:
            strand = "-"
        return pybedtools.create_interval_from_list([
            region.chrom, str(region.start + start),
            str(region.start + stop), seq, '0', strand])

    def motifs_in_region(self, region):
        seq = self.current_seq
        for hit in self.pwm.find(seq, threshold=self.jaspar_thresh):
            yield helpers.normalize(self._hit_to_interval(hit, region), region)


class DummyExample(BasePrinter):
    def __init__(self, regions, genome_fasta, letter='G'):
        super(DummyExample, self).__init__(regions, genome_fasta)
        self.letter = letter
        self.trackfuncs.append(self.g_dots)

    def g_dots(self, region):
        seq = self.current_seq
        s = []
        for i in seq:
            if i == self.letter:
                s.append('.')
            else:
                s.append(' ')
        return ''.join(s)


if __name__ == "__main__":
    from helpers import data_file
    bedfile = data_file('regions.bed')
    fasta = data_file('chr11_subset.fa')
    annotations = data_file('CTCF_ENCODE_subset.bed')
    jaspar = data_file('ctcf.jaspar')
    p = BasePrinter(bedfile, fasta)
    p = MotifPrinter(bedfile, fasta, jaspar_file=jaspar, jaspar_thresh=1.5,
            annotations={"=": annotations},
            motif_positions=[3, 4, 5, 8, 9, 10, 12, 13, 14])
    print p.printseq()
