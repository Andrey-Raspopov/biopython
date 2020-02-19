# Copyright 2007-2016 by Peter Cock.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.SeqIO support for the "genbank" and "embl" file formats.

You are expected to use this module via the Bio.SeqIO functions.
Note that internally this module calls Bio.GenBank to do the actual
parsing of GenBank, EMBL and IMGT files.

See Also:
International Nucleotide Sequence Database Collaboration
http://www.insdc.org/

GenBank
http://www.ncbi.nlm.nih.gov/Genbank/

EMBL Nucleotide Sequence Database
http://www.ebi.ac.uk/embl/

DDBJ (DNA Data Bank of Japan)
http://www.ddbj.nig.ac.jp/

IMGT (use a variant of EMBL format with longer feature indents)
http://imgt.cines.fr/download/LIGM-DB/userman_doc.html
http://imgt.cines.fr/download/LIGM-DB/ftable_doc.html
http://www.ebi.ac.uk/imgt/hla/docs/manual.html

"""
import logging
import re

from Bio.SeqIO.InsdcIO import parse_feature
from Bio.SeqIO.EmblIO import EmblWriter, EmblScanner


class ImgtWriter(EmblWriter):
    """IMGT writer (EMBL format variant)."""

    HEADER_WIDTH = 5
    QUALIFIER_INDENT = 25  # Not 21 as in EMBL
    QUALIFIER_INDENT_STR = "FT" + " " * (QUALIFIER_INDENT - 2)
    QUALIFIER_INDENT_TMP = "FT   %s                    "  # 25 if %s is empty
    FEATURE_HEADER = "FH   Key                 Location/Qualifiers\nFH\n"


def ImgtIterator(source):
    """Break up an IMGT file into SeqRecord objects.

    Argument source is a file-like object opened in text mode or a path to a file.
    Every section from the LOCUS line to the terminating // becomes
    a single SeqRecord with associated annotation and features.

    Note that for genomes or chromosomes, there is typically only
    one record.
    """
    try:
        handle = open(source)
    except TypeError:
        handle = source
        if handle.read(0) != "":
            raise ValueError("IMGT files must be opened in text mode.") from None

    try:
        records = ImgtScanner().parse_records(handle)
        yield from records
    finally:
        if handle is not source:
            handle.close()


class ImgtScanner(EmblScanner):
    """For extracting chunks of information in IMGT (EMBL like) files (PRIVATE).

    IMGT files are like EMBL files but in order to allow longer feature types
    the features should be indented by 25 characters not 21 characters. In
    practice the IMGT flat files tend to use either 21 or 25 characters, so we
    must cope with both.

    This is private to encourage use of Bio.SeqIO rather than Bio.GenBank.
    """

    FEATURE_START_MARKERS = [
        "FH   Key             Location/Qualifiers",
        "FH   Key             Location/Qualifiers (from EMBL)",
        "FH   Key                 Location/Qualifiers",
        "FH",
    ]

    def _feed_first_line(self, consumer, line):
        assert line[: self.HEADER_WIDTH].rstrip() == "ID"
        if line[self.HEADER_WIDTH :].count(";") != 5:
            # Assume its an older EMBL-like line,
            return EmblScanner._feed_first_line(self, consumer, line)
        # Otherwise assume its the new (circa 2016) IMGT style
        # as used in the IPD-IMGT/HLA Database
        #
        # https://github.com/ANHIG/IMGTHLA/
        #
        # The key changes post 3.16 are the addition of an SV value
        # to the ID line, these additions should make the format more
        # similar to the ENA style.
        #
        # ID   HLA00001   standard; DNA; HUM; 3503 BP.
        #
        # becomes
        #
        # ID   HLA00001; SV 1; standard; DNA; HUM; 3503 BP.
        fields = [data.strip() for data in line[self.HEADER_WIDTH :].strip().split(";")]
        assert len(fields) == 6
        """
        The tokens represent:

           0. Primary accession number (eg 'HLA00001')
           1. Sequence version number (eg 'SV 1')
           2. ??? eg 'standard'
           3. Molecule type (e.g. 'DNA')
           4. Taxonomic division (e.g. 'HUM')
           5. Sequence length (e.g. '3503 BP.')
        """
        consumer.locus(fields[0])

        # See TODO on the EMBL _feed_first_line_new about version field
        version_parts = fields[1].split()
        if (
            len(version_parts) == 2
            and version_parts[0] == "SV"
            and version_parts[1].isdigit()
        ):
            consumer.version_suffix(version_parts[1])

        consumer.residue_type(fields[3])
        if "circular" in fields[3]:
            consumer.topology("circular")
            consumer.molecule_type(fields[3].replace("circular", "").strip())
        elif "linear" in fields[3]:
            consumer.topology("linear")
            consumer.molecule_type(fields[3].replace("linear", "").strip())
        else:
            consumer.molecule_type(fields[3].strip())
        consumer.data_file_division(fields[4])
        self._feed_seq_length(consumer, fields[5])

    def parse_features(self, skip=False):
        """Return list of tuples for the features (if present).

        Each feature is returned as a tuple (key, location, qualifiers)
        where key and location are strings (e.g. "CDS" and
        "complement(join(490883..490885,1..879))") while qualifiers
        is a list of two string tuples (feature qualifier keys and values).

        Assumes you have already read to the start of the features table.
        """
        if self.line.rstrip() not in self.FEATURE_START_MARKERS:
            logging.debug("Didn't find any feature table")
            return []

        while self.line.rstrip() in self.FEATURE_START_MARKERS:
            self.line = self.handle.readline()

        bad_position_re = re.compile(r"([0-9]+)>")

        features = []
        line = self.line
        while True:
            assert line, "Premature end of line during features table"
            if line[: self.HEADER_WIDTH].rstrip() in self.SEQUENCE_HEADERS:
                logging.debug("Found start of sequence")
                break
            line = line.rstrip()
            assert line != "//", "Premature end of features table, marker '//' found"
            if line in self.FEATURE_END_MARKERS:
                logging.debug("Found end of features")
                line = self.handle.readline()
                break
            if line[2 : self.FEATURE_QUALIFIER_INDENT].strip() == "":
                # This is an empty feature line between qualifiers. Empty
                # feature lines within qualifiers are handled below (ignored).
                line = self.handle.readline()
                continue

            if skip:
                line = self.handle.readline()
                while (
                    line[: self.FEATURE_QUALIFIER_INDENT]
                    == self.FEATURE_QUALIFIER_SPACER
                ):
                    line = self.handle.readline()
            else:
                assert line[:2] == "FT"
                try:
                    feature_key, location_start = line[2:].strip().split()
                except ValueError:
                    feature_key = line[2:25].strip()
                    location_start = line[25:].strip()
                feature_lines = [location_start]
                line = self.handle.readline()
                while (
                    line[: self.FEATURE_QUALIFIER_INDENT]
                    == self.FEATURE_QUALIFIER_SPACER
                    or line.rstrip() == ""
                ):  # cope with blank lines in the midst of a feature
                    # Use strip to remove any harmless trailing white space AND and leading
                    # white space (copes with 21 or 26 indents and other variants)
                    assert line[:2] == "FT"
                    feature_lines.append(line[self.FEATURE_QUALIFIER_INDENT :].strip())
                    line = self.handle.readline()
                feature_key, location, qualifiers = parse_feature(
                    feature_key, feature_lines
                )
                # Try to handle known problems with IMGT locations here:
                if ">" in location:
                    # Nasty hack for common IMGT bug, should be >123 not 123>
                    # in a location string. At least here the meaning is clear,
                    # and since it is so common I don't want to issue a warning
                    # warnings.warn("Feature location %s is invalid, "
                    #              "moving greater than sign before position"
                    #              % location, BiopythonParserWarning)
                    location = bad_position_re.sub(r">\1", location)
                features.append((feature_key, location, qualifiers))
        self.line = line
        return features
