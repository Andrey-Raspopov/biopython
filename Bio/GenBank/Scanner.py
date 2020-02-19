# Copyright 2007-2017 by Peter Cock.  All rights reserved.
# Revisions copyright 2010 by Uri Laserson.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Internal code for parsing GenBank and EMBL files (PRIVATE).

This code is NOT intended for direct use.  It provides a basic scanner
(for use with a event consumer such as Bio.GenBank._FeatureConsumer)
to parse a GenBank or EMBL file (with their shared INSDC feature table).

It is used by Bio.GenBank to parse GenBank files
It is also used by Bio.SeqIO to parse GenBank and EMBL files

Feature Table Documentation:

- http://www.insdc.org/files/feature_table.html
- http://www.ncbi.nlm.nih.gov/projects/collab/FT/index.html
- ftp://ftp.ncbi.nih.gov/genbank/docs/
"""
# 17-MAR-2009: added wgs, wgs_scafld for GenBank whole genome shotgun master records.
# These are GenBank files that summarize the content of a project, and provide lists of
# scaffold and contig files in the project. These will be in annotations['wgs'] and
# annotations['wgs_scafld']. These GenBank files do not have sequences. See
# http://groups.google.com/group/bionet.molbio.genbank/browse_thread/thread/51fb88bf39e7dc36
# http://is.gd/nNgk
# for more details of this format, and an example.
# Added by Ying Huang & Iddo Friedberg
import logging
import warnings
import re
import sys
from collections import OrderedDict

from Bio.File import as_handle
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
from Bio import BiopythonParserWarning


class InsdcScanner:
    """Basic functions for breaking up a GenBank/EMBL file into sub sections.

    The International Nucleotide Sequence Database Collaboration (INSDC)
    between the DDBJ, EMBL, and GenBank.  These organisations all use the
    same "Feature Table" layout in their plain text flat file formats.

    However, the header and sequence sections of an EMBL file are very
    different in layout to those produced by GenBank/DDBJ.
    """

    # These constants get redefined with sensible values in the sub classes:
    RECORD_START = "XXX"  # "LOCUS       " or "ID   "
    HEADER_WIDTH = 3  # 12 or 5
    FEATURE_START_MARKERS = ["XXX***FEATURES***XXX"]
    FEATURE_END_MARKERS = ["XXX***END FEATURES***XXX"]
    FEATURE_QUALIFIER_INDENT = 0
    FEATURE_QUALIFIER_SPACER = ""
    SEQUENCE_HEADERS = ["XXX"]  # with right hand side spaces removed

    def __init__(self):
        """Initialize."""
        assert len(self.RECORD_START) == self.HEADER_WIDTH
        for marker in self.SEQUENCE_HEADERS:
            assert marker == marker.rstrip()
        assert len(self.FEATURE_QUALIFIER_SPACER) == self.FEATURE_QUALIFIER_INDENT
        self.handle = None
        self.line = None

    def set_handle(self, handle):
        """Set the handle attribute."""
        self.handle = handle
        self.line = ""

    def find_start(self):
        """Read in lines until find the ID/LOCUS line, which is returned.

        Any preamble (such as the header used by the NCBI on ``*.seq.gz`` archives)
        will we ignored.
        """
        while True:
            if self.line:
                line = self.line
                self.line = ""
            else:
                line = self.handle.readline()
            if not line:
                logging.debug("End of file")
                return None
            assert not isinstance(
                line[0], int
            ), "Is this handle in binary mode not text mode?"
            if line[: self.HEADER_WIDTH] == self.RECORD_START:
                logging.debug("Found the start of a record:\n" + line)
                break
            line = line.rstrip()
            if line == "//":
                logging.debug("Skipping // marking end of last record")
            elif line == "":
                logging.debug("Skipping blank line before record")
            else:
                logging.debug("Skipping header line before record:\n" + line)
        self.line = line
        return line

    def parse_header(self):
        """Return list of strings making up the header.

        New line characters are removed.

        Assumes you have just read in the ID/LOCUS line.
        """
        assert (
            self.line[: self.HEADER_WIDTH] == self.RECORD_START
        ), "Not at start of record"

        header_lines = []
        while True:
            line = self.handle.readline()
            assert line, "Premature end of line during sequence data"
            line = line.rstrip()
            if line in self.FEATURE_START_MARKERS:
                logging.debug("Found feature table")
                break
            if line[: self.HEADER_WIDTH].rstrip() in self.SEQUENCE_HEADERS:
                logging.debug("Found start of sequence")
                break
            assert line != "//", "Premature end of sequence data marker '//' found"
            header_lines.append(line)
        self.line = line
        return header_lines

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
            if len(line) < self.FEATURE_QUALIFIER_INDENT:
                warnings.warn(
                    "line too short to contain a feature: %r" % line,
                    BiopythonParserWarning,
                )
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
                # Build up a list of the lines making up this feature:
                if (
                    line[self.FEATURE_QUALIFIER_INDENT] != " "
                    and " " in line[self.FEATURE_QUALIFIER_INDENT :]
                ):
                    # The feature table design enforces a length limit on the feature keys.
                    # Some third party files (e.g. IGMT's EMBL like files) solve this by
                    # over indenting the location and qualifiers.
                    feature_key, line = line[2:].strip().split(None, 1)
                    feature_lines = [line]
                    warnings.warn(
                        "Over indented %s feature?" % feature_key,
                        BiopythonParserWarning,
                    )
                else:
                    feature_key = line[2 : self.FEATURE_QUALIFIER_INDENT].strip()
                    feature_lines = [line[self.FEATURE_QUALIFIER_INDENT :]]
                line = self.handle.readline()
                while line[
                    : self.FEATURE_QUALIFIER_INDENT
                ] == self.FEATURE_QUALIFIER_SPACER or (
                    line != "" and line.rstrip() == ""
                ):  # cope with blank lines in the midst of a feature
                    # Use strip to remove any harmless trailing white space AND and leading
                    # white space (e.g. out of spec files with too much indentation)
                    feature_lines.append(line[self.FEATURE_QUALIFIER_INDENT :].strip())
                    line = self.handle.readline()
                features.append(self.parse_feature(feature_key, feature_lines))
        self.line = line
        return features

    def parse_feature(self, feature_key, lines):
        r"""Parse a feature given as a list of strings into a tuple.

        Expects a feature as a list of strings, returns a tuple (key, location,
        qualifiers)

        For example given this GenBank feature::

             CDS             complement(join(490883..490885,1..879))
                             /locus_tag="NEQ001"
                             /note="conserved hypothetical [Methanococcus jannaschii];
                             COG1583:Uncharacterized ACR; IPR001472:Bipartite nuclear
                             localization signal; IPR002743: Protein of unknown
                             function DUF57"
                             /codon_start=1
                             /transl_table=11
                             /product="hypothetical protein"
                             /protein_id="NP_963295.1"
                             /db_xref="GI:41614797"
                             /db_xref="GeneID:2732620"
                             /translation="MRLLLELKALNSIDKKQLSNYLIQGFIYNILKNTEYSWLHNWKK
                             EKYFNFTLIPKKDIIENKRYYLIISSPDKRFIEVLHNKIKDLDIITIGLAQFQLRKTK
                             KFDPKLRFPWVTITPIVLREGKIVILKGDKYYKVFVKRLEELKKYNLIKKKEPILEEP
                             IEISLNQIKDGWKIIDVKDRYYDFRNKSFSAFSNWLRDLKEQSLRKYNNFCGKNFYFE
                             EAIFEGFTFYKTVSIRIRINRGEAVYIGTLWKELNVYRKLDKEEREFYKFLYDCGLGS
                             LNSMGFGFVNTKKNSAR"

        Then should give input key="CDS" and the rest of the data as a list of strings
        lines=["complement(join(490883..490885,1..879))", ..., "LNSMGFGFVNTKKNSAR"]
        where the leading spaces and trailing newlines have been removed.

        Returns tuple containing: (key as string, location string, qualifiers as list)
        as follows for this example:

        key = "CDS", string
        location = "complement(join(490883..490885,1..879))", string
        qualifiers = list of string tuples:

        [('locus_tag', '"NEQ001"'),
         ('note', '"conserved hypothetical [Methanococcus jannaschii];\nCOG1583:..."'),
         ('codon_start', '1'),
         ('transl_table', '11'),
         ('product', '"hypothetical protein"'),
         ('protein_id', '"NP_963295.1"'),
         ('db_xref', '"GI:41614797"'),
         ('db_xref', '"GeneID:2732620"'),
         ('translation', '"MRLLLELKALNSIDKKQLSNYLIQGFIYNILKNTEYSWLHNWKK\nEKYFNFT..."')]

        In the above example, the "note" and "translation" were edited for compactness,
        and they would contain multiple new line characters (displayed above as \n)

        If a qualifier is quoted (in this case, everything except codon_start and
        transl_table) then the quotes are NOT removed.

        Note that no whitespace is removed.
        """
        iterator = (x for x in lines if x)
        try:
            line = next(iterator)

            feature_location = line.strip()
            while feature_location[-1:] == ",":
                line = next(iterator)
                feature_location += line.strip()
            if feature_location.count("(") > feature_location.count(")"):
                warnings.warn(
                    "Non-standard feature line wrapping (didn't break on comma)?",
                    BiopythonParserWarning,
                )
                while feature_location[-1:] == "," or feature_location.count(
                    "("
                ) > feature_location.count(")"):
                    line = next(iterator)
                    feature_location += line.strip()

            qualifiers = []

            for line_number, line in enumerate(iterator):
                if line_number == 0 and line.startswith(")"):
                    feature_location += line.strip()
                elif line[0] == "/":
                    i = line.find("=")
                    key = line[1:i]
                    value = line[i + 1 :]
                    if i and value.startswith(" ") and value.lstrip().startswith('"'):
                        warnings.warn(
                            "White space after equals in qualifier",
                            BiopythonParserWarning,
                        )
                        value = value.lstrip()
                    if i == -1:
                        # Qualifier with no key, e.g. /pseudo
                        key = line[1:]
                        qualifiers.append((key, None))
                    elif not value:
                        # ApE can output /note=
                        qualifiers.append((key, ""))
                    elif value == '"':
                        logging.debug("Single quote %s:%s" % (key, value))
                        qualifiers.append((key, value))
                    elif value[0] == '"':
                        value_list = [value]
                        while value_list[-1][-1] != '"':
                            value_list.append(next(iterator))
                        value = "\n".join(value_list)
                        qualifiers.append((key, value))
                    else:
                        qualifiers.append((key, value))
                else:
                    assert len(qualifiers) > 0
                    assert key == qualifiers[-1][0]
                    assert qualifiers[-1][1], StopIteration
                    qualifiers[-1] = (key, qualifiers[-1][1] + "\n" + line)
            return feature_key, feature_location, qualifiers
        except StopIteration:
            raise ValueError(
                "Problem with '%s' feature:\n%s" % (feature_key, "\n".join(lines))
            ) from None

    def parse_footer(self):
        """Return a tuple containing a list of any misc strings, and the sequence."""
        if self.line in self.FEATURE_END_MARKERS:
            while self.line[: self.HEADER_WIDTH].rstrip() not in self.SEQUENCE_HEADERS:
                self.line = self.handle.readline()
                assert self.line, "Premature end of file"
                self.line = self.line.rstrip()

        assert (
            self.line[: self.HEADER_WIDTH].rstrip() in self.SEQUENCE_HEADERS
        ), "Not at start of sequence"
        while True:
            line = self.handle.readline()
            assert line, "Premature end of line during sequence data"
            line = line.rstrip()
            if line == "//":
                break
        self.line = line
        return [], ""  # Dummy values!

    def _feed_first_line(self, consumer, line):
        """Handle the LOCUS/ID line, passing data to the consumer (PRIVATE).

        This should be implemented by the EMBL / GenBank specific subclass

        Used by the parse_records() and parse() methods.
        """
        pass

    def _feed_header_lines(self, consumer, lines):
        """Handle the header lines (list of strings), passing data to the consumer (PRIVATE).

        This should be implemented by the EMBL / GenBank specific subclass

        Used by the parse_records() and parse() methods.
        """
        pass

    @staticmethod
    def _feed_feature_table(consumer, feature_tuples):
        """Handle the feature table (list of tuples), passing data to the consumer (PRIVATE).

        Used by the parse_records() and parse() methods.
        """
        consumer.start_feature_table()
        for feature_key, location_string, qualifiers in feature_tuples:
            consumer.feature_key(feature_key)
            consumer.location(location_string)
            for q_key, q_value in qualifiers:
                if q_value is None:
                    consumer.feature_qualifier(q_key, q_value)
                else:
                    consumer.feature_qualifier(q_key, q_value.replace("\n", " "))

    def _feed_misc_lines(self, consumer, lines):
        """Handle any lines between features and sequence (list of strings), passing data to the consumer (PRIVATE).

        This should be implemented by the EMBL / GenBank specific subclass

        Used by the parse_records() and parse() methods.
        """
        pass

    def feed(self, handle, consumer, do_features=True):
        """Feed a set of data into the consumer.

        This method is intended for use with the "old" code in Bio.GenBank

        Arguments:
         - handle - A handle with the information to parse.
         - consumer - The consumer that should be informed of events.
         - do_features - Boolean, should the features be parsed?
           Skipping the features can be much faster.

        Return values:
         - true  - Passed a record
         - false - Did not find a record

        """
        # Should work with both EMBL and GenBank files provided the
        # equivalent Bio.GenBank._FeatureConsumer methods are called...
        self.set_handle(handle)
        if not self.find_start():
            # Could not find (another) record
            consumer.data = None
            return False

        # We use the above class methods to parse the file into a simplified format.
        # The first line, header lines and any misc lines after the features will be
        # dealt with by GenBank / EMBL specific derived classes.

        # First line and header:
        self._feed_first_line(consumer, self.line)
        self._feed_header_lines(consumer, self.parse_header())

        # Features (common to both EMBL and GenBank):
        if do_features:
            self._feed_feature_table(consumer, self.parse_features(skip=False))
        else:
            self.parse_features(skip=True)  # ignore the data

        # Footer and sequence
        misc_lines, sequence_string = self.parse_footer()
        self._feed_misc_lines(consumer, misc_lines)

        consumer.sequence(sequence_string)
        # Calls to consumer.base_number() do nothing anyway
        consumer.record_end("//")

        assert self.line == "//"

        # And we are done
        return True

    def parse(self, handle, do_features=True):
        """Return a SeqRecord (with SeqFeatures if do_features=True).

        See also the method parse_records() for use on multi-record files.
        """
        from Bio.GenBank import _FeatureConsumer
        from Bio.GenBank.utils import FeatureValueCleaner

        consumer = _FeatureConsumer(
            use_fuzziness=1, feature_cleaner=FeatureValueCleaner()
        )

        if self.feed(handle, consumer, do_features):
            return consumer.data
        else:
            return None

    def parse_records(self, handle, do_features=True):
        """Parse records, return a SeqRecord object iterator.

        Each record (from the ID/LOCUS line to the // line) becomes a SeqRecord

        The SeqRecord objects include SeqFeatures if do_features=True

        This method is intended for use in Bio.SeqIO
        """
        # This is a generator function
        with as_handle(handle) as handle:
            while True:
                record = self.parse(handle, do_features)
                if record is None:
                    break
                assert record.id, "Failed to parse the record's ID. Invalid ID line?"
                assert (
                    record.name != "<unknown name>"
                ), "Failed to parse the record's name. Invalid ID line?"
                assert (
                    record.description != "<unknown description>"
                ), "Failed to parse the record's description"
                yield record

    def parse_cds_features(
        self,
        handle,
        alphabet=generic_protein,
        tags2id=("protein_id", "locus_tag", "product"),
    ):
        """Parse CDS features, return SeqRecord object iterator.

        Each CDS feature becomes a SeqRecord.

        Arguments:
         - alphabet - Used for any sequence found in a translation field.
         - tags2id  - Tuple of three strings, the feature keys to use
           for the record id, name and description,

        This method is intended for use in Bio.SeqIO

        """
        with as_handle(handle) as handle:
            self.set_handle(handle)
            while self.find_start():
                # Got an EMBL or GenBank record...
                self.parse_header()  # ignore header lines!
                feature_tuples = self.parse_features()
                # self.parse_footer() # ignore footer lines!
                while True:
                    line = self.handle.readline()
                    if not line:
                        break
                    if line[:2] == "//":
                        break
                self.line = line.rstrip()

                # Now go though those features...
                for key, location_string, qualifiers in feature_tuples:
                    if key == "CDS":
                        # Create SeqRecord
                        # ================
                        # SeqRecord objects cannot be created with annotations, they
                        # must be added afterwards.  So create an empty record and
                        # then populate it:
                        record = SeqRecord(seq=None)
                        annotations = record.annotations

                        # Should we add a location object to the annotations?
                        # I *think* that only makes sense for SeqFeatures with their
                        # sub features...
                        annotations["raw_location"] = location_string.replace(" ", "")

                        for (qualifier_name, qualifier_data) in qualifiers:
                            if (
                                qualifier_data is not None
                                and qualifier_data[0] == '"'
                                and qualifier_data[-1] == '"'
                            ):
                                # Remove quotes
                                qualifier_data = qualifier_data[1:-1]
                            # Append the data to the annotation qualifier...
                            if qualifier_name == "translation":
                                assert record.seq is None, "Multiple translations!"
                                record.seq = Seq(
                                    qualifier_data.replace("\n", ""), alphabet
                                )
                            elif qualifier_name == "db_xref":
                                # its a list, possibly empty.  Its safe to extend
                                record.dbxrefs.append(qualifier_data)
                            else:
                                if qualifier_data is not None:
                                    qualifier_data = qualifier_data.replace(
                                        "\n", " "
                                    ).replace("  ", " ")
                                try:
                                    annotations[qualifier_name] += " " + qualifier_data
                                except KeyError:
                                    # Not an addition to existing data, its the first bit
                                    annotations[qualifier_name] = qualifier_data

                        # Fill in the ID, Name, Description
                        # =================================
                        try:
                            record.id = annotations[tags2id[0]]
                        except KeyError:
                            pass
                        try:
                            record.name = annotations[tags2id[1]]
                        except KeyError:
                            pass
                        try:
                            record.description = annotations[tags2id[2]]
                        except KeyError:
                            pass

                        yield record


class EmblScanner(InsdcScanner):
    """For extracting chunks of information in EMBL files."""

    RECORD_START = "ID   "
    HEADER_WIDTH = 5
    FEATURE_START_MARKERS = ["FH   Key             Location/Qualifiers", "FH"]
    FEATURE_END_MARKERS = ["XX"]  # XX can also mark the end of many things!
    FEATURE_QUALIFIER_INDENT = 21
    FEATURE_QUALIFIER_SPACER = "FT" + " " * (FEATURE_QUALIFIER_INDENT - 2)
    SEQUENCE_HEADERS = ["SQ", "CO"]  # Remove trailing spaces

    EMBL_INDENT = HEADER_WIDTH
    EMBL_SPACER = " " * EMBL_INDENT

    def parse_footer(self):
        """Return a tuple containing a list of any misc strings, and the sequence."""
        assert self.line[: self.HEADER_WIDTH].rstrip() in self.SEQUENCE_HEADERS, (
            "Footer format unexpected: '%s'" % self.line
        )

        misc_lines = []
        while self.line[: self.HEADER_WIDTH].rstrip() in self.SEQUENCE_HEADERS:
            misc_lines.append(self.line)
            self.line = self.handle.readline()
            assert self.line, "Premature end of file"
            self.line = self.line.rstrip()

        assert (
            self.line[: self.HEADER_WIDTH] == " " * self.HEADER_WIDTH
            or self.line.strip() == "//"
        ), ("Unexpected content after SQ or CO line: %r" % self.line)

        seq_lines = []
        line = self.line
        while True:
            assert line, "Premature end of file in sequence data"
            line = line.strip()
            assert line, "Blank line in sequence data"
            if line == "//":
                break
            assert self.line[: self.HEADER_WIDTH] == (" " * self.HEADER_WIDTH), (
                "Problem with characters in header line, or incorrect header width: "
                + self.line
            )
            linersplit = line.rsplit(None, 1)
            if len(linersplit) == 2 and linersplit[1].isdigit():
                seq_lines.append(linersplit[0])
            elif line.isdigit():
                # Special case of final blank line with no bases
                # just the sequence coordinate
                pass
            else:
                warnings.warn(
                    "EMBL sequence line missing coordinates", BiopythonParserWarning
                )
                seq_lines.append(line)
            line = self.handle.readline()
        self.line = line
        return misc_lines, "".join(seq_lines).replace(" ", "")

    def _feed_first_line(self, consumer, line):
        assert line[: self.HEADER_WIDTH].rstrip() == "ID"
        if line[self.HEADER_WIDTH :].count(";") == 6:
            # Looks like the semi colon separated style introduced in 2006
            self._feed_first_line_new(consumer, line)
        elif line[self.HEADER_WIDTH :].count(";") == 3:
            if line.rstrip().endswith(" SQ"):
                self._feed_first_line_patents(consumer, line)
            else:
                self._feed_first_line_old(consumer, line)
        elif line[self.HEADER_WIDTH :].count(";") == 2:
            self._feed_first_line_patents_kipo(consumer, line)
        else:
            raise ValueError("Did not recognise the ID line layout:\n" + line)

    def _feed_first_line_patents(self, consumer, line):
        # Old style EMBL patent records where ID line ended SQ
        # Not 100% sure that PRT here is really molecule type and
        # not the data file division...
        #
        # Either Non-Redundant Level 1 database records,
        # ID <accession>; <molecule type>; <non-redundant level 1>; <cluster size L1>
        # e.g. ID   NRP_AX000635; PRT; NR1; 15 SQ
        #
        # Or, Non-Redundant Level 2 database records:
        # ID <L2-accession>; <molecule type>; <non-redundant level 2>; <cluster size L2>
        # e.g. ID   NRP0000016E; PRT; NR2; 5 SQ
        # e.g. ID   NRP_AX000635; PRT; NR1; 15 SQ
        fields = [
            data.strip() for data in line[self.HEADER_WIDTH :].strip()[:-3].split(";")
        ]
        assert len(fields) == 4
        consumer.locus(fields[0])
        consumer.residue_type(fields[1])  # semi-redundant
        consumer.data_file_division(fields[2])
        # TODO - Record cluster size?

    def _feed_first_line_patents_kipo(self, consumer, line):
        # EMBL format patent sequence from KIPO, e.g.
        # ftp://ftp.ebi.ac.uk/pub/databases/patentdata/kipo_prt.dat.gz
        #
        # e.g. ID   DI500001       STANDARD;      PRT;   111 AA.
        #
        # This follows the style of _feed_first_line_old
        assert line[: self.HEADER_WIDTH].rstrip() == "ID"
        fields = [line[self.HEADER_WIDTH :].split(None, 1)[0]]
        fields.extend(line[self.HEADER_WIDTH :].split(None, 1)[1].split(";"))
        fields = [entry.strip() for entry in fields]
        """
        The tokens represent:

           0. Primary accession number
           (space sep)
           1. ??? (e.g. standard)
           (semi-colon)
           2. Molecule type (protein)? Division? Always 'PRT'
           3. Sequence length (e.g. '111 AA.')
        """
        consumer.locus(fields[0])  # Should we also call the accession consumer?
        # consumer.molecule_type(fields[2])
        self._feed_seq_length(consumer, fields[3])

    def _feed_first_line_old(self, consumer, line):
        # Expects an ID line in the style before 2006, e.g.
        # ID   SC10H5 standard; DNA; PRO; 4870 BP.
        # ID   BSUB9999   standard; circular DNA; PRO; 4214630 BP.
        assert line[: self.HEADER_WIDTH].rstrip() == "ID"
        fields = [line[self.HEADER_WIDTH :].split(None, 1)[0]]
        fields.extend(line[self.HEADER_WIDTH :].split(None, 1)[1].split(";"))
        fields = [entry.strip() for entry in fields]
        """
        The tokens represent:

           0. Primary accession number
           (space sep)
           1. ??? (e.g. standard)
           (semi-colon)
           2. Topology and/or Molecule type (e.g. 'circular DNA' or 'DNA')
           3. Taxonomic division (e.g. 'PRO')
           4. Sequence length (e.g. '4639675 BP.')

        """
        consumer.locus(fields[0])  # Should we also call the accession consumer?
        consumer.residue_type(fields[2])
        if "circular" in fields[2]:
            consumer.topology("circular")
            consumer.molecule_type(fields[2].replace("circular", "").strip())
        elif "linear" in fields[2]:
            consumer.topology("linear")
            consumer.molecule_type(fields[2].replace("linear", "").strip())
        else:
            consumer.molecule_type(fields[2].strip())
        consumer.data_file_division(fields[3])
        self._feed_seq_length(consumer, fields[4])

    def _feed_first_line_new(self, consumer, line):
        # Expects an ID line in the style introduced in 2006, e.g.
        # ID   X56734; SV 1; linear; mRNA; STD; PLN; 1859 BP.
        # ID   CD789012; SV 4; linear; genomic DNA; HTG; MAM; 500 BP.
        assert line[: self.HEADER_WIDTH].rstrip() == "ID"
        fields = [data.strip() for data in line[self.HEADER_WIDTH :].strip().split(";")]
        assert len(fields) == 7
        """
        The tokens represent:

           0. Primary accession number
           1. Sequence version number
           2. Topology: 'circular' or 'linear'
           3. Molecule type (e.g. 'genomic DNA')
           4. Data class (e.g. 'STD')
           5. Taxonomic division (e.g. 'PRO')
           6. Sequence length (e.g. '4639675 BP.')

        """

        consumer.locus(fields[0])

        # Call the accession consumer now, to make sure we record
        # something as the record.id, in case there is no AC line
        consumer.accession(fields[0])

        # TODO - How to deal with the version field?  At the moment the consumer
        # will try and use this for the ID which isn't ideal for EMBL files.
        version_parts = fields[1].split()
        if (
            len(version_parts) == 2
            and version_parts[0] == "SV"
            and version_parts[1].isdigit()
        ):
            consumer.version_suffix(version_parts[1])

        # Based on how the old GenBank parser worked, merge these two:
        consumer.residue_type(" ".join(fields[2:4]))  # Semi-obsolete

        consumer.topology(fields[2])
        consumer.molecule_type(fields[3])

        # consumer.xxx(fields[4]) # TODO - What should we do with the data class?

        consumer.data_file_division(fields[5])

        self._feed_seq_length(consumer, fields[6])

    @staticmethod
    def _feed_seq_length(consumer, text):
        length_parts = text.split()
        assert len(length_parts) == 2, "Invalid sequence length string %r" % text
        assert length_parts[1].upper() in ["BP", "BP.", "AA", "AA."]
        consumer.size(length_parts[0])

    def _feed_header_lines(self, consumer, lines):
        consumer_dict = {
            "AC": "accession",
            "SV": "version",  # SV line removed in June 2006, now part of ID line
            "DE": "definition",
            # 'RN' : 'reference_num',
            # 'RC' : reference comment... TODO
            # 'RP' : 'reference_bases',
            # 'RX' : reference cross reference... DOI or Pubmed
            "RG": "consrtm",  # optional consortium
            # 'RA' : 'authors',
            # 'RT' : 'title',
            "RL": "journal",
            "OS": "organism",
            "OC": "taxonomy",
            # 'DR' : data reference
            "CC": "comment",
            # 'XX' : splitter
        }
        # We have to handle the following specially:
        # RX (depending on reference type...)
        for line in lines:
            line_type = line[: self.EMBL_INDENT].strip()
            data = line[self.EMBL_INDENT :].strip()
            if line_type == "XX":
                pass
            elif line_type == "RN":
                # Reformat reference numbers for the GenBank based consumer
                # e.g. '[1]' becomes '1'
                if data[0] == "[" and data[-1] == "]":
                    data = data[1:-1]
                consumer.reference_num(data)
            elif line_type == "RP":
                if data.strip() == "[-]":
                    # Patent EMBL files from KIPO just use: RN  [-]
                    pass
                else:
                    # Reformat reference numbers for the GenBank based consumer
                    # e.g. '1-4639675' becomes '(bases 1 to 4639675)'
                    # and '160-550, 904-1055' becomes '(bases 160 to 550; 904 to 1055)'
                    # Note could be multi-line, and end with a comma
                    parts = [
                        bases.replace("-", " to ").strip()
                        for bases in data.split(",")
                        if bases.strip()
                    ]
                    consumer.reference_bases("(bases %s)" % "; ".join(parts))
            elif line_type == "RT":
                # Remove the enclosing quotes and trailing semi colon.
                # Note the title can be split over multiple lines.
                if data.startswith('"'):
                    data = data[1:]
                if data.endswith('";'):
                    data = data[:-2]
                consumer.title(data)
            elif line_type == "RX":
                # EMBL support three reference types at the moment:
                # - PUBMED    PUBMED bibliographic database (NLM)
                # - DOI       Digital Object Identifier (International DOI Foundation)
                # - AGRICOLA  US National Agriculture Library (NAL) of the US Department
                #             of Agriculture (USDA)
                #
                # Format:
                # RX  resource_identifier; identifier.
                #
                # e.g.
                # RX   DOI; 10.1016/0024-3205(83)90010-3.
                # RX   PUBMED; 264242.
                #
                # Currently our reference object only supports PUBMED and MEDLINE
                # (as these were in GenBank files?).
                key, value = data.split(";", 1)
                if value.endswith("."):
                    value = value[:-1]
                value = value.strip()
                if key == "PUBMED":
                    consumer.pubmed_id(value)
                # TODO - Handle other reference types (here and in BioSQL bindings)
            elif line_type == "CC":
                # Have to pass a list of strings for this one (not just a string)
                consumer.comment([data])
            elif line_type == "DR":
                # Database Cross-reference, format:
                # DR   database_identifier; primary_identifier; secondary_identifier.
                #
                # e.g.
                # DR   MGI; 98599; Tcrb-V4.
                #
                # TODO - How should we store any secondary identifier?
                parts = data.rstrip(".").split(";")
                # Turn it into "database_identifier:primary_identifier" to
                # mimic the GenBank parser. e.g. "MGI:98599"
                if len(parts) == 1:
                    warnings.warn(
                        "Malformed DR line in EMBL file.", BiopythonParserWarning
                    )
                else:
                    consumer.dblink("%s:%s" % (parts[0].strip(), parts[1].strip()))
            elif line_type == "RA":
                # Remove trailing ; at end of authors list
                consumer.authors(data.rstrip(";"))
            elif line_type == "PR":
                # In the EMBL patent files, this is a PR (PRiority) line which
                # provides the earliest active priority within the family.
                # The priority  number comes first, followed by the priority date.
                #
                # e.g.
                # PR   JP19990377484 16-DEC-1999
                #
                # However, in most EMBL files this is a PR (PRoject) line which
                # gives the BioProject reference number.
                #
                # e.g.
                # PR   Project:PRJNA60715;
                #
                # In GenBank files this corresponds to the old PROJECT line
                # which was later replaced with the DBLINK line.
                if data.startswith("Project:"):
                    # Remove trailing ; at end of the project reference
                    consumer.project(data.rstrip(";"))
            elif line_type == "KW":
                consumer.keywords(data.rstrip(";"))
            elif line_type in consumer_dict:
                # Its a semi-automatic entry!
                getattr(consumer, consumer_dict[line_type])(data)
            else:
                logging.debug("Ignoring EMBL header line:\n%s" % line)

    def _feed_misc_lines(self, consumer, lines):
        # TODO - Should we do something with the information on the SQ line(s)?
        lines.append("")
        line_iter = iter(lines)
        try:
            for line in line_iter:
                if line.startswith("CO   "):
                    line = line[5:].strip()
                    contig_location = line
                    while True:
                        line = next(line_iter)
                        if not line:
                            break
                        elif line.startswith("CO   "):
                            contig_location += line[5:].strip()
                        else:
                            raise "Expected CO (contig) continuation line, got:\n" + line
                    consumer.contig_location(contig_location)
                if line.startswith("SQ   Sequence "):
                    self._feed_seq_length(
                        consumer, line[14:].rstrip().rstrip(";").split(";", 1)[0]
                    )
                    # TODO - Record the checksum etc?
            return
        except StopIteration:
            raise ValueError("Problem in misc lines before sequence") from None


class _ImgtScanner(EmblScanner):
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
                feature_key, location, qualifiers = self.parse_feature(
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


def _feed_header_organism(consumer, data):
    organism_data = ""
    lineage_data = ""
    for line in data:
        if lineage_data or ";" in line:
            lineage_data += f" {line}"
        elif line != ".":
            organism_data += f" {line}"
    if organism_data == "":
        organism_data = "."
    consumer.organism(organism_data.strip())
    if lineage_data.strip() == "":
        logging.debug("Taxonomy line(s) missing or blank")
    consumer.taxonomy(lineage_data.strip())
    del organism_data, lineage_data


def _feed_header_reference(consumer, data):
    logging.debug("Found reference [" + data[0] + "]")
    for line in data:
        logging.debug("Extended reference text [" + line + "]")
    data = " ".join(data)
    while "  " in data:
        data = data.replace("  ", " ")
    if " " not in data:
        logging.debug('Reference number "' + data + '"')
        consumer.reference_num(data)
    else:
        logging.debug(
            'Reference number "'
            + data[: data.find(" ")]
            + '", "'
            + data[data.find(" ") + 1 :]
            + '"'
        )
        consumer.reference_num(data[: data.find(" ")])
        consumer.reference_bases(data[data.find(" ") + 1 :])


def _feed_header_dblink(consumer, data):
    for line in data:
        consumer.dblink(line.strip())


def _feed_header_version(consumer, data):
    while "  " in data:
        data = data.replace("  ", " ")
    if " GI:" not in data:
        consumer.version(data)
    else:
        logging.debug(
            f"Version [{data.split(' GI:')[0]}], gi [{data.split(' GI:')[1]}]"
        )
        consumer.version(data.split(" GI:")[0])
        consumer.gi(data.split(" GI:")[1])


def _feed_header_common(consumer, data, line_type, consumer_dict):
    data = " ".join(data)
    if line_type == "DEFINITION" and data.endswith("."):
        data = data[:-1]
    getattr(consumer, consumer_dict[line_type])(data)


def _feed_header_comment(consumer, data, comment_start, comment_delim, comment_end):
    logging.debug("Found comment")
    comment_list = []
    structured_comment_key = None
    structured_comment_dict = OrderedDict()
    for line in data:
        if comment_start in line:
            pattern = r"([^#]+){0}$".format(comment_start)
            structured_comment_key = re.search(pattern, line)
            if structured_comment_key is not None:
                structured_comment_key = structured_comment_key.group(1)
                logging.debug("Found Structured Comment")
            else:
                comment_list.append(line)
        elif structured_comment_key is not None and comment_delim in line:
            pattern = r"(.+?)\s*{0}\s*(.+)".format(comment_delim)
            match = re.search(pattern, line)
            structured_comment_dict.setdefault(structured_comment_key, OrderedDict())
            structured_comment_dict[structured_comment_key][
                match.group(1)
            ] = match.group(2)
            logging.debug(f"Structured Comment continuation [{line}]")
        elif structured_comment_key is not None and comment_end not in line:
            previous_value_line = structured_comment_dict[structured_comment_key][
                match.group(1)
            ]
            structured_comment_dict[structured_comment_key][match.group(1)] = (
                previous_value_line + " " + line.strip()
            )
        elif comment_end in line:
            structured_comment_key = None
        else:
            comment_list.append(line)
            logging.debug("Comment continuation [" + line + "]")
    if comment_list:
        consumer.comment(comment_list)
    if structured_comment_dict:
        consumer.structured_comment(structured_comment_dict)
    del comment_list, structured_comment_key, structured_comment_dict


def _feed_misc_common(consumer, consumer_dict, key, line):
    if line.startswith(key):
        line = line[len(key)].strip()
        if line:
            getattr(consumer, consumer_dict[key])(line)
    return line


class GenBankScanner(InsdcScanner):
    """For extracting chunks of information in GenBank files."""

    CONSUMER_DICT = {
        "DEFINITION": "definition",
        "ACCESSION": "accession",
        "NID": "nid",
        "PID": "pid",
        "DBSOURCE": "db_source",
        "KEYWORDS": "keywords",
        "SEGMENT": "segment",
        "SOURCE": "source",
        "AUTHORS": "authors",
        "CONSRTM": "consrtm",
        "PROJECT": "project",
        "TITLE": "title",
        "JOURNAL": "journal",
        "MEDLINE": "medline_id",
        "PUBMED": "pubmed_id",
        "REMARK": "remark",
    }
    RECORD_START = "LOCUS       "
    HEADER_WIDTH = 12
    FEATURE_START_MARKERS = ["FEATURES             Location/Qualifiers", "FEATURES"]
    FEATURE_END_MARKERS = []
    FEATURE_QUALIFIER_INDENT = 21
    FEATURE_QUALIFIER_SPACER = " " * FEATURE_QUALIFIER_INDENT
    SEQUENCE_HEADERS = [
        "CONTIG",
        "ORIGIN",
        "BASE COUNT",
        "WGS",
        "TSA",
        "TLS",
    ]

    GENBANK_INDENT = HEADER_WIDTH
    GENBANK_SPACER = " " * GENBANK_INDENT

    STRUCTURED_COMMENT_START = "-START##"
    STRUCTURED_COMMENT_END = "-END##"
    STRUCTURED_COMMENT_DELIM = " :: "

    def parse_footer(self):
        """Return a tuple containing a list of any misc strings, and the sequence."""
        assert self.line[: self.HEADER_WIDTH].rstrip() in self.SEQUENCE_HEADERS, (
            "Footer format unexpected:  '%s'" % self.line
        )

        misc_lines = []
        while (
            self.line[: self.HEADER_WIDTH].rstrip() in self.SEQUENCE_HEADERS
            or self.line[: self.HEADER_WIDTH] == " " * self.HEADER_WIDTH
            or "WGS" == self.line[:3]
        ):
            misc_lines.append(self.line.rstrip())
            self.line = self.handle.readline()
            assert self.line, "Premature end of file"
            self.line = self.line

        assert self.line[: self.HEADER_WIDTH].rstrip() not in self.SEQUENCE_HEADERS, (
            "Eh? '%s'" % self.line
        )

        seq_lines = []
        line = self.line
        while True:
            if not line:
                warnings.warn(
                    "Premature end of file in sequence data", BiopythonParserWarning
                )
                line = "//"
                break
            line = line.rstrip()
            if not line:
                warnings.warn("Blank line in sequence data", BiopythonParserWarning)
                line = self.handle.readline()
                continue
            if line == "//":
                break
            if line.startswith("CONTIG"):
                break
            if len(line) > 9 and line[9:10] != " ":
                # Some broken programs indent the sequence by one space too many
                # so try to get rid of that and test again.
                warnings.warn(
                    "Invalid indentation for sequence line", BiopythonParserWarning
                )
                line = line[1:]
                assert len(line) <= 9 or line[9:10] == " ", (
                    "Sequence line mal-formed, '%s'" % line
                )
            seq_lines.append(line[10:])  # remove spaces later
            line = self.handle.readline()

        self.line = line
        # Seq("".join(seq_lines), self.alphabet)
        return misc_lines, "".join(seq_lines).replace(" ", "")

    def _feed_first_line(self, consumer, line):
        assert line.startswith(self.RECORD_START), (
            "LOCUS line does not start correctly:\n" + line
        )
        pattern = (
            r"LOCUS"
            r" +([^\s]+)"
            r" *([0-9]+)?"
            r" *(bp|aa|rc)?"
            r" *(.*DNA|.*RNA|.*dna|.*rna)?"
            r" *(linear|circular)?"
            r" *(?!.*DNA|.*RNA)([A-Z]{3})?"
            r" *([0-9]{2}-[A-Z]{3}-[0-9]{4})?"
        )
        matches = re.match(pattern, line)
        res = dict(
            zip(
                [
                    "locus_name",
                    "size",
                    "unit",
                    "mol_type",
                    "topology",
                    "division",
                    "date",
                ],
                matches.groups(),
            )
        )

        residue_type = ""
        if len(res["locus_name"]) > 16:
            warnings.warn(
                "GenBank LOCUS line identifier over 16 characters",
                BiopythonParserWarning,
            )
        consumer.locus(res["locus_name"])
        if res["size"]:
            assert int(res["size"]) <= sys.maxsize, (
                "Tried to load a sequence with a length %s, your installation of python can only load sequences of length %s"
                % (res["size"], sys.maxsize)
            )
            consumer.size(res["size"])
        if res["mol_type"]:
            consumer.molecule_type(res["mol_type"])
            residue_type += res["mol_type"]
        if res["topology"]:
            consumer.topology(res["topology"])
            residue_type += " " + res["topology"]
        consumer.residue_type(residue_type.strip())
        if res["unit"] == "aa":
            consumer.residue_type("PROTEIN")
        if res["division"]:
            consumer.data_file_division(res["division"])
        if res["date"]:
            consumer.date(res["date"])

    def _feed_header_lines(self, consumer, lines):
        lines = [_f for _f in lines if _f]
        data_acc = []
        working_line_type = None
        for line in lines:
            line_type = line[: self.GENBANK_INDENT].strip()
            data = line[self.GENBANK_INDENT :].strip()
            if working_line_type:
                if line_type == "":
                    data_acc.append(data)
                else:
                    self._process_header_data(consumer, working_line_type, data_acc)
                    working_line_type = line_type
                    data_acc = [data]
            else:
                data_acc.append(data)
                working_line_type = line_type
        self._process_header_data(consumer, working_line_type, data_acc)

    def _process_header_data(self, consumer, line_type, data):
        if line_type == "VERSION":
            _feed_header_version(consumer, data[0])
        elif line_type == "DBLINK":
            _feed_header_dblink(consumer, data)
        elif line_type == "REFERENCE":
            _feed_header_reference(consumer, data)
        elif line_type == "ORGANISM":
            _feed_header_organism(consumer, data)
        elif line_type == "COMMENT":
            _feed_header_comment(
                consumer,
                data,
                self.STRUCTURED_COMMENT_START,
                self.STRUCTURED_COMMENT_DELIM,
                self.STRUCTURED_COMMENT_END,
            )
        elif line_type in self.CONSUMER_DICT:
            _feed_header_common(consumer, data, line_type, self.CONSUMER_DICT)
        else:
            logging.debug("Ignoring GenBank header line:\n" % data)

    def _feed_misc_lines(self, consumer, lines):
        consumer_dict = {
            "BASE COUNT": "base_count",
            "ORIGIN": "origin_name",
            "TLS": "tls",
            "TSA": "tsa",
            "WGS": "wgs",
            "WGS_SCAFLD": "add_wgs_scafld",
        }
        contig = False
        contig_location = ""
        for line in lines:
            if contig:
                if line[: self.GENBANK_INDENT] == self.GENBANK_SPACER:
                    contig_location += line[self.GENBANK_INDENT :].rstrip()
                elif line.startswith("ORIGIN"):
                    line = line[6:].strip()
                    if line:
                        consumer.origin_name(line)
                else:
                    raise ValueError("Expected CONTIG continuation line, got:\n" + line)
            else:
                if line.startswith("CONTIG"):
                    line = line[6:].strip()
                    contig_location = line
                    contig = True
                else:
                    for key in consumer_dict:
                        if line.startswith(key):
                            line = line[len(key) :].strip()
                            if line:
                                getattr(consumer, consumer_dict[key])(line)
                            break
        if contig_location:
            consumer.contig_location(contig_location)
        return
