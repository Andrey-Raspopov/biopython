# Copyright 2007-2016 by Peter Cock.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.SeqIO support for the "genbank" and "embl" file formats.

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
import warnings
from Bio import BiopythonWarning, BiopythonParserWarning
from Bio.Alphabet import generic_protein
from Bio.File import as_handle
from Bio.Seq import Seq

from Bio.SeqIO.Interfaces import SequentialSequenceWriter
from Bio import SeqFeature

from Bio.SeqRecord import SeqRecord


def _insdc_feature_position_string(pos, offset=0):
    """Build a GenBank/EMBL position string (PRIVATE).

    Use offset=1 to add one to convert a start position from python counting.
    """
    if isinstance(pos, SeqFeature.ExactPosition):
        return "%i" % (pos.position + offset)
    elif isinstance(pos, SeqFeature.WithinPosition):
        return "(%i.%i)" % (
            pos.position + offset,
            pos.position + pos.extension + offset,
        )
    elif isinstance(pos, SeqFeature.BetweenPosition):
        return "(%i^%i)" % (
            pos.position + offset,
            pos.position + pos.extension + offset,
        )
    elif isinstance(pos, SeqFeature.BeforePosition):
        return "<%i" % (pos.position + offset)
    elif isinstance(pos, SeqFeature.AfterPosition):
        return ">%i" % (pos.position + offset)
    elif isinstance(pos, SeqFeature.OneOfPosition):
        return "one-of(%s)" % ",".join(
            _insdc_feature_position_string(p, offset) for p in pos.position_choices
        )
    elif isinstance(pos, SeqFeature.AbstractPosition):
        raise NotImplementedError("Please report this as a bug in Biopython.")
    else:
        raise ValueError("Expected a SeqFeature position object.")


def _insdc_location_string_ignoring_strand_and_subfeatures(location, rec_length):
    if location.ref:
        ref = "%s:" % location.ref
    else:
        ref = ""
    assert not location.ref_db
    if (
        isinstance(location.start, SeqFeature.ExactPosition)
        and isinstance(location.end, SeqFeature.ExactPosition)
        and location.start.position == location.end.position
    ):
        if location.end.position == rec_length:
            return "%s%i^1" % (ref, rec_length)
        else:
            return "%s%i^%i" % (ref, location.end.position, location.end.position + 1)
    if (
        isinstance(location.start, SeqFeature.ExactPosition)
        and isinstance(location.end, SeqFeature.ExactPosition)
        and location.start.position + 1 == location.end.position
    ):
        return "%s%i" % (ref, location.end.position)
    elif isinstance(location.start, SeqFeature.UnknownPosition) or isinstance(
        location.end, SeqFeature.UnknownPosition
    ):
        if isinstance(location.start, SeqFeature.UnknownPosition) and isinstance(
            location.end, SeqFeature.UnknownPosition
        ):
            raise ValueError("Feature with unknown location")
        elif isinstance(location.start, SeqFeature.UnknownPosition):
            return "%s<%i..%s" % (
                ref,
                location.nofuzzy_end,
                _insdc_feature_position_string(location.end),
            )
        else:
            return "%s%s..>%i" % (
                ref,
                _insdc_feature_position_string(location.start, +1),
                location.nofuzzy_start + 1,
            )
    else:
        return (
            ref
            + _insdc_feature_position_string(location.start, +1)
            + ".."
            + _insdc_feature_position_string(location.end)
        )


def _insdc_location_string(location, rec_length):
    """Build a GenBank/EMBL location from a (Compound) FeatureLocation (PRIVATE).

    There is a choice of how to show joins on the reverse complement strand,
    GenBank used "complement(join(1,10),(20,100))" while EMBL used to use
    "join(complement(20,100),complement(1,10))" instead (but appears to have
    now adopted the GenBank convention). Notice that the order of the entries
    is reversed! This function therefore uses the first form. In this situation
    we expect the CompoundFeatureLocation and its parts to all be marked as
    strand == -1, and to be in the order 19:100 then 0:10.
    """
    try:
        parts = location.parts
        if location.strand == -1:
            return "complement(%s(%s))" % (
                location.operator,
                ",".join(
                    _insdc_location_string_ignoring_strand_and_subfeatures(
                        p, rec_length
                    )
                    for p in parts[::-1]
                ),
            )
        else:
            return "%s(%s)" % (
                location.operator,
                ",".join(_insdc_location_string(p, rec_length) for p in parts),
            )
    except AttributeError:
        loc = _insdc_location_string_ignoring_strand_and_subfeatures(
            location, rec_length
        )
        if location.strand == -1:
            return "complement(%s)" % loc
        else:
            return loc


def parse_feature(feature_key, lines):
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
        feature_location = _extract_feature_location(iterator)

        qualifiers = []
        key = None

        for line_number, line in enumerate(iterator):
            if line[0] == "/":
                line = line[1:]
                try:
                    key, value = line.split("=", 1)
                except ValueError:
                    key = line
                    value = None
                if value:
                    if value.startswith(" ") and value.lstrip().startswith('"'):
                        warnings.warn(
                            "White space after equals in qualifier",
                            BiopythonParserWarning,
                        )
                    elif value == '"':
                        logging.debug("Single quote %s:%s" % (key, value))
                    elif value[0] == '"':
                        while value[-1] != '"':
                            value += f"\n{next(iterator)}"
                    value = value.strip()
                qualifiers.append((key, value))
            else:
                assert len(qualifiers) > 0
                assert key == qualifiers[-1][0]
                assert qualifiers[-1][1], StopIteration
                qualifiers[-1] = (key, f"{qualifiers[-1][1]}\n{line}")
        return feature_key, feature_location, qualifiers
    except StopIteration:
        raise ValueError(
            "Problem with '%s' feature:\n%s" % (feature_key, "\n".join(lines))
        ) from None


def _extract_feature_location(iterator):
    feature_location = next(iterator).strip()
    while True:
        if feature_location[-1:] != ",":
            if feature_location.count("(") > feature_location.count(")"):
                warnings.warn(
                    "Non-standard feature line wrapping (didn't break on comma)?",
                    BiopythonParserWarning,
                )
                feature_location += next(iterator).strip()
            else:
                break
        else:
            feature_location += next(iterator).strip()
    return feature_location


class _InsdcWriter(SequentialSequenceWriter):
    """Base class for GenBank and EMBL writers (PRIVATE)."""

    MAX_WIDTH = 80
    QUALIFIER_INDENT = 21
    QUALIFIER_INDENT_STR = " " * QUALIFIER_INDENT
    QUALIFIER_INDENT_TMP = "     %s                "
    FTQUAL_NO_QUOTE = (
        "anticodon",
        "citation",
        "codon_start",
        "compare",
        "direction",
        "estimated_length",
        "mod_base",
        "number",
        "rpt_type",
        "rpt_unit_range",
        "tag_peptide",
        "transl_except",
        "transl_table",
    )

    def _write_feature_qualifier(self, key, value=None, quote=None):
        if value is None:
            self.handle.write("%s/%s\n" % (self.QUALIFIER_INDENT_STR, key))
            return

        if type(value) == str:
            value = value.replace(
                '"', '""'
            )

        if quote is None:
            if isinstance(value, int) or key in self.FTQUAL_NO_QUOTE:
                quote = False
            else:
                quote = True
        if quote:
            line = '%s/%s="%s"' % (self.QUALIFIER_INDENT_STR, key, value)
        else:
            line = "%s/%s=%s" % (self.QUALIFIER_INDENT_STR, key, value)
        if len(line) <= self.MAX_WIDTH:
            self.handle.write(line + "\n")
            return
        while line.lstrip():
            if len(line) <= self.MAX_WIDTH:
                self.handle.write(line + "\n")
                return
            for index in range(
                min(len(line) - 1, self.MAX_WIDTH), self.QUALIFIER_INDENT + 1, -1
            ):
                if line[index] == " ":
                    break
            if line[index] != " ":
                index = self.MAX_WIDTH
            assert index <= self.MAX_WIDTH
            self.handle.write(line[:index] + "\n")
            line = self.QUALIFIER_INDENT_STR + line[index:].lstrip()

    def _wrap_location(self, location):
        # TODO - Rewrite this not to recurse!
        length = self.MAX_WIDTH - self.QUALIFIER_INDENT
        if len(location) <= length:
            return location
        index = location[:length].rfind(",")
        if index == -1:
            warnings.warn("Couldn't split location:\n%s" % location, BiopythonWarning)
            return location
        return (
            location[: index + 1]
            + "\n"
            + self.QUALIFIER_INDENT_STR
            + self._wrap_location(location[index + 1 :])
        )

    def _write_feature(self, feature, record_length):
        assert feature.type, feature
        location = _insdc_location_string(feature.location, record_length)
        f_type = feature.type.replace(" ", "_")
        line = (
            (self.QUALIFIER_INDENT_TMP % f_type)[: self.QUALIFIER_INDENT]
            + self._wrap_location(location)
            + "\n"
        )
        self.handle.write(line)
        for key, values in feature.qualifiers.items():
            if isinstance(values, (list, tuple)):
                for value in values:
                    self._write_feature_qualifier(key, value)
            else:
                self._write_feature_qualifier(key, values)

    @staticmethod
    def _get_annotation_str(record, key, default=".", just_first=False):
        """Get an annotation dictionary entry (as a string) (PRIVATE).

        Some entries are lists, in which case if just_first=True the first entry
        is returned.  If just_first=False (default) this verifies there is only
        one entry before returning it.
        """
        try:
            answer = record.annotations[key]
        except KeyError:
            return default
        if isinstance(answer, list):
            if not just_first:
                assert len(answer) == 1
            return str(answer[0])
        else:
            return str(answer)

    @staticmethod
    def _split_multi_line(text, max_len):
        """Return a list of strings (PRIVATE).

        Any single words which are too long get returned as a whole line
        (e.g. URLs) without an exception or warning.
        """
        # TODO - Do the line spliting while preserving white space?
        text = text.strip()
        if len(text) <= max_len:
            return [text]

        words = text.split()
        text = ""
        while words and len(text) + 1 + len(words[0]) <= max_len:
            text += " " + words.pop(0)
            text = text.strip()
        answer = [text]
        while words:
            text = words.pop(0)
            while words and len(text) + 1 + len(words[0]) <= max_len:
                text += " " + words.pop(0)
                text = text.strip()
            answer.append(text)
        assert not words
        return answer

    def _split_contig(self, record, max_len):
        """Return a list of strings, splits on commas (PRIVATE)."""
        # TODO - Merge this with _write_multi_line method?
        # It would need the addition of the comma splitting logic...
        # are there any other cases where that would be sensible?
        contig = record.annotations.get("contig", "")
        if isinstance(contig, (list, tuple)):
            contig = "".join(contig)
        contig = self.clean(contig)
        answer = []
        while contig:
            if len(contig) > max_len:
                pos = contig[: max_len - 1].rfind(",")
                if pos == -1:
                    raise ValueError("Could not break up CONTIG")
                text, contig = contig[: pos + 1], contig[pos + 1 :]
            else:
                text, contig = contig, ""
            answer.append(text)
        return answer


class _InsdcScanner:
    """Basic functions for breaking up a GenBank/EMBL file into sub sections.

    The International Nucleotide Sequence Database Collaboration (INSDC)
    between the DDBJ, EMBL, and GenBank.  These organisations all use the
    same "Feature Table" layout in their plain text flat file formats.

    However, the header and sequence sections of an EMBL file are very
    different in layout to those produced by GenBank/DDBJ.
    """

    RECORD_START = "XXX"
    HEADER_WIDTH = 3
    FEATURE_START_MARKERS = ["XXX***FEATURES***XXX"]
    FEATURE_END_MARKERS = ["XXX***END FEATURES***XXX"]
    FEATURE_QUALIFIER_INDENT = 0
    FEATURE_QUALIFIER_SPACER = ""
    SEQUENCE_HEADERS = ["XXX"]

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
                if (
                    line[self.FEATURE_QUALIFIER_INDENT] != " "
                    and " " in line[self.FEATURE_QUALIFIER_INDENT :]
                ):
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
                ):
                    feature_lines.append(line[self.FEATURE_QUALIFIER_INDENT :].strip())
                    line = self.handle.readline()
                features.append(parse_feature(feature_key, feature_lines))
        self.line = line
        return features

    def parse_footer(self):
        """Return a tuple containing a list of any misc strings, and the sequence."""
        pass

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
        pass

    def feed(self, handle, consumer, do_features=True):
        """Feed a set of data into the consumer.

        Arguments:
         - handle - A handle with the information to parse.
         - consumer - The consumer that should be informed of events.
         - do_features - Boolean, should the features be parsed?
           Skipping the features can be much faster.

        Return values:
         - true  - Passed a record
         - false - Did not find a record

        """
        self.set_handle(handle)
        if not self.find_start():
            consumer.data = None
            return False
        self._feed_first_line(consumer, self.line)
        self._feed_header_lines(consumer, self.parse_header())
        if do_features:
            self._feed_feature_table(consumer, self.parse_features(skip=False))
        else:
            self.parse_features(skip=True)
        misc_lines, sequence_string = self.parse_footer()
        self._feed_misc_lines(consumer, misc_lines)
        consumer.sequence(sequence_string)
        consumer.record_end("//")
        assert self.line == "//"
        return True

    def parse(self, handle, do_features=True):
        """Return a SeqRecord (with SeqFeatures if do_features=True).

        See also the method parse_records() for use on multi-record files.
        """
        from Bio.SeqIO.GenBankIO import _FeatureConsumer
        from Bio.SeqIO.utils import FeatureValueCleaner

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
                self.parse_header()
                feature_tuples = self.parse_features()
                while True:
                    line = self.handle.readline()
                    if not line:
                        break
                    if line[:2] == "//":
                        break
                self.line = line.rstrip()

                for key, location_string, qualifiers in feature_tuples:
                    if key == "CDS":
                        record = SeqRecord(seq=None)
                        annotations = record.annotations
                        annotations["raw_location"] = location_string.replace(" ", "")

                        for (qualifier_name, qualifier_data) in qualifiers:
                            if (
                                qualifier_data is not None
                                and qualifier_data[0] == '"'
                                and qualifier_data[-1] == '"'
                            ):
                                qualifier_data = qualifier_data[1:-1]
                            if qualifier_name == "translation":
                                assert record.seq is None, "Multiple translations!"
                                record.seq = Seq(
                                    qualifier_data.replace("\n", ""), alphabet
                                )
                            elif qualifier_name == "db_xref":
                                record.dbxrefs.append(qualifier_data)
                            else:
                                if qualifier_data is not None:
                                    qualifier_data = qualifier_data.replace(
                                        "\n", " "
                                    ).replace("  ", " ")
                                try:
                                    annotations[qualifier_name] += " " + qualifier_data
                                except KeyError:
                                    annotations[qualifier_name] = qualifier_data

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


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest(verbose=0)
