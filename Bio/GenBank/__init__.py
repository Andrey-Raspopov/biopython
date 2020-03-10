# Copyright 2000 by Jeffrey Chang, Brad Chapman.  All rights reserved.
# Copyright 2006-2017 by Peter Cock.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Code to work with GenBank formatted files.

Rather than using Bio.GenBank, you are now encouraged to use Bio.SeqIO with
the "genbank" or "embl" format names to parse GenBank or EMBL files into
SeqRecord and SeqFeature objects (see the Biopython tutorial for details).

Using Bio.GenBank directly to parse GenBank files is only useful if you want
to obtain GenBank-specific Record objects, which is a much closer
representation to the raw file contents than the SeqRecord alternative from
the FeatureParser (used in Bio.SeqIO).

To use the Bio.GenBank parser, there are two helper functions:

    - read                  Parse a handle containing a single GenBank record
      as Bio.GenBank specific Record objects.
    - parse                 Iterate over a handle containing multiple GenBank
      records as Bio.GenBank specific Record objects.

The following internal classes are not intended for direct use and may
be deprecated in a future release.

Classes:
 - Iterator              Iterate through a file of GenBank entries
 - ErrorFeatureParser    Catch errors caused during parsing.
 - FeatureParser         Parse GenBank data in SeqRecord and SeqFeature objects.
 - RecordParser          Parse GenBank data into a Record object.

Exceptions:
 - ParserFailureError    Exception indicating a failure in the parser (ie.
   scanner or consumer)
 - LocationParserError   Exception indicating a problem with the spark based
   location parser.

"""

import re

# other Biopython stuff
import warnings
from Bio import BiopythonParserWarning

# other Bio.GenBank stuff
from Bio.SeqIO.utils import FeatureValueCleaner
from Bio.SeqIO.GenBankIO import GenBankScanner, _FeatureConsumer, _BaseGenBankConsumer, _re_simple_location, _solo_bond, \
    _re_simple_compound, _re_complex_location, _re_complex_compound, ParserFailureError, _within_location, \
    _oneof_location, _re_within_position, _re_oneof_position

# Constants used to parse GenBank header lines
GENBANK_INDENT = 12
GENBANK_SPACER = " " * GENBANK_INDENT

# Constants for parsing GenBank feature lines
FEATURE_KEY_INDENT = 5
FEATURE_QUALIFIER_INDENT = 21
FEATURE_KEY_SPACER = " " * FEATURE_KEY_INDENT
FEATURE_QUALIFIER_SPACER = " " * FEATURE_QUALIFIER_INDENT

# Regular expressions for location parsing

assert _re_within_position.match("(3.9)")
assert re.compile(_within_location).match("(3.9)..10")
assert re.compile(_within_location).match("26..(30.33)")
assert re.compile(_within_location).match("(13.19)..(20.28)")

assert _re_oneof_position.match("one-of(6,9)")
assert re.compile(_oneof_location).match("one-of(6,9)..101")
assert re.compile(_oneof_location).match("one-of(6,9)..one-of(101,104)")
assert re.compile(_oneof_location).match("6..one-of(101,104)")

assert not _re_oneof_position.match("one-of(3)")
assert _re_oneof_position.match("one-of(3,6)")
assert _re_oneof_position.match("one-of(3,6,9)")

assert _re_simple_location.match("104..160")
assert not _re_simple_location.match("68451760..68452073^68452074")
assert not _re_simple_location.match("<104..>160")
assert not _re_simple_location.match("104")
assert not _re_simple_location.match("<1")
assert not _re_simple_location.match(">99999")
assert not _re_simple_location.match("join(104..160,320..390,504..579)")
assert not _re_simple_compound.match("bond(12,63)")
assert _re_simple_compound.match("join(104..160,320..390,504..579)")
assert _re_simple_compound.match("order(1..69,1308..1465)")
assert not _re_simple_compound.match("order(1..69,1308..1465,1524)")
assert not _re_simple_compound.match("join(<1..442,992..1228,1524..>1983)")
assert not _re_simple_compound.match("join(<1..181,254..336,422..497,574..>590)")
assert not _re_simple_compound.match(
    "join(1475..1577,2841..2986,3074..3193,3314..3481,4126..>4215)"
)
assert not _re_simple_compound.match("test(1..69,1308..1465)")
assert not _re_simple_compound.match("complement(1..69)")
assert not _re_simple_compound.match("(1..69)")
assert _re_complex_location.match("(3.9)..10")
assert _re_complex_location.match("26..(30.33)")
assert _re_complex_location.match("(13.19)..(20.28)")
assert _re_complex_location.match("41^42")  # between
assert _re_complex_location.match("AL121804:41^42")
assert _re_complex_location.match("AL121804:41..610")
assert _re_complex_location.match("AL121804.2:41..610")
assert _re_complex_location.match(
    "AL358792.24.1.166931:3274..3461"
)  # lots of dots in external reference
assert _re_complex_location.match("one-of(3,6)..101")
assert _re_complex_compound.match(
    "join(153490..154269,AL121804.2:41..610,AL121804.2:672..1487)"
)
assert not _re_simple_compound.match(
    "join(153490..154269,AL121804.2:41..610,AL121804.2:672..1487)"
)
assert _re_complex_compound.match("join(complement(69611..69724),139856..140650)")
assert _re_complex_compound.match(
    "join(complement(AL354868.10.1.164018:80837..81016),complement(AL354868.10.1.164018:80539..80835))"
)

# Trans-spliced example from NC_016406, note underscore in reference name:
assert _re_complex_location.match("NC_016402.1:6618..6676")
assert _re_complex_location.match("181647..181905")
assert _re_complex_compound.match(
    "join(complement(149815..150200),complement(293787..295573),NC_016402.1:6618..6676,181647..181905)"
)
assert not _re_complex_location.match(
    "join(complement(149815..150200),complement(293787..295573),NC_016402.1:6618..6676,181647..181905)"
)
assert not _re_simple_compound.match(
    "join(complement(149815..150200),complement(293787..295573),NC_016402.1:6618..6676,181647..181905)"
)
assert not _re_complex_location.match(
    "join(complement(149815..150200),complement(293787..295573),NC_016402.1:6618..6676,181647..181905)"
)
assert not _re_simple_location.match(
    "join(complement(149815..150200),complement(293787..295573),NC_016402.1:6618..6676,181647..181905)"
)

assert _solo_bond.match("bond(196)")
assert _solo_bond.search("bond(196)")
assert _solo_bond.search("join(bond(284),bond(305),bond(309),bond(305))")


class Iterator:
    """Iterator interface to move over a file of GenBank entries one at a time (OBSOLETE).

    This class is likely to be deprecated in a future release of Biopython.
    Please use Bio.SeqIO.parse(..., format="gb") or Bio.GenBank.parse(...)
    for SeqRecord and GenBank specific Record objects respectively instead.
    """

    def __init__(self, handle, parser=None):
        """Initialize the iterator.

        Arguments:
         - handle - A handle with GenBank entries to iterate through.
         - parser - An optional parser to pass the entries through before
           returning them. If None, then the raw entry will be returned.

        """
        self.handle = handle
        self._parser = parser

    def __next__(self):
        """Return the next GenBank record from the handle.

        Will return None if we ran out of records.
        """
        if self._parser is None:
            lines = []
            while True:
                line = self.handle.readline()
                if not line:
                    return None  # Premature end of file?
                lines.append(line)
                if line.rstrip() == "//":
                    break
            return "".join(lines)
        try:
            return self._parser.parse(self.handle)
        except StopIteration:
            return None

    def __iter__(self):
        """Iterate over the records."""
        return iter(self.__next__, None)


_cleaner = FeatureValueCleaner()


class FeatureParser:
    """Parse GenBank files into Seq + Feature objects (OBSOLETE).

    Direct use of this class is discouraged, and may be deprecated in
    a future release of Biopython.

    Please use Bio.SeqIO.parse(...) or Bio.SeqIO.read(...) instead.
    """

    def __init__(self, use_fuzziness=1, feature_cleaner=_cleaner):
        """Initialize a GenBank parser and Feature consumer.

        Arguments:
         - use_fuzziness - Specify whether or not to use fuzzy representations.
           The default is 1 (use fuzziness).
         - feature_cleaner - A class which will be used to clean out the
           values of features. This class must implement the function
           clean_value. GenBank.utils has a "standard" cleaner class, which
           is used by default.

        """
        self._scanner = GenBankScanner()
        self.use_fuzziness = use_fuzziness
        self._cleaner = feature_cleaner

    def parse(self, handle):
        """Parse the specified handle."""
        _consumer = _FeatureConsumer(self.use_fuzziness, self._cleaner)
        self._scanner.feed(handle, _consumer)
        return _consumer.data


class RecordParser:
    """Parse GenBank files into Record objects (OBSOLETE).

    Direct use of this class is discouraged, and may be deprecated in
    a future release of Biopython.

    Please use the Bio.GenBank.parse(...) or Bio.GenBank.read(...) functions
    instead.
    """

    def __init__(self):
        """Initialize the parser."""
        self._scanner = GenBankScanner()

    def parse(self, handle):
        """Parse the specified handle into a GenBank record."""
        _consumer = _RecordConsumer()

        self._scanner.feed(handle, _consumer)
        return _consumer.data


class _RecordConsumer(_BaseGenBankConsumer):
    """Create a GenBank Record object from scanner generated information (PRIVATE)."""

    def __init__(self):
        _BaseGenBankConsumer.__init__(self)
        from . import Record

        self.data = Record.Record()

        self._seq_data = []
        self._cur_reference = None
        self._cur_feature = None
        self._cur_qualifier = None

    def tls(self, content):
        self.data.tls = content.split("-")

    def tsa(self, content):
        self.data.tsa = content.split("-")

    def wgs(self, content):
        self.data.wgs = content.split("-")

    def add_wgs_scafld(self, content):
        self.data.wgs_scafld.append(content.split("-"))

    def locus(self, content):
        self.data.locus = content

    def size(self, content):
        self.data.size = content

    def residue_type(self, content):
        # Be lenient about parsing, but technically lowercase residue types are malformed.
        if "dna" in content or "rna" in content:
            warnings.warn(
                "Invalid seq_type (%s): DNA/RNA should be uppercase." % content,
                BiopythonParserWarning,
            )
        self.data.residue_type = content

    def data_file_division(self, content):
        self.data.data_file_division = content

    def date(self, content):
        self.data.date = content

    def definition(self, content):
        self.data.definition = content

    def accession(self, content):
        for acc in self._split_accessions(content):
            if acc not in self.data.accession:
                self.data.accession.append(acc)

    def molecule_type(self, mol_type):
        """Validate and record the molecule type (for round-trip etc)."""
        if mol_type:
            if "circular" in mol_type or "linear" in mol_type:
                raise ParserFailureError(
                    "Molecule type %r should not include topology" % mol_type
                )

            # Writing out records will fail if we have a lower case DNA
            # or RNA string in here, so upper case it.
            # This is a bit ugly, but we don't want to upper case e.g.
            # the m in mRNA, but thanks to the strip we lost the spaces
            # so we need to index from the back
            if mol_type[-3:].upper() in ("DNA", "RNA") and not mol_type[-3:].isupper():
                warnings.warn(
                    "Non-upper case molecule type in LOCUS line: %s" % mol_type,
                    BiopythonParserWarning,
                )

            self.data.molecule_type = mol_type

    def topology(
        self, topology
    ):  # noqa: D402  # flake8 thinks this line is a function signature. It ain't
        """Validate and record sequence topology (linear or circular as strings)."""
        if topology:
            if topology not in ["linear", "circular"]:
                raise ParserFailureError(
                    "Unexpected topology %r should be linear or circular" % topology
                )
            self.data.topology = topology

    def nid(self, content):
        self.data.nid = content

    def pid(self, content):
        self.data.pid = content

    def version(self, content):
        self.data.version = content

    def db_source(self, content):
        self.data.db_source = content.rstrip()

    def gi(self, content):
        self.data.gi = content

    def keywords(self, content):
        self.data.keywords = self._split_keywords(content)

    def project(self, content):
        self.data.projects.extend(p for p in content.split() if p)

    def dblink(self, content):
        self.data.dblinks.append(content)

    def segment(self, content):
        self.data.segment = content

    def source(self, content):
        self.data.source = content

    def organism(self, content):
        self.data.organism = content

    def taxonomy(self, content):
        self.data.taxonomy = self._split_taxonomy(content)

    def reference_num(self, content):
        """Grab the reference number and signal the start of a new reference."""
        # check if we have a reference to add
        if self._cur_reference is not None:
            self.data.references.append(self._cur_reference)

        from . import Record

        self._cur_reference = Record.Reference()
        self._cur_reference.number = content

    def reference_bases(self, content):
        self._cur_reference.bases = content

    def authors(self, content):
        self._cur_reference.authors = content

    def consrtm(self, content):
        self._cur_reference.consrtm = content

    def title(self, content):
        if self._cur_reference is None:
            warnings.warn(
                "GenBank TITLE line without REFERENCE line.", BiopythonParserWarning
            )
            return
        self._cur_reference.title = content

    def journal(self, content):
        self._cur_reference.journal = content

    def medline_id(self, content):
        self._cur_reference.medline_id = content

    def pubmed_id(self, content):
        self._cur_reference.pubmed_id = content

    def remark(self, content):
        self._cur_reference.remark = content

    def comment(self, content):
        self.data.comment += "\n".join(content)

    def structured_comment(self, content):
        self.data.structured_comment = content

    def primary_ref_line(self, content):
        """Save reference data for the PRIMARY line."""
        self.data.primary.append(content)

    def primary(self, content):
        pass

    def features_line(self, content):
        """Get ready for the feature table when we reach the FEATURE line."""
        self.start_feature_table()

    def start_feature_table(self):
        """Signal the start of the feature table."""
        # we need to add on the last reference
        if self._cur_reference is not None:
            self.data.references.append(self._cur_reference)

    def feature_key(self, content):
        """Grab the key of the feature and signal the start of a new feature."""
        # first add on feature information if we've got any
        self._add_feature()

        from . import Record

        self._cur_feature = Record.Feature()
        self._cur_feature.key = content

    def _add_feature(self):
        """Add a feature to the record, with relevant checks (PRIVATE).

        This does all of the appropriate checking to make sure we haven't
        left any info behind, and that we are only adding info if it
        exists.
        """
        if self._cur_feature is not None:
            # if we have a left over qualifier, add it to the qualifiers
            # on the current feature
            if self._cur_qualifier is not None:
                self._cur_feature.qualifiers.append(self._cur_qualifier)

            self._cur_qualifier = None
            self.data.features.append(self._cur_feature)

    def location(self, content):
        self._cur_feature.location = self._clean_location(content)

    def feature_qualifier(self, key, value):
        self.feature_qualifier_name([key])
        if value is not None:
            self.feature_qualifier_description(value)

    def feature_qualifier_name(self, content_list):
        """Deal with qualifier names.

        We receive a list of keys, since you can have valueless keys such as
        /pseudo which would be passed in with the next key (since no other
        tags separate them in the file)
        """
        from . import Record

        for content in content_list:
            # the record parser keeps the /s -- add them if we don't have 'em
            if not content.startswith("/"):
                content = "/%s" % content
            # add on a qualifier if we've got one
            if self._cur_qualifier is not None:
                self._cur_feature.qualifiers.append(self._cur_qualifier)

            self._cur_qualifier = Record.Qualifier()
            self._cur_qualifier.key = content

    def feature_qualifier_description(self, content):
        # if we have info then the qualifier key should have a ='s
        if "=" not in self._cur_qualifier.key:
            self._cur_qualifier.key = "%s=" % self._cur_qualifier.key
        cur_content = self._remove_newlines(content)
        # remove all spaces from the value if it is a type where spaces
        # are not important
        for remove_space_key in self.__class__.remove_space_keys:
            if remove_space_key in self._cur_qualifier.key:
                cur_content = self._remove_spaces(cur_content)
        self._cur_qualifier.value = self._normalize_spaces(cur_content)

    def base_count(self, content):
        self.data.base_counts = content

    def origin_name(self, content):
        self.data.origin = content

    def contig_location(self, content):
        """Signal that we have contig information to add to the record."""
        self.data.contig = self._clean_location(content)

    def sequence(self, content):
        """Add sequence information to a list of sequence strings.

        This removes spaces in the data and uppercases the sequence, and
        then adds it to a list of sequences. Later on we'll join this
        list together to make the final sequence. This is faster than
        adding on the new string every time.
        """
        assert " " not in content
        self._seq_data.append(content.upper())

    def record_end(self, content):
        """Signal the end of the record and do any necessary clean-up."""
        # add together all of the sequence parts to create the
        # final sequence string
        self.data.sequence = "".join(self._seq_data)
        # add on the last feature
        self._add_feature()


def parse(handle):
    """Iterate over GenBank formatted entries as Record objects.

    >>> from Bio import GenBank
    >>> with open("GenBank/NC_000932.gb") as handle:
    ...     for record in GenBank.parse(handle):
    ...         print(record.accession)
    ['NC_000932']

    To get SeqRecord objects use Bio.SeqIO.parse(..., format="gb")
    instead.
    """
    return iter(Iterator(handle, RecordParser()))


def read(handle):
    """Read a handle containing a single GenBank entry as a Record object.

    >>> from Bio import GenBank
    >>> with open("GenBank/NC_000932.gb") as handle:
    ...     record = GenBank.read(handle)
    ...     print(record.accession)
    ['NC_000932']

    To get a SeqRecord object use Bio.SeqIO.read(..., format="gb")
    instead.
    """
    iterator = parse(handle)
    try:
        record = next(iterator)
    except StopIteration:
        raise ValueError("No records found in handle") from None
    try:
        next(iterator)
        raise ValueError("More than one record found in handle")
    except StopIteration:
        pass
    return record


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
