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

from Bio import Alphabet, BiopythonWarning, SeqFeature, BiopythonParserWarning
from Bio.Seq import UnknownSeq
from Bio.SeqIO.InsdcIO import _InsdcWriter, _InsdcScanner


class EmblWriter(_InsdcWriter):
    """EMBL writer."""

    HEADER_WIDTH = 5
    QUALIFIER_INDENT = 21
    QUALIFIER_INDENT_STR = "FT" + " " * (QUALIFIER_INDENT - 2)
    QUALIFIER_INDENT_TMP = "FT   %s                "  # 21 if %s is empty
    # Note second spacer line of just FH is expected:
    FEATURE_HEADER = "FH   Key             Location/Qualifiers\nFH\n"

    LETTERS_PER_BLOCK = 10
    BLOCKS_PER_LINE = 6
    LETTERS_PER_LINE = LETTERS_PER_BLOCK * BLOCKS_PER_LINE
    POSITION_PADDING = 10

    def _write_contig(self, record):
        max_len = self.MAX_WIDTH - self.HEADER_WIDTH
        lines = self._split_contig(record, max_len)
        for text in lines:
            self._write_single_line("CO", text)

    def _write_sequence(self, record):
        handle = self.handle  # save looking up this multiple times

        if isinstance(record.seq, UnknownSeq):
            # We have already recorded the length, and there is no need
            # to record a long sequence of NNNNNNN...NNN or whatever.
            if "contig" in record.annotations:
                self._write_contig(record)
            else:
                # TODO - Can the sequence just be left out as in GenBank files?
                handle.write("SQ   \n")
            return

        # Catches sequence being None
        data = self._get_seq_string(record).lower()
        seq_len = len(data)

        # Get the base alphabet (underneath any Gapped or StopCodon encoding)
        a = Alphabet._get_base_alphabet(record.seq.alphabet)
        if isinstance(a, Alphabet.DNAAlphabet):
            # TODO - What if we have RNA?
            a_count = data.count("A") + data.count("a")
            c_count = data.count("C") + data.count("c")
            g_count = data.count("G") + data.count("g")
            t_count = data.count("T") + data.count("t")
            other = seq_len - (a_count + c_count + g_count + t_count)
            handle.write(
                "SQ   Sequence %i BP; %i A; %i C; %i G; %i T; %i other;\n"
                % (seq_len, a_count, c_count, g_count, t_count, other)
            )
        else:
            handle.write("SQ   \n")

        for line_number in range(0, seq_len // self.LETTERS_PER_LINE):
            handle.write("    ")  # Just four, not five
            for block in range(self.BLOCKS_PER_LINE):
                index = (
                    self.LETTERS_PER_LINE * line_number + self.LETTERS_PER_BLOCK * block
                )
                handle.write(" %s" % data[index : index + self.LETTERS_PER_BLOCK])
            handle.write(
                str((line_number + 1) * self.LETTERS_PER_LINE).rjust(
                    self.POSITION_PADDING
                )
            )
            handle.write("\n")
        if seq_len % self.LETTERS_PER_LINE:
            # Final (partial) line
            line_number = seq_len // self.LETTERS_PER_LINE
            handle.write("    ")  # Just four, not five
            for block in range(self.BLOCKS_PER_LINE):
                index = (
                    self.LETTERS_PER_LINE * line_number + self.LETTERS_PER_BLOCK * block
                )
                handle.write(
                    (" %s" % data[index : index + self.LETTERS_PER_BLOCK]).ljust(11)
                )
            handle.write(str(seq_len).rjust(self.POSITION_PADDING))
            handle.write("\n")

    def _write_single_line(self, tag, text):
        assert len(tag) == 2
        line = tag + "   " + text
        if len(text) > self.MAX_WIDTH:
            warnings.warn("Line %r too long" % line, BiopythonWarning)
        self.handle.write(line + "\n")

    def _write_multi_line(self, tag, text):
        max_len = self.MAX_WIDTH - self.HEADER_WIDTH
        lines = self._split_multi_line(text, max_len)
        for line in lines:
            self._write_single_line(tag, line)

    def _write_the_first_lines(self, record):
        """Write the ID and AC lines (PRIVATE)."""
        if "." in record.id and record.id.rsplit(".", 1)[1].isdigit():
            version = "SV " + record.id.rsplit(".", 1)[1]
            accession = self._get_annotation_str(
                record, "accession", record.id.rsplit(".", 1)[0], just_first=True
            )
        else:
            version = ""
            accession = self._get_annotation_str(
                record, "accession", record.id, just_first=True
            )

        if ";" in accession:
            raise ValueError(
                "Cannot have semi-colon in EMBL accession, %s" % repr(str(accession))
            )
        if " " in accession:
            # This is out of practicallity... might it be allowed?
            raise ValueError(
                "Cannot have spaces in EMBL accession, %s" % repr(str(accession))
            )

        topology = self._get_annotation_str(record, "topology", default="")

        # Get the molecule type
        # TODO - record this explicitly in the parser?
        # Get the base alphabet (underneath any Gapped or StopCodon encoding)
        a = Alphabet._get_base_alphabet(record.seq.alphabet)
        if not isinstance(a, Alphabet.Alphabet):
            raise TypeError("Invalid alphabet")
        elif isinstance(a, Alphabet.DNAAlphabet):
            mol_type = "DNA"
            units = "BP"
        elif isinstance(a, Alphabet.RNAAlphabet):
            mol_type = "RNA"
            units = "BP"
        elif isinstance(a, Alphabet.ProteinAlphabet):
            mol_type = "PROTEIN"
            units = "AA"
        else:
            # Must be something like NucleotideAlphabet
            raise ValueError("Need a DNA, RNA or Protein alphabet")

        if record.annotations.get("molecule_type", None):
            # Note often get RNA vs DNA discrepancy in real EMBL/NCBI files
            mol_type = record.annotations["molecule_type"]
            if mol_type in ["protein"]:
                mol_type = "PROTEIN"

        # Get the taxonomy division
        division = self._get_data_division(record)

        # TODO - Full ID line
        handle = self.handle
        # ID   <1>; SV <2>; <3>; <4>; <5>; <6>; <7> BP.
        # 1. Primary accession number
        # 2. Sequence version number
        # 3. Topology: 'circular' or 'linear'
        # 4. Molecule type
        # 5. Data class
        # 6. Taxonomic division
        # 7. Sequence length
        self._write_single_line(
            "ID",
            "%s; %s; %s; %s; ; %s; %i %s."
            % (accession, version, topology, mol_type, division, len(record), units),
        )
        handle.write("XX\n")
        self._write_single_line("AC", accession + ";")
        handle.write("XX\n")

    @staticmethod
    def _get_data_division(record):
        try:
            division = record.annotations["data_file_division"]
        except KeyError:
            division = "UNC"
        if division in [
            "PHG",
            "ENV",
            "FUN",
            "HUM",
            "INV",
            "MAM",
            "VRT",
            "MUS",
            "PLN",
            "PRO",
            "ROD",
            "SYN",
            "TGN",
            "UNC",
            "VRL",
            "XXX",
        ]:
            # Good, already EMBL style
            #    Division                 Code
            #    -----------------        ----
            #    Bacteriophage            PHG
            #    Environmental Sample     ENV
            #    Fungal                   FUN
            #    Human                    HUM
            #    Invertebrate             INV
            #    Other Mammal             MAM
            #    Other Vertebrate         VRT
            #    Mus musculus             MUS
            #    Plant                    PLN
            #    Prokaryote               PRO
            #    Other Rodent             ROD
            #    Synthetic                SYN
            #    Transgenic               TGN
            #    Unclassified             UNC (i.e. unknown)
            #    Viral                    VRL
            #
            # (plus XXX used for submiting data to EMBL)
            pass
        else:
            # See if this is in GenBank style & can be converted.
            # Generally a problem as the GenBank groups are wider
            # than those of EMBL. Note that GenBank use "BCT" for
            # both bacteria and acherea thus this maps to EMBL's
            # "PRO" nicely.
            gbk_to_embl = {"BCT": "PRO", "UNK": "UNC"}
            try:
                division = gbk_to_embl[division]
            except KeyError:
                division = "UNC"
        assert len(division) == 3
        return division

    def _write_keywords(self, record):
        # Put the keywords right after DE line.
        # Each 'keyword' can have multiple words and spaces, but we
        # must not split any 'keyword' between lines.
        # TODO - Combine short keywords onto one line
        for keyword in record.annotations["keywords"]:
            self._write_single_line("KW", keyword)
        self.handle.write("XX\n")

    def _write_references(self, record):
        # The order should be RN, RC, RP, RX, RG, RA, RT, RL
        number = 0
        for ref in record.annotations["references"]:
            if not isinstance(ref, SeqFeature.Reference):
                continue
            number += 1
            self._write_single_line("RN", "[%i]" % number)
            # TODO - support for RC line (needed in parser too)
            # TODO - support more complex record reference locations?
            if ref.location and len(ref.location) == 1:
                self._write_single_line(
                    "RP",
                    "%i-%i"
                    % (ref.location[0].nofuzzy_start + 1, ref.location[0].nofuzzy_end),
                )
            # TODO - record any DOI or AGRICOLA identifier in the reference object?
            if ref.pubmed_id:
                self._write_single_line("RX", "PUBMED; %s." % ref.pubmed_id)
            if ref.consrtm:
                self._write_single_line("RG", "%s" % ref.consrtm)
            if ref.authors:
                # We store the AUTHORS data as a single string
                self._write_multi_line("RA", ref.authors + ";")
            if ref.title:
                # We store the title as a single string
                self._write_multi_line("RT", '"%s";' % ref.title)
            if ref.journal:
                # We store this as a single string - holds the journal name,
                # volume, year, and page numbers of the citation
                self._write_multi_line("RL", ref.journal)
            self.handle.write("XX\n")

    def _write_comment(self, record):
        # This is a bit complicated due to the range of possible
        # ways people might have done their annotation...
        # Currently the parser uses a single string with newlines.
        # A list of lines is also reasonable.
        # A single (long) string is perhaps the most natural of all.
        # This means we may need to deal with line wrapping.
        comment = record.annotations["comment"]
        if isinstance(comment, str):
            lines = comment.split("\n")
        elif isinstance(comment, (list, tuple)):
            lines = comment
        else:
            raise ValueError("Could not understand comment annotation")
        # TODO - Merge this with the GenBank comment code?
        if not lines:
            return
        for line in lines:
            self._write_multi_line("CC", line)
        self.handle.write("XX\n")

    def write_record(self, record):
        """Write a single record to the output file."""
        handle = self.handle
        self._write_the_first_lines(record)

        # PR line (0 or 1 lines only), project identifier
        #
        # Assuming can't use 2 lines, we should prefer newer GenBank
        # DBLINK BioProject:... entries over the older GenBank DBLINK
        # Project:... lines.
        #
        # In either case, seems EMBL usess just "PR    Project:..."
        # regardless of the type of ID (old numeric only, or new
        # with alpha prefix), e.g. for CP002497 NCBI now uses:
        #
        # DBLINK      BioProject: PRJNA60715
        #             BioSample: SAMN03081426
        #
        # While EMBL uses:
        #
        # XX
        # PR   Project:PRJNA60715;
        # XX
        #
        # Sorting ensures (new) BioProject:... is before old Project:...
        for xref in sorted(record.dbxrefs):
            if xref.startswith("BioProject:"):
                self._write_single_line("PR", xref[3:] + ";")
                handle.write("XX\n")
                break
            if xref.startswith("Project:"):
                self._write_single_line("PR", xref + ";")
                handle.write("XX\n")
                break

        # TODO - DT lines (date)

        descr = record.description
        if descr == "<unknown description>":
            descr = "."
        self._write_multi_line("DE", descr)
        handle.write("XX\n")

        if "keywords" in record.annotations:
            self._write_keywords(record)

        # Should this be "source" or "organism"?
        self._write_multi_line("OS", self._get_annotation_str(record, "organism"))
        try:
            # List of strings
            taxonomy = "; ".join(record.annotations["taxonomy"]) + "."
        except KeyError:
            taxonomy = "."
        self._write_multi_line("OC", taxonomy)
        handle.write("XX\n")

        if "references" in record.annotations:
            self._write_references(record)

        if "comment" in record.annotations:
            self._write_comment(record)

        handle.write(self.FEATURE_HEADER)
        rec_length = len(record)
        for feature in record.features:
            self._write_feature(feature, rec_length)
        handle.write("XX\n")

        self._write_sequence(record)
        handle.write("//\n")


def EmblIterator(source):
    """Break up an EMBL file into SeqRecord objects.

    Argument source is a file-like object opened in text mode or a path to a file.
    Every section from the LOCUS line to the terminating // becomes
    a single SeqRecord with associated annotation and features.

    Note that for genomes or chromosomes, there is typically only
    one record.

    This gets called internally by Bio.SeqIO for the EMBL file format:

    >>> from Bio import SeqIO
    >>> for record in SeqIO.parse("EMBL/epo_prt_selection.embl", "embl"):
    ...     print(record.id)
    ...
    A00022.1
    A00028.1
    A00031.1
    A00034.1
    A00060.1
    A00071.1
    A00072.1
    A00078.1
    CQ797900.1

    Equivalently,

    >>> with open("EMBL/epo_prt_selection.embl") as handle:
    ...     for record in EmblIterator(handle):
    ...         print(record.id)
    ...
    A00022.1
    A00028.1
    A00031.1
    A00034.1
    A00060.1
    A00071.1
    A00072.1
    A00078.1
    CQ797900.1

    """
    try:
        handle = open(source)
    except TypeError:
        handle = source
        if handle.read(0) != "":
            raise ValueError("EMBL files must be opened in text mode.") from None

    try:
        records = EmblScanner().parse_records(handle)
        yield from records
    finally:
        if handle is not source:
            handle.close()


def EmblCdsFeatureIterator(source, alphabet=Alphabet.generic_protein):
    """Break up a EMBL file into SeqRecord objects for each CDS feature.

    Argument source is a file-like object opened in text mode or a path to a file.

    Every section from the LOCUS line to the terminating // can contain
    many CDS features.  These are returned as with the stated amino acid
    translation sequence (if given).
    """
    try:
        handle = open(source)
    except TypeError:
        handle = source
        if handle.read(0) != "":
            raise ValueError("EMBL files must be opened in text mode.") from None

    try:
        records = EmblScanner().parse_cds_features(handle, alphabet)
        yield from records
    finally:
        if handle is not source:
            handle.close()


class EmblScanner(_InsdcScanner):
    """For extracting chunks of information in EMBL files."""

    RECORD_START = "ID   "
    HEADER_WIDTH = 5
    FEATURE_START_MARKERS = ["FH   Key             Location/Qualifiers", "FH"]
    FEATURE_END_MARKERS = ["XX"]
    FEATURE_QUALIFIER_INDENT = 21
    FEATURE_QUALIFIER_SPACER = "FT" + " " * (FEATURE_QUALIFIER_INDENT - 2)
    SEQUENCE_HEADERS = ["SQ", "CO"]

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
