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
from datetime import datetime

import sys
import warnings
from collections import OrderedDict

from Bio import BiopythonParserWarning, BiopythonWarning, Alphabet, SeqFeature
from Bio.Seq import UnknownSeq
from Bio.SeqIO.InsdcIO import _InsdcWriter, _InsdcScanner


class GenBankScanner(_InsdcScanner):
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

    @staticmethod
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

    @staticmethod
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

    @staticmethod
    def _feed_header_dblink(consumer, data):
        for line in data:
            consumer.dblink(line.strip())

    @staticmethod
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

    @staticmethod
    def _feed_header_common(consumer, data, line_type, consumer_dict):
        data = " ".join(data)
        if line_type == "DEFINITION" and data.endswith("."):
            data = data[:-1]
        getattr(consumer, consumer_dict[line_type])(data)

    @staticmethod
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
                structured_comment_dict.setdefault(
                    structured_comment_key, OrderedDict()
                )
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

    @staticmethod
    def _feed_misc_common(consumer, consumer_dict, key, line):
        if line.startswith(key):
            line = line[len(key)].strip()
            if line:
                getattr(consumer, consumer_dict[key])(line)
        return line

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
            self._feed_header_version(consumer, data[0])
        elif line_type == "DBLINK":
            self._feed_header_dblink(consumer, data)
        elif line_type == "REFERENCE":
            self._feed_header_reference(consumer, data)
        elif line_type == "ORGANISM":
            self._feed_header_organism(consumer, data)
        elif line_type == "COMMENT":
            self._feed_header_comment(
                consumer,
                data,
                self.STRUCTURED_COMMENT_START,
                self.STRUCTURED_COMMENT_DELIM,
                self.STRUCTURED_COMMENT_END,
            )
        elif line_type in self.CONSUMER_DICT:
            self._feed_header_common(consumer, data, line_type, self.CONSUMER_DICT)
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


class GenBankWriter(_InsdcWriter):
    """GenBank writer."""

    HEADER_WIDTH = 12
    QUALIFIER_INDENT = 21
    STRUCTURED_COMMENT_START = "-START##"
    STRUCTURED_COMMENT_END = "-END##"
    STRUCTURED_COMMENT_DELIM = " :: "
    LETTERS_PER_LINE = 60
    SEQUENCE_INDENT = 9

    def _write_single_line(self, tag, text):
        """Write single line in each GenBank record (PRIVATE).

        Used in the 'header' of each GenBank record.
        """
        assert len(tag) < self.HEADER_WIDTH
        if len(text) > self.MAX_WIDTH - self.HEADER_WIDTH:
            if tag:
                warnings.warn(
                    "Annotation %r too long for %r line" % (text, tag), BiopythonWarning
                )
            else:
                # Can't give such a precise warning
                warnings.warn("Annotation %r too long" % text, BiopythonWarning)
        self.handle.write(
            "%s%s\n" % (tag.ljust(self.HEADER_WIDTH), text.replace("\n", " "))
        )

    def _write_multi_line(self, tag, text):
        """Write multiple lines in each GenBank record (PRIVATE).

        Used in the 'header' of each GenBank record.
        """
        # TODO - Do the line spliting while preserving white space?
        max_len = self.MAX_WIDTH - self.HEADER_WIDTH
        lines = self._split_multi_line(text, max_len)
        self._write_single_line(tag, lines[0])
        for line in lines[1:]:
            self._write_single_line("", line)

    def _write_multi_entries(self, tag, text_list):
        # used for DBLINK and any similar later line types.
        # If the list of strings is empty, nothing is written.
        for i, text in enumerate(text_list):
            if i == 0:
                self._write_single_line(tag, text)
            else:
                self._write_single_line("", text)

    @staticmethod
    def _get_date(record):
        default = "01-JAN-1980"
        try:
            date = record.annotations["date"]
        except KeyError:
            return default
        # Cope with a list of one string:
        if isinstance(date, list) and len(date) == 1:
            date = date[0]
        if isinstance(date, datetime):
            date = date.strftime("%d-%b-%Y").upper()

        months = [
            "JAN",
            "FEB",
            "MAR",
            "APR",
            "MAY",
            "JUN",
            "JUL",
            "AUG",
            "SEP",
            "OCT",
            "NOV",
            "DEC",
        ]
        if not isinstance(date, str) or len(date) != 11:
            return default
        try:
            datetime(int(date[-4:]), months.index(date[3:6]) + 1, int(date[0:2]))
        except ValueError:
            date = default
        return date

    @staticmethod
    def _get_data_division(record):
        try:
            division = record.annotations["data_file_division"]
        except KeyError:
            division = "UNK"
        if division in [
            "PRI",
            "ROD",
            "MAM",
            "VRT",
            "INV",
            "PLN",
            "BCT",
            "VRL",
            "PHG",
            "SYN",
            "UNA",
            "EST",
            "PAT",
            "STS",
            "GSS",
            "HTG",
            "HTC",
            "ENV",
            "CON",
            "TSA",
        ]:
            # Good, already GenBank style
            #    PRI - primate sequences
            #    ROD - rodent sequences
            #    MAM - other mammalian sequences
            #    VRT - other vertebrate sequences
            #    INV - invertebrate sequences
            #    PLN - plant, fungal, and algal sequences
            #    BCT - bacterial sequences [plus archea]
            #    VRL - viral sequences
            #    PHG - bacteriophage sequences
            #    SYN - synthetic sequences
            #    UNA - unannotated sequences
            #    EST - EST sequences (expressed sequence tags)
            #    PAT - patent sequences
            #    STS - STS sequences (sequence tagged sites)
            #    GSS - GSS sequences (genome survey sequences)
            #    HTG - HTGS sequences (high throughput genomic sequences)
            #    HTC - HTC sequences (high throughput cDNA sequences)
            #    ENV - Environmental sampling sequences
            #    CON - Constructed sequences
            #    TSA - Transcriptome Shotgun Assembly
            #
            # (plus UNK for unknown)
            pass
        else:
            # See if this is in EMBL style:
            #    Division                 Code
            #    -----------------        ----
            #    Bacteriophage            PHG - common
            #    Environmental Sample     ENV - common
            #    Fungal                   FUN - map to PLN (plants + fungal)
            #    Human                    HUM - map to PRI (primates)
            #    Invertebrate             INV - common
            #    Other Mammal             MAM - common
            #    Other Vertebrate         VRT - common
            #    Mus musculus             MUS - map to ROD (rodent)
            #    Plant                    PLN - common
            #    Prokaryote               PRO - map to BCT (poor name)
            #    Other Rodent             ROD - common
            #    Synthetic                SYN - common
            #    Transgenic               TGN - ??? map to SYN ???
            #    Unclassified             UNC - map to UNK
            #    Viral                    VRL - common
            #
            # (plus XXX for submiting which we can map to UNK)
            embl_to_gbk = {
                "FUN": "PLN",
                "HUM": "PRI",
                "MUS": "ROD",
                "PRO": "BCT",
                "UNC": "UNK",
                "XXX": "UNK",
            }
            try:
                division = embl_to_gbk[division]
            except KeyError:
                division = "UNK"
        assert len(division) == 3
        return division

    def _get_topology(self, record):
        """Set the topology to 'circular', 'linear' if defined (PRIVATE)."""
        max_topology_len = len("circular")

        topology = self._get_annotation_str(record, "topology", default="")
        if topology and len(topology) <= max_topology_len:
            return topology.ljust(max_topology_len)
        else:
            return " " * max_topology_len

    def _write_the_first_line(self, record):
        """Write the LOCUS line (PRIVATE)."""
        locus = record.name
        if not locus or locus == "<unknown name>":
            locus = record.id
        if not locus or locus == "<unknown id>":
            locus = self._get_annotation_str(record, "accession", just_first=True)
        if len(locus) > 16:
            if len(locus) + 1 + len(str(len(record))) > 28:
                # Locus name and record length to long to squeeze in.
                # Per updated GenBank standard (Dec 15, 2018) 229.0
                # the Locus identifier can be any length, and a space
                # is added after the identifier to keep the identifier
                # and length fields separated
                warnings.warn(
                    "Increasing length of locus line to allow "
                    "long name. This will result in fields that "
                    "are not in usual positions.",
                    BiopythonWarning,
                )

        if len(locus.split()) > 1:
            raise ValueError("Invalid whitespace in %r for LOCUS line" % locus)
        if len(record) > 99999999999:
            # As of the GenBank release notes 229.0, the locus line can be
            # any length. However, long locus lines may not be compatible
            # with all software.
            warnings.warn(
                "The sequence length is very long. The LOCUS "
                "line will be increased in length to compensate. "
                "This may cause unexpected behavior.",
                BiopythonWarning,
            )

        # Get the base alphabet (underneath any Gapped or StopCodon encoding)
        a = Alphabet._get_base_alphabet(record.seq.alphabet)
        if not isinstance(a, Alphabet.Alphabet):
            raise TypeError("Invalid alphabet")
        elif isinstance(a, Alphabet.ProteinAlphabet):
            units = "aa"
        elif isinstance(a, Alphabet.NucleotideAlphabet):
            units = "bp"
        else:
            # Must be something like NucleotideAlphabet or
            # just the generic Alphabet (default for fasta files)
            raise ValueError("Need a Nucleotide or Protein alphabet")

        # Get the molecule type
        mol_type = self._get_annotation_str(record, "molecule_type", default=None)
        if mol_type and len(mol_type) > 7:
            # Deal with common cases from EMBL to GenBank
            mol_type = mol_type.replace("unassigned ", "").replace("genomic ", "")
            if len(mol_type) > 7:
                warnings.warn("Molecule type %r too long" % mol_type, BiopythonWarning)
                mol_type = None
        if mol_type in ["protein", "PROTEIN"]:
            mol_type = ""

        if mol_type:
            pass
        elif isinstance(a, Alphabet.ProteinAlphabet):
            mol_type = ""
        elif isinstance(a, Alphabet.DNAAlphabet):
            mol_type = "DNA"
        elif isinstance(a, Alphabet.RNAAlphabet):
            mol_type = "RNA"
        else:
            # Must be something like NucleotideAlphabet or
            # just the generic Alphabet (default for fasta files)
            raise ValueError("Need a DNA, RNA or Protein alphabet")

        topology = self._get_topology(record)

        division = self._get_data_division(record)

        # Accommodate longer header, with long accessions and lengths
        if len(locus) > 16 and len(str(len(record))) > (11 - (len(locus) - 16)):
            name_length = locus + " " + str(len(record))

        # This is the older, standard 80 position header
        else:
            name_length = str(len(record)).rjust(28)
            name_length = locus + name_length[len(locus) :]
            assert len(name_length) == 28, name_length
            assert " " in name_length, name_length

        assert len(units) == 2
        assert len(division) == 3
        line = "LOCUS       %s %s    %s %s %s %s\n" % (
            name_length,
            units,
            mol_type.ljust(7),
            topology,
            division,
            self._get_date(record),
        )
        # Extra long header
        if len(line) > 80:
            splitline = line.split()
            if splitline[3] not in ["bp", "aa"]:
                raise ValueError(
                    "LOCUS line does not contain size units at "
                    "expected position:\n" + line
                )

            if not (
                splitline[4].strip() == ""
                or "DNA" in splitline[4].strip().upper()
                or "RNA" in splitline[4].strip().upper()
            ):
                raise ValueError(
                    "LOCUS line does not contain valid "
                    "sequence type (DNA, RNA, ...):\n" + line
                )

            self.handle.write(line)

        # 80 position header
        else:
            assert len(line) == 79 + 1, repr(line)  # plus one for new line

            # We're bending the rules to allow an identifier over 16 characters
            # if we can steal spaces from the length field:
            # assert line[12:28].rstrip() == locus, \
            #     'LOCUS line does not contain the locus at the expected position:\n' + line
            # assert line[28:29] == " "
            # assert line[29:40].lstrip() == str(len(record)), \
            #     'LOCUS line does not contain the length at the expected position:\n' + line
            assert line[12:40].split() == [locus, str(len(record))], line

            # Tests copied from Bio.GenBank.Scanner
            if line[40:44] not in [" bp ", " aa "]:
                raise ValueError(
                    "LOCUS line does not contain size units at "
                    "expected position:\n" + line
                )
            if line[44:47] not in ["   ", "ss-", "ds-", "ms-"]:
                raise ValueError(
                    "LOCUS line does not have valid strand "
                    "type (Single stranded, ...):\n" + line
                )
            if not (
                line[47:54].strip() == ""
                or "DNA" in line[47:54].strip().upper()
                or "RNA" in line[47:54].strip().upper()
            ):
                raise ValueError(
                    "LOCUS line does not contain valid "
                    "sequence type (DNA, RNA, ...):\n" + line
                )
            if line[54:55] != " ":
                raise ValueError(
                    "LOCUS line does not contain space at position 55:\n" + line
                )
            if line[55:63].strip() not in ["", "linear", "circular"]:
                raise ValueError(
                    "LOCUS line does not contain valid "
                    "entry (linear, circular, ...):\n" + line
                )
            if line[63:64] != " ":
                raise ValueError(
                    "LOCUS line does not contain space at position 64:\n" + line
                )
            if line[67:68] != " ":
                raise ValueError(
                    "LOCUS line does not contain space at position 68:\n" + line
                )
            if line[70:71] != "-":
                raise ValueError(
                    "LOCUS line does not contain - at position 71 in date:\n" + line
                )
            if line[74:75] != "-":
                raise ValueError(
                    "LOCUS line does not contain - at position 75 in date:\n" + line
                )

            self.handle.write(line)

    def _write_references(self, record):
        number = 0
        for ref in record.annotations["references"]:
            if not isinstance(ref, SeqFeature.Reference):
                continue
            number += 1
            data = str(number)
            # TODO - support more complex record reference locations?
            if ref.location and len(ref.location) == 1:
                a = Alphabet._get_base_alphabet(record.seq.alphabet)
                if isinstance(a, Alphabet.ProteinAlphabet):
                    units = "residues"
                else:
                    units = "bases"
                data += "  (%s %i to %i)" % (
                    units,
                    ref.location[0].nofuzzy_start + 1,
                    ref.location[0].nofuzzy_end,
                )
            self._write_single_line("REFERENCE", data)
            if ref.authors:
                # We store the AUTHORS data as a single string
                self._write_multi_line("  AUTHORS", ref.authors)
            if ref.consrtm:
                # We store the consortium as a single string
                self._write_multi_line("  CONSRTM", ref.consrtm)
            if ref.title:
                # We store the title as a single string
                self._write_multi_line("  TITLE", ref.title)
            if ref.journal:
                # We store this as a single string - holds the journal name,
                # volume, year, and page numbers of the citation
                self._write_multi_line("  JOURNAL", ref.journal)
            if ref.medline_id:
                # This line type is obsolete and was removed from the GenBank
                # flatfile format in April 2005. Should we write it?
                # Note this has a two space indent:
                self._write_multi_line("  MEDLINE", ref.medline_id)
            if ref.pubmed_id:
                # Note this has a THREE space indent:
                self._write_multi_line("   PUBMED", ref.pubmed_id)
            if ref.comment:
                self._write_multi_line("  REMARK", ref.comment)

    def _write_comment(self, record):
        # This is a bit complicated due to the range of possible
        # ways people might have done their annotation...
        # Currently the parser uses a single string with newlines.
        # A list of lines is also reasonable.
        # A single (long) string is perhaps the most natural of all.
        # This means we may need to deal with line wrapping.
        lines = []
        if "structured_comment" in record.annotations:
            comment = record.annotations["structured_comment"]
            # Find max length of keys for equal padded printing
            padding = 0
            for key, data in comment.items():
                for subkey, subdata in data.items():
                    padding = len(subkey) if len(subkey) > padding else padding
            # Construct output
            for key, data in comment.items():
                lines.append(f"##{key}{self.STRUCTURED_COMMENT_START}")
                for subkey, subdata in data.items():
                    spaces = " " * (padding - len(subkey))
                    lines.append(
                        f"{subkey}{spaces}{self.STRUCTURED_COMMENT_DELIM}{subdata}"
                    )
                lines.append(f"##{key}{self.STRUCTURED_COMMENT_END}")
        if "comment" in record.annotations:
            comment = record.annotations["comment"]
            if isinstance(comment, str):
                lines += comment.split("\n")
            elif isinstance(comment, (list, tuple)):
                lines += list(comment)
            else:
                raise ValueError("Could not understand comment annotation")
        self._write_multi_line("COMMENT", lines[0])
        for line in lines[1:]:
            self._write_multi_line("", line)

    def _write_contig(self, record):
        max_len = self.MAX_WIDTH - self.HEADER_WIDTH
        lines = self._split_contig(record, max_len)
        self._write_single_line("CONTIG", lines[0])
        for text in lines[1:]:
            self._write_single_line("", text)

    def _write_sequence(self, record):
        # Loosely based on code from Howard Salis
        # TODO - Force lower case?

        if isinstance(record.seq, UnknownSeq):
            # We have already recorded the length, and there is no need
            # to record a long sequence of NNNNNNN...NNN or whatever.
            if "contig" in record.annotations:
                self._write_contig(record)
            else:
                self.handle.write("ORIGIN\n")
            return

        # Catches sequence being None:
        data = self._get_seq_string(record).lower()
        seq_len = len(data)
        self.handle.write("ORIGIN\n")
        for line_number in range(0, seq_len, self.LETTERS_PER_LINE):
            self.handle.write(str(line_number + 1).rjust(self.SEQUENCE_INDENT))
            for words in range(
                line_number, min(line_number + self.LETTERS_PER_LINE, seq_len), 10
            ):
                self.handle.write(" %s" % data[words : words + 10])
            self.handle.write("\n")

    def write_record(self, record):
        """Write a single record to the output file."""
        handle = self.handle
        self._write_the_first_line(record)

        default = record.id
        if default.count(".") == 1 and default[default.index(".") + 1 :].isdigit():
            # Good, looks like accesion.version and not something
            # else like identifier.start-end
            default = record.id.split(".", 1)[0]
        accession = self._get_annotation_str(
            record, "accession", default, just_first=True
        )
        acc_with_version = accession
        if record.id.startswith(accession + "."):
            try:
                acc_with_version = "%s.%i" % (
                    accession,
                    int(record.id.split(".", 1)[1]),
                )
            except ValueError:
                pass
        gi = self._get_annotation_str(record, "gi", just_first=True)

        descr = record.description
        if descr == "<unknown description>":
            descr = ""  # Trailing dot will be added later

        # The DEFINITION field must end with a period
        # see ftp://ftp.ncbi.nih.gov/genbank/gbrel.txt [3.4.5]
        # and discussion https://github.com/biopython/biopython/pull/616
        # So let's add a period
        descr += "."
        self._write_multi_line("DEFINITION", descr)

        self._write_single_line("ACCESSION", accession)
        if gi != ".":
            self._write_single_line("VERSION", "%s  GI:%s" % (acc_with_version, gi))
        else:
            self._write_single_line("VERSION", "%s" % acc_with_version)

        # The NCBI initially expected two types of link,
        # e.g. "Project:28471" and "Trace Assembly Archive:123456"
        #
        # This changed and at some point the formatting switched to
        # include a space after the colon, e.g.
        #
        # LOCUS       NC_000011               1606 bp    DNA     linear   CON 06-JUN-2016
        # DEFINITION  Homo sapiens chromosome 11, GRCh38.p7 Primary Assembly.
        # ACCESSION   NC_000011 REGION: complement(5225466..5227071) GPC_000001303
        # VERSION     NC_000011.10  GI:568815587
        # DBLINK      BioProject: PRJNA168
        #             Assembly: GCF_000001405.33
        # ...
        #
        # Or,
        #
        # LOCUS       JU120277                1044 bp    mRNA    linear   TSA 27-NOV-2012
        # DEFINITION  TSA: Tupaia chinensis tbc000002.Tuchadli mRNA sequence.
        # ACCESSION   JU120277
        # VERSION     JU120277.1  GI:379775257
        # DBLINK      BioProject: PRJNA87013
        #             Sequence Read Archive: SRR433859
        # ...
        dbxrefs_with_space = []
        for x in record.dbxrefs:
            if ": " not in x:
                x = x.replace(":", ": ")
            dbxrefs_with_space.append(x)
        self._write_multi_entries("DBLINK", dbxrefs_with_space)
        del dbxrefs_with_space

        try:
            # List of strings
            # Keywords should be given separated with semi colons,
            keywords = "; ".join(record.annotations["keywords"])
            # with a trailing period:
            if not keywords.endswith("."):
                keywords += "."
        except KeyError:
            # If no keywords, there should be just a period:
            keywords = "."
        self._write_multi_line("KEYWORDS", keywords)

        if "segment" in record.annotations:
            # Deal with SEGMENT line found only in segmented records,
            # e.g. AH000819
            segment = record.annotations["segment"]
            if isinstance(segment, list):
                assert len(segment) == 1, segment
                segment = segment[0]
            self._write_single_line("SEGMENT", segment)

        self._write_multi_line("SOURCE", self._get_annotation_str(record, "source"))
        # The ORGANISM line MUST be a single line, as any continuation is the taxonomy
        org = self._get_annotation_str(record, "organism")
        if len(org) > self.MAX_WIDTH - self.HEADER_WIDTH:
            org = org[: self.MAX_WIDTH - self.HEADER_WIDTH - 4] + "..."
        self._write_single_line("  ORGANISM", org)
        try:
            # List of strings
            # Taxonomy should be given separated with semi colons,
            taxonomy = "; ".join(record.annotations["taxonomy"])
            # with a trailing period:
            if not taxonomy.endswith("."):
                taxonomy += "."
        except KeyError:
            taxonomy = "."
        self._write_multi_line("", taxonomy)

        if "references" in record.annotations:
            self._write_references(record)

        if (
            "comment" in record.annotations
            or "structured_comment" in record.annotations
        ):
            self._write_comment(record)

        handle.write("FEATURES             Location/Qualifiers\n")
        rec_length = len(record)
        for feature in record.features:
            self._write_feature(feature, rec_length)
        self._write_sequence(record)
        handle.write("//\n")


def GenBankIterator(source):
    """Break up a Genbank file into SeqRecord objects.

    Argument source is a file-like object opened in text mode or a path to a file.

    Every section from the LOCUS line to the terminating // becomes
    a single SeqRecord with associated annotation and features.

    Note that for genomes or chromosomes, there is typically only
    one record.

    This gets called internally by Bio.SeqIO for the GenBank file format:

    >>> from Bio import SeqIO
    >>> for record in SeqIO.parse("GenBank/cor6_6.gb", "gb"):
    ...     print(record.id)
    ...
    X55053.1
    X62281.1
    M81224.1
    AJ237582.1
    L31939.1
    AF297471.1

    Equivalently,

    >>> with open("GenBank/cor6_6.gb") as handle:
    ...     for record in GenBankIterator(handle):
    ...         print(record.id)
    ...
    X55053.1
    X62281.1
    M81224.1
    AJ237582.1
    L31939.1
    AF297471.1

    """
    try:
        handle = open(source)
    except TypeError:
        handle = source
        if handle.read(0) != "":
            raise ValueError("GenBank files must be opened in text mode.") from None

    try:
        records = GenBankScanner().parse_records(handle)
        yield from records
    finally:
        if handle is not source:
            handle.close()


def GenBankCdsFeatureIterator(source, alphabet=Alphabet.generic_protein):
    """Break up a Genbank file into SeqRecord objects for each CDS feature.

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
            raise ValueError("GenBank files must be opened in text mode.") from None

    try:
        records = GenBankScanner().parse_cds_features(handle, alphabet)
        yield from records
    finally:
        if handle is not source:
            handle.close()
