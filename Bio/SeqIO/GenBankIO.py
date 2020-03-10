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


class _BaseGenBankConsumer:
    """Abstract GenBank consumer providing useful general functions (PRIVATE).

    This just helps to eliminate some duplication in things that most
    GenBank consumers want to do.
    """

    # Special keys in GenBank records that we should remove spaces from
    # For instance, \translation keys have values which are proteins and
    # should have spaces and newlines removed from them. This class
    # attribute gives us more control over specific formatting problems.
    remove_space_keys = ["translation"]

    def __init__(self):
        pass

    @staticmethod
    def _split_keywords(keyword_string):
        """Split a string of keywords into a nice clean list (PRIVATE)."""
        # process the keywords into a python list
        if keyword_string == "" or keyword_string == ".":
            keywords = ""
        elif keyword_string[-1] == ".":
            keywords = keyword_string[:-1]
        else:
            keywords = keyword_string
        keyword_list = keywords.split(";")
        return [x.strip() for x in keyword_list]

    @staticmethod
    def _split_accessions(accession_string):
        """Split a string of accession numbers into a list (PRIVATE)."""
        # first replace all line feeds with spaces
        # Also, EMBL style accessions are split with ';'
        accession = accession_string.replace("\n", " ").replace(";", " ")

        return [x.strip() for x in accession.split() if x.strip()]

    @staticmethod
    def _split_taxonomy(taxonomy_string):
        """Split a string with taxonomy info into a list (PRIVATE)."""
        if not taxonomy_string or taxonomy_string == ".":
            # Missing data, no taxonomy
            return []

        if taxonomy_string[-1] == ".":
            tax_info = taxonomy_string[:-1]
        else:
            tax_info = taxonomy_string
        tax_list = tax_info.split(";")
        new_tax_list = []
        for tax_item in tax_list:
            new_items = tax_item.split("\n")
            new_tax_list.extend(new_items)
        while "" in new_tax_list:
            new_tax_list.remove("")
        return [x.strip() for x in new_tax_list]

    @staticmethod
    def _clean_location(location_string):
        """Clean whitespace out of a location string (PRIVATE).

        The location parser isn't a fan of whitespace, so we clean it out
        before feeding it into the parser.
        """
        # Originally this imported string.whitespace and did a replace
        # via a loop.  It's simpler to just split on whitespace and rejoin
        # the string - and this avoids importing string too.  See Bug 2684.
        return "".join(location_string.split())

    @staticmethod
    def _remove_newlines(text):
        """Remove any newlines in the passed text, returning the new string (PRIVATE)."""
        # get rid of newlines in the qualifier value
        newlines = ["\n", "\r"]
        for ws in newlines:
            text = text.replace(ws, "")

        return text

    @staticmethod
    def _normalize_spaces(text):
        """Replace multiple spaces in the passed text with single spaces (PRIVATE)."""
        # get rid of excessive spaces
        return " ".join(x for x in text.split(" ") if x)

    @staticmethod
    def _remove_spaces(text):
        """Remove all spaces from the passed text (PRIVATE)."""
        return text.replace(" ", "")

    @staticmethod
    def _convert_to_python_numbers(start, end):
        """Convert a start and end range to python notation (PRIVATE).

        In GenBank, starts and ends are defined in "biological" coordinates,
        where 1 is the first base and [i, j] means to include both i and j.

        In python, 0 is the first base and [i, j] means to include i, but
        not j.

        So, to convert "biological" to python coordinates, we need to
        subtract 1 from the start, and leave the end and things should
        be converted happily.
        """
        new_start = start - 1
        new_end = end

        return new_start, new_end


class _FeatureConsumer(_BaseGenBankConsumer):
    """Create a SeqRecord object with Features to return (PRIVATE).

    Attributes:
     - use_fuzziness - specify whether or not to parse with fuzziness in
       feature locations.
     - feature_cleaner - a class that will be used to provide specialized
       cleaning-up of feature values.

    """

    def __init__(self, use_fuzziness, feature_cleaner=None):
        from Bio.SeqRecord import SeqRecord

        _BaseGenBankConsumer.__init__(self)
        self.data = SeqRecord(None, id=None)
        self.data.id = None
        self.data.description = ""

        self._use_fuzziness = use_fuzziness
        self._feature_cleaner = feature_cleaner

        self._seq_type = ""
        self._seq_data = []
        self._cur_reference = None
        self._cur_feature = None
        self._expected_size = None

    def locus(self, locus_name):
        """Set the locus name is set as the name of the Sequence."""
        self.data.name = locus_name

    def size(self, content):
        """Record the sequence length."""
        self._expected_size = int(content)

    def residue_type(self, type):
        """Record the sequence type (SEMI-OBSOLETE).

        This reflects the fact that the topology (linear/circular) and
        molecule type (e.g. DNA vs RNA) were a single field in early
        files. Current GenBank/EMBL files have two fields.
        """
        self._seq_type = type.strip()

    def topology(self, topology):  # noqa: D402
        """Validate and record sequence topology (linear or circular as strings)."""
        if topology:
            if topology not in ["linear", "circular"]:
                raise ParserFailureError(
                    "Unexpected topology %r should be linear or circular" % topology
                )
            self.data.annotations["topology"] = topology

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

            self.data.annotations["molecule_type"] = mol_type

    def data_file_division(self, division):
        self.data.annotations["data_file_division"] = division

    def date(self, submit_date):
        self.data.annotations["date"] = submit_date

    def definition(self, definition):
        """Set the definition as the description of the sequence."""
        if self.data.description:
            # Append to any existing description
            # e.g. EMBL files with two DE lines.
            self.data.description += " " + definition
        else:
            self.data.description = definition

    def accession(self, acc_num):
        """Set the accession number as the id of the sequence.

        If we have multiple accession numbers, the first one passed is
        used.
        """
        new_acc_nums = self._split_accessions(acc_num)

        # Also record them ALL in the annotations
        try:
            # On the off chance there was more than one accession line:
            for acc in new_acc_nums:
                # Prevent repeat entries
                if acc not in self.data.annotations["accessions"]:
                    self.data.annotations["accessions"].append(acc)
        except KeyError:
            self.data.annotations["accessions"] = new_acc_nums

        # if we haven't set the id information yet, add the first acc num
        if not self.data.id:
            if len(new_acc_nums) > 0:
                # self.data.id = new_acc_nums[0]
                # Use the FIRST accession as the ID, not the first on this line!
                self.data.id = self.data.annotations["accessions"][0]

    def tls(self, content):
        self.data.annotations["tls"] = content.split("-")

    def tsa(self, content):
        self.data.annotations["tsa"] = content.split("-")

    def wgs(self, content):
        self.data.annotations["wgs"] = content.split("-")

    def add_wgs_scafld(self, content):
        self.data.annotations.setdefault("wgs_scafld", []).append(content.split("-"))

    def nid(self, content):
        self.data.annotations["nid"] = content

    def pid(self, content):
        self.data.annotations["pid"] = content

    def version(self, version_id):
        # Want to use the versioned accession as the record.id
        # This comes from the VERSION line in GenBank files, or the
        # obsolete SV line in EMBL.  For the new EMBL files we need
        # both the version suffix from the ID line and the accession
        # from the AC line.
        if version_id.count(".") == 1 and version_id.split(".")[1].isdigit():
            self.accession(version_id.split(".")[0])
            self.version_suffix(version_id.split(".")[1])
        elif version_id:
            # For backwards compatibility...
            self.data.id = version_id

    def project(self, content):
        """Handle the information from the PROJECT line as a list of projects.

        e.g.::

            PROJECT     GenomeProject:28471

        or::

            PROJECT     GenomeProject:13543  GenomeProject:99999

        This is stored as dbxrefs in the SeqRecord to be consistent with the
        projected switch of this line to DBLINK in future GenBank versions.
        Note the NCBI plan to replace "GenomeProject:28471" with the shorter
        "Project:28471" as part of this transition.
        """
        content = content.replace("GenomeProject:", "Project:")
        self.data.dbxrefs.extend(p for p in content.split() if p)

    def dblink(self, content):
        """Store DBLINK cross references as dbxrefs in our record object.

        This line type is expected to replace the PROJECT line in 2009. e.g.

        During transition::

            PROJECT     GenomeProject:28471
            DBLINK      Project:28471
                        Trace Assembly Archive:123456

        Once the project line is dropped::

            DBLINK      Project:28471
                        Trace Assembly Archive:123456

        Note GenomeProject -> Project.

        We'll have to see some real examples to be sure, but based on the
        above example we can expect one reference per line.

        Note that at some point the NCBI have included an extra space, e.g.::

            DBLINK      Project: 28471

        """
        # During the transition period with both PROJECT and DBLINK lines,
        # we don't want to add the same cross reference twice.
        while ": " in content:
            content = content.replace(": ", ":")
        if content.strip() not in self.data.dbxrefs:
            self.data.dbxrefs.append(content.strip())

    def version_suffix(self, version):
        """Set the version to overwrite the id.

        Since the version provides the same information as the accession
        number, plus some extra info, we set this as the id if we have
        a version.
        """
        # e.g. GenBank line:
        # VERSION     U49845.1  GI:1293613
        # or the obsolete EMBL line:
        # SV   U49845.1
        # Scanner calls consumer.version("U49845.1")
        # which then calls consumer.version_suffix(1)
        #
        # e.g. EMBL new line:
        # ID   X56734; SV 1; linear; mRNA; STD; PLN; 1859 BP.
        # Scanner calls consumer.version_suffix(1)
        assert version.isdigit()
        self.data.annotations["sequence_version"] = int(version)

    def db_source(self, content):
        self.data.annotations["db_source"] = content.rstrip()

    def gi(self, content):
        self.data.annotations["gi"] = content

    def keywords(self, content):
        if "keywords" in self.data.annotations:
            # Multi-line keywords, append to list
            # Note EMBL states "A keyword is never split between lines."
            self.data.annotations["keywords"].extend(self._split_keywords(content))
        else:
            self.data.annotations["keywords"] = self._split_keywords(content)

    def segment(self, content):
        self.data.annotations["segment"] = content

    def source(self, content):
        # Note that some software (e.g. VectorNTI) may produce an empty
        # source (rather than using a dot/period as might be expected).
        if content == "":
            source_info = ""
        elif content[-1] == ".":
            source_info = content[:-1]
        else:
            source_info = content
        self.data.annotations["source"] = source_info

    def organism(self, content):
        self.data.annotations["organism"] = content

    def taxonomy(self, content):
        """Record (another line of) the taxonomy lineage."""
        lineage = self._split_taxonomy(content)
        try:
            self.data.annotations["taxonomy"].extend(lineage)
        except KeyError:
            self.data.annotations["taxonomy"] = lineage

    def reference_num(self, content):
        """Signal the beginning of a new reference object."""
        # if we have a current reference that hasn't been added to
        # the list of references, add it.
        if self._cur_reference is not None:
            self.data.annotations["references"].append(self._cur_reference)
        else:
            self.data.annotations["references"] = []

        self._cur_reference = SeqFeature.Reference()

    def reference_bases(self, content):
        """Attempt to determine the sequence region the reference entails.

        Possible types of information we may have to deal with:

        (bases 1 to 86436)
        (sites)
        (bases 1 to 105654; 110423 to 111122)
        1  (residues 1 to 182)
        """
        # first remove the parentheses
        assert content.endswith(")"), content
        ref_base_info = content[1:-1]

        all_locations = []
        # parse if we've got 'bases' and 'to'
        if "bases" in ref_base_info and "to" in ref_base_info:
            # get rid of the beginning 'bases'
            ref_base_info = ref_base_info[5:]
            locations = self._split_reference_locations(ref_base_info)
            all_locations.extend(locations)
        elif "residues" in ref_base_info and "to" in ref_base_info:
            residues_start = ref_base_info.find("residues")
            # get only the information after "residues"
            ref_base_info = ref_base_info[(residues_start + len("residues ")) :]
            locations = self._split_reference_locations(ref_base_info)
            all_locations.extend(locations)

        # make sure if we are not finding information then we have
        # the string 'sites' or the string 'bases'
        elif ref_base_info == "sites" or ref_base_info.strip() == "bases":
            pass
        # otherwise raise an error
        else:
            raise ValueError(
                "Could not parse base info %s in record %s"
                % (ref_base_info, self.data.id)
            )

        self._cur_reference.location = all_locations

    def _split_reference_locations(self, location_string):
        """Get reference locations out of a string of reference information (PRIVATE).

        The passed string should be of the form::

            1 to 20; 20 to 100

        This splits the information out and returns a list of location objects
        based on the reference locations.
        """
        # split possibly multiple locations using the ';'
        all_base_info = location_string.split(";")

        new_locations = []
        for base_info in all_base_info:
            start, end = base_info.split("to")
            new_start, new_end = self._convert_to_python_numbers(
                int(start.strip()), int(end.strip())
            )
            this_location = SeqFeature.FeatureLocation(new_start, new_end)
            new_locations.append(this_location)
        return new_locations

    def authors(self, content):
        if self._cur_reference.authors:
            self._cur_reference.authors += " " + content
        else:
            self._cur_reference.authors = content

    def consrtm(self, content):
        if self._cur_reference.consrtm:
            self._cur_reference.consrtm += " " + content
        else:
            self._cur_reference.consrtm = content

    def title(self, content):
        if self._cur_reference is None:
            warnings.warn(
                "GenBank TITLE line without REFERENCE line.", BiopythonParserWarning
            )
        elif self._cur_reference.title:
            self._cur_reference.title += " " + content
        else:
            self._cur_reference.title = content

    def journal(self, content):
        if self._cur_reference.journal:
            self._cur_reference.journal += " " + content
        else:
            self._cur_reference.journal = content

    def medline_id(self, content):
        self._cur_reference.medline_id = content

    def pubmed_id(self, content):
        self._cur_reference.pubmed_id = content

    def remark(self, content):
        """Deal with a reference comment."""
        if self._cur_reference.comment:
            self._cur_reference.comment += " " + content
        else:
            self._cur_reference.comment = content

    def comment(self, content):
        try:
            self.data.annotations["comment"] += "\n" + "\n".join(content)
        except KeyError:
            self.data.annotations["comment"] = "\n".join(content)

    def structured_comment(self, content):
        self.data.annotations["structured_comment"] = content

    def features_line(self, content):
        """Get ready for the feature table when we reach the FEATURE line."""
        self.start_feature_table()

    def start_feature_table(self):
        """Indicate we've got to the start of the feature table."""
        # make sure we've added on our last reference object
        if self._cur_reference is not None:
            self.data.annotations["references"].append(self._cur_reference)
            self._cur_reference = None

    def feature_key(self, content):
        # start a new feature
        self._cur_feature = SeqFeature.SeqFeature()
        self._cur_feature.type = content
        self.data.features.append(self._cur_feature)

    def location(self, content):
        """Parse out location information from the location string.

        This uses simple Python code with some regular expressions to do the
        parsing, and then translates the results into appropriate objects.
        """
        # clean up newlines and other whitespace inside the location before
        # parsing - locations should have no whitespace whatsoever
        location_line = self._clean_location(content)

        # Older records have junk like replace(266,"c") in the
        # location line. Newer records just replace this with
        # the number 266 and have the information in a more reasonable
        # place. So we'll just grab out the number and feed this to the
        # parser. We shouldn't really be losing any info this way.
        if "replace" in location_line:
            comma_pos = location_line.find(",")
            location_line = location_line[8:comma_pos]

        cur_feature = self._cur_feature

        # Handle top level complement here for speed
        if location_line.startswith("complement("):
            assert location_line.endswith(")")
            location_line = location_line[11:-1]
            strand = -1
        elif "PROTEIN" in self._seq_type.upper():
            strand = None
        else:
            # Assume nucleotide otherwise feature strand for
            # GenBank files with bad LOCUS lines set to None
            strand = 1

        # Special case handling of the most common cases for speed
        if _re_simple_location.match(location_line):
            # e.g. "123..456"
            s, e = location_line.split("..")
            try:
                cur_feature.location = SeqFeature.FeatureLocation(
                    int(s) - 1, int(e), strand
                )
            except ValueError:
                # Could be non-integers, more likely bad origin wrapping
                cur_feature.location = _loc(
                    location_line,
                    self._expected_size,
                    strand,
                    seq_type=self._seq_type.lower(),
                )
            return

        if ",)" in location_line:
            warnings.warn(
                "Dropping trailing comma in malformed feature location",
                BiopythonParserWarning,
            )
            location_line = location_line.replace(",)", ")")

        if _solo_bond.search(location_line):
            # e.g. bond(196)
            # e.g. join(bond(284),bond(305),bond(309),bond(305))
            warnings.warn(
                "Dropping bond qualifier in feature location", BiopythonParserWarning
            )
            # There ought to be a better way to do this...
            for x in _solo_bond.finditer(location_line):
                x = x.group()
                location_line = location_line.replace(x, x[5:-1])

        if _re_simple_compound.match(location_line):
            # e.g. join(<123..456,480..>500)
            i = location_line.find("(")
            # cur_feature.location_operator = location_line[:i]
            # we can split on the comma because these are simple locations
            locs = []
            for part in location_line[i + 1 : -1].split(","):
                s, e = part.split("..")

                try:
                    locs.append(SeqFeature.FeatureLocation(int(s) - 1, int(e), strand))
                except ValueError:
                    # Could be non-integers, more likely bad origin wrapping

                    # In the case of bad origin wrapping, _loc will return
                    # a CompoundLocation. CompoundLocation.parts returns a
                    # list of the FeatureLocation objects inside the
                    # CompoundLocation.
                    locs.extend(
                        _loc(
                            part, self._expected_size, strand, self._seq_type.lower()
                        ).parts
                    )

            if len(locs) < 2:
                # The CompoundLocation will raise a ValueError here!
                warnings.warn(
                    "Should have at least 2 parts for compound location",
                    BiopythonParserWarning,
                )
                cur_feature.location = None
                return
            if strand == -1:
                cur_feature.location = SeqFeature.CompoundLocation(
                    locs[::-1], operator=location_line[:i]
                )
            else:
                cur_feature.location = SeqFeature.CompoundLocation(
                    locs, operator=location_line[:i]
                )
            return

        # Handle the general case with more complex regular expressions
        if _re_complex_location.match(location_line):
            # e.g. "AL121804.2:41..610"
            cur_feature.location = _loc(
                location_line,
                self._expected_size,
                strand,
                seq_type=self._seq_type.lower(),
            )
            return

        if _re_complex_compound.match(location_line):
            i = location_line.find("(")
            # cur_feature.location_operator = location_line[:i]
            # Can't split on the comma because of positions like one-of(1,2,3)
            locs = []
            for part in _split_compound_loc(location_line[i + 1 : -1]):
                if part.startswith("complement("):
                    assert part[-1] == ")"
                    part = part[11:-1]
                    assert strand != -1, "Double complement?"
                    part_strand = -1
                else:
                    part_strand = strand
                try:
                    # There is likely a problem with origin wrapping.
                    # Using _loc to return a CompoundLocation of the
                    # wrapped feature and returning the two FeatureLocation
                    # objects to extend to the list of feature locations.
                    loc = _loc(
                        part,
                        self._expected_size,
                        part_strand,
                        seq_type=self._seq_type.lower(),
                    ).parts

                except ValueError:
                    print(location_line)
                    print(part)
                    raise
                # loc will be a list of one or two FeatureLocation items.
                locs.extend(loc)
            # Historically a join on the reverse strand has been represented
            # in Biopython with both the parent SeqFeature and its children
            # (the exons for a CDS) all given a strand of -1.  Likewise, for
            # a join feature on the forward strand they all have strand +1.
            # However, we must also consider evil mixed strand examples like
            # this, join(complement(69611..69724),139856..140087,140625..140650)
            if strand == -1:
                # Whole thing was wrapped in complement(...)
                for l in locs:
                    assert l.strand == -1
                # Reverse the backwards order used in GenBank files
                # with complement(join(...))
                cur_feature.location = SeqFeature.CompoundLocation(
                    locs[::-1], operator=location_line[:i]
                )
            else:
                cur_feature.location = SeqFeature.CompoundLocation(
                    locs, operator=location_line[:i]
                )
            return
        # Not recognised
        if "order" in location_line and "join" in location_line:
            # See Bug 3197
            msg = (
                'Combinations of "join" and "order" within the same '
                "location (nested operators) are illegal:\n" + location_line
            )
            raise LocationParserError(msg)
        # This used to be an error....
        cur_feature.location = None
        warnings.warn(
            BiopythonParserWarning(
                "Couldn't parse feature location: %r" % location_line
            )
        )

    def feature_qualifier(self, key, value):
        """When we get a qualifier key and its value.

        Can receive None, since you can have valueless keys such as /pseudo
        """
        # Hack to try to preserve historical behaviour of /pseudo etc
        if value is None:
            # if the key doesn't exist yet, add an empty string
            if key not in self._cur_feature.qualifiers:
                self._cur_feature.qualifiers[key] = [""]
                return
            # otherwise just skip this key
            return

        # Remove enclosing quotation marks
        value = re.sub('^"|"$', "", value)

        # Handle NCBI escaping
        # Warn if escaping is not according to standard
        if re.search(r'[^"]"[^"]|^"[^"]|[^"]"$', value):
            warnings.warn(
                'The NCBI states double-quote characters like " should be escaped as "" '
                "(two double - quotes), but here it was not: %r" % value,
                BiopythonParserWarning,
            )
        # Undo escaping, repeated double quotes -> one double quote
        value = value.replace('""', '"')

        if self._feature_cleaner is not None:
            value = self._feature_cleaner.clean_value(key, value)

        # if the qualifier name exists, append the value
        if key in self._cur_feature.qualifiers:
            self._cur_feature.qualifiers[key].append(value)
        # otherwise start a new list of the key with its values
        else:
            self._cur_feature.qualifiers[key] = [value]

    def feature_qualifier_name(self, content_list):
        """Use feature_qualifier instead (OBSOLETE)."""
        raise NotImplementedError("Use the feature_qualifier method instead.")

    def feature_qualifier_description(self, content):
        """Use feature_qualifier instead (OBSOLETE)."""
        raise NotImplementedError("Use the feature_qualifier method instead.")

    def contig_location(self, content):
        """Deal with CONTIG information."""
        # Historically this was stored as a SeqFeature object, but it was
        # stored under record.annotations["contig"] and not under
        # record.features with the other SeqFeature objects.
        #
        # The CONTIG location line can include additional tokens like
        # Gap(), Gap(100) or Gap(unk100) which are not used in the feature
        # location lines, so storing it using SeqFeature based location
        # objects is difficult.
        #
        # We now store this a string, which means for BioSQL we are now in
        # much better agreement with how BioPerl records the CONTIG line
        # in the database.
        #
        # NOTE - This code assumes the scanner will return all the CONTIG
        # lines already combined into one long string!
        self.data.annotations["contig"] = content

    def origin_name(self, content):
        pass

    def base_count(self, content):
        pass

    def base_number(self, content):
        pass

    def sequence(self, content):
        """Add up sequence information as we get it.

        To try and make things speedier, this puts all of the strings
        into a list of strings, and then uses string.join later to put
        them together. Supposedly, this is a big time savings
        """
        assert " " not in content
        self._seq_data.append(content.upper())

    def record_end(self, content):
        """Clean up when we've finished the record."""
        from Bio import Alphabet
        from Bio.Alphabet import IUPAC
        from Bio.Seq import Seq, UnknownSeq

        # Try and append the version number to the accession for the full id
        if not self.data.id:
            if "accessions" in self.data.annotations:
                raise ValueError(
                    "Problem adding version number to accession: "
                    + str(self.data.annotations["accessions"])
                )
            self.data.id = self.data.name  # Good fall back?
        elif self.data.id.count(".") == 0:
            try:
                self.data.id += ".%i" % self.data.annotations["sequence_version"]
            except KeyError:
                pass

        # add the sequence information
        # first, determine the alphabet
        # we default to an generic alphabet if we don't have a
        # seq type or have strange sequence information.
        seq_alphabet = Alphabet.generic_alphabet

        # now set the sequence
        sequence = "".join(self._seq_data)

        if (
            self._expected_size is not None
            and len(sequence) != 0
            and self._expected_size != len(sequence)
        ):
            warnings.warn(
                "Expected sequence length %i, found %i (%s)."
                % (self._expected_size, len(sequence), self.data.id),
                BiopythonParserWarning,
            )

        if self._seq_type:
            # mRNA is really also DNA, since it is actually cDNA
            if "DNA" in self._seq_type.upper() or "MRNA" in self._seq_type.upper():
                seq_alphabet = IUPAC.ambiguous_dna
            # are there ever really RNA sequences in GenBank?
            elif "RNA" in self._seq_type.upper():
                # Even for data which was from RNA, the sequence string
                # is usually given as DNA (T not U).  Bug 2408
                if "T" in sequence and "U" not in sequence:
                    seq_alphabet = IUPAC.ambiguous_dna
                else:
                    seq_alphabet = IUPAC.ambiguous_rna
            elif (
                "PROTEIN" in self._seq_type.upper() or self._seq_type == "PRT"
            ):  # PRT is used in EMBL-bank for patents
                seq_alphabet = IUPAC.protein  # or extended protein?
            # work around ugly GenBank records which have circular or
            # linear but no indication of sequence type
            elif self._seq_type in ["circular", "linear", "unspecified"]:
                pass
            # we have a bug if we get here
            else:
                raise ValueError(
                    "Could not determine alphabet for seq_type %s" % self._seq_type
                )

        if not sequence and self._expected_size:
            self.data.seq = UnknownSeq(self._expected_size, seq_alphabet)
        else:
            self.data.seq = Seq(sequence, seq_alphabet)


_simple_location = r"\d+\.\.\d+"
_re_simple_location = re.compile(r"^%s$" % _simple_location)
_solo_location = r"[<>]?\d+"
_solo_bond = re.compile(r"bond\(%s\)" % _solo_location)
_re_simple_compound = re.compile(
    r"^(join|order|bond)\(%s(,%s)*\)$" % (_simple_location, _simple_location)
)
_pair_location = r"[<>]?\d+\.\.[<>]?\d+"
_between_location = r"\d+\^\d+"
_within_position = r"\(\d+\.\d+\)"
_within_location = r"([<>]?\d+|%s)\.\.([<>]?\d+|%s)" % (
    _within_position,
    _within_position,
)
_oneof_position = r"one\-of\(\d+(,\d+)+\)"
_oneof_location = r"([<>]?\d+|%s)\.\.([<>]?\d+|%s)" % (_oneof_position, _oneof_position)
_complex_location = r"([a-zA-Z][a-zA-Z0-9_\.\|]*[a-zA-Z0-9]?\:)?(%s|%s|%s|%s|%s)" % (
    _pair_location,
    _solo_location,
    _between_location,
    _within_location,
    _oneof_location,
)
_re_complex_location = re.compile(r"^%s$" % _complex_location)
_possibly_complemented_complex_location = r"(%s|complement\(%s\))" % (
    _complex_location,
    _complex_location,
)
_re_complex_compound = re.compile(
    r"^(join|order|bond)\(%s(,%s)*\)$"
    % (_possibly_complemented_complex_location, _possibly_complemented_complex_location)
)


def _split_compound_loc(compound_loc):
    """Split a tricky compound location string (PRIVATE).

    >>> list(_split_compound_loc("123..145"))
    ['123..145']
    >>> list(_split_compound_loc("123..145,200..209"))
    ['123..145', '200..209']
    >>> list(_split_compound_loc("one-of(200,203)..300"))
    ['one-of(200,203)..300']
    >>> list(_split_compound_loc("complement(123..145),200..209"))
    ['complement(123..145)', '200..209']
    >>> list(_split_compound_loc("123..145,one-of(200,203)..209"))
    ['123..145', 'one-of(200,203)..209']
    >>> list(_split_compound_loc("123..145,one-of(200,203)..one-of(209,211),300"))
    ['123..145', 'one-of(200,203)..one-of(209,211)', '300']
    >>> list(_split_compound_loc("123..145,complement(one-of(200,203)..one-of(209,211)),300"))
    ['123..145', 'complement(one-of(200,203)..one-of(209,211))', '300']
    >>> list(_split_compound_loc("123..145,200..one-of(209,211),300"))
    ['123..145', '200..one-of(209,211)', '300']
    >>> list(_split_compound_loc("123..145,200..one-of(209,211)"))
    ['123..145', '200..one-of(209,211)']
    >>> list(_split_compound_loc("complement(149815..150200),complement(293787..295573),NC_016402.1:6618..6676,181647..181905"))
    ['complement(149815..150200)', 'complement(293787..295573)', 'NC_016402.1:6618..6676', '181647..181905']
    """
    if "one-of(" in compound_loc:
        # Hard case
        while "," in compound_loc:
            assert compound_loc[0] != ","
            assert compound_loc[0:2] != ".."
            i = compound_loc.find(",")
            part = compound_loc[:i]
            compound_loc = compound_loc[i:]  # includes the comma
            while part.count("(") > part.count(")"):
                assert "one-of(" in part, (part, compound_loc)
                i = compound_loc.find(")")
                part += compound_loc[: i + 1]
                compound_loc = compound_loc[i + 1 :]
            if compound_loc.startswith(".."):
                i = compound_loc.find(",")
                if i == -1:
                    part += compound_loc
                    compound_loc = ""
                else:
                    part += compound_loc[:i]
                    compound_loc = compound_loc[i:]  # includes the comma
            while part.count("(") > part.count(")"):
                assert part.count("one-of(") == 2
                i = compound_loc.find(")")
                part += compound_loc[: i + 1]
                compound_loc = compound_loc[i + 1 :]
            if compound_loc.startswith(","):
                compound_loc = compound_loc[1:]
            assert part
            yield part
        if compound_loc:
            yield compound_loc
    else:
        # Easy case
        yield from compound_loc.split(",")


class ParserFailureError(Exception):
    """Failure caused by some kind of problem in the parser."""

    pass


class LocationParserError(Exception):
    """Could not Properly parse out a location from a GenBank file."""

    pass


def _loc(loc_str, expected_seq_length, strand, seq_type=None):
    """Make FeatureLocation from non-compound non-complement location (PRIVATE).

    This is also invoked to 'automatically' fix ambiguous formatting of features
    that span the origin of a circular sequence.

    Simple examples,

    >>> _loc("123..456", 1000, +1)
    FeatureLocation(ExactPosition(122), ExactPosition(456), strand=1)
    >>> _loc("<123..>456", 1000, strand = -1)
    FeatureLocation(BeforePosition(122), AfterPosition(456), strand=-1)

    A more complex location using within positions,

    >>> _loc("(9.10)..(20.25)", 1000, 1)
    FeatureLocation(WithinPosition(8, left=8, right=9), WithinPosition(25, left=20, right=25), strand=1)

    Notice how that will act as though it has overall start 8 and end 25.

    Zero length between feature,

    >>> _loc("123^124", 1000, 0)
    FeatureLocation(ExactPosition(123), ExactPosition(123), strand=0)

    The expected sequence length is needed for a special case, a between
    position at the start/end of a circular genome:

    >>> _loc("1000^1", 1000, 1)
    FeatureLocation(ExactPosition(1000), ExactPosition(1000), strand=1)

    Apart from this special case, between positions P^Q must have P+1==Q,

    >>> _loc("123^456", 1000, 1)
    Traceback (most recent call last):
       ...
    ValueError: Invalid between location '123^456'

    You can optionally provide a reference name:

    >>> _loc("AL391218.9:105173..108462", 2000000, 1)
    FeatureLocation(ExactPosition(105172), ExactPosition(108462), strand=1, ref='AL391218.9')

    >>> _loc("<2644..159", 2868, 1, "circular")
    CompoundLocation([FeatureLocation(BeforePosition(2643), ExactPosition(2868), strand=1), FeatureLocation(ExactPosition(0), ExactPosition(159), strand=1)], 'join')
    """
    if ":" in loc_str:
        ref, loc_str = loc_str.split(":")
    else:
        ref = None
    try:
        s, e = loc_str.split("..")
    except ValueError:
        assert ".." not in loc_str
        if "^" in loc_str:
            # A between location like "67^68" (one based counting) is a
            # special case (note it has zero length). In python slice
            # notation this is 67:67, a zero length slice.  See Bug 2622
            # Further more, on a circular genome of length N you can have
            # a location N^1 meaning the junction at the origin. See Bug 3098.
            # NOTE - We can imagine between locations like "2^4", but this
            # is just "3".  Similarly, "2^5" is just "3..4"
            s, e = loc_str.split("^")
            if int(s) + 1 == int(e):
                pos = _pos(s)
            elif int(s) == expected_seq_length and e == "1":
                pos = _pos(s)
            else:
                raise ValueError(
                    "Invalid between location %s" % repr(loc_str)
                ) from None
            return SeqFeature.FeatureLocation(pos, pos, strand, ref=ref)
        else:
            # e.g. "123"
            s = loc_str
            e = loc_str

    # Attempt to fix features that span the origin
    s_pos = _pos(s, -1)
    e_pos = _pos(e)
    if int(s_pos) > int(e_pos):
        if seq_type is None or "circular" not in seq_type.lower():
            warnings.warn(
                "It appears that %r is a feature that spans "
                "the origin, but the sequence topology is "
                "undefined. Skipping feature." % loc_str,
                BiopythonParserWarning,
            )
            return None
        warnings.warn(
            "Attempting to fix invalid location %r as "
            "it looks like incorrect origin wrapping. "
            "Please fix input file, this could have "
            "unintended behavior." % loc_str,
            BiopythonParserWarning,
        )

        f1 = SeqFeature.FeatureLocation(s_pos, expected_seq_length, strand)
        f2 = SeqFeature.FeatureLocation(0, int(e_pos), strand)

        if strand == -1:
            # For complementary features spanning the origin
            return f2 + f1
        else:
            return f1 + f2

    return SeqFeature.FeatureLocation(_pos(s, -1), _pos(e), strand, ref=ref)


def _pos(pos_str, offset=0):
    """Build a Position object (PRIVATE).

    For an end position, leave offset as zero (default):

    >>> _pos("5")
    ExactPosition(5)

    For a start position, set offset to minus one (for Python counting):

    >>> _pos("5", -1)
    ExactPosition(4)

    This also covers fuzzy positions:

    >>> p = _pos("<5")
    >>> p
    BeforePosition(5)
    >>> print(p)
    <5
    >>> int(p)
    5

    >>> _pos(">5")
    AfterPosition(5)

    By default assumes an end position, so note the integer behaviour:

    >>> p = _pos("one-of(5,8,11)")
    >>> p
    OneOfPosition(11, choices=[ExactPosition(5), ExactPosition(8), ExactPosition(11)])
    >>> print(p)
    one-of(5,8,11)
    >>> int(p)
    11

    >>> _pos("(8.10)")
    WithinPosition(10, left=8, right=10)

    Fuzzy start positions:

    >>> p = _pos("<5", -1)
    >>> p
    BeforePosition(4)
    >>> print(p)
    <4
    >>> int(p)
    4

    Notice how the integer behaviour changes too!

    >>> p = _pos("one-of(5,8,11)", -1)
    >>> p
    OneOfPosition(4, choices=[ExactPosition(4), ExactPosition(7), ExactPosition(10)])
    >>> print(p)
    one-of(4,7,10)
    >>> int(p)
    4

    """
    if pos_str.startswith("<"):
        return SeqFeature.BeforePosition(int(pos_str[1:]) + offset)
    elif pos_str.startswith(">"):
        return SeqFeature.AfterPosition(int(pos_str[1:]) + offset)
    elif _re_within_position.match(pos_str):
        s, e = pos_str[1:-1].split(".")
        s = int(s) + offset
        e = int(e) + offset
        if offset == -1:
            default = s
        else:
            default = e
        return SeqFeature.WithinPosition(default, left=s, right=e)
    elif _re_oneof_position.match(pos_str):
        assert pos_str.startswith("one-of(")
        assert pos_str[-1] == ")"
        parts = [
            SeqFeature.ExactPosition(int(pos) + offset)
            for pos in pos_str[7:-1].split(",")
        ]
        if offset == -1:
            default = min(int(pos) for pos in parts)
        else:
            default = max(int(pos) for pos in parts)
        return SeqFeature.OneOfPosition(default, choices=parts)
    else:
        return SeqFeature.ExactPosition(int(pos_str) + offset)


_re_within_position = re.compile(_within_position)
_re_oneof_position = re.compile(_oneof_position)
