import requests
import requests_cache
from pyteomics.parser import cleave, expasy_rules
from Bio import SeqIO
from io import StringIO

requests_cache.install_cache()


def get_sequence_from_uniprot(uniprot_id):
    params = {'query': uniprot_id, 'format': 'fasta'}
    r = requests.get('http://www.uniprot.org/uniprot/', params)
    # grabbing first element and returning as string
    return str(next(SeqIO.parse(StringIO(r.text), 'fasta')).seq)


def predict_tryptic_peptides(sequence):
    return cleave(sequence, expasy_rules['trypsin'])


class Peptide:
    delchars = {ord(c): None for c in map(chr, range(256)) if not c.isalpha()}

    def __init__(self, sequence, mass=None, ratio=None, run=None):
        self.sequence = sequence
        self.clean_sequence = sequence.split('.')[1].translate(Peptide.delchars)
        self.mass = mass
        self.ratio = ratio
        self.run = run

    def __repr__(self):
        return 'Peptide(sequence={}, clean_sequence={}, mass={}, ratio={}, run={})'.format(
            self.sequence,
            self.clean_sequence,
            self.mass,
            self.ratio,
            self.run
        )


class Protein:
    def __init__(self, name, uniprot_id, experimental_peptides):
        self.name = name
        self.uniprot_id = uniprot_id
        self.experimental_peptides = experimental_peptides
        self.full_sequence = get_sequence_from_uniprot(uniprot_id)
        self.predicted_tryptic_sequences = predict_tryptic_peptides(self.full_sequence)
        self.unique_experimental_sequences = self.get_unique_experimental_sequences()

    def get_unique_experimental_sequences(self):
        return set(p.clean_sequence for p in self.experimental_peptides)

    def __repr__(self):
        return 'Protein(name={}, uniprot_id={})'.format(
            self.name,
            self.uniprot_id
        )


class ProteinGroup(object):
    def __init__(self, proteins=None):
        if proteins:
            self.proteins = proteins
        else:
            self.proteins = []

    def add(self, protein):
        self.proteins.append(protein)

    def get_unique(self, ensure_full_unique=True):
        unique = {}

        for index, protein in enumerate(self.proteins):
            other_proteins = self.proteins[:index] + self.proteins[(index + 1):]
            other_sequences = (p.unique_experimental_sequences for p in other_proteins)

            unique_experimental = protein.unique_experimental_sequences.difference(*other_sequences)

            # filter unique experimental sequences by those which are not found in
            # the full sequences of other comparison group members
            if ensure_full_unique:
                other_full_sequences = (p.full_sequence for p in other_proteins)
                unique_experimental = (x for x in unique_experimental if not any(x in s for s in other_full_sequences))

            unique[protein.name] = list(unique_experimental)

        return unique

    def get_unique_predicted(self):
        unique = {}
        for index, protein in enumerate(self.proteins):
            other_proteins = self.proteins[:index] + self.proteins[(index + 1):]
            other_sequences = (p.predicted_tryptic_sequences for p in other_proteins)
            other_full_sequences = (p.full_sequence for p in other_proteins)

            unique_predicted = protein.predicted_tryptic_sequences.difference(*other_sequences)
            unique_predicted = (x for x in unique_predicted if not any(x in s for s in other_full_sequences))
            unique[protein.name] = list(unique_predicted)

        return unique

    def __repr__(self):
        return 'ProteinGroup(proteins={})'.format(
            self.proteins
        )
