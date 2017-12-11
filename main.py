import pathlib
from unique_tryptic import Peptide, Protein, ProteinGroup

DATA_PATH = pathlib.Path('data')

def main():
    # nras IPI00315334.5    Nras neuroblastoma ras oncogene ENSMUSG00000027852
    # hras OPI00403929.3    Hras1 Hras1 protein (Fragment) ENSMUSG00000025499
    # kras IPI00113248.1    Kras Isoform 2A of GTPase KRas precursor ENSMUSG00000030265
    # proteins = {'nras': 'ENSMUSG00000027852', 'hras': 'ENSMUSG00000025499', 'kras': 'ENSMUSG00000030265'}
    proteins = {'nras': 'P08556', 'hras': ' Q61411', 'kras': 'P32883'}

    group = ProteinGroup()

    for name, _id in proteins.items():
        raw = DATA_PATH.joinpath(name).read_text().splitlines()
        peptides = (Peptide(*line.split()) for line in raw)
        protein = Protein(name=name, uniprot_id=_id, experimental_peptides=list(peptides))
        group.add(protein)

    for name, items in group.get_unique(ensure_full_unique=True).items():
        print(name)
        print('\n'.join(sorted(items)), '\n')

    return group

if __name__ == '__main__':
    main()
