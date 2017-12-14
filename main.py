import pathlib
import requests
import requests_cache
from io import StringIO
from unique_tryptic import Peptide, Protein, ProteinGroup

import sys
sys.path.insert(0, '../cimagex')
from cimagex import Dataset
from cimagex import make_mudpit_dataset as make_dataset


requests_cache.install_cache()


def main():
    # note, not working properly when instantiating more than one ProteinGroup...
    print('human')
    process(
        pathlib.Path('data/human'),
        {'nras': 'P01111', 'hras': 'P01112', 'kras': 'P01116'}
    )
    # nras IPI00315334.5    Nras neuroblastoma ras oncogene ENSMUSG00000027852
    # hras OPI00403929.3    Hras1 Hras1 protein (Fragment) ENSMUSG00000025499
    # kras IPI00113248.1    Kras Isoform 2A of GTPase KRas precursor ENSMUSG00000030265

    # print('mouse')
    # process(
    #     pathlib.Path('data/mouse'),
    #     # {'nras': 'P08556', 'hras': ' Q61411', 'kras': 'P32883'}
    #     {'nras': 'ENSMUSG00000027852', 'hras': 'ENSMUSG00000025499', 'kras': 'ENSMUSG00000030265'}
    # )


def process(data_path, proteins):
    group = ProteinGroup()

    for name, _id in proteins.items():
        raw = data_path.joinpath(name).read_text().splitlines()
        peptides = (Peptide(*line.split()) for line in raw)
        protein = Protein(name=name, uniprot_id=_id, experimental_peptides=list(peptides))
        group.add(protein)

    # for name, items in group.get_unique(ensure_full_unique=True).items():
    for name, items in group.get_unique_predicted().items():
        print(name)
        print('\n'.join(sorted(items)), '\n')

    return group

def get_peptides(proteins):
    BASE_URL = 'http://bfclabcomp4.scripps.edu/~remsberg/data/{}/combined_dta.txt'

    datasets = [
        '20170915_hydroxylamine_rep1',
        '20171005_hydroxylamine_rep2',
        '20171110_hydroxylamine_rep3',
        '20171110_hydroxylamine_rep4',
        '20171111_timecourse_1hr_rep1',
        '20171111_timecourse_1hr_rep2',
        '20171111_timecourse_30min_rep1',
        '20171111_timecourse_30min_rep2'
    ]

    master = Dataset()

    for d in datasets:
        r = requests.get(BASE_URL.format(d))
        dataset = make_dataset(StringIO(r.text), parse_as_file=False)
        master += dataset

    for name, uniprot_id in proteins:
        print(name)
        print('\n'.join('\t'.join((p.sequence, str(p.mass), str(p.ratio))) for p in master.get(uniprot_id).peptides))

    return [master.get(_id) for _id in uniprot_ids]


if __name__ == '__main__':
    main()
