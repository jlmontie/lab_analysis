import csv
from ncbi_taxonomy_utils import ncbi_taxonomy
import numpy as np


def bin_reads(sample_composition_path, ncbi_class=None,
              quantification='relative', ctrl_taxids=None):
    with open(sample_composition_path) as file:
        reader = csv.DictReader(file, delimiter='\t',
                                fieldnames=['taxid', 'count'])
        taxid_counts = []
        for entry in reader:
            taxid_counts.append(entry)

    if ncbi_class is None:
        ncbi = ncbi_taxonomy()
    else:
        ncbi = ncbi_class

    human_taxid = 9606
    bacteria_taxid = 2
    virus_taxid = 10239
    eukaryota_taxid = 2759
    parasite_taxids = {
        7563,
        188941,
        6029,
        5653,
        6935,
        6178,
        5794,
        6308,
        31277,
        119088,
        6199,
        85819,
        33083,
        33084,
        75966,
        41165,
        7509,
        6236,
        198624,
        33634,
        5988,
        6249,
        5738,
        1489900,
        740972,
        1485168,
        37104,
        10232
    }
    fungus_taxid = 4751
    human_counts = []
    bacteria_counts = []
    virus_counts = []
    parasite_counts = []
    fungus_counts = []
    unclassified_counts = []

    for taxid in taxid_counts:
        if taxid in ctrl_taxids:
            continue
        taxid_path = ncbi.get_path(int(taxid['taxid']))
        count_value = taxid['count']
        if isinstance(count_value, int):
            count = count_value
        elif isinstance(count_value, str):
            if len(count_value) > 0:
                count = np.int64(taxid['count'].split('.')[0])
            else:
                continue
        if human_taxid in taxid_path:
            human_counts.append(count)
        elif fungus_taxid in taxid_path:
            fungus_counts.append(count)
        elif parasite_taxids.intersection(set(taxid_path)) != set():
            parasite_counts.append(count)
        elif bacteria_taxid in taxid_path:
            bacteria_counts.append(count)
        elif virus_taxid in taxid_path:
            virus_counts.append(count)
        else:
            unclassified_counts.append(count)

    human_sum = np.array(human_counts).sum()
    bacteria_sum = np.array(bacteria_counts).sum()
    virus_sum = np.array(virus_counts).sum()
    parasite_sum = np.array(parasite_counts).sum()
    fungus_sum = np.array(fungus_counts).sum()
    unclassified_sum = np.array(unclassified_counts).sum()
    cumalitve_sum = human_sum + bacteria_sum + virus_sum + parasite_sum + \
        fungus_sum + unclassified_sum

    if quantification == 'relative':
        count_dict = {
            "Human": 100 * human_sum / cumalitve_sum,
            "Bacteria": 100 * bacteria_sum / cumalitve_sum,
            "Virus": 100 * virus_sum / cumalitve_sum,
            "Parasite": 100 * parasite_sum / cumalitve_sum,
            "Fungus": 100 * fungus_sum / cumalitve_sum,
            "Unclassified": 100 * unclassified_sum / cumalitve_sum
        }
    elif quantification == 'absolute':
        count_dict = {
            "Human": human_sum,
            "Bacteria": bacteria_sum,
            "Virus": virus_sum,
            "Parasite": parasite_sum,
            "Fungus": fungus_sum,
            "Unclassified": unclassified_sum
        }
    elif quantification == 'both':
        count_dict = {
            "Human": {
                "relative": 100 * human_sum / cumalitve_sum,
                "absolute": human_sum
            },
            "Bacteria": {
                "relative": 100 * bacteria_sum / cumalitve_sum,
                "absolute": bacteria_sum
            },
            "Virus": {
                "relative": 100 * virus_sum / cumalitve_sum,
                "absolute": virus_sum
            },
            "Parasite": {
                "relative": 100 * parasite_sum / cumalitve_sum,
                "absolute": parasite_sum
            },
            "Fungus": {
                "relative": 100 * fungus_sum / cumalitve_sum,
                "absolute": fungus_sum
            },
            "Unclassified": {
                "relative": 100 * unclassified_sum / cumalitve_sum,
                "absolute": unclassified_sum
            },
        }
    return count_dict


if __name__ == "__main__":
    import argparse
    import matplotlib.pyplot as plt
    import seaborn as sns
    import json


    def circle_plot(count_dict):
        names = [f"{name} ({count:0.2f})%" for name, count in
                 count_dict.items()]
        counts = count_dict.values()
        sorted_zip = sorted(zip(counts, names), reverse=True)
        counts_sorted = [a for a, b in sorted_zip]
        names_sorted = [b for a, b in sorted_zip]
        _, ax = plt.subplots()
        patches, _ = ax.pie(counts_sorted)
        circle = plt.Circle((0, 0), 0.7, color='white')
        ax.add_artist(circle)
        ax.set_title('Sample Composition')
        plt.legend(patches, names_sorted, loc='center',
                   bbox_to_anchor=(1.2, 0.5))
        plt.show()


    parser = argparse.ArgumentParser(description='Bin read kmers into human, '
                                     'bacterial, viral, parasitic, fungal, '
                                     'and unclassified categories.')
    parser.add_argument('composition_path',
                        help='Path to the sample composition file.')
    parser.add_argument('-o', '--output', nargs=1,
                        help='Path to save the JSON of binned reads.')
    parser.add_argument('-pc', '--plot_circle', action='store_true',
                        help='Produce a circle plot of sample composition.')
    args = parser.parse_args()

    binned_reads = bin_reads(args.composition_path)
    if args.output is not None:
        with open(args.output[0], 'w') as file:
            json.dump(binned_reads, file)
    else:
        print(binned_reads)
    if args.plot_circle:
        circle_plot(binned_reads)
