"""
I don't really know what this does... I kinda just found it from my previous pipeline and using it I guess??
"""

from ete3 import Tree
import csv
import sys

RESULT_TREE = sys.argv[1]
REF_TREE = sys.argv[2]
DISTANCES_INFO_CSV = sys.argv[3]
PRUNED_RESULT_TREE = sys.argv[4]
PRUNED_REF_TREE = sys.argv[5]


def main():
    tree = read_in_tree(RESULT_TREE)
    r_tree = read_in_tree(REF_TREE)

    prune_tree(tree, DISTANCES_INFO_CSV, PRUNED_RESULT_TREE, "result")

    prune_tree(r_tree, DISTANCES_INFO_CSV, PRUNED_REF_TREE, "ref")


def read_in_tree(file):
    with open(file, mode='r') as f:
        t_file = f.read()
    return Tree(t_file, format=1)


def prune_tree(tree, query_info, output_file, t_type):
    prune_list = get_prune_list(query_info, t_type)

    pruned_tree = tree.prune(prune_list, preserve_branch_length=True)

    pruned_tree.write(format=1, outfile=output_file)


def get_prune_list(query_info, t_type):
    q_prune_list = []
    r_prune_list = []

    reader = csv.DictReader(open(query_info, mode="r"))

    special_char_list = ['Bilophila_sp._4_1_30__BilSpec3',
                         'Anaerobic_bacterium_sk.prop8__Bv2Prop2',
                         'Desulfovibrio_sp._3_1_syn3__DsvSp230',
                         'Desulfovibrio_sp._6_1_46AFAA__DsvSp231',
                         'Desulfovibrio_desulfuricans_subsp._desulfuricans_str._ATCC_27774__CP001358',
                         'Desulfovibrio_sp._Dsv1__DsvSp233',
                         'Desulfitobacterium_sp._PCE-1__DstSpe13']

    for row in reader:
        query = row["query_names"]
        q_prune_list.append(query)

        ref = row["closest_ref"].strip("[]\"\'")
        ref = ref.strip("\'")

        if ref not in r_prune_list:
            if ref not in special_char_list:
                r_prune_list.append(ref)

    if t_type == "result":
        return q_prune_list + r_prune_list

    elif t_type == "ref":
        return r_prune_list


if __name__ == '__main__':
    main()
