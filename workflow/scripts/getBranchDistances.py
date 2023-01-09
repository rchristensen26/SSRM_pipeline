"""
This code returns the distance between every query sequence leaf and the closest reference leaf.

INPUT:
    INPUT_TREE: unrooted tree containing reference and query sequences (newick format)
    QUERY_INFO: csv file with query information, including query name by sample # and number of duplicates

OUTPUT:
    DISTANCES_INFO_CSV: csv file with shortest branch length between every query sequence and the closest
                    reference leaf
"""

from ete3 import Tree
import csv
import sys
import json

RESULT_TREE = sys.argv[1]
QUERY_INFO = sys.argv[2]
DISTANCES_INFO_CSV = sys.argv[3]
DISTANCES_INFO_JSON = sys.argv[4]


def main():
    tree = read_in_tree(RESULT_TREE)

    dist_dict = get_distances(tree)
    with open(DISTANCES_INFO_JSON, 'w') as f:
        json.dump(dist_dict, f)

    add_dist_info(dist_dict, QUERY_INFO, DISTANCES_INFO_CSV)


def read_in_tree(file):
    with open(file, mode='r') as f:
        t_file = f.read()
    return Tree(t_file, format=1)


def get_distances(tree):
    # make list of query names in tree
    query_list = []
    for node in tree:
        if node.is_leaf():
            if "QUERY" in node.name:
                query_list.append(node.name)

    # make list of ref names in tree
    ref_list = []
    for node in tree:
        if node.is_leaf():
            if "QUERY" not in node.name:
                ref_list.append(node.name)

    closest_rleaf_dict = {}

    for query in query_list:
        q_leaf = tree.search_nodes(name=query)[0]
        closest_rleaf_dict[query] = {}

        find_closest_rleaf(q_leaf, tree, closest_rleaf_dict, query, ref_list)

    return closest_rleaf_dict


def find_closest_rleaf(q_leaf, tree, closest_rleaf_dict, query, ref_list):
    for ref in ref_list:
        r_leaf = tree.search_nodes(name=ref)[0]
        dist = q_leaf.get_distance(r_leaf)

        if closest_rleaf_dict[query] == {}:
            closest_rleaf_dict[query]["closest_ref"] = [ref]
            closest_rleaf_dict[query]["dist"] = dist

        elif closest_rleaf_dict[query]["dist"] == dist:
            closest_rleaf_dict[query]["closest_ref"].append(ref)

        elif closest_rleaf_dict[query]["dist"] > dist:
            closest_rleaf_dict[query]["closest_ref"] = [ref]
            closest_rleaf_dict[query]["dist"] = dist


def add_dist_info(dist_dict, query_info, output_f):
    with open(query_info, mode='r') as f_read:
        with open(output_f, mode='w') as f_write:
            reader = csv.DictReader(f_read)
            o_fieldnames = reader.fieldnames
            o_fieldnames.extend(["closest_ref", "dist", "query_names"])
            writer = csv.DictWriter(f_write, fieldnames=o_fieldnames)

            writer.writeheader()

            for row in reader:
                query_name = "QUERY___" + row["sample_name"]

                if query_name in dist_dict:
                    row["closest_ref"] = dist_dict[query_name]["closest_ref"]

                    row["dist"] = dist_dict[query_name]["dist"]

                    row["query_names"] = query_name

                    writer.writerow(row)


if __name__ == '__main__':
    main()
