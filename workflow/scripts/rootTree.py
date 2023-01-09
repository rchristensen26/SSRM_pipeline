""""
This code roots a phylogenetic tree between archaeal and bacterial dsrAB clusters.
"""

from ete3 import Tree
import sys
import os

UNROOTED_TREE_F = sys.argv[1]
rooted_tree_fpath = os.path.splitext(os.path.basename(UNROOTED_TREE_F))[0] + "_rooted.newick"

with open(UNROOTED_TREE_F, 'r') as f:
    t = f.read()
    tree = Tree(t, format=1)

rootset = "Pyrobaculum_aerophilum_str._IM2_copy_1__PrbAero3", \
          "Pyrobaculum_aerophilum_str._IM2_copy_2__PrbAero4",\
          "Pyrobaculum_islandicum__PrbIslan",\
          "Pyrobaculum_arsenaticum_copy_1__PrbArse3",\
          "Pyrobaculum_arsenaticum_copy_2__PrbArse4",\
          "Pyrobaculum_arsenaticum_copy_3__PrbArse5",\
          "Pyrobaculum_calidifontis_copy_1__PrbCali3",\
          "Pyrobaculum_calidifontis_copy_2__PrbCali4",\
          "Caldivirga_maquilingensis__CvgMaqu3",\
          "Pyrobaculum_oguniense_copy_3__PyrOgun5",\
          "Pyrobaculum_oguniense_copy_1__PyrOgun6",\
          "Pyrobaculum_oguniense_copy_2__PyrOgun7",\
          "Pyrobaculum_neutrophilum__TptNeut2",\
          "Pyrobaculum_sp._1860__CP003098",\
          "Microbial_hot_spring_community__entry547",\
          "Vulcanisaeta_distributa_IC-017T__CP002100",\
          "Vulcanisaeta_moutnovskia_768-28__CP002529"

common_ancestor = tree.get_common_ancestor(rootset)
# print(common_ancestor)

if common_ancestor is tree:  # check if common_ancestor is the current root node
    print("node is root!")
    tree.unroot()
    common_ancestor = tree.get_common_ancestor(rootset)
    tree.set_outgroup(common_ancestor)
else:
    tree.set_outgroup(common_ancestor)

tree.write(format=1, outfile=rooted_tree_fpath)
