rule removeBootstrapValues: # for formatting reasons
    input:
        join(config["raxmlOutputDir"],"RAxML_originalLabelledTree.raxmlfinal")
    output:
        join(config["raxmlOutputDir"],"RAxML_originalLabelledTree_noBootstrap.newick")
    shell:
        """
        python3 workflow/scripts/removeBootstrapValues.py {input} {output}
        """

rule getBranchDistances:
    input:
        tree=join(config["raxmlOutputDir"],"RAxML_originalLabelledTree_noBootstrap.newick"),
        query_info=join(config["cleanHitsDir"],"compiled_dsrAB_scoreThreshold_noDups_gapPercentageInfo.csv")
    output:
        csv="workflow/out/treeInfo/queryDistanceInfo.csv",
        json="workflow/out/treeInfo/queryDistanceInfo.json"
    shell:
        """
        python3 workflow/scripts/getBranchDistances.py {input.tree} {input.query_info} {output.csv} {output.json}
        """

rule pruneTree:
    input:
        result_tree=join(config["raxmlOutputDir"],"RAxML_originalLabelledTree_noBootstrap.newick"),
        ref_tree=join(config["raxmlOutputDir"], "RAxML_originalLabelledTree_noBootstrap.newick"),
        distance_info="workflow/out/treeInfo/queryDistanceInfo.csv",
    output:
        pruned_result_tree="workflow/out/treeInfo/prunedResultTree.nw",
        pruned_ref_tree="workflow/out/treeInfo/prunedRefTree.nw"
    shell:
        """
        python3 workflow/scripts/pruneTree.py {input.result_tree} {input.ref_tree} {input.distance_info} {output.pruned_result_tree} {output.pruned_ref_tree}
        """

rule rootTree:
    input:
        join(config["raxmlOutputDir"],"RAxML_originalLabelledTree_noBootstrap.newick")
    output:
        join(config["raxmlOutputDir"],"RAxML_originalLabelledTree_noBootstrap_rooted.newick")
    shell:
        """
        python3 workflow/scripts/rootTree.py {input} 
        """
