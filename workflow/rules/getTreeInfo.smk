rule getAnantharamanQueriesIDs:
    input: join(config["cleanHitsDir"], "Anantharaman2018_dsrA_dsrB_noDups.faa")
    output: "workflow/out/treeInfo/Anantharaman2018_queryIDs_inTree.txt"
    shell:
        "python3 Ananthtaraman2018_querySeqsInTree.py {input} {output}"


rule getBranchDistances:
    input:
        tree=join(config["raxmlOutputDir"],"RAxML_originalLabelledTree_noBootstrap.newick"),
        query_info=join(config["cleanHitsDir"],"compiled_dsrAB_scoreThreshold_noDups_gapPercentageInfo.csv"),
        Anantharaman2018_seq_list="workflow/out/treeInfo/Anantharaman2018_queryIDs_inTree.txt"
    output:
        csv="workflow/out/treeInfo/queryDistanceInfo.csv",
        json="workflow/out/treeInfo/queryDistanceInfo.json"
    shell:
        """
        python3 workflow/scripts/getBranchDistances.py {input.tree} {input.query_info} {input.Anantharaman2018_seq_list} {output.csv} {output.json}
        """

