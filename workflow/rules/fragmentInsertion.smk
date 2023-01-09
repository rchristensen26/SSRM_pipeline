rule estimateModelParams:
    input:
        refMSA=config["refMSA"],
        refTree="config/dsrAB_consensus_phylogeny.newick"
    params:
        outputDir=config["modelParamsDir"],
        name=config["paramName"]
    shell:
        """
        raxmlHPC -f e -m PROTGAMMADAYHOFF -s {input.refMSA} -t {input.refTree} -n {params.name} -w {params.outputDir}
        """


rule runRAXML:
    input:
        seqAlignment=join(config["cleanHitsDir"],"compiled_dsrAB_scoreThreshold_noDups_msa_withRef_trimmedGaps.faa"),
        modelParams=join(config["binaryModelParams"], config["paramName"]), # will be output from rule: estimateModelParams
        refTree=(config["refTree"], config["paramName"]) # will be output from rule: estimateModelParams
    params:
        outputDir=config["raxmlOutputDir"],
        fileExtension=config["raxmlOutputFileExtension"]
    output:
        join(config["raxmlOutputDir"], "RAxML_info.raxml"),
        join(config["raxmlOutputDir"], "RAxML_classification.raxml"),
        join(config["raxmlOutputDir"], "RAxML_classificationLikelihoodWeights.raxml"),
        join(config["raxmlOutputDir"], "RAxML_entropy.raxml"),
        join(config["raxmlOutputDir"], "RAxML_labelledTree.raxml"),
        join(config["raxmlOutputDir"], "RAxML_originalLabelledTree.raxml"),
        join(config["raxmlOutputDir"], "RAxML_portableTree.raxml.jplace")
    shell:
        """
        raxmlHPC -f v -R {input.modelParams} -r {input.refTree} -s {input.seqAlignment} -m PROTGAMMADAYHOFF -G 0.1 -n {params.fileExtension} -w {params.outputDir}
        """
