import os


repertoire = os.getcwd()


fichier = open("config.yaml", "w")

fichier.write("projectTitle: Detection of RNA Outlier Pipeline\n")
fichier.write("root: " + repertoire + "/output\n")
fichier.write("htmlOutputPath: " + repertoire + "/output/html\n")
fichier.write("indexWithFolderName: true\n")

fichier.write("hpoFile: /media/jbogoin/Data1/Annotations/HPO/genes_to_disease.txt\n")
fichier.write("sampleAnnotation: sample_annotation.tsv\n")

fichier.write("geneAnnotation:\n")
fichier.write("    v43: /media/jbogoin/Data1/References/RNA-seq/hg38/gencode.v43.primary_assembly.basic.annotation.gtf\n")
fichier.write("genomeAssembly: hg38\n")
fichier.write("genome: /media/jbogoin/Data1/References/fa_hg38/hg38_rnaseq/GRCh38.v43.primary_assembly.genome.fa\n")

###############
fichier.write("\n")
###############

fichier.write("exportCounts:\n")
fichier.write("    geneAnnotations:\n")
fichier.write("        - v43\n")
fichier.write("    excludeGroups:\n")
fichier.write("        - group1\n")

###############
fichier.write("\n")
###############

fichier.write("aberrantExpression:\n")
fichier.write("    run: true\n")
fichier.write("    groups:\n")
fichier.write("        - outrider\n")
# fichier.write("        - outrider_external\n")
fichier.write("    fpkmCutoff: 1\n")
fichier.write("    implementation: autoencoder\n")
fichier.write("    padjCutoff: 1\n")
fichier.write("    zScoreCutoff: 0\n")
fichier.write("    genesToTest: null\n")
fichier.write("    maxTestedDimensionProportion: 3\n")
fichier.write("    dassie:\n")
fichier.write("        tssWindow: 500\n")
fichier.write("        pasWindow: 1000\n")

###############
fichier.write("\n")
###############

fichier.write("aberrantSplicing:\n")
fichier.write("    run: true\n")
fichier.write("    groups:\n")
fichier.write("        - fraser\n")
# fichier.write("        - fraser_external\n")
fichier.write("    recount: true\n")
fichier.write("    longRead: false\n")
fichier.write("    keepNonStandardChrs: false\n")
fichier.write("    filter: false\n")
fichier.write("    minExpressionInOneSample: 10\n")
fichier.write("    minDeltaPsi: 0.05\n")
fichier.write("    implementation: PCA\n")
fichier.write("    padjCutoff: 1\n")
fichier.write("    maxTestedDimensionProportion: 6\n")
    
    ### FRASER1 configuration
fichier.write('    FRASER_version: "FRASER"\n')
fichier.write('    deltaPsiCutoff : 0.3\n') 
fichier.write('    quantileForFiltering: 0.95\n') 
    
    ### For FRASER2, use the follwing parameters instead of the 3 lines above:
    # FRASER_version: "FRASER2"
    # deltaPsiCutoff : 0.1
    # quantileForFiltering: 0.75

###############
fichier.write("\n")
###############

fichier.write("mae:\n")
fichier.write("    run: false\n")
fichier.write("    groups:\n")
fichier.write("        - group1\n")
fichier.write("    gatkIgnoreHeaderCheck: true\n")
fichier.write("    padjCutoff: 0.05\n")
fichier.write("    allelicRatioCutoff: 0.8\n")
fichier.write("    addAF: true\n")
fichier.write("    maxAF: 0.001\n")
fichier.write("    maxVarFreqCohort: 0.05\n")
    # VCF-BAM matching
fichier.write("    qcVcf: /media/jbogoin/Data1/References/RNA-seq/hg38/DROP/qc_vcf_1000G_hg38.vcf.gz\n")
fichier.write("    qcGroups:\n")
fichier.write("        - mae\n")
fichier.write("    dnaRnaMatchCutoff: 0.85\n")

###############
fichier.write("\n")
###############

fichier.write("rnaVariantCalling:\n")
fichier.write("    run: false\n")
fichier.write("    groups:\n")
fichier.write("        - batch_0\n")
fichier.write("    highQualityVCFs:\n")
fichier.write("        - /media/jbogoin/Data1/References/RNA-seq/hg38/DROP/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz\n")
fichier.write("        - /media/jbogoin/Data1/References/RNA-seq/hg38/DROP/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz\n")
fichier.write("    dbSNP: /media/jbogoin/Data1/References/RNA-seq/hg38/DROP/dbsnp132_20101103.vcf.gz\n")
fichier.write("    repeat_mask: /media/jbogoin/Data1/References/RNA-seq/hg38/DROP/hg38_repeatMasker_sorted.bed\n")
fichier.write("    createSingleVCF: true\n")
fichier.write("    addAF: true\n")
fichier.write("    maxAF: 0.001\n")
fichier.write("    maxVarFreqCohort: 0.05\n")
fichier.write('    hcArgs: ""\n')
fichier.write("    minAlt: 3\n")
fichier.write("    yieldSize: 100000\n")

###############
fichier.write("\n")
###############
fichier.write("tools:\n")
fichier.write("    gatkCmd: gatk\n")
fichier.write("    bcftoolsCmd: bcftools\n")
fichier.write("    samtoolsCmd: samtools\n")

fichier.close()


print('#####\nFichier config.yaml genere !\n')