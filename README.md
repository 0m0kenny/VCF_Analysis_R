# VCF_Analysis_R
An R markdown tutorial for analysing a vcf file.

## Background
- The use of VCF files is important and plays a critical role in
describing variants.
* It stores many genetic information allowing identification of
potential mutations and variants present in the family.
+ Using these variants, we can analyse and reveal important
insights in inheritance patterns, potential-disease causing
mutation that could effect the individual.

## Aims:
- Read, Filter, re-structure and process VCF file using R studio.
- Summarise variants types and function
- Identify the parent and child sample
- Identify the De novo mutations.
- Annotate different variants (De novo mutation variants).
- To find potential diseases or a disorder that might affect the
child.

[!NOTE]
> The 'DellyVariation.vcf' which contains multiple samples but this tutorial focused
on NA19238, NA19239, NA19240 samples as they are from the same family. 
> After filtering, the VCF file is uploaded to Ensembl's Variant Effect Predictor and the resulting output is the the 'VEPannotated.txt.zip' file is the compressed version of the VEPannotated.txt file which is too large to upload on git. 
> Unzip the 'VEPannotated.txt.zip' before loading into R.

[!TIP]
> Explore different filtering parameters and different families.
> Explore [Ensembl's Variant Effect Predictor](https://www.ensembl.org/info/docs/tools/vep/index.html) after filtering the vcf using your own parameters. 
> Can also use different variant annotation tools such as [ANNOVAR](https://annovar.openbioinformatics.org/en/latest/).
> Visualise variants via [Integrative Genomics Viewer (IGV)](https://igv.org/app/)


### I also included the R SCript to run the full analysis.

