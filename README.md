# GPRuler

## Execution of the pipeline
When GPRuler is executed, the user is given the option to choose between two alternative inputs:
* the name of the target organism
* an existing metabolic model

Both kind of inputs are firstly processed to obtain the list of metabolic genes associated with each metabolic reaction in the target organism/model. This intermediate output is then used as input for the core pipeline, which returns as ultimate output the GPR rule of each metabolic reactions.
The sequential execution steps of the proposed pipeline are detailed in the following.  

When an existing metabolic model is chosen:  
**Step 1. *metaboliteIdentification.py***: identification of the metabolites involved in the model reactions.
 * Inputs:
   * modelXml: the SBML model, which needs to be saved into the rawData directory
   * dfmetsInfo: a string to name the output file name
 * Outputs saved in the outputs directory:
   * dfmetsInfo + \'\_\' + timeStamp + \'.csv\': a file including all the information stored in the model for each metabolite. In particular: Id column includes the identifier associated to each metabolite, Name column includes the name associated to each metabolite, KeggId column includes the KEGG identifier associated to each metabolite, ChebiId  column includes the ChEBI identifier associated to each metabolite, PubchemId column includes the PubChem identifier associated to each metabolite, boundaryCondition column includes a boolean value indicating if the metabolite is in the boundary compartment, chemicalFormula column includes the chemical formula associated to each metabolite, Inchi column includes the InChi associated to each metabolite.
   * dfmetsInfo + \'\_wInferredIds\_\' + timeStamp + \'.csv\': a file storing in each row the name of the metabolite (Name column) and a list of the inferred identifiers(Identifiers column).

**Step 2. *metabolitesIdentification_FuzzyWuzzy.py***: identification of the metabolites involved in the model reactions through the FuzzyWuzzy package. 
* Inputs:
   * modelXml: the SBML model, which needs to be saved into the rawData directory
   * outputFileName: a string to name the output files name
 * Outputs saved in the outputs directory:
   * outputFileName + \'\_mappingMetaCyc\_allResults\_\' + timeStamp + \'.tsv\': the output of the FuzzyWuzzy package execution on MetaCyc database reporting for each metabolite the first ten tuple returned from the exploited package where the first element is the proposed found match and the second one is the associated score.
   * outputFileName + \'\_mappingKeggC\_allResults\_\' + timeStamp + \'.tsv\': the output of the FuzzyWuzzy package execution on KEGG COMPOUND database reporting for each metabolite the first ten tuple returned from the exploited package where the first element is the proposed found match and the second one is the associated score.
   * outputFileName + \'\_mappingKeggG\_allResults\_\' + timeStamp + \'.tsv\': the output of the FuzzyWuzzy package execution on KEGG GLYCAN database reporting for each metabolite the first ten tuple returned from the exploited package where the first element is the proposed found match and the second one is the associated score.
   * outputFileName + \'\_mappingChebi\_allResults\_\' + timeStamp + \'.tsv\': the output of the FuzzyWuzzy package execution on ChEBI database reporting for each metabolite the first ten tuple returned from the exploited package where the first element is the proposed found match and the second one is the associated score.  

**Step 3. *metabolitesIdentification_FuzzyWuzzy_part2***: refinement of the outcomes of Step 2.
* Inputs:
   * the 4 output files from the metabolitesIdentification_FuzzyWuzzy.py script
   * outputFileName: a string to name the output file name
 * Outputs saved in the outputs directory:
   * outputFileName + '_mappingMetaCyc_100_' + timeStamp + '.tsv': the output of the FuzzyWuzzy package execution on MetaCyc database reporting for each metabolite name the list of retrieved matches with 100 score.
   * outputFileName + '_mappingMetaCyc_91_99_' + timeStamp + '.tsv': the output of the FuzzyWuzzy package execution on MetaCyc database reporting for each metabolite name the list of retrieved matches with score between 91 and 99 needing of a manual curation.
   * outputFileName + '_mappingMetaCyc_empty_' + timeStamp + '.tsv': the output of the FuzzyWuzzy package execution on MetaCyc database reporting for each metabolite name the list of retrieved matches with score below 91, thus associated to an empty list.
   * outputFileName + '_mappingKeggC_100_' + timeStamp + '.tsv': the output of the FuzzyWuzzy package execution on KEGG COMPOUND database reporting for each metabolite name the list of retrieved matches with 100 score.
   * outputFileName + '_mappingKeggC_91_99_' + timeStamp + '.tsv': the output of the FuzzyWuzzy package execution on KEGG COMPOUND database reporting for each metabolite name the list of retrieved matches with score between 91 and 99 needing of a manual curation.
   * outputFileName + '_mappingKeggC_empty_' + timeStamp + '.tsv': the output of the FuzzyWuzzy package execution on KEGG COMPOUND database reporting for each metabolite name the list of retrieved matches with score below 91, thus associated to an empty list.
   * outputFileName + '_mappingKeggG_100_' + timeStamp + '.tsv': the output of the FuzzyWuzzy package execution on KEGG GLYCAN database reporting for each metabolite name the list of retrieved matches with 100 score.
   * outputFileName + '_mappingKeggG_91_99_' + timeStamp + '.tsv': the output of the FuzzyWuzzy package execution on KEGG GLYCAN database reporting for each metabolite name the list of retrieved matches with score between 91 and 99 needing of a manual curation.
   * outputFileName + '_mappingKeggG_empty_' + timeStamp + '.tsv': the output of the FuzzyWuzzy package execution on KEGG GLYCAN database reporting for each metabolite name the list of retrieved matches with score below 91, thus associated to an empty list.
   * outputFileName + '_mappingChebi_100_' + timeStamp + '.tsv': the output of the FuzzyWuzzy package execution on ChEBI database reporting for each metabolite name the list of retrieved matches with 100 score.
   * outputFileName + '_mappingChebi_91_99_' + timeStamp + '.tsv': the output of the FuzzyWuzzy package execution on ChEBI database reporting for each metabolite name the list of retrieved matches with score between 91 and 99 needing of a manual curation.
   * outputFileName + '_mappingChebi_empty_' + timeStamp + '.tsv': the output of the FuzzyWuzzy package execution on ChEBI database reporting for each metabolite name the list of retrieved matches with score below 91, thus associated to an empty list.
  
**Step 4. *metabolitesIdentification_joiningData***

## Input data
All input data are saved or required to be saved in a folder name `inputData`.  

When GPRuler is executed, the user is given the option to choose between the use of test models (option 1) and the use of your own data (option 2) - **Do you want to use test models (1) or your own data (2)?**:
* *If option 1 is chosen*, the user is also given the option to choose among HMRcore (option 1), Recon 3 (option 2), Yeast 7 (option 3) and Yeast8 (option 4) models - **Which test model you wanto to use? HMRcore (1), Recon 3 (2), Yeast 7 (3), Yeast 8 (4)**. According to the choice, the inputData directory includes for each model two files: the complete list of genes included in the model, which is a tab-separated file reporting the KEGG (`keggId` column) and the Uniprot (`uniprotId` column) identifiers of each gene (file name ends with `*_Kegg2UniprotGenes`); the complete list of reactions included in the model, which is a tab-separated file reporting for each reaction of the *Rxn* column the list of catalysing genes in the `Genes` column (file name ends with `*_Rxns2Genes`).
* *If option 2 is chosen*, the user is given the option to use own data starting from the name of the organism of interest (option 1) or from the list of reactions included in his/her own model and the associated metabolic genes (option 2) - **Which type of input data you have? Organism Name (1) or Reactions-Genes associations (2)**.
  * *If option 1 is chosen*, the user is firstly prompted for the name of the model he/she wants to generate - **Which is the model name?**. Then, the user is prompted for the target organism - **Do you have the organism name (1) or the KEGG code (2) of the organism under investigation?** - starting from the name (option 1) or the KEGG code (option 2). In option 1 - **Insert the organism name:** - the name is inserted and putative KEGG codes to choose from are returned - **Type the correct KEGG code among the returned ones:**. *In option 2* - **Insert the KEGG organism code:** - the KEGG code is known and inserted from the user. Starting from this information, KEGG database is queried to search from the list of metabolic genes annotated for the target organism and the corresponding list of catalysed reactions. In this regard, three additional tab-separated files are saved in the inputData folder:
    * file name ending with `*_GeneId2Rxns`: for each metabolic gene in the `GeneId` column the corresponding list of catalysed reaction is reported in `Rxns` column;
    * file name ending with `*_RxnId2Equation`: for each retrieved reaction in the `RxnId` column the corresponding equation is reported where the involved compounds are expressed through their KEGG identifier (`Equation` column) or their complete name (`Definition` column);
    * file name ending with `*_RxnId2ECs`: for each retrieved reaction in the `RxnId` column the corresponding list of EC number is reported in the `EC number` column;
    * file name ending with `*_Rxns2Genes`: for each retrieved reaction in the `Rxn` column the list of corresponding metabolic catalysing genes is reported in `Genes` column;
    * file name ending with `*_Kegg2UniprotGenes`: for each KEGG identifier in the `keggId` column of the retrieved list of metabolic genes the corresponding UNIPROT identifier is reported in `uniprotId` column.
  Finally, the user is given the option to choose between the possibility to insert the right regular expression used in order to identify the strings underlying the list of KEGG genes identifiers (option 1), or to infer the correct regular expression from the retrieved previous data (option 2) - **Do you want to manually insert your regex (1) or infer it from your data (2)?**.
  * *If option 2 is chosen*, the user is firstly prompted for the name of the model he/she wants to generate - **Which is the model name?**. Secondly, the user is prompted for the file name of the reactions to genes associations: the required input is a tab-separated file with two columns named `Rxn` and `Genes`, which respectively report the list of target reactions and the corresponding catalysing genes. Then, the user is prompted for the target organism - **Do you have the organism name (1) or the KEGG code (2) of the organism under investigation?** - starting from the name (option 1) or the KEGG code (option 2): see above for the details. Finally, the user is given the option to choose between the possibility to insert the right regular expression used in order to identify the strings underlying the list of KEGG genes identifiers (option 1), or to infer the correct regular expression from the retrieved previous data (option 2) - **Do you want to manually insert your regex (1) or infer it from your data (2)?**.

## Output data
When all the pipeline is executed, three output files are generated:
* a tab-separated file named as follow: `model + '_GenesData_' + timeStamp`, where model is the model name inserted as input data and timeStamp reports the date and time of the pipeline execution. This file includes 27 columns where for each gene all the retrieved relative information are returned.
* a tab-separated file named as follow: `model + '_GenesRelationships_' + timeStamp`, where model is the model name inserted as input data and timeStamp reports the date and time of the pipeline execution. This file includes 4 columns where for each gene all the genes having an AND or an OR relationship are reported.
* a tab-separated file named as follow: `model + '_gprRules_' + timeStamp`, where model is the model name inserted as input data and timeStamp reports the date and time of the pipeline execution. This file, which includes 3 columns, is an enrichment of the input file reporting for each input reaction the list of catalysing gene due to the inclusion in the `GPR rule` column of the corresponding reconstructed GPR rule.
