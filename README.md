# GPRuler

When GPRuler is executed, the user is given the option to choose between two alternative inputs:
* the name of the target organism
* an existing metabolic model

Both kind of inputs are firstly processed to obtain the list of metabolic genes associated with each metabolic reaction in the target organism/model. This intermediate output is then used as input for the core pipeline, which returns as ultimate output the GPR rule of each metabolic reactions.
The sequential execution steps of the proposed pipeline are detailed in the following.  

## Execution of the pipeline from an existing metabolic model
When an existing metabolic model is chosen:  
**Step 1. *metaboliteIdentification.py***: identification of the metabolites involved in the model reactions.
 * Inputs:
   * modelXml: the SBML model, which needs to be saved into the rawData directory
   * dfmetsInfo: a string to name the output files
   * includeCompartment: a boolean variable to indicate if metabolites name also includes the compartment information between square brackets.
 * Outputs saved in the outputs directory:
   * dfmetsInfo + \'.csv\': a file including all the information stored in the model for each metabolite. In particular: Id column includes the identifier associated to each metabolite, Name column includes the name associated to each metabolite, KeggId column includes the KEGG identifier associated to each metabolite, ChebiId  column includes the ChEBI identifier associated to each metabolite, PubchemId column includes the PubChem identifier associated to each metabolite, boundaryCondition column includes a boolean value indicating if the metabolite is in the boundary compartment, chemicalFormula column includes the chemical formula associated to each metabolite, Inchi column includes the InChi associated to each metabolite.
   * dfmetsInfo + \'\_wInferredIds\_.csv\': a file storing in each row the name of the metabolite (Name column) and a list of the inferred identifiers(Identifiers column).

**Step 2. *metabolitesIdentification_FuzzyWuzzy.py***: identification of the metabolites involved in the model reactions through the FuzzyWuzzy package.
* Inputs:
   * modelXml: the SBML model, which needs to be saved into the rawData directory
   * outputFileName: a string to name the output files
   * dfmetsInfo: the first output of Step 1
 * Outputs saved in the outputs directory:
   * outputFileName + \'\_mappingMetaCyc\_allResults.tsv\': the output of the FuzzyWuzzy package execution on MetaCyc database reporting for each metabolite the first ten tuple returned from the exploited package where the first element is the proposed found match and the second one is the associated score.
   * outputFileName + \'\_mappingKeggC\_allResults.tsv\': the output of the FuzzyWuzzy package execution on KEGG COMPOUND database reporting for each metabolite the first ten tuple returned from the exploited package where the first element is the proposed found match and the second one is the associated score.
   * outputFileName + \'\_mappingKeggG\_allResults.tsv\': the output of the FuzzyWuzzy package execution on KEGG GLYCAN database reporting for each metabolite the first ten tuple returned from the exploited package where the first element is the proposed found match and the second one is the associated score.
   * outputFileName + \'\_mappingChebi\_allResults.tsv\': the output of the FuzzyWuzzy package execution on ChEBI database reporting for each metabolite the first ten tuple returned from the exploited package where the first element is the proposed found match and the second one is the associated score.  

**Step 3. *metabolitesIdentification_FuzzyWuzzy_part2***: refinement of the outcomes of Step 2.
  * Inputs:
    * the 4 output files from the Step 2
    * outputFileName: a string to name the output files
  * Outputs saved in the outputs directory:
    * outputFileName + \'\_mappingMetaCyc\_100.tsv\': the output of the FuzzyWuzzy package execution on MetaCyc database reporting for each metabolite name the list of retrieved matches with 100 score.
    * outputFileName + \'\_mappingMetaCyc_91_99.tsv\': the output of the FuzzyWuzzy package execution on MetaCyc database reporting for each metabolite name the list of retrieved matches with score between 91 and 99 needing of a manual curation.
    * outputFileName + \'\_mappingMetaCyc_empty.tsv\': the output of the FuzzyWuzzy package execution on MetaCyc database reporting for each metabolite name the list of retrieved matches with score below 91, thus associated to an empty list.
    * outputFileName + \'\_mappingKeggC\_100.tsv\': the output of the FuzzyWuzzy package execution on KEGG COMPOUND database reporting for each metabolite name the list of retrieved matches with 100 score.
    * outputFileName + \'\_mappingKeggC\_91\_99.tsv\': the output of the FuzzyWuzzy package execution on KEGG COMPOUND database reporting for each metabolite name the list of retrieved matches with score between 91 and 99 needing of a manual curation.
    * outputFileName + \'\_mappingKeggC\_empty.tsv\': the output of the FuzzyWuzzy package execution on KEGG COMPOUND database reporting for each metabolite name the list of retrieved matches with score below 91, thus associated to an empty list.
    * outputFileName + \'\_mappingKeggG\_100.tsv\': the output of the FuzzyWuzzy package execution on KEGG GLYCAN database reporting for each metabolite name the list of retrieved matches with 100 score.
    * outputFileName + \'\_mappingKeggG\_91\_99.tsv\': the output of the FuzzyWuzzy package execution on KEGG GLYCAN database reporting for each metabolite name the list of retrieved matches with score between 91 and 99 needing of a manual curation.
    * outputFileName + \'\_mappingKeggG\_empty.tsv\': the output of the FuzzyWuzzy package execution on KEGG GLYCAN database reporting for each metabolite name the list of retrieved matches with score below 91, thus associated to an empty list.
    * outputFileName + \'\_mappingChebi\_100.tsv\': the output of the FuzzyWuzzy package execution on ChEBI database reporting for each metabolite name the list of retrieved matches with 100 score.
    * outputFileName + \'\_mappingChebi\_91\_99.tsv\': the output of the FuzzyWuzzy package execution on ChEBI database reporting for each metabolite name the list of retrieved matches with score between 91 and 99 needing of a manual curation.
    * outputFileName + \'\_mappingChebi\_empty.tsv\': the output of the FuzzyWuzzy package execution on ChEBI database reporting for each metabolite name the list of retrieved matches with score below 91, thus associated to an empty list.

**Step 4. *metabolitesIdentification_joiningData***: join the results obtained from both the Steps 1 and 2.
  * Input:
    * dfmetsInfo + \'\_wInferredIds\_\' + timeStamp + \'.csv\': the second output of Step 1
    * the outputs of Step 3
    * outputFileName: a string to name the output files
  * Output saved in the outputs directory:
    * outputFileName + \'\_mappingFuzzyAndClassic\_.tsv\': a file storing in each row the name of the metabolite (Name column), a list of the inferred identifiers from the Step 1 (Identifiers_classic), a list of the inferred identifiers from the Step 2 (Identifiers_fuzzy), the joining of the two columns Identifiers_classic and Identifiers_fuzzy (Identifiers column).

**Step 5. *reactionsIdentification***: identification of the reactions metabolites involved in the model.
  * Input:
    * modelXml: the SBML model, which needs to be saved into the rawData directory
    * dfmetsInfo: the first output of Step 1
    * dfmetsIds: the output of Step 4
    * dfrxnsInfo: a string to name the output files
  * Output saved in the outputs directory:
    * dfmetsInfo + \'\_enriched.csv: the first output of Step 1 enriched with the metabolites identifiers retrieved from GPRuler
    * dfrxnsInfo + \'.csv\': a file including all the information stored in the model for each reaction
    * dfrxnsInfo + \'\_wIds.csv\': a file storing in each row the name of the reaction (RxnId column) and a list of the inferred identifiers (PutativeIdentifiers column)
    * dfrxnsInfo + \'\_enriched.csv\': the second output file enriched with the inferred reactions identifiers

**Step 6. *reactionsIdentification_TCDB***: identification of the reactions metabolites involved in the model by querying the TCDB database.
* Input:
  * modelXml: the SBML model, which needs to be saved into the rawData directory
  * dfmetsInfo: the first output of Step 1
  * dfmetsIds: the output of Step 4
  * dfrxnsInfo: a string to name the output files
* Output saved in the outputs directory:
  * dfrxnsInfo + \'\_enriched\_tcdb.csv\': the second output of Step 5 is enriched with a dictionary for each reaction (Identifiers_fromTCDB column) including the associated TC numbers and the corresponding Uniprot identifiers

**Step 7. *fromReactions2Genes***: identification of the genes list associated to each model reaction.
* Input:
  * modelXml: the SBML model, which needs to be saved into the rawData directory
  * rxns: the last output of Step 5
  * transportRxns: the output of Step 6
  * outputName : a string to name the output file
  * metModelFile: the first output of Step 5
  * orgCode: the KEGG organism code of the target organism
  * taxId: a list of the NCBI Taxonomy ID of the target organism
* Output saved in the outputs directory:
  * outputName + \'.csv\': the second output of Step 5 is enriched with the retrieved identifiers from Step 5 (PutativeIdentifiers column) and 6 (Identifiers_fromTCDB column), and with the identified list of catalysing genes (lGenes). Each element of this last column is a list of lists where the genes identified in the internal lists are joined by the AND operator.


**Step 8. *genesLocationFilter***: apply a filter on the identified list of genes of each reaction to be in line with the compartment associated to the reaction.
* Input:
  * dfrxnsInfo: a string to name the output file
  * dfrxns2Genes: the output of Step 7
  * orgCode: the KEGG organism code of the target organism
  * lCompartmentsOrganization: according to the compartmentalization of the model, this dictionary includes as keys the cell compartment of the model and as values a list of compartments all pointing to the corresponding key. For example: mitochondrion is associated to the following list: ['mitochondrion', 'mitochondrial matrix', 'mitochondrial membrane','mitochondrion inner membrane', 'mitochondrion inner membrane', 'mitochondrial inner membrane', 'mitochondrial intermembrane space', 'mitochondrial outer membrane']. The list includes both alternative names with which the compartment is annotated in Uniprot and Gene Ontology database. Moreover, depending on the compartments included in the model if internal compartment to the target one are not included in the model as separate compartments, it is necessary to include them within the list to avoid that the gene is excluded from the filter, despite included in this cellular localization.
* Output saved in the outputs directory:
  * dfrxnsInfo + \'.csv\': a file storing for each gene (Gene column) the list of corresponding annotated compartments (lCompartments column)


**Step 9. *fromReactions2Genes_wFilteredData***: filter the list of genes retrieved for each reaction in accordance with its compartment into the model.
* Input:
  * modelXml: the SBML model, which needs to be saved into the rawData directory
  * dfGenes2Comp: the output of Step 8
  * dfrxns2Genes: the output of Step 7
  * outputFileName: a string to name the output file
  * lCompartmentsOrganization: according to the compartmentalization of the model, this dictionary includes as keys the cell compartment of the model and as values a list of compartments all pointing to the corresponding key. For example: mitochondrion is associated to the following list: ['mitochondrion', 'mitochondrial matrix', 'mitochondrial membrane','mitochondrion inner membrane', 'mitochondrion inner membrane', 'mitochondrial inner membrane', 'mitochondrial intermembrane space', 'mitochondrial outer membrane']. The list includes both alternative names with which the compartment is annotated in Uniprot and Gene Ontology database. Moreover, depending on the compartments included in the model if internal compartment to the target one are not included in the model as separate compartments, it is necessary to include them within the list to avoid that the gene is excluded from the filter, despite included in this cellular localization.
* Output saved in the outputs directory:
  * outputFileName + \'.csv\':  a file storing for each reaction (Rxn column) the corresponding list of catalysing genes filtered according to the compartment associated to the reaction (lGenes_filtered column).

**Step 10. *prepareGPRulerInput***: prepare the two input file of Step 11
* Input:
  * rxnswGenesFileName: the output of Step 9
  * outputFileName: a string to name the output file
  * organismCode: the KEGG organism code of the target organism
* Output saved in the outputs directory:
  * outputFileName + \'\_Rxns2Genes.csv\': a file storing for each reaction (Rxn column) the corresponding list of catalysing genes retrieved in Step 9 (Genes column).
  * outputFileName + \'\_Kegg2UniprotGenes.csv\': a file storing for each KEGG gene identifier of the investigated organism (keggId column) the corresponding Uniprot identifier (uniprotId column).


**Step 11. *GPRULER**: reconstruct the GPR rules
* Input:
  * model: a string to get the input files and name the output files
  * organismCode: the KEGG organism code of the target organism. The user has two options: insert the organism name (option 1) or insert the KEGG organism code (option 2). Choosing option 1, the user is asked to enter the organism name. The most putative KEGG organism codes will be proposed among which the user will choose the most correct one. Choosing option 2, the user will directly enter the KEGG code fof the target organism.
  * regexOrgSpecific: define the regular expression to extract the KEGG genes of the target organism. The user has two options: insert the correct regular expression (option 1) or infer this expression from the KEGG genes list (option 2)
* Output saved in the outputs directory:
  * model + \'\_GenesData.csv\': this file includes for each gene the corresponding information retrieved from Uniprot, Complex Portal, STRING and KEGG database about the established interactions with other proteins.
  * model + \'\_GenesRelationships.csv\': this file includes for each gene the list of the AND (AND column) and OR (OR column) established interactions
  * model + \'\_gprRules.csv\': for each reaction, this file includes the reaction name (Rxn column), the corresponding list of catalysing genes (Genes column), and the reconstructed GPR rule (GPR rule column).

## Execution of the pipeline from the organism name
When the organism name is inserted:
**Step 1. *fromOrganismName2ReactionsList**: generate the list of metabolic reactions of the target organism
* Input:
  * model = input('Which is the model name? ')
  * organismCode: the KEGG organism code of the target organism. The user has two options: insert the organism name (option 1) or insert the KEGG organism code (option 2). Choosing option 1, the user is asked to enter the organism name. The most putative KEGG organism codes will be proposed among which the user will choose the most correct one. Choosing option 2, the user will directly enter the KEGG code fof the target organism.
* Output saved in the outputs directory:
  * model + \'\_GeneId2Rxns.csv\': for each metabolic gene (GeneId column) the corresponding list of catalysed reaction is reported (Rxns column);
  * model + \'\_RxnId2Equation.csv\': for each retrieved metabolic reaction (RxnId column) the corresponding equation is reported expressing the involved compounds through their KEGG identifier (Equation column) or their complete name (Definition column);
  * model + \'\_RxnId2ECs.csv\': for each retrieved reaction (RxnId column), the corresponding list of EC numbers is reported (EC number) column;
  * model + \'\_Rxns2Genes.csv\': for each retrieved reaction (RxnId column), the list of corresponding metabolic catalysing genes is reported (Genes column);
  * model + \'\_Kegg2UniprotGenes.csv\': for each KEGG identifier (keggId column) of the retrieved list of metabolic genes the corresponding Uniprot (uniprotId column).

**Step 2. *fromReactions2Genes_fromOrgName**:  identification of the genes list associated to each model reaction according to the macrodatabase annotation.
* Input:
  * dfRxns2GenesFile: the fourth output from Step 1
  * dfRxnId2Equation: the second output from Step 1
  * orgCode: the KEGG organism code of the target organism
  * taxId: a list of the NCBI Taxonomy ID of the target organism
  * outputName: a string to name the output file
* Output saved in the outputs directory:
  * outputName + \'\_Rxns2Genes.csv\':  for each retrieved reaction (RxnId column), the list of corresponding metabolic catalysing genes is reported (Genes_fromKEGG column), and the list of catalysing genes annotated in the macrodatabase is reported (Genes_fromMacroDb column), and the sum of Genes_fromKEGG and Genes_fromMacroDb columns (Genes column).

**Step 3. *GPRULER**: reconstruct the GPR rules
* Input:
  * model: a string to get the input files and name the output files
  * organismCode: the KEGG organism code of the target organism. The user has two options: insert the organism name (option 1) or insert the KEGG organism code (option 2). Choosing option 1, the user is asked to enter the organism name. The most putative KEGG organism codes will be proposed among which the user will choose the most correct one. Choosing option 2, the user will directly enter the KEGG code fof the target organism.
  * regexOrgSpecific: define the regular expression to extract the KEGG genes of the target organism. The user has two options: insert the correct regular expression (option 1) or infer this expression from the KEGG genes list (option 2)
* Output saved in the outputs directory:
  * model + \'\_GenesData.csv\': this file includes for each gene the corresponding information retrieved from Uniprot, Complex Portal, STRING and KEGG database about the established interactions with other proteins.
  * model + \'\_GenesRelationships.csv\': this file includes for each gene the list of the AND (AND column) and OR (OR column) established interactions
  * model + \'\_gprRules.csv\': for each reaction, this file includes the reaction name (Rxn column), the corresponding list of catalysing genes (Genes column), and the reconstructed GPR rule (GPR rule column).
