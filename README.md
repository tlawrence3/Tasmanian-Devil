# Tasmanian Devil
Each module's help page contains documentation about usage. The following are test examples and general advice for using each module.
**model module**: <br />
Simplest command for yeast metabolism:
```bash
tas model TASMANIAN_DEVIL/data/testing/iMM904_NADcorrected_1127_MTHFDi.xml _e -s > test_iMM904_1.txt
```
Adjusting for lower boundaries, gene2rxn mapings, and 1st attempt at adding adaptations to test fixing the biomass reaction for yeast metabolism:
```bash
tas model TASMANIAN_DEVIL/data/testing/iMM904_NADcorrected_1127_MTHFDi.xml _e -l TASMANIAN_DEVIL/data/testing/YPD_lb.csv -g TASMANIAN_DEVIL/data/testing/iMM904_NADcorrected_1127_MTHFDi_genes_genes2rxns.csv -a TASMANIAN_DEVIL/data/testing/iMM904_NADcorrected_1127_MTHFDi_adaptations_V1.csv -d TASMANIAN_DEVIL/data/testing/iMM904_NADcorrected_1127_MTHFDi_metabolite_dict.csv -s > test_iMM904_2.txt
```
After inspecting output for which reactions are unbalanced, try a 2nd adaptation to correct for these reactions. <br />
Adjusting for lower boundaries, gene2rxn mapings, and 2nd attempt at adding adaptations to test fixing unbalanced reactions for yeast metabolism:
```bash
tas model TASMANIAN_DEVIL/data/testing/iMM904_NADcorrected_1127_MTHFDi.xml _e -l TASMANIAN_DEVIL/data/testing/YPD_lb.csv -g TASMANIAN_DEVIL/data/testing/iMM904_NADcorrected_1127_MTHFDi_genes_genes2rxns.csv -a TASMANIAN_DEVIL/data/testing/iMM904_NADcorrected_1127_MTHFDi_adaptations_V2.csv -d TASMANIAN_DEVIL/data/testing/iMM904_NADcorrected_1127_MTHFDi_metabolite_dict.csv -s > test_iMM904_3.txt
```
After inspecting output for which reactions are unbalanced still after the adapations, force the model to be carbon balanced. The algorithm will identify any metabolites to be removed from the biomass reaction as needed. <br />
Adjusting for lower boundaries, gene2rxn mapings, adding adaptations to fix unbalanced reactions, accounting for nucleoside conversions and metabolite mapping complexes, removing metabolites without carbons, and remove inactive reactions and carbon balance the model for yeast metabolism:
```bash
tas model TASMANIAN_DEVIL/data/testing/iMM904_NADcorrected_1127_MTHFDi.xml _e -l TASMANIAN_DEVIL/data/testing/YPD_lb.csv -g TASMANIAN_DEVIL/data/testing/iMM904_NADcorrected_1127_MTHFDi_genes_genes2rxns.csv -a TASMANIAN_DEVIL/data/testing/iMM904_NADcorrected_1127_MTHFDi_adaptations_V2.csv -m TASMANIAN_DEVIL/data/testing/150723_iMM904_NADcorrected_1127_MTHFDi_metabolite_mappings.csv -n TASMANIAN_DEVIL/data/testing/150723_iMM904_NADcorrected_1127_MTHFDi_nucleotide_conversions.csv -d TASMANIAN_DEVIL/data/testing/iMM904_NADcorrected_1127_MTHFDi_metabolite_dict.csv -s -z -r > test_iMM904_4.txt
```
Simplest command for human metabolism:
```bash
tas model TASMANIAN_DEVIL/data/testing/Recon2_v04.xml _e -s > test_Recon2_1.txt
```
Model for human metabolism accounting for nucleoside conversions and metabolite mapping complexes, removing metabolites without carbons, and removing inactive reactions and carbon balancing the model:
```bash
tas model TASMANIAN_DEVIL/data/testing/Recon2_v04.xml _e -m TASMANIAN_DEVIL/data/testing/150722_Recon2.v04_metabolite_mappings.csv -n TASMANIAN_DEVIL/data/testing/150721_Recon2.v04_nucleotide_conversions.csv -d TASMANIAN_DEVIL/data/testing/Recon2_metabolite_carbon_dict4.csv -s -z -r > test_Recon2_2.txt
```
**gene module**: <br />
```bash
tas gene TASMANIAN_DEVIL/data/testing/151012_Gasch_glucose.txt -m TASMANIAN_DEVIL/data/testing/lgmncmodiMM904_NADcorrected_1127_MTHFDi.mat -o TASMANIAN_DEVIL/data/testing/151012_glucose_0.25.csv -c
```
**flux module**: <br />
simplest command:
```bash
tas flux TASMANIAN_DEVIL/data/testing/lgmncmodiMM904_NADcorrected_1127_MTHFDi.mat TASMANIAN_DEVIL/data/testing/151012_ethanol_0.25.csv _e 1 1 1 -c
```
For multiple concurrent processes, repetitions of concurrent processes and repetitions with pruning reactions in a desired order by compartment:
```bash
tas flux TASMANIAN_DEVIL/data/testing/lgmncmodiMM904_NADcorrected_1127_MTHFDi.mat TASMANIAN_DEVIL/data/testing/151012_ethanol_0.25.csv _e 2 2 2 -c -b 0.2879 -EXrxns TASMANIAN_DEVIL/data/testing/EXrxns.csv -EXtrrxns TASMANIAN_DEVIL/data/testing/EXtrrxns.csv -Othertrrxns TASMANIAN_DEVIL/data/testing/Othertrrxns.csv
```
**visualization module**: <br />
Visualizing 1 condition
```bash
tas visualization TASMANIAN_DEVIL/data/testing/lgmncmodiMM904_NADcorrected_1127_MTHFDi.mat TASMANIAN_DEVIL/data/testing/151012_glucose_0.25.csv TASMANIAN_DEVIL/data/testing/metabolicState_151012_glucose_0.25_lgmncmodiMM904_NADcorrected_1127_MTHFDi Glycolysis_PPP_Serine_Alanine_shortened 1 _e -c -c1 TASMANIAN_DEVIL/data/testing/RxnsClassifiedByExpression_151012_glucose_0.25_lgmncmodiMM904_NADcorrected_1127_MTHFDi.pkl -b1 TASMANIAN_DEVIL/data/testing/freqBasedRxns_151012_glucose_0.25_lgmncmodiMM904_NADcorrected_1127_MTHFDi.pkl
```
Comparing 2 conditions
```bash
tas visualization TASMANIAN_DEVIL/data/testing/lgmncmodiMM904_NADcorrected_1127_MTHFDi.mat TASMANIAN_DEVIL/data/testing/151012_glucose_0.25.csv TASMANIAN_DEVIL/data/testing/metabolicState_151012_glucose_0.25_lgmncmodiMM904_NADcorrected_1127_MTHFDi Glycolysis_PPP_Serine_Alanine_shortened 1 _e -c -c1 TASMANIAN_DEVIL/data/testing/RxnsClassifiedByExpression_151012_glucose_0.25_lgmncmodiMM904_NADcorrected_1127_MTHFDi.pkl -b1 TASMANIAN_DEVIL/data/testing/freqBasedRxns_151012_glucose_0.25_lgmncmodiMM904_NADcorrected_1127_MTHFDi.pkl -c2 TASMANIAN_DEVIL/data/testing/RxnsClassifiedByExpression_151012_ethanol_0.25_lgmncmodiMM904_NADcorrected_1127_MTHFDi.pkl -b2 TASMANIAN_DEVIL/data/testing/freqBasedRxns_151012_ethanol_0.25_lgmncmodiMM904_NADcorrected_1127_MTHFDi.pkl -m2 TASMANIAN_DEVIL/data/testing/lgmncmodiMM904_NADcorrected_1127_MTHFDi.mat -g2 TASMANIAN_DEVIL/data/testing/151012_ethanol_0.25.csv -f2 TASMANIAN_DEVIL/data/testing/metabolicState_151012_ethanol_0.25_lgmncmodiMM904_NADcorrected_1127_MTHFDi
```