# Tasmanian Devil
gene module: <br />
```bash
tas gene TASMANIAN_DEVIL/data/testing/151012_Gasch_glucose.txt -m TASMANIAN_DEVIL/data/testing/lgmncmodiMM904_NADcorrected_1127_MTHFDi.mat -o TASMANIAN_DEVIL/data/testing/151012_glucose_0.25.csv -c
```
flux module: <br />
simplest command
```bash
tas flux TASMANIAN_DEVIL/data/testing/lgmncmodiMM904_NADcorrected_1127_MTHFDi.mat TASMANIAN_DEVIL/data/testing/151012_ethanol_0.25.csv _e 1 1 1 -c
```
For multiple concurrent processes, repetitions of concurrent processes and repetitions with pruning reactions in a desired order by compartment
```bash
tas flux TASMANIAN_DEVIL/data/testing/lgmncmodiMM904_NADcorrected_1127_MTHFDi.mat TASMANIAN_DEVIL/data/testing/151012_ethanol_0.25.csv _e 2 2 2 -c -b 0.2879 -EXrxns TASMANIAN_DEVIL/data/testing/EXrxns.csv -EXtrrxns TASMANIAN_DEVIL/data/testing/EXtrrxns.csv -Othertrrxns TASMANIAN_DEVIL/data/testing/Othertrrxns.csv
```
visualization module: <br />
1 condition
```bash
tas visualization TASMANIAN_DEVIL/data/testing/lgmncmodiMM904_NADcorrected_1127_MTHFDi.mat TASMANIAN_DEVIL/data/testing/151012_glucose_0.25.csv TASMANIAN_DEVIL/data/testing/metabolicState_151012_glucose_0.25_lgmncmodiMM904_NADcorrected_1127_MTHFDi Glycolysis_PPP_Serine_Alanine_shortened 1 _e -c -c1 TASMANIAN_DEVIL/data/testing/RxnsClassifiedByExpression_151012_glucose_0.25_lgmncmodiMM904_NADcorrected_1127_MTHFDi.pkl -b1 TASMANIAN_DEVIL/data/testing/freqBasedRxns_151012_glucose_0.25_lgmncmodiMM904_NADcorrected_1127_MTHFDi.pkl
```
Comparing 2 conditions
```bash
tas visualization TASMANIAN_DEVIL/data/testing/lgmncmodiMM904_NADcorrected_1127_MTHFDi.mat TASMANIAN_DEVIL/data/testing/151012_glucose_0.25.csv TASMANIAN_DEVIL/data/testing/metabolicState_151012_glucose_0.25_lgmncmodiMM904_NADcorrected_1127_MTHFDi Glycolysis_PPP_Serine_Alanine_shortened 1 _e -c -c1 TASMANIAN_DEVIL/data/testing/RxnsClassifiedByExpression_151012_glucose_0.25_lgmncmodiMM904_NADcorrected_1127_MTHFDi.pkl -b1 TASMANIAN_DEVIL/data/testing/freqBasedRxns_151012_glucose_0.25_lgmncmodiMM904_NADcorrected_1127_MTHFDi.pkl -c2 TASMANIAN_DEVIL/data/testing/RxnsClassifiedByExpression_151012_ethanol_0.25_lgmncmodiMM904_NADcorrected_1127_MTHFDi.pkl -b2 TASMANIAN_DEVIL/data/testing/freqBasedRxns_151012_ethanol_0.25_lgmncmodiMM904_NADcorrected_1127_MTHFDi.pkl -m2 TASMANIAN_DEVIL/data/testing/lgmncmodiMM904_NADcorrected_1127_MTHFDi.mat -g2 TASMANIAN_DEVIL/data/testing/151012_ethanol_0.25.csv -f2 TASMANIAN_DEVIL/data/testing/metabolicState_151012_ethanol_0.25_lgmncmodiMM904_NADcorrected_1127_MTHFDi
```