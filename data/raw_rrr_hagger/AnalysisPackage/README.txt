BDEC Analysis Package
Developed at Chandra Sripada's (sripada@umich.edu) lab at University of Michigan Dept of Psychiatry.
Coded by Daniel Kessler (kesslerd@umich.edu).

NOTE: This README will look best read in a monospaced font (e.g. courier)

Unpack these scripts into the same directory, and edit BDEC_AnalysisDriver.R as follows:

1. Edit the second line so that it is the path where you have unpacked the scripts.

2. If you are running the task using a modified set of experiment scripts (e.g. not in English),
you may need to change lines 3-5 to reflect this. For each, you need to modify the entry so that it corresponds
to what you find in your edat files for each task.

Additionally, in this directory, place the following files

- SubjectStatus.csv
    Needs to have the following columns
    1 - SubjectID (consistent with subject IDs in all exported ePrime data)
    2 - Task ("E":Easy | "H":Hard)
    Any additional columns can contain between-subject phenotypics (e.g. gender).

- LetterE_EData.txt
    An exported e-merge file (using SPSS/Statview Unicode settings)

- MSIT200_EData.txt
    An exported e-merge file (using SPSS/Statview Unicode settings)

Run BDEC_AnalysisDriver.R using R. The first time you run it, it may ask you
about downloading and installing some R packages necessary to complete the processing.

It will produce the following output files (in the same directory as the script)
- BDEC_Results_Full.csv
    A CSV containing a wealth of subject-level statistics
- BDEC_Results_Lite.csv
    A CSV containing only the most critical variables (Subject, Task, Inclusion, and RTVar)

Included in Distribution:

- BDEC_AnalysisDriver.R
- func.R
- Cond_LoadParse.R
- LetE_LoadParse.R
- MSIT_LoadParse.R

BDEC_AnalysisDriver.R
This script sets things up, calls the other scripts, and ultimately produces
output in the form of "BDEC_Results_Full.csv" and "BDEC_Results_Lite.csv"

func.R
Contains some small functions that may be called in other scripts

Cond_LoadParse.R
Loads a "SubjectStatus.csv" file. Needs to have the following columns
1 - SubjectID (consistent with subject IDs in all exported ePrime data)
2 - Task ("E":Easy | "H":Hard)

LetE_LoadParse.R
Expects an exported e-merge file (using SPSS/Statview Unicode settings)
named "LetterE_EData.txt".

MSIT_LoadParse.R
Expects an exported e-merge file (using SPSS/Statview Unicode settings)
named "MSIT200_EData.txt".


Variable Key

-BDEC_Results_Lite.csv

Subject: Subject ID
Task: Assigned condition. E for easy, H for hard
Include.Overall: TRUE/FALSE: does this subject pass all inclusion criteria
ExGauss.I.RTVar.MSIT: RTVar (sum of sigma and tau from ExGauss fit) for accurate Incongruent Trials
ExGauss.C.RTVar.MSIT: RTVar (sum of sigma and tau from ExGauss fit) for accurate Congruent Trials

- BDEC_Results_Full.csv

| Variable Name                 | Description                                                                                                         |
|-------------------------------+---------------------------------------------------------------------------------------------------------------------|
| Subject                       | Subject ID                                                                                                          |
| Task                          | Assigned condition. E for Easy; H for hard                                                                          |
|-------------------------------+---------------------------------------------------------------------------------------------------------------------|
| Letter E Stats                | All accuracies expressed as proportion of 1. Time counted in seconds.                                               |
|-------------------------------+---------------------------------------------------------------------------------------------------------------------|
| Acc.Overall.LetE              | Overall accuracy on all trials                                                                                      |
| Acc.No.E.LetE                 | Accuracy on No E Trials                                                                                             |
| Acc.Non.Lonely.E.LetE         | Accuracy on Non-Lonely E Trials (note: this stat may be irrelevant for subjects in "easy" condition)                |
| Acc.Lonely.E.LetE             | Accuracy on Lonely E Trials (note: this stat may be irrelevant for subjects in "easy" condition)                    |
| Go_Acc.LetE                   | Accuracy on "Go" Trials (requirement for "Go" trial varies based on easy/hard condition)                            |
| NoGo_Acc.LetE                 | Accuracy on "NoGo" Trials (requirement for "NoGo" trial varies based on easy/hard condition)                        |
| SigDet_Correct.Rejection.LetE | Count of correct rejections                                                                                         |
| SigDet_False.Alarm.LetE       | Count of false alarms                                                                                               |
| SigDet_Hit.LetE               | Count of hits                                                                                                       |
| SigDet_Miss.LetE              | Count of misses                                                                                                     |
| Include.LetE                  | Does subject pass Letter E inclusion criteria (accurate > 80%)?                                                     |
| LetEPresent.LetE              | Does subject have Letter E data available?                                                                          |
|-------------------------------+---------------------------------------------------------------------------------------------------------------------|
| MSIT Stats                    | All accuracies expressed as proportion of 1. Time counted in seconds.                                               |
|-------------------------------+---------------------------------------------------------------------------------------------------------------------|
| C_1_MeanRT.MSIT               | Mean RT on correct congruent trials                                                                                 |
| I_1_MeanRT.MSIT               | Mean RT on correct incongruent trials                                                                               |
| Overall_1_RT.MSIT             | Mean RT on all correct trials                                                                                       |
| C_1_SDRT.MSIT                 | SD of RT on correct congruent trials                                                                                |
| I_1_SDRT.MSIT                 | SD of RT on correct incongruent trials                                                                              |
| Overall_1_SD.MSIT             | SD of RT on all correct trials                                                                                      |
| Acc.C.MSIT                    | Proportion of correct congruent trials                                                                              |
| Acc.I.MSIT                    | Proportion of correct incongruent trials                                                                            |
| Acc.Overall.MSIT              | Proportion of correct trials (overall)                                                                              |
| C_NonResponse.MSIT            | Proportion of nonresponses on congruent trials                                                                      |
| I_NonResponse.MSIT            | Proportion of nonresponses on incongruent trials                                                                    |
| TrialTypeEffect_Mean.MSIT     | Difference in mean RT on accurate trials (incongruent minus congruent)                                              |
| TrialTypeEffect_SD.MSIT       | Different in SD of RT on accurate trials (incongruent minus congruent)                                              |
| ExGauss.C.mu.MSIT             | Mu parameter from Ex Gaussian fit of accurate congruent trials                                                      |
| ExGauss.C.sigma.MSIT          | Sigma parameter from Ex Gaussian fit of accurate congruent trials                                                   |
| ExGauss.C.tau.MSIT            | Tau parameter from Ex Gaussian fit of accurate congruent trials                                                     |
| ExGauss.I.mu.MSIT             | Mu parameter from Ex Gaussian fit of accurate incongruent trials                                                    |
| ExGauss.I.sigma.MSIT          | Sigma parameter from Ex Gaussian fit of accurate incongruent trials                                                 |
| ExGauss.I.tau.MSIT            | Tau parameter from Ex Gaussian fit of accurate incongruent trials                                                   |
| ExGauss.C.RTVar.MSIT          | RTVar (sum of sigma and tau parameters) for accurate congruent trials                                               |
| ExGauss.I.RTVar.MSIT          | RTVar (sum of sigma and tau parameters) for accurate incongruent trials                                             |
| Include.Errors.MSIT           | Does the subject pass inclusion criteria based on accuracy (> 80%)?                                                 |
| Include.I2SD.MSIT             | Does the subject pass inclusion criteria based on incongruent mean RT being no more than two SDs above sample mean? |
| Include.IRTV2SD               | Does the subject pass inclusion criteria based on incongruent RTV being no more than two SDs above sample mean?     |
| Include.Task.MSIT             | Does the subject pass inclusion criteria based on having an assigned task condition (Easy/Hard)?                    |
| Include.Overall.MSIT          | Does the subject pass all MSIT-based inclusion criteria?                                                            |
| Include.MSITPresent.MSIT      | Does the subject pass inclusion criteria based on having MSIT data available?                                       |
|-------------------------------+---------------------------------------------------------------------------------------------------------------------|
| Include.Overall               | Does subject pass all study inclusion criteria?                                                                     |



**Developer Info**
VersionHash: 3fdbcafee691c0022f6b50cfbc1c9add66f51d9d
