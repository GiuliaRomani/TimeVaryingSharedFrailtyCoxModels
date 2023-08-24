## General description of the Project of the PACS course
The Project of the PACS course is associated to my work thesis, developed in the statistical context of survival analysis and based on the paper "Centre-Effect on Survival After Bone Marrow Transplantation: Application of Time-Dependent Frailty Models", written by C.M. Wintrebert, H. Putter, A.H. Zwinderman and J.C. van Houwelingen. 

Briefly, this paper describes three different "Time-Varying Shared Frailty Cox Models", that are elaborated and complex Cox regression models, at which a random term (called frailty) that varies in time is added for studying the temporal behaviour of a portion of data heterogeneity that cannot be explained by the sole regression coefficients. More details are provided in the project report and in my thesis. We apply these models to the PoliMi dataset.

However, no numerical codes fitting these time-varying models are available and the ones implemented by us in R are quite slow. Therefore, this project re-implements in C++ a consistent portion of our codes and suggest how to further speed up them through the introduction of parallel computing (precisely, OpenMP).

## Structure of the .txt files
The implemented codes need two input .txt files to work and they are called "DataIndividualsFile.txt" and "DataToolFile.txt". They contain the dataset on which the models will be applied, the optimal results and some other variables.

In detail:

DataIndividualsFile.txt: The first line of the dataset must contain the number of individual constituting the dataset and the number of regressors.
All the other lines must have the following structure:
- First element: a "string" corresponding to the cluster the individual belongs to. Actually, it must be composed of four letters and it must be enclosed by quotes.
- Block of regressors: categorical variables must be converted into dummy variables and numerical variables must be standardized.
- Last element: the time-to-event variable (i.e. when the event happens)

DataToolFile.txt: It is composed of different blocks read using GetPot and they are:
- TimeDomain
- Parameters
- DiscretizationStep
- Model
- ParallelVersion

## What the user needs to decide
Other than the variables related to the dataset, the user has to decide two important blocks of the previous list: Model and ParallelVersion.

Model: The user has to provide the numeric id of one of the three Time-Varying Shared Frailty Cox Model he/she wants to apply.
The only possibilities are: 
- (1) for the "Adapted Paik et al's Model"
- (2) for the "Centre-Specific Frailty Model with Power Parameter"
- (3) for the "Stochastic Time-Dependent Centre-Specific Frailty Model".

ParallelVersion: The C++ implementation of these models provides also a parallel version of some code sections, to speed up the computations.
The user has to provide three quantities:
- (n_threads) for the number of threads to be used for the parallel computing. If he/she does not want to use it, just set value (1).
- (chunk size) for the number of iterations each thread need to evaluate. It depends also on the following variables.
- (schedule_type) for the type of scheduling strategy to be used inside the parallel for loop. The only possibilities are:
    - (1) for "static"
    - (2) for "dynamic"
    - (3) for "guided"
    - (4) for "auto".
If any of the previous values is omitted, the project provides them a default value equal to (1). 

## Structure of the Project_PACS github folder
- BashScript: Several bash scripts are contained and each one refers to a precide academic year, except "bash_test.sh" that tests the models.
- Data: The .txt input files are divided in other two folders: DataTool and DataIndividuals.
- Doc: The slides and the report of the Project, plus the doxygen configuration file.
- Src: The implemented codes

## To execute a model:
- Enter the BashScript folder
- Choose a bash script (it's commented).
- Execute it on terminal. That's it.
- If you want to change the academic year, change the bash script. 
- If you want to change some input variables, you must change only the files in Data/DataTool and only the indicated variables.






