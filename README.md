# General description of the Project of the PACS course
The Project of the PACS course is associated to my work thesis, developed in the statistical context of "survival analysis" and based on the paper "Centre-Effect on Survival After Bone Marrow Transplantation: Application of Time-Dependent Frailty Models", written by C.M. Wintrebert, H. Putter, A.H. Zwinderman and J.C. van Houwelingen. 

Briefly, this paper describes three different "Time-Varying Shared Frailty Cox Models", that are elaborated and complex Cox regression models, at which a time-varying random term (called frailty) is added for studying the temporal behaviour of a portion of data heterogeneity that cannot be solely explained by the regression coefficients. More details are provided in the project report and in my thesis. 
We apply these models to the PoliMi dataset.

No numerical codes fitting these time-varying models are available and the ones implemented by me and my advisors in R are quite slow. Therefore, this project re-implements in C++ a consistent portion of the codes and suggests how to further speed up them through the introduction of parallel computing (precisely, OpenMP).

### Structure of the .txt files
The implemented codes require two input .txt files, that are called "DataIndividualsFile.txt" and "DataToolFile.txt". They contain the dataset on which the models will be applied, the optimal results and some other variables.
In detail:

"DataIndividualsFile.txt": The first line of the dataset must contain the number of individual constituting the dataset and the number of regressors.
All the other lines must have the following structure:
- First element: a "string" corresponding to the cluster the individual belongs to. Actually, it must be composed of four letters and it must be enclosed by quotes.
- Block of regressors: categorical variables must be converted into dummy variables and numerical variables must be standardized.
- Last element: the time-to-event variable 

"DataToolFile.txt": It is composed of different blocks, that are:
- TimeDomain
- Parameters
- DiscretizationStep
- Model
- ParallelVersion

Four couples of files are present: a couple refers to the test case, to test the codes, and each one of the others refers to a dataset that contains students enrolled only in a precise academic year (i.e. 2010, 2018, 2017-2018). These files differ in the number of students, TimeDomain and Parameters variables.

### What the user can do
While the variables related to the dataset are fixed, the user can change three important blocks of the previous list: DiscretizationStep, Model and ParallelVersion.

DiscretizationStep: It is the discretization step used for the numerical approximation of the second derivative of a precise function.

Model: The user has to provide the numeric id of one of the three Time-Varying Shared Frailty Cox Models he/she wants to apply.
The only possibilities are: 
- (1) for the "Adapted Paik et al's Model"
- (2) for the "Centre-Specific Frailty Model with Power Parameter"
- (3) for the "Stochastic Time-Dependent Centre-Specific Frailty Model".

ParallelVersion: The C++ implementation of these models provides also a parallel version of a code section, to speed up the computations.
The user has to provide three quantities:
- (n_threads) for the number of threads to be used for the parallel computing. If he/she does not want to use it, just set value (1).
- (chunk size) for the number of iterations each thread is going to evaluate. It depends also on the chosen scheduling strategy.
- (schedule_type) for the scheduling strategy to be used inside the parallel for loop. The only possibilities are:
    - (1) for "static"
    - (2) for "dynamic"
    - (3) for "guided"
    - (4) for "auto".
If any of the previous values is omitted, a default value equal to (1) or (0) is provided.

### Structure of the TimeVaryingSharedFrailtyCoxModels folder
- BashScript: several bash scripts are contained and each one refers to a precide academic year, except "bash_test.sh" that tests the models. They contain all the commands to compile and run the desired model, with the input and data provided. 
- Data: the .txt input files are grouped in other two folders: DataTool and DataIndividuals.
- Doc: the slides and the report of the Project, plus the doxygen configuration file.
- Src: the implemented codes.

### To execute a model:
- Change the variable "PACS_PATH" in the "Src/Makefile".
- Enter the "Data/DataTool" folder and select the file related to an academic year. If you want, change its variables as previously indicated and as also indicated in the file. 
- Enter the "BashScript" folder, choose the bash script associated to the selected academic year and execute it on terminal. That's it.

### Recommendations:
- If you want to change the academic year, change the bash script.
- If you want to change some input variables and play with the models, you must change only the files in "Data/DataTool" and only the indicated variables. Then, follow the instructions in the bash script to run the model, without compiling again the codes.
- If you want to measure the execution time required by sole the evaluation of the log-likelihood function, go into "Src/ModelDerived... .cpp" and follow the instructions of the method "evaluate_loglikelihood(...)".
- Any bash script contains the command for generating the doxygen documentation (make docs) and all the files are stored in the subfolder "Doc/doxygen". This command open the index.html doxygen page. There is also a command for removing the folder (make docsclean). Actually, they are both commented. If the user wants to use them, he/she has to change both the "OUTPUT_DIRECTORY" tag (line 61) and "INPUT" (line 867) tag of the doxygen configuration file.







