## General description of the Project
The Project of the PACS course is associated to my work thesis, developed in the statistical context of survival analysis. Briefly it describes and applies to a PoliMi dataset three different "Time-Varying Shared Frailty Cox Models", that are elaborated and complex Cox regression models, at which we add a random term (called frailty) that varies in time. 

The entire analysis is performed through the help of the software R. Due to the lack of the numerical codes fitting the time-varying models previously cited, we implement them in R, starting from scratch. This project re-implements in the C++ language, a consistent portion of our R codes but it lacks of a foundamental part, that gives us some problems in R. No pre-existsent C++ codes are available.

We suggest to the user to read the report of the project to know how the models are built and work. In this way, it will be easier to understand the next two paragraphs.

## Structure of the .txt files
A dataset, some optimal results and some necessary variables are provided in two .txt files to test the functioning of the models. 

DataIndividualsFile.txt: The first line of the dataset must contain the number of individual constituting the dataset and the number of regressors.
All the other lines must have the following structure:
- First element: a "string" corresponding to the cluster the individual belongs to. Actually, it must be composed of four letters and it must be enclosed by quotes.
- Block of regressors: categorical variables must be converted into dummy variables and numerical variables must be standardized.
- Last element: the time-to-event variable (i.e. when the dropout event happens)

DataToolFile.txt: It is composed of different blocks read using GetPot and they are:
- TimeDomain
- Parameters
- DiscretizationStep
- Model
- ParallelVersion

## What the user needs to decide
Other than the variables related to the dataset, the user has to decide two important blocks of the previous list: Model and ParallelVersion.

Model: The user has to provide the numeric id of one of the three Time-Varying Shared Frailty Cox Model he/she wants to apply.
The possibilities are: 
- (1) for the "Adapted Paik et al's Model"
- (2) for the "Centre-Specific Frailty Model with Power Parameter"
- (3) for the "Stochastic Time-Dependent Centre-Specific Frailty Model".
If a different value is chosen, the project throws an exception.

ParallelVersion: The C++ implementation of these models provides also a parallel version of some code sections, to speed up the computations.
The user has to provide three quantities:
- (n_threads) for the number of threads to be used for the parallel application. If he/she does not want to use it, just set value (1).
- (chunk size) for the ...
- (schedule_type) for the type of schedule to be used inside the for loop. The possibilities are:
    - (1) for "static"
    - (2) for "dynamic"
    - (3) for "guided"
    - (4) for "auto".
If any of the previous values is omitted, the project provides them a default value equal to (1). 
Only for the schedule_type, if a different value is chosen, the project throws an exception.





