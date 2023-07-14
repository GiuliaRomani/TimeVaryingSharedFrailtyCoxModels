The ProjectPACS is associated to the work Thesis of its author and it re-implements in the "C++" language some parts of the "R" codes for three different "Time-Varying Shared Frailty Cox Models", in the context of "Survival Analysis". 

With respect to the "R" codes, this project simply implements the models and do not optimize them. 
This implementation starts from scratch. No pre-existent model structure were available.

A dataset and some optimal results are already provided in a .txt file to test the correct functioning of the models. 
The first line of the dataset must contain the number of individual constituting the dataset and the number of regressors.
All the other lines must have the following structure:
- First element: a "string" corresponding to the cluster the individual belongs to. Actually, it must be composed of four letters and it must be enclosed by quotes.
- Block of regressors: categorical variables must be converted into dummy variables and numerical variables must be scaled to have null mean and unitary standard deviation
- Last element: the time-to-event variable (that is, when the event happens)

The user has only to decide wich method he/she wants to apply: 
- (1) for the "Adapted Paik et al's model"
- (2) for the "Centre-Specific Frailty Model" with Power Parameter
- (3) for the "Stochastic Time-Dependent Centre-Specific Frailty Model"
All models are able to plot also the standard deviation of the frailty, the basleine hazard function.
Only model (1) can compute ad plot the posterior frailty estimate. 

A suggestion: all models work perfectly with at most 3000 individuals. For greater number, model (1) may not work.



