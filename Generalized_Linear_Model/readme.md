## Generalized Linear Model Fitting for Student Performance Dataset:
The project for Generalized Linear Model with Instructor Prof.Hans-Georg Mueller </br>
Dataset: https://archive.ics.uci.edu/ml/datasets/student+performance
1. Binarize the students' final grade into Pass/Failure. Fit a logistic regression model with Group-Lasso Penalty to extract important factors related to exam failure. Choose the optimal amount of penalty via LOOCV(leave-one-out cross validation) and AUPR(Area Under Precision-recall Curve). 
2. Categorize the final grades into division(A, B, C, D, F). Fit a Multinomial Logistic Regression model(proportional hazard model) and perform the model diagnostics and goodness-of-fit test.
3. Apply PCA on students' grade in multiple stages(2 midterms, 1 finals), apply box-cox transformation on the first principle component and fit linear regression model to identify important features influencing students' progressive performance.
