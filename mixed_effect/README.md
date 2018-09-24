# An example of mixed effects regression model selection in global health research
### Methods overview as in steps:
1. lasso regression to eliminate independent variables with zero coefficients
  + we use `glmnet` for this stage: a package that fits a generalized linear model via penalized maximum likelihood. The regularization path is computed for the lasso at a grid of values for the regularization parameter _lambda_.
  + a brief introduction of `glmnet` here: https://www.analyticsvidhya.com/blog/2017/06/a-comprehensive-guide-for-linear-ridge-and-lasso-regression/

2. forward and backward model selection on mixed effects models
  + I found this package called `AICcmodavg` that works efficiently with mixed effects regression models. I am using 
`aictab` and `bictab` to construct model selection tables based on my selected parameters, and calculates _AIC_ and _BIC_ values for each model matrix.
  + convenient link: https://www.rdocumentation.org/packages/AICcmodavg/versions/2.1-1/topics/AICcmodavg-package

3. cross validation on 1000 best models from step 2 based on their RMSE values
  + Because what we measure is a continuous variable ($ expenditure), i use _root mean square error (RMSE) - square root of the average of squared differences between prediction and actual observation_ as the final metric to measure model performace. Another reason i choose _RMSE_ is because we model spend of a certain disease as the log of odds of proportion 
(logit). Therefore RMSE it has the benefit of penalizing large errors more since our response variable (log of odds) is expected to be small. So in our case, we think an error of 6 is more than twice as bad as an error of 3.
