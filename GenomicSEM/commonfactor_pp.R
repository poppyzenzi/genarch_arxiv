# common factor function model from step 3 of GSEM (R script on Eddie)


#To run using DWLS estimation#
CommonFactor_DWLS<- commonfactor(covstruc = LDSCoutput, estimation="DWLS")

#print CommonFactor_DWLs output#
CommonFactor_DWLS

$modelfit

$results
# use these output vals to put into path diagram


# Watch for Heywood case.
# Instances when the standardized factor loading exceeds 1 - the indicator has a negative residual variance.
# Not possible and  l produce non-interpretable estimates of model fit
# if so:  user-specified function below, should be used to impose a model constraint to keep the residual variance above 0.

# commonfactor.model<-'F1=~ NA*MDD + PTSD + ANX + ALCH
# F1~~1*F1