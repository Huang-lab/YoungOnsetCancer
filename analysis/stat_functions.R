# some_stat_functions.R

# some examples of common scenarios encountered in implementing statisitical tests

##### Fisher's exact test given a universe, and setA, setB #####
two_set_fisher_test = function(universe, setA, setB, alternative = "greater"){
  # p = NA; OR = NA
  
  fisher_elements = c(sum(!universe %in% c(setA,setB)),sum((universe %in% setA) & !(universe %in% setB)),
                      sum(!(universe %in% setA) & (universe %in% setB)),sum((universe %in% setA) & (universe %in% setB)))
  if (fisher_elements[1] > 0 && fisher_elements[3] > 0 && fisher_elements[2] >= 0 && fisher_elements[4] >= 0){
    test.table = matrix(as.numeric(fisher_elements), nrow=2)
    f.test = fisher.test(test.table, alternative = alternative) # only test for enrichment for now
    OR = f.test$estimate
    p = f.test$p.value
    
    count00 = test.table[1,1]
    count10 = test.table[2,1]
    count01 = test.table[1,2]
    count11 = test.table[2,2]
    
    return(list("p"=p, "OR"=OR, "notAnotB" = count00
                , "AnotB" = count10, "BnotA" = count01, "AB" = count11))
    
    # result_row = c(fisher_elements,p,OR)
    # return(result_row)
  }
}

##### multivariate regression modeling #####

run_glm = function(data=NULL, xi = "", yi = "", covi=NULL, ytype="Continuous") {
  row_stat=NULL
  
  ## Check data input
  cat(paste("Processing: yi =", yi, " xi =", xi, " covi =", covi, "\n") )
  
  ## determine the formula, ie covariates to be included
  covars = ""
  # only add covariates if they have more than one level in the data subset
  for (covar in covi){ 
    if (sum(!is.na(unique(data[,covar])))>1){
      covars = paste(covars,covar,sep="+")
    } 
  }
  if (covars != ""){
    model=formula(paste(yi,"~",xi,covars))
  } else{
    model=formula(paste(yi,"~",xi))
  }
  
  test = ""
  ## fit the model
  if (ytype=="Continuous") {
    glmfit= try(glm(formula=model,data=data,family=gaussian(link = "identity")))
    test = "F" # F test
  } else if (ytype=="Binary") {
    glmfit= try(glm(formula=model,data=data,family=binomial(link = "logit")))
    test = "Chisq" # Chi Square test
  } else {
    stop("Unknown model ytype ", ytype)
  }
  
  if(class(glmfit)[1] == "try-error") {
    cat(paste("    Error caught fitting glm, continuing.  yi =", yi, " xi =", xi, " covi =", covi, "\n") )
    next
  }
  
  # retrieve coefficient
  if (length(names(coefficients(glmfit)))>1 & xi %in% names(coefficients(glmfit))){ 
    coeff = coefficients(glmfit)[[xi]]
  } else if (length(names(coefficients(glmfit)))>1) { # for binary the names become [xi][level1ofxi]
    coeff = coefficients(glmfit)[[2]]
  } else {
    coeff = NA
  }
  
  ## ANOVA
  fit = try(anova(glmfit,test=test))
  if(class(fit)[1] == "try-error") {
    cat(paste("    Error caught conducting ANOVA, continuing.  yi =", yi, " xi =", xi, " covi =", covi, "\n") )
    next
  } else {
    fit=as.matrix(fit)
    if (xi %in% rownames(fit)) (row_stat = cbind(yi,ytype,xi,as.data.frame(t(fit[xi,])),coeff,covars))
  }
  return(row_stat)
  
}