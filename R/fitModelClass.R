# S4 Helper classes to allow NULL values
setClassUnion("numericOrNULL", c("numeric", "NULL"))
setClassUnion("dfOrNULL", c("data.frame", "NULL"))
setClassUnion("listOrNULL", c("list", "NULL"))
setClassUnion("characterOrNULL", c("character", "NULL"))
setClassUnion("glmOrNULL", c("glm", "lm", "NULL"))
setClassUnion("POSIXctOrNULL", c("POSIXct", "NULL"))

#' A S4 class representing an analysis
#' 
#' @name fitModel-class
#' @exportClass fitModel
fitModel <- setClass("fitModel",
                    slots=c(dateOfCreation="POSIXct",
                            lastUpdated="POSIXct",
			    data="data.frame",
			    meta="data.frame",
			    var="characterOrNULL",
			    type="characterOrNULL",
			    signTest="listOrNULL",
			    signP="numericOrNULL",
			    signP2="numericOrNULL",
			    signPAdj="numericOrNULL",
			    signID="characterOrNULL",
			    dataSign="dfOrNULL",
			    cv="dfOrNULL",
			    model="glmOrNULL",
			    pred="dfOrNULL"
                            ))


#' @title Constructor fitModel
#' 
#' @description Constructor for new fitModel instance
#' 
#' @param data feature data frame (features in rows, samples in cols)
#' @param meta metadata, with covariates as columns and samples in rows
#' @param type lm for linear model, lr for logistic regression
#' @param var variable of interest from the metadata dataframe
#' 
#' @param .Object fitModel object
setMethod("initialize", "fitModel",
          function(.Object,
		   data=data.frame,
		   meta=data.frame,
		   type=character,
		   var=character,
                   dateOfCreation=POSIXct) {

              .Object@dateOfCreation=Sys.time()
	      .Object@data = data
	      .Object@meta = meta
	      .Object@type = type
	      .Object@var = var
              .Object
          })

