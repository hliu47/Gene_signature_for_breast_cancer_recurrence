perma_generic_function = function(x, ...) {
  # message('perma_generic_function is on ...')
  UseMethod("perma_generic_name")
}
perma_generic_name.default = function(x, ...) {
  # message('perma_generic_function.default is on ...')
  if(!is.any(class(x) %in% "matrix")) stop('perma_generic_function is only implemented for matrix class.')
}

perma_generic_name.matrix = function(x, y, alpha = 1, tau = 1, ...) {
  # message('perma_generic_name.matrix is on ...')
  if(!is.matrix(x)) x <- as.matrix(x)
  if(!is.factor(y)) stop("y must be a factor")
  RET <- list(learn = list(y = y, X = x))
  # RET$k <- k
  RET$alpha = alpha
  RET$tau = tau
  RET$terms <- NULL
  RET$contrasts <- NULL
  RET$xlevels <- NULL
  RET$theDots <- list(...)
  class(RET) <- "perma_class"
  RET
}

print.perma_class <- function (x, ...)
{
  # cat(x$k, "-nearest neighbor model\n", sep = "")
  cat("Permanental approach with alpha", x$alpha, "and tau", x$tau, "\n", sep= " ")
  cat("Training set outcome distribution:\n")
  if(is.factor(x$learn$y)) {
    print(table(x$learn$y))  # classificiation give table.
  } else print(summary(x$learn$y))  # for regression give summary.
  
  cat("\n")
  invisible(x)
}

predict.perma_class <- function (object, newdata, type = c("prob", "class"), ...) {
  # message('predict.perma_class is on ...')
  type <- match.arg(type)
  if (!inherits(object, "perma_class"))
    stop("object not of class perma_class")
  if (!is.null(Terms <- object$terms)) {
    if (missing(newdata))
      newdata <- model.frame(object)
    else {
      newdata <- model.frame(as.formula(delete.response(Terms)),
                             newdata, na.action = function(x) x, xlev = object$xlevels)
    }
    x <- model.matrix(delete.response(Terms), newdata, contrasts = object$contrasts)
    xint <- match("(Intercept)", colnames(x), nomatch = 0)
    if (xint > 0)
      x <- x[, -xint, drop = FALSE]
  }
  else {
    x <- as.matrix(newdata)
  }
  
  argList <- list(
    train = object$learn$X,     # for train, objects$learn$X from fit above. 
    test = x,                   # for test, x is as.matrix(newdata)
    cl = object$learn$y,        # train response. 
    # k = object$k
    alpha = object$alpha,
    tau = object$tau
    )
  
  if(length(object$theDots) == 0) object$theDots <- list(prob = TRUE)
  if(any(names(object$theDots) == "prob")) object$theDots$prob <- TRUE
  
  argList <- c(argList, object$theDots)
  
  RET <- do.call(
    "permaTrain",
    argList)
  # saveRDS(RET, "RET.RDS")
  # print('RET is')
  # print(RET)
  if (type == "prob") {
    return(attr(RET, "prob"))    # give predicted probability
  }  else {
    # message('RET before levels to y')
    # print(RET)
    # print(class(RET))
    
    
    
    RET <- factor(RET, levels = levels(object$learn$y))   # give predicted class
    
    # message('predict.perma_class is done.')
    # message('RET after levels to y is')
    # print(RET)
    # message('object learn y levels are')
    # print(levels(object$learn$y))
    
    return(RET)
  }
  
}

permaTrain = function(train, test, cl, alpha, tau, prob = TRUE) {
  # message('permaTrain is on ...')
  # message(paste('train dim is', dim(train)[1], 'and', dim(train)[2]))
  perma_predict = perma_predict_give_prob(
    train_X = train,
    test_X = test,
    train_y = cl,
    perma_alpha = alpha,
    perma_tau = tau
  )
  # message('perma_predict is done.')
  predClassPerma = perma_predict$predClassPerma
  
  predProbClass1 = perma_predict$predProbClass1K2
  
  class_prob = data.frame(
    prob_class1 = predProbClass1,
    prob_class2 = 1 - predProbClass1
  )
  
  class_prob = as.matrix(class_prob)
  colnames(class_prob) = sort(levels(cl))     # sort the character levels of response.
  
  predClassPerma[predClassPerma == 1] = sort(levels(cl))[1]   # convert 1 to level 1
  predClassPerma[predClassPerma == 2] = sort(levels(cl))[2]   # convert 2 to level 2.
  
  res = predClassPerma   # here class is 1, 2 with respect to sort(levels(train_y))
  
  
  if (prob) attr(res, 'prob') = class_prob
  
  # message('res is done.')
  # print('res is')
  # print(res)
  res
}
