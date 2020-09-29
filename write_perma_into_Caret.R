perma_caret = list(
  label = 'Permanental Approach',
  library = NULL,
  loop = NULL,
  type = 'Classification',
  parameters = data.frame(
    parameter = c('alpha', 'tau'),
    class = c('numeric', 'numeric'),
    label = c('Alpha for permanent calculation', 'Tau for covariance calculation') #
  ),
  
  grid = function (x, y, len = NULL, search = "grid") {
    # message('grid is on ...')
    if (is.numeric(len)) {
      out <- expand.grid(alpha = len, tau = len)                      
    } else if (length(strsplit(len, "_")[[1]]) == 3) {
      # len = start_end_by
      len_start = as.numeric(strsplit(len, "_")[[1]][1])
      len_end = as.numeric(strsplit(len, "_")[[1]][2])
      len_by = as.numeric(strsplit(len, "_")[[1]][3])

      out <- expand.grid(
        alpha = seq(len_start, len_end, len_by), 
        tau = seq(len_start, len_end, len_by))
      
    } else if (length(strsplit(len, "_")[[1]]) == 2) {
      # This mean first is alpha, second is tau.
      out = expand.grid(
        alpha=as.numeric(strsplit(len, "_")[[1]][1]),
        tau=as.numeric(strsplit(len, "_")[[1]][2])
      )
    }
    out
  },
  
  fit = function (x, y, wts, param, lev, last, classProbs, ...) {
    # message('fit is on ...')
    # message(paste('dim of x is', dim(x)[1], 'and', dim(x)[2]))
    # message(paste('length of y is', length(y)))
    
    out = perma_generic_function(
      x = as.matrix(x), 
      y = y, 
      # k = param$k,
      alpha = param$alpha,
      tau = param$tau,
      ...)
    # message('fit is done.')
    # message('fit out is')
    # print(out)
    out
  },
  predict = function (modelFit, newdata, submodels = NULL) {
    # message('predict is on ...')
    # message(paste('newdata dim is', dim(newdata)[1], 'and', dim(newdata)[2]))
    
    out <- predict(modelFit, newdata, type = "class")  # here give class prediction.
    # message('predict is done')
    # message('predict out is')
    # print(out)
    out
  },
  predictors = NULL,
  tags = 'Prototype models',
  prob = function (modelFit, newdata, submodels = NULL) {
    predict(modelFit, newdata, type = "prob")  # here give probability prediciton.
  },
  levels = NULL,
  sort = NULL
)