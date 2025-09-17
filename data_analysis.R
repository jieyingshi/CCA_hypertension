# feature reduction -----
load("data/ACEi_1.RData")
load("data/beta_1.RData")
load("data/CCB_1.RData")
load("data/diuretic_1.RData")
total <- bind_rows(ACEi_1, beta_1, CCB_1, diuretic_1);nrow(total) #138479

## 1. 缺失值处理 ----
library(VIM)
library(mice)
pdf("缺失值.pdf", width = 22, height = 14)  
aggr(total[,c("sex","age","dwell","ethnicity","work","marriage","education","alco",
              "smk","sport","dNa","mental","Fhis_hype","hype_duration_5",
              "dbp_arv","dbp_sd","dbp_cv","sbp_arv","sbp_sd","sbp_cv",
              "bmi","waist","HR","AF","cvd","hf","stroke_w_o_TIA","DM","ckdepi","CKD_risk","b_K",
              "b_Na","hcy","ua","dyslipidemia","drug_time")],
     col=c("navyblue",'yellow'), border=NA,
     combined=F,
     numbers=F, 
     varheightT=T,
     sortVars=T,
     bars=F,
     prop = T,
     only.miss=F,
     labels=c("sex","age","dwell","ethnicity","work","marriage","education","alco",
              "smk","sport","dNa","mental","Fhis_hype","hype_duration_5",
              "dbp_arv","dbp_sd","dbp_cv","sbp_arv","sbp_sd","sbp_cv",
              "bmi","waist","HR","AF","cvd","hf","stroke_w_o_TIA","DM","ckdepi","CKD_risk","b_K",
              "b_Na","hcy","ua","dyslipidemia","drug_time"), 
     cex.axis=0.6, gap=0, 
     ylab=c("Histogram  of missing data"," "))
dev.off()
### 1.1 列 ----
pMiss <- function(x){round(sum(is.na(x))/length(x),3)}
cols_to_delete <- which(apply(total[,c("sex","age","dwell","ethnicity","work","marriage","education","alco",
                                       "smk","sport","dNa","mental","Fhis_hype","hype_duration_5",
                                       "dbp_arv","dbp_sd","dbp_cv","sbp_arv","sbp_sd","sbp_cv",
                                       "bmi","waist","HR","AF","cvd","hf","stroke_w_o_TIA","DM","ckdepi","CKD_risk","b_K",
                                       "b_Na","hcy","ua","dyslipidemia","drug_time")], 2, pMiss) > 0.2);cols_to_delete
# dwell waist   b_K  b_Na   hcy 
### 1.2 行 ----
total <- total[, !(colnames(total) %in% cols_to_delete)]
rows_to_delete <- which(apply(total[,c("sex","age","ethnicity","work","marriage","education","alco",
                                       "smk","sport","dNa","mental","Fhis_hype","hype_duration_5",
                                       "dbp_arv","dbp_sd","dbp_cv","sbp_arv","sbp_sd","sbp_cv",
                                       "bmi","HR","AF","cvd","hf","stroke_w_o_TIA","DM","ckdepi","CKD_risk",
                                       "ua","dyslipidemia","drug_time")], 1, pMiss) > 0.2);length(rows_to_delete)
# 294
uids_to_delete <- total$patient_uid[rows_to_delete]
ACEi     <- ACEi_1[!ACEi_1$patient_uid %in% uids_to_delete, ]
beta     <- beta_1[!beta_1$patient_uid %in% uids_to_delete, ]
CCB      <- CCB_1[!CCB_1$patient_uid %in% uids_to_delete, ]
diuretic <- diuretic_1[!diuretic_1$patient_uid %in% uids_to_delete, ]
total      <- total[!total$patient_uid %in% uids_to_delete, ]

### 插补 ----
library(parallel)
library(dplyr)
library(mice)
responses <- c("sbpTTR1", "SBPload1", "meanSBP")

data_list <- list(
  ACEi = ACEi,
  beta = beta,
  CCB = CCB,
  diuretic = diuretic
)

predictors <- c(c("sex","age","ethnicity","work","marriage","education","alco",
                  "smk","sport","dNa","mental","Fhis_hype","hype_duration_5","hype_duration_10",
                  "dbp_arv","dbp_sd","dbp_cv","sbp_arv","sbp_sd","sbp_cv",
                  "bmi","HR","AF","cvd","hf","stroke_w_o_TIA","DM","ckdepi","CKD_risk",
                  "ua","dyslipidemia","drug_time"))

mice_impute <- function(df) {
  uid <- df$patient_uid
  vars_needed <- c(predictors, responses)
  df_sub <- df[, vars_needed]
  df_sub <- df_sub %>% mutate(across(where(is.character), as.factor))
  df_sub[df_sub == ""] <- NA
  meth <- make.method(df_sub)
  
  for(v in names(df_sub)) {
    if(is.numeric(df_sub[[v]])) {
      meth[v] <- "pmm"
    } else if(nlevels(as.factor(df_sub[[v]])) == 2) {
      meth[v] <- "logreg"  
    } else {
      meth[v] <- "polyreg"  
    }
  }
  
  imp <- mice(df_sub, m = 1, method = meth, maxit = 10, seed = 123)
  df_completed <- complete(imp, 1)
  df_completed <- dplyr::bind_cols(patient_uid = uid, df_completed)
  return(df_completed)}

ncores <- 2
cl <- makeCluster(ncores)
clusterExport(cl, varlist = c("mice_impute", "predictors", "responses"))
clusterEvalQ(cl, { library(dplyr); library(mice) })
data_list_imp <- parLapply(cl, data_list, mice_impute)
stopCluster(cl)

save(data_list_imp, file = "data_list_imp.Rdata")

## 2. 单变量分析 ----
library(dplyr)
library(broom)
library(parallel)
load("data_list_imp.Rdata")

data_list_imp <- lapply(data_list_imp, function(df) {
  ttr_cut <- quantile(df$sbpTTR1, probs = 1/3, na.rm = TRUE)   
  load_cut <- quantile(df$SBPload1, probs = 2/3, na.rm = TRUE) 
  
  df <- df %>%
    mutate(sbp_control_bad =case_when(
      sbpTTR1 <= ttr_cut & SBPload1 >= load_cut & meanSBP >= 130 ~ 1,  
      TRUE ~ 0
    ))
  return(df)
})

data_list_imp <- lapply(data_list_imp, function(df) {
  df$mental <- factor(df$mental,
                      levels = c("Good Mental Health", "Anxiety or Depression"))
  df$DM <- factor(df$DM,
                  levels = c("no","preDM","DM"))
  return(df)
})

responses <- c("sbpTTR1", "SBPload1", "meanSBP","sbp_control_bad")
predictors <- c(c("sex","age","ethnicity","work","marriage","education","alco",
                  "smk","sport","dNa","mental","Fhis_hype","hype_duration_5","hype_duration_10",
                  "dbp_arv","dbp_sd","dbp_cv","sbp_arv","sbp_sd","sbp_cv",
                  "bmi","HR","AF","cvd","hf","stroke_w_o_TIA","DM","ckdepi","CKD_risk",
                  "ua","dyslipidemia","drug_time"))


batch_univar_parallel <- function(data_list_imp, predictors, responses, ncores = 4) {
  cl <- makeCluster(ncores)
  clusterExport(cl, varlist = c("data_list_imp", "predictors", "responses"), envir = environment())
  clusterEvalQ(cl, {
    library(dplyr)
    library(broom)
  })
  
  results_list <- parLapply(cl, names(data_list_imp), function(drug) {
    df <- data_list_imp[[drug]]
    drug_results <- list()
    
    for(resp in responses) {
      for(pred in predictors) {
        formula_str <- paste(resp, "~", pred)
        is_binary <- all(df[[resp]] %in% c(0,1))
        numeric_vars <- setdiff(names(df)[sapply(df, is.numeric)], resp)
        data_std <- df %>% mutate(across(all_of(numeric_vars), scale))
        
        if(is_binary) {
          model <- glm(as.formula(formula_str), data = df, family = binomial)
          res <- broom::tidy(model) %>%
            rename(B = estimate, Std_Error = std.error) %>% 
            select(term, B, Std_Error)
          res$Model <- "Logistic"
          
          model_std <- glm(formula(formula_str), data = data_std, family = binomial)
          res_std <- broom::tidy(model_std) %>%
            mutate(x = statistic^2) %>%
            rename(beta = estimate) %>% 
            select(term, beta, x, p.value)
        } else {
          model <- lm(as.formula(formula_str), data = df)
          res <- broom::tidy(model) %>%
            rename(B = estimate, Std_Error = std.error) %>% 
            select(term, B, Std_Error)
          res$Model <- "Linear"
          model_std <- lm(formula(formula_str), data = data_std)
          res_std <- broom::tidy(model_std) %>%
            mutate(x = statistic^2) %>%
            rename(beta = estimate) %>% 
            select(term, beta, x, p.value)
        }
        
        res <- left_join(res, res_std, by = "term") 
        res <- res %>%
          mutate(Response = resp,
                 Predictor = pred,
                 Drug = drug)
        
        drug_results[[paste(resp, pred, sep = "_")]] <- res
      }
    }
    
    bind_rows(drug_results)
  })
  
  stopCluster(cl)
  
  final_df <- bind_rows(results_list)
  expected_cols <- c("Drug","Response", "Predictor", "term", "B", "Std_Error", "beta", "x", "p.value")
  if (!all(expected_cols %in% names(final_df))) {
    stop("Missing columns in the final data frame. Check model outputs.")
  }
  
  final_df %>%
    dplyr::select(Drug,Response, Predictor, term, B, Std_Error, beta, x, p.value)
}

final_results <- batch_univar_parallel(data_list_imp, predictors, responses)
final_results <- final_results %>%
  filter(term != "(Intercept)")

library(openxlsx)
wb <- createWorkbook()

for(d in unique(final_results$Drug)) {
  df_drug <- final_results %>% filter(Drug == d)
  addWorksheet(wb, d)
  
  start_row <- 1
  for(r in unique(df_drug$Response)) {
    df_sub <- df_drug %>% filter(Response == r)
    writeData(wb, sheet = d, x = r, startRow = start_row, startCol = 1)
    writeData(wb, sheet = d, x = df_sub, startRow = start_row + 1, startCol = 1)
    
    start_row <- start_row + nrow(df_sub) + 3
  }
}

saveWorkbook(wb, "final_results_by_drug.xlsx", overwrite = TRUE)

## 3. 共线性 ----
install.packages("C:/Users/zhangy/Downloads/zip_2.3.3.zip", 
                 repos = NULL, 
                 type = "win.binary")

library(corrplot)
total <- bind_rows(data_list_imp$ACEi, data_list_imp$beta, data_list_imp$CCB, data_list_imp$diuretic);nrow(total) #138184
predictors <- c(c("sex","age","ethnicity","work","marriage","education","alco",
                  "smk","sport","dNa","mental","Fhis_hype","hype_duration_5","hype_duration_10",
                  "dbp_arv","dbp_sd","dbp_cv","sbp_arv","sbp_sd","sbp_cv",
                  "bmi","HR","AF","cvd","hf","stroke_w_o_TIA","DM","ckdepi","CKD_risk",
                  "ua","dyslipidemia","drug_time"))

df_sub <- total %>%
  dplyr::select(dplyr::all_of(predictors)) %>%
  mutate(across(where(is.factor), as.numeric))

pdf("热图.pdf", width = 15, height = 15)
cor_mat <- cor(df_sub, use = "pairwise.complete.obs")
corrplot(cor_mat, method = "color", type = "upper", 
         tl.col = "black", tl.srt = 45, 
         title = "Correlation Heatmap", mar=c(0,0,1,0))
dev.off()

## 4. VIF ----
library(car)
responses <- c("sbpTTR1", "SBPload1", "meanSBP","sbp_control_bad")
predictors <- c(c("sex","age","ethnicity","work","marriage","education","alco",
                  "smk","sport","dNa","mental","Fhis_hype","hype_duration_5",
                  "dbp_arv","sbp_arv",
                  "bmi","HR","AF","cvd","hf","stroke_w_o_TIA","DM","ckdepi",
                  "ua","dyslipidemia","drug_time"))


vif_results <- list()

for(df_name in names(data_list_imp)){
  df <- data_list_imp[[df_name]]
  df[] <- lapply(df, function(x) if(is.character(x)) as.factor(x) else x)
  
  vif_list <- list()
  
  for(resp in responses){
    all_vars <- c(resp, predictors)
    df_sub <- df[, intersect(all_vars, names(df))]
    form <- as.formula(paste(resp, "~", paste(predictors, collapse = "+")))
    
    if(resp == "sbp_control_bad"){
      fit <- glm(form, data = df_sub, family = binomial)
    } else {
      fit <- lm(form, data = df_sub)
    }
    
    vif_values <- vif(fit)[,'GVIF^(1/(2*Df))']
    
    vif_list[[resp]] <- data.frame(
      Drug = df_name,
      Response = resp,
      Variable = names(vif_values),
      VIF = round(vif_values, 2)
    )
  }
  
  if(length(vif_list) > 0){
    vif_results[[df_name]] <- bind_rows(vif_list)
  }
}


# 合并所有药物结果
all_vif <- bind_rows(vif_results)

rownames(all_vif) <- NULL
vif_split <- list()
for(drug in unique(all_vif$Drug)){
  df_drug <- all_vif[all_vif$Drug == drug, ]
  vif_split[[drug]] <- split(df_drug, df_drug$Response)
}

library(openxlsx)
wb <- createWorkbook()

for(drug in names(vif_split)){
  df_list <- vif_split[[drug]]
  sheet_data <- do.call(rbind, lapply(names(df_list), function(resp){
    df <- df_list[[resp]]
    df$Response <- resp  
    df
  }))
  addWorksheet(wb, drug)
  writeData(wb, sheet = drug, sheet_data, rowNames = FALSE)
}
saveWorkbook(wb, "VIF_by_drug.xlsx", overwrite = TRUE)


# boruta  ----
library(Boruta)
predictors <- c(c("sex","age","ethnicity","work","marriage","education","alco",
                  "smk","sport","dNa","mental","Fhis_hype","hype_duration_5",
                  "dbp_arv","sbp_arv",
                  "bmi","HR","AF","cvd","hf","stroke_w_o_TIA","DM","ckdepi",
                  "ua","dyslipidemia","drug_time"))

# formula <- as.formula(paste("sbp_control_bad ~", paste(predictors, collapse = " + ")))
# 
# KKK <- Boruta(formula ,
#               data = data_list_imp$diuretic,
#               ntree = 10)
# attStats(KKK)[,1]
# plot(KKK)

boruta_results <- list()
boruta_plots <- list()

for(df_name in names(data_list_imp)){
  df <- data_list_imp[[df_name]]
  formula <- as.formula(paste("sbp_control_bad ~", paste(predictors, collapse = " + ")))
  
  set.seed(123)  
  KKK <- Boruta(formula, data = df, ntree = 1000, maxRuns = 50)  
  stats <- attStats(KKK)
  importance_df <- data.frame(
    Variable = rownames(stats),
    importance = stats[, "meanImp"],
    stringsAsFactors = FALSE
  )
  
  boruta_results[[df_name]] <- importance_df
  
  # 保存 plot 对象
  pdf_filename <- paste0("Boruta_", df_name, ".pdf")
  pdf(pdf_filename, width = 16, height = 10)
  plot(KKK, main = paste("Boruta Importance -", df_name), cex.axis = 1.5, las = 2)
  dev.off()
}
save(boruta_results,file = "boruta.Rdata")

# grid search + XGBoost ----
library(xgboost)
library(dplyr)
library(Ckmeans.1d.dp)
library(Matrix)
library(caret)

myCl <- makeCluster(detectCores() - 1)
registerDoParallel(myCl)

predictors <- c("sex","age","ethnicity","work","marriage","education","alco",
                "smk","sport","dNa","mental","Fhis_hype","hype_duration_5",
                "dbp_arv","sbp_arv","bmi","HR","AF","cvd","hf",
                "stroke_w_o_TIA","DM","ckdepi","ua","dyslipidemia","drug_time")

xgb.grid <- expand.grid(
  nrounds = c(25, 50, 100),
  eta = c(0.01, 0.05),
  max_depth = 3,
  gamma = 0,
  subsample = c(0.7, 0.8, 0.9),
  min_child_weight = c(1, 2),
  colsample_bytree = 1
)

xgb.control <- trainControl(
  method = "cv",
  number = 5,
  verboseIter = TRUE,
  returnData = FALSE,
  returnResamp = "none",
  classProbs = TRUE,
  allowParallel = TRUE
)

best_params_list <- list()

for(name in names(data_list_imp)){ 
  cat("Processing:", name, "\n") 
  df <- data_list_imp[[name]] 
  y_train <- factor(df$sbp_control_bad, levels = c(0,1), labels = c("class0", "class1")) 
  X_train <- df[, predictors] 
  X_train_mat <- model.matrix(~ . -1, data = X_train) 
  
  # Grid search + CV 
  xgb.train.model <- train( 
    x = X_train_mat, 
    y = y_train, 
    method = "xgbTree", 
    trControl = xgb.control, 
    tuneGrid = xgb.grid) 
  best_params <- xgb.train.model$bestTune 
  dtrain <- xgb.DMatrix(data = X_train_mat, label = as.numeric(y_train) - 1)
  
  params <- list(
    eta = best_params$eta,
    max_depth = best_params$max_depth,
    gamma = best_params$gamma,
    min_child_weight = best_params$min_child_weight,
    subsample = best_params$subsample,
    colsample_bytree = best_params$colsample_bytree,
    objective = "binary:logistic",
    eval_metric = "auc",
    nthread = detectCores() - 1
  )
  
  xgb.cv.model <- xgb.cv(
    params = params,
    data = dtrain,
    nrounds = 5000,
    nfold = 5,
    showsd = TRUE,
    metrics = "auc",
    stratified = TRUE,
    verbose = TRUE,
    print_every_n = 10L,
    early_stopping_rounds = 100
  )
  
  best_nrounds <- xgb.cv.model$best_iteration
  cat("Best nrounds for", name, ":", best_nrounds, "\n")
  best_params_list[[name]] <- list(
    bestTune = best_params,
    best_nrounds = best_nrounds
  )
  
  
  final_model <- xgboost(
    data = X_train_mat,
    label = as.numeric(y_train) - 1,
    max.depth = best_params$max_depth,
    eta = best_params$eta,
    nrounds = best_nrounds,
    gamma = best_params$gamma,
    min_child_weight = best_params$min_child_weight,
    subsample = best_params$subsample,
    colsample_bytree = best_params$colsample_bytree,
    objective = "binary:logistic",
    eval_metric = "auc",
    verbose = 0,
    nthread = detectCores() - 1
  )
  
  importance_matrix <- xgb.importance(feature_names = colnames(X_train_mat), model = final_model)
  
  importance_matrix$OriginalFeature <- sapply(importance_matrix$Feature, function(x) {
    matched <- predictors[sapply(predictors, function(p) grepl(p, x))]
    if(length(matched) > 0) matched[1] else x
  })
  
  importance_summary <- importance_matrix %>%
    group_by(OriginalFeature) %>%
    summarise(TotalGain = sum(Gain)) %>%
    arrange(desc(TotalGain))
  
  write.csv(importance_summary, file = paste0("XGB_importance_", name, ".csv"), row.names = FALSE)
  write.csv(importance_matrix, file = paste0("XGB_importance_c_", name, ".csv"), row.names = FALSE)
  
  pdf(file = paste0("XGB_importance_", name, ".pdf"), width = 10, height = 8)
  importance_plot <- importance_matrix[order(-importance_matrix$Gain), ]
  xgb.plot.importance(importance_plot, rel_to_first = TRUE, xlab = "Relative Gain",
                      main = paste("XGBoost Feature Importance -", name))
  dev.off()
}

stopCluster(myCl)
saveRDS(best_params_list, file = "best_params_list.rds")

# xgb_params <- list(
#   booster = "gbtree",
#   eta = 0.05,
#   max_depth = 4,
#   gamma = 4,
#   subsample = 0.75,
#   colsample_bytree = 0.75,
#   objective = "binary:logistic",
#   eval_metric = "logloss"
# )
# 
# # 循环data_list_imp
# for(name in names(data_list_imp)){
#   df <- data_list_imp[[name]]
#   y_train <- df$sbp_control_bad
#   X_train <- df[, predictors]
#   # one-hot 编码
#   X_train_mat <- model.matrix(~ . -1, data = X_train)
#   
#   xgb_train <- xgb.DMatrix(data = X_train_mat, label = y_train)
#   xgb_model <- xgb.train(
#     params = xgb_params,
#     data = xgb_train,
#     nrounds = 1000,
#     verbose = 0
#   )
#   
#   # 特征重要性
#   importance_matrix <- xgb.importance(feature_names = colnames(X_train_mat), model = xgb_model)
#   importance_matrix$OriginalFeature <- sapply(importance_matrix$Feature, function(x) {
#     matched <- predictors[sapply(predictors, function(p) grepl(p, x))]
#     if(length(matched) > 0) matched[1] else x
#   })
#   
#   importance_summary <- importance_matrix %>%
#     group_by(OriginalFeature) %>%
#     summarise(TotalGain = sum(Gain)) %>%
#     arrange(desc(TotalGain))
#   
#     write.csv(importance_summary, file = paste0("XGB_importance_", name, ".csv"), row.names = FALSE)
#   
#   pdf(file = paste0("XGB_importance_", name, ".pdf"), width = 10, height = 8)
#   importance_plot <- importance_matrix[order(-importance_matrix$Gain), ]
#   xgb.plot.importance(importance_plot, rel_to_first = TRUE, xlab = "Relative Gain",
#                       main = paste("XGBoost Feature Importance -", name))
#   dev.off()
# }


# LightGBM ----
library(lightgbm)
library(dplyr)
library(Matrix)   
library(Ckmeans.1d.dp)  
library(ggplot2)

predictors <- c("sex","age","ethnicity","work","marriage","education","alco",
                "smk","sport","dNa","mental","Fhis_hype","hype_duration_5",
                "dbp_arv","sbp_arv","bmi","HR","AF","cvd","hf","stroke_w_o_TIA",
                "DM","ckdepi","ua","dyslipidemia","drug_time")

lgb_params <- list(
  objective = "binary",
  metric = "binary_logloss",
  boosting = "gbdt",
  learning_rate = 0.05,
  num_leaves = 20,
  max_depth = 4,
  feature_fraction = 0.75,
  bagging_fraction = 0.75,
  bagging_freq = 1,
  min_data_in_leaf = 15,
  nthread = 8
)

for(name in names(data_list_imp)){
  df <- data_list_imp[[name]]
  y_train <- df$sbp_control_bad
  X_train <- df[, predictors]
  
  # one-hot 编码
  X_train_mat <- model.matrix(~ . -1, data = X_train)
  
  dtrain <- lgb.Dataset(data = X_train_mat, label = y_train)
  
  lgb_model <- lgb.train(
    params = lgb_params,
    data = dtrain,
    nrounds = 1000,
    verbose = 0,
  )
  
  importance_df <- lgb.importance(lgb_model)
  
  importance_df$OriginalFeature <- sapply(importance_df$Feature, function(x) {
    matched <- predictors[sapply(predictors, function(p) grepl(p, x))]
    if(length(matched) > 0) matched[1] else x
  })
  
  importance_summary <- importance_df %>%
    group_by(OriginalFeature) %>%
    summarise(TotalGain = sum(Gain)) %>%
    arrange(desc(TotalGain))
  write.csv(importance_summary, file = paste0("LGB_importance_", name, ".csv"), row.names = FALSE)
  
  pdf(file = paste0("LGB_importance_", name, ".pdf"), width = 10, height = 8)
  lgb.importance(lgb_model) %>% lgb.plot.importance(top_n = 30)
  title(main = paste("             LightGBM                                    -", name))
  dev.off()
}

# multi logistic ----
library(broom)
library(dplyr)

for(drug_name in names(data_list_imp)) {
  
  df <- data_list_imp[[drug_name]]
  formula_str <- paste("sbp_control_bad ~", paste(predictors, collapse = " + "))
  model_std <- glm(as.formula(formula_str), data = df, family = binomial)
  
  res_std <- broom::tidy(model_std) %>%
    mutate(x = statistic^2) %>%
    filter(term != "(Intercept)")
  res_std$var_base <- sapply(res_std$term, function(t) {
    matched <- predictors[sapply(predictors, function(p) grepl(p, t))]
    if(length(matched) > 0) matched[1] else t
  })
  
  importance_summary <- res_std %>%
    group_by(var_base) %>%
    summarise(TotalImportance = sum(x, na.rm = TRUE)) %>%
    arrange(desc(TotalImportance))
  
  write.csv(importance_summary, file = paste0("logi_", drug_name, ".csv"), row.names = FALSE)
  
}



# XGBoost - RFE ----
library(xgboost)
library(parallel)
library(doParallel)

myCl <- makeCluster(detectCores() - 1)
registerDoParallel(myCl)
set.seed(123)

predictors <- c("sex","age","ethnicity","work","marriage","education","alco",
                "smk","sport","dNa","mental","Fhis_hype","hype_duration_5",
                "dbp_arv","sbp_arv","bmi","HR","AF","cvd","hf",
                "stroke_w_o_TIA","DM","ckdepi","ua","dyslipidemia","drug_time")
best_params_list <- readRDS("best_params_list.rds")

rfe_results <- list()
for(drug_name in names(data_list_imp)){
  cat("Running RFE for:", drug_name, "\n")
  
  df <- data_list_imp[[drug_name]]
  y_train <- df$sbp_control_bad
  X_train <- df[, predictors]
  X_train_mat <- model.matrix(~ . -1, data = X_train)
  
  bestTune <- best_params_list[[drug_name]]$bestTune
  best_nrounds <- best_params_list[[drug_name]]$best_nrounds
  
  remaining_features <- colnames(X_train_mat)
  auc_history <- c()
  removed_order <- c()
  
  repeat {
    # DMatrix
    dtrain <- xgb.DMatrix(data = X_train_mat[, remaining_features, drop = FALSE],
                          label = y_train)
    cv <- xgb.cv(
      params = list(
        eta = bestTune$eta,
        max_depth = bestTune$max_depth,
        gamma = bestTune$gamma,
        min_child_weight = bestTune$min_child_weight,
        subsample = bestTune$subsample,
        colsample_bytree = bestTune$colsample_bytree,
        objective = "binary:logistic",
        eval_metric = "auc",
        nthread = detectCores() - 1
      ),
      data = dtrain,
      nrounds = best_nrounds,
      nfold = 10,
      stratified = TRUE,
      verbose = 0,
      early_stopping_rounds = 50
    )
    
    auc_history <- c(auc_history, max(cv$evaluation_log$test_auc_mean))
    
    model <- xgboost(
      data = dtrain,
      max.depth = bestTune$max_depth,
      eta = bestTune$eta,
      nrounds = best_nrounds,
      gamma = bestTune$gamma,
      min_child_weight = bestTune$min_child_weight,
      subsample = bestTune$subsample,
      colsample_bytree = bestTune$colsample_bytree,
      objective = "binary:logistic",
      eval_metric = "auc",
      verbose = 0,
      nthread = detectCores() - 1
    )
    
    importance_matrix <- xgb.importance(feature_names = remaining_features, model = model)
    min_gain_feat <- importance_matrix$Feature[which.min(importance_matrix$Gain)]
    
    if(length(remaining_features) <= 2) break
    remaining_features <- setdiff(remaining_features, min_gain_feat)
    removed_order <- c(removed_order, min_gain_feat)
    
    cat("Removed:", min_gain_feat, "Remaining features:", length(remaining_features),
        "Current CV AUC:", round(tail(auc_history,1),4), "\n")
  }
  removed_order <- c(removed_order, remaining_features)
  
  rfe_results[[drug_name]] <- list(
    final_features = remaining_features,
    auc_history = auc_history,
    removed_order = removed_order
  )
}

stopCluster(myCl)
saveRDS(rfe_results, file = "xgb_rfe_results_dynamicTune.rds")

## RFE+XGB 绘图 ----
library(xgboost)
library(dplyr)
library(Ckmeans.1d.dp)
drug_names <- c("ACEi", "beta", "CCB", "diuretic")

importance_summary_list <- list()

for(drug in drug_names){
  imp_df <- read_csv(paste0("XGB_importance_c_", drug, ".csv"))
  importance_summary_list[[drug]] <- imp_df
}


for(drug_name in drug_names){
  cat("Plotting RFE for:", drug_name, "\n")
  feature_order <- switch(drug_name,
                          "ACEi" = gain_ACEi,
                          "beta" = gain_beta,
                          "CCB" = gain_CCB,
                          "diuretic" = gain_diuretic)
  
  
  tmp_gain <- importance_summary_list[[drug_name]]$Gain
  names(tmp_gain) <- importance_summary_list[[drug_name]]$Feature
  barplot_height <- tmp_gain[feature_order]
  barplot_height <- barplot_height[!is.na(barplot_height)]
  feature_order <- names(barplot_height) 
  
  auc_history <- rfe_results[[drug_name]]$auc_history
  x_vals <- seq(from = length(barplot_height), to = length(barplot_height) - length(auc_history) + 1, by = -1)
  
  pdf(file = paste0("RFE_XGB_", drug_name, ".pdf"), width = 10, height = 8)
  par(mar = c(10, 4, 4, 4))  
  bar_positions <- barplot(barplot_height,
                           col = "steelblue1",
                           names.arg = NA,
                           ylim = c(0, max(barplot_height)*1.2),
                           main = paste("Feature Elimination -", drug_name),
                           ylab = "Feature Gain")
  text(x = bar_positions, 
       y = par("usr")[3] - max(barplot_height)*0.05, 
       labels = feature_order, 
       srt = 45,  
       adj = 1,   # 右对齐
       xpd = TRUE, 
       cex = 0.8)
  
  top_line <- max(barplot_height)*1.15
  text(x = bar_positions, 
       y = rep(top_line, length(bar_positions)), 
       labels = x_vals,
       col = "black", 
       cex = 1)
  
  par(new = TRUE)
  plot(x_vals, auc_history,
       type = "p", pch = 16, col = "black",
       axes = FALSE, xlab = "", ylab = "",
       ylim = c(min(auc_history)*0.95, max(auc_history)*1.05))
  axis(side = 4, col.axis = "black", cex.axis = 0.8)
  mtext("RFE AUROC", side = 4, line = 2.5, col = "black")
  
  suggestion_strictness <- 0.05
  if(!is.null(suggestion_strictness)){
    tmp.median.AUROCs <- auc_history
    for(i in suggestion_strictness){
      AUC_cutoff <- ((1 - i) * (max(tmp.median.AUROCs) - min(tmp.median.AUROCs))) + min(tmp.median.AUROCs)
      feature_cutoff <- max(which(tmp.median.AUROCs >= AUC_cutoff)) + 0.5
      
      abline(v = feature_cutoff, col = "red", lty = 2)
      tmp.cutoff.msg <- paste("Strictness:", round(i, 2))
      
      tmp.strictness.bg <- paste(rep("-", nchar(tmp.cutoff.msg) - 2), collapse = "")
      mtext(text = tmp.strictness.bg,
            side = 3,
            at = feature_cutoff,
            line = -12,
            cex = 0.8,
            col = "white",
            padj = 0.56,
            las = 2,
            family = "mono",
            font = 2)
      mtext(text = tmp.cutoff.msg,
            side = 3,
            at = feature_cutoff,
            line = -12,
            cex = 0.8,
            col = "red",
            las = 2,
            family = "mono",
            font = 2)
    }
  }
  dev.off()
}


# SIVS (LASSO+glmnet) ----
# 依赖交叉验证的嵌入式特征选择
library(sivs)
library(varhandle)
library(glmnet)
library(foreach)
library(doParallel)
normalize_data <- function(df, predictors){
  df_num <- df[, predictors]
  
  for(col in colnames(df_num)){
    if(is.factor(df_num[[col]]) | is.character(df_num[[col]])){
      df_num[[col]] <- as.numeric(as.factor(df_num[[col]]))
    }
    rng <- range(df_num[[col]], na.rm = TRUE)
    if(diff(rng) > 0){
      df_num[[col]] <- (df_num[[col]] - rng[1]) / diff(rng)
    } else {
      df_num[[col]] <- 0  
    }
  }
  
  return(df_num)
}

run_sivs_for_all <- function(data_list, predictors, iter.count = 100, nfolds = 10){
  results <- list()
  
  for(drug in names(data_list)){
    cat("Running SIVS for:", drug, "\n")
    
    df_sub <- data_list[[drug]][, c("sbp_control_bad", predictors)]
    RESP <- as.factor(df_sub[["sbp_control_bad"]])
    DATA <- as.matrix(normalize_data(df_sub, predictors))
    
    sivs_obj <- sivs(x = DATA, y = RESP,
                     iter.count = iter.count,
                     nfolds = nfolds,
                     parallel.cores = "grace")
    
    suggested <- sivs::suggest(sivs_obj)
    results[[drug]] <- list(
      sivs_vimp = sivs_obj$vimp,
      suggested_features = suggested
    )
    
    pdf(paste0("SIVS_plot_", drug, ".pdf"), width = 8, height = 6)
    layout(mat = matrix(c(1,2,
                          3,3),
                        nrow = 2, byrow = TRUE))
    plot(sivs_obj)
    layout(1)
    dev.off()
  }
  
  return(results)
}

sivs_results_all <- run_sivs_for_all(data_list_imp, predictors)

merge_sivs_results <- function(sivs_results_all){
  merged_list <- list()
  
  for(drug in names(sivs_results_all)){
    vimp <- sivs_results_all[[drug]]$sivs_vimp
    suggested <- sivs_results_all[[drug]]$suggested_features
    
    if(is.vector(vimp)) {
      df_vimp <- data.frame(Feature = names(vimp),
                            VIMP = as.numeric(vimp))
    } else {
      df_vimp <- vimp
      if(!all(c("Feature","VIMP") %in% colnames(df_vimp))){
        df_vimp <- data.frame(Feature = rownames(df_vimp),
                              VIMP = df_vimp[,1])
      }
    }
    
    df_vimp$Suggested <- ifelse(df_vimp$Feature %in% suggested, TRUE, FALSE)
    df_vimp$Drug <- drug
    
    merged_list[[drug]] <- df_vimp
  }
  
  merged_df <- do.call(rbind, merged_list)
  rownames(merged_df) <- NULL
  merged_df <- merged_df %>% arrange(Drug, desc(VIMP))
  return(merged_df)
}

sivs_summary_df <- merge_sivs_results(sivs_results_all)
write.csv(sivs_summary_df,"sivs.csv")


## 预处理-归一化 
# normalize_data <- function(df, predictors){
#   df_num <- df[, predictors]
#   
#   for(col in colnames(df_num)){
#     if(is.factor(df_num[[col]]) | is.character(df_num[[col]])){
#       df_num[[col]] <- as.numeric(as.factor(df_num[[col]]))
#     }
#     rng <- range(df_num[[col]], na.rm = TRUE)
#     if(diff(rng) > 0){
#       df_num[[col]] <- (df_num[[col]] - rng[1]) / diff(rng)
#     } else {
#       df_num[[col]] <- 0  
#     }
#   }
#   
#   return(df_num)
# }
# 
# ## VIMP
# calculate_vimp <- function(coef.df){
#   # I(ci)：符号一致性
#   tmp.I <- apply(coef.df, 1, function(r){
#     r <- r[r != 0]
#     if(length(r) == 0) return(0)
#     if(all(r > 0) | all(r < 0)) return(1)
#     return(0)
#   })
#   
#   # mabs(ci)：非零系数绝对中位数
#   tmp.mabs <- apply(coef.df, 1, function(r){
#     r <- r[r != 0]
#     if(length(r) == 0) return(0)
#     return(abs(median(r)))
#   })
#   
#   # |ci|：非零系数绝对均值
#   tmp.abs <- apply(coef.df, 1, function(r){
#     r <- r[r != 0]
#     if(length(r) == 0) return(0)
#     return(mean(abs(r)))
#   })
#   
#   # IQR(ci)：四分位差
#   tmp.iqr <- apply(coef.df, 1, function(r){
#     r <- r[r != 0]
#     if(length(r) == 0) return(0)
#     return(IQR(r))
#   })
#   
#   # VIMP 公式
#   tmp.vimp <- (tmp.I * tmp.mabs * tmp.abs) / (1 + tmp.iqr)
#   tmp.vimp[is.nan(tmp.vimp)] <- 0
#   
#   # 排序
#   tmp.vimp <- sort(tmp.vimp, decreasing = TRUE)
#   tmp.vimp <- tmp.vimp[!is.element(names(tmp.vimp), "(Intercept)")]
#   
#   return(tmp.vimp)
# }
# 
# ## LASSO + glmnet 
# run_iterative_LASSO <- function(df_sub, response_name, predictors,
#                                 iter.count = 100, test.ratio = 1/3, nfolds = 10,
#                                 ncores = 11, pdf_filename = NULL, main_title = NULL){
#   
#   RESP <- as.factor(df_sub[[response_name]])
#   DATA <- as.matrix(normalize_data(df_sub, predictors))
#   features <- colnames(DATA)
#   feature_selected_count <- setNames(rep(0, length(features)), features)
#   
#   cl <- makeCluster(ncores)
#   registerDoParallel(cl)
#   
#   iterative_res <- foreach(i = 1:iter.count, .packages = "glmnet") %dopar% {
#     train_idx <- sample(1:nrow(DATA), size = floor((1 - test.ratio) * nrow(DATA)))
#     x_train <- DATA[train_idx, ]
#     y_train <- RESP[train_idx]
#     
#     cv_fit <- tryCatch({
#       cv.glmnet(x=x_train, y=y_train, family="binomial", nfolds=nfolds,alpha=1)
#     }, error=function(e) NULL)
#     
#     if(!is.null(cv_fit)){
#       coef_vec <- coef(cv_fit, s="lambda.min")
#       data.frame(names=rownames(coef_vec), coef=coef_vec[,1], stringsAsFactors=FALSE)
#     } else {
#       NULL
#     }
#   }
#   
#   stopCluster(cl)  
#   
#   ## 结果
#   clean_iter <- iterative_res[!sapply(iterative_res, is.null)]
#   coef.df <- Reduce(function(...){ merge(..., by="names", all=TRUE) },
#                     lapply(seq_along(clean_iter), function(i){
#                       temp <- clean_iter[[i]]
#                       colnames(temp)[2] <- paste0("coef.iter", i)
#                       temp
#                     }))
#   rownames(coef.df) <- coef.df$names
#   coef.df <- coef.df[, -match("names", colnames(coef.df))]
#   
#   ### 选中次数
#   for(i in seq_along(clean_iter)){
#     coef_vec <- clean_iter[[i]]
#     selected <- coef_vec$names[coef_vec$coef != 0]
#     selected <- setdiff(selected, "(Intercept)")
#     feature_selected_count[selected] <- feature_selected_count[selected] + 1
#   }
#   
#   ### VIMP
#   vimp_res <- calculate_vimp(coef.df)
#   
#   ### 图
#   if(!is.null(pdf_filename)){
#     pdf(pdf_filename, width=16, height=10)
#     par(mfrow=c(1,2))
#     
#     # 1. 系数
#     coef_df_plot <- coef.df[apply(coef.df,1,function(x) any(x!=0)), ]
#     coef_df_plot <- t(coef_df_plot)
#     coef_df_plot <- coef_df_plot[, order(apply(coef_df_plot,2,median), decreasing=TRUE)]
#     if(is.null(main_title)) main_title <- "Iterative LASSO Coefficients"
#     boxplot(coef_df_plot, col="darkolivegreen3",
#             main=paste(main_title,"- Coefficients"),
#             ylab="Coefficient", las=2, cex.axis=0.7)
#     abline(h=0,col="gray",lty=2)
#     
#     # 2. VIMP
#     barplot(vimp_res, col="skyblue", las=2,
#             main=paste(main_title,"- VIMP"), ylab="VIMP")
#     
#     par(mfrow=c(1,1))
#     dev.off()
#   }
#   
#   return(list(
#     iterative_res = iterative_res,
#     coef.df = coef.df,
#     selection_count = feature_selected_count,
#     vimp = vimp_res
#   ))
# }
# 
# 
# predictors <- c("sex","age","ethnicity","work","marriage","education","alco",
#                 "smk","sport","dNa","mental","Fhis_hype","hype_duration_5",
#                 "dbp_arv","sbp_arv",
#                 "bmi","HR","AF","cvd","hf","stroke_w_o_TIA","DM","ckdepi",
#                 "ua","dyslipidemia","drug_time")
# 
# sivs_results <- list()
# for(df_name in names(data_list_imp)){
#   df <- data_list_imp[[df_name]]
#   df_sub <- df[, c("sbp_control_bad", predictors)]
#   
#   pdf_filename <- paste0("SIVS_", df_name, ".pdf")
#   sivs_results[[df_name]] <- run_iterative_LASSO(
#     df_sub=df_sub,
#     response_name="sbp_control_bad",
#     predictors=predictors,
#     main_title=df_name,
#     pdf_filename=pdf_filename
#   )
#   
#   vimp_df <- data.frame(
#     Feature = names(sivs_results[[df_name]]$vimp),
#     VIMP = as.numeric(sivs_results[[df_name]]$vimp),
#     Selection_Count = sivs_results[[df_name]]$selection_count[names(sivs_results[[df_name]]$vimp)]
#   )
#   write.csv(vimp_df, paste0("VIMP_", df_name, ".csv"), row.names=FALSE)
# }

## importance 排序 ----
load("boruta.Rdata")
XGB_importance_CCB <- read_csv("XGB_importance_CCB.csv")
XGB_importance_ACEi <- read_csv("XGB_importance_ACEi.csv")
XGB_importance_beta <- read_csv("XGB_importance_beta.csv")
XGB_importance_diuretic <- read_csv("XGB_importance_diuretic.csv")
# LGB_importance_ACEi <- read_csv("LGB_importance_ACEi.csv")
# LGB_importance_CCB <- read_csv("LGB_importance_CCB.csv")
# LGB_importance_beta <- read_csv("LGB_importance_beta.csv")
# LGB_importance_diuretic <- read_csv("LGB_importance_diuretic.csv")
logi_ACEi <- read_csv("logi_ACEi.csv")
logi_CCB <- read_csv("logi_CCB.csv")
logi_diuretic <- read_csv("logi_diuretic.csv")
logi_beta <- read_csv("logi_beta.csv")
str(boruta_results$ACEi)
str(XGB_importance_ACEi)
# str(LGB_importance_ACEi)
str(logi_ACEi)


merge_four_importance <- function(boruta_df, logi_df, xgb_df) {
  
  # Boruta
  boruta_clean <- dplyr::select(boruta_df, Variable, importance) %>%
    dplyr::rename(Feature = Variable, Boruta_importance = importance)
  
  # 逻辑回归
  logi_clean <- dplyr::select(logi_df, var_base, TotalImportance) %>%
    dplyr::rename(Feature = var_base, Logi_importance = TotalImportance)
  
  # XGBoost
  xgb_clean <- dplyr::select(xgb_df, OriginalFeature, TotalGain) %>%
    dplyr::rename(Feature = OriginalFeature, XGB_importance = TotalGain)
  
  # # LightGBM
  # lgb_clean <- dplyr::select(lgb_df, OriginalFeature, TotalGain) %>%
  #   dplyr::rename(Feature = OriginalFeature, LGB_importance = TotalGain)
  
  # 合并
  merged <- boruta_clean %>%
    dplyr::full_join(logi_clean, by = "Feature") %>%
    dplyr::full_join(xgb_clean, by = "Feature") 
  # %>% dplyr::full_join(lgb_clean, by = "Feature")
  
  return(merged)
}


merged_all_list <- list(
  ACEi     = merge_four_importance(boruta_results$ACEi, logi_ACEi, XGB_importance_ACEi),
  beta     = merge_four_importance(boruta_results$beta, logi_beta, XGB_importance_beta),
  CCB      = merge_four_importance(boruta_results$CCB, logi_CCB, XGB_importance_CCB),
  diuretic = merge_four_importance(boruta_results$diuretic, logi_diuretic, XGB_importance_diuretic)
)

merged_all <- dplyr::bind_rows(
  lapply(names(merged_all_list), function(drug) {
    merged_all_list[[drug]] %>% dplyr::mutate(Drug = drug)
  })
)%>%
  mutate(across(where(is.numeric), ~ tidyr::replace_na(.x, 0)))

write.csv(merged_all,"importance_total.csv")



merged_long <- merged_all %>%
  pivot_longer(
    cols = c(Boruta_importance, Logi_importance , XGB_importance),
    names_to = "Importance_Type",
    values_to = "Importance"
  ) 

merged_long <- merged_long %>%
  group_by(Drug, Importance_Type) %>%
  mutate(
    minv = min(Importance, na.rm = TRUE),
    maxv = max(Importance, na.rm = TRUE),
    Importance_norm = if_else(
      maxv > minv,
      (Importance - minv) / (maxv - minv),
      0
    )
  ) %>%
  ungroup() %>%
  dplyr::select(-minv, -maxv)

merged_long <- merged_long %>%
  group_by(Drug, Feature) %>%
  mutate(Total_Imp_norm = sum(Importance_norm, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(Drug, desc(Total_Imp_norm)) %>%
  mutate(Feature = factor(Feature, levels = unique(Feature)))

ggplot(merged_long, aes(x = Importance_Type, y = Feature, fill = Importance_norm)) +
  geom_tile(color = "white") +
  facet_wrap(~Drug, scales = "free_y") +
  scale_fill_gradient(low = "white", high = "steelblue", na.value = "grey90") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 8),
    strip.text = element_text(face = "bold", size = 12),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Normalized Feature Importance Heatmap Across Drugs",
    x = "Importance Type",
    y = "Feature",
    fill = "Normalized Importance (0-1)"
  )

ggsave("feature_importance_heatmap.pdf", width = 12, height = 8)
