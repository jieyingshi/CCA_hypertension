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
              "smk","sport","still","dNa","mental","Fhis_hype","hype_duration_5",
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
              "smk","sport","still","dNa","mental","Fhis_hype","hype_duration_5",
              "dbp_arv","dbp_sd","dbp_cv","sbp_arv","sbp_sd","sbp_cv",
              "bmi","waist","HR","AF","cvd","hf","stroke_w_o_TIA","DM","ckdepi","CKD_risk","b_K",
              "b_Na","hcy","ua","dyslipidemia","drug_time"), 
     cex.axis=0.6, gap=0, 
     ylab=c("Histogram  of missing data"," "))
dev.off()
### 1.1 列 ----
pMiss <- function(x){round(sum(is.na(x))/length(x),3)}
cols_to_delete <- which(apply(total[,c("sex","age","dwell","ethnicity","work","marriage","education","alco",
                                       "smk","sport","still","dNa","mental","Fhis_hype","hype_duration_5",
                                       "dbp_arv","dbp_sd","dbp_cv","sbp_arv","sbp_sd","sbp_cv",
                                       "bmi","waist","HR","AF","cvd","hf","stroke_w_o_TIA","DM","ckdepi","CKD_risk","b_K",
                                       "b_Na","hcy","ua","dyslipidemia","drug_time")], 2, pMiss) > 0.2);cols_to_delete
# dwell waist   b_K  b_Na   hcy 
### 1.2 行 ----
total <- total[, !(colnames(total) %in% cols_to_delete)]
rows_to_delete <- which(apply(total[,c("sex","age","ethnicity","work","marriage","education","alco",
                                       "smk","sport","still","dNa","mental","Fhis_hype","hype_duration_5",
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
                  "smk","sport","still","dNa","mental","Fhis_hype","hype_duration_5","hype_duration_10",
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
                  "smk","sport","still","dNa","mental","Fhis_hype","hype_duration_5","hype_duration_10",
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
library(corrplot)
total <- bind_rows(data_list_imp$ACEi, data_list_imp$beta, data_list_imp$CCB, data_list_imp$diuretic);nrow(total) #138184
predictors <- c(c("sex","age","ethnicity","work","marriage","education","alco",
                  "smk","sport","still","dNa","mental","Fhis_hype","hype_duration_5","hype_duration_10",
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
                  "smk","sport","still","dNa","mental","Fhis_hype","hype_duration_5",
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


# boruta - RF ----
library(Boruta)
predictors <- c(c("sex","age","ethnicity","work","marriage","education","alco",
                  "smk","sport","still","dNa","mental","Fhis_hype","hype_duration_5",
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

# SIVS (lasso+glmnet) ----
# 依赖交叉验证的嵌入式特征选择
library(sivs)
library(varhandle)
library(glmnet)
library(foreach)
library(doParallel)
# df_sub <- data_list_imp$diuretic[, c("sbp_control_bad", predictors)]
# DATA <- data.matrix(~ . - 1, data = df_sub[, predictors])
# RESP <- as.factor(df_sub$sbp_control_bad)
# sivs_object <- sivs(x = DATA, y = RESP,
#                     iter.count = 20,
#                     nfolds = 10,
#                     parallel.cores ='grace')
# sivs_object$selection.freq

## 预处理-归一化 
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

## VIMP
calculate_vimp <- function(coef.df){
  # I(ci)：符号一致性
  tmp.I <- apply(coef.df, 1, function(r){
    r <- r[r != 0]
    if(length(r) == 0) return(0)
    if(all(r > 0) | all(r < 0)) return(1)
    return(0)
  })
  
  # mabs(ci)：非零系数绝对中位数
  tmp.mabs <- apply(coef.df, 1, function(r){
    r <- r[r != 0]
    if(length(r) == 0) return(0)
    return(abs(median(r)))
  })
  
  # |ci|：非零系数绝对均值
  tmp.abs <- apply(coef.df, 1, function(r){
    r <- r[r != 0]
    if(length(r) == 0) return(0)
    return(mean(abs(r)))
  })
  
  # IQR(ci)：四分位差
  tmp.iqr <- apply(coef.df, 1, function(r){
    r <- r[r != 0]
    if(length(r) == 0) return(0)
    return(IQR(r))
  })
  
  # VIMP 公式
  tmp.vimp <- (tmp.I * tmp.mabs * tmp.abs) / (1 + tmp.iqr)
  tmp.vimp[is.nan(tmp.vimp)] <- 0
  
  # 排序
  tmp.vimp <- sort(tmp.vimp, decreasing = TRUE)
  tmp.vimp <- tmp.vimp[!is.element(names(tmp.vimp), "(Intercept)")]
  
  return(tmp.vimp)
}

## Elastic Net + glmnet 
run_iterative_lasso <- function(df_sub, response_name, predictors,
                                iter.count = 100, test.ratio = 1/3, nfolds = 10,
                                ncores = 11, pdf_filename = NULL, main_title = NULL){
  
  RESP <- as.factor(df_sub[[response_name]])
  DATA <- as.matrix(normalize_data(df_sub, predictors))
  features <- colnames(DATA)
  feature_selected_count <- setNames(rep(0, length(features)), features)
  
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  
  iterative_res <- foreach(i = 1:iter.count, .packages = "glmnet") %dopar% {
    train_idx <- sample(1:nrow(DATA), size = floor((1 - test.ratio) * nrow(DATA)))
    x_train <- DATA[train_idx, ]
    y_train <- RESP[train_idx]
    
    cv_fit <- tryCatch({
      cv.glmnet(x=x_train, y=y_train, family="binomial", nfolds=nfolds,alpha=0.2)
    }, error=function(e) NULL)
    
    if(!is.null(cv_fit)){
      coef_vec <- coef(cv_fit, s="lambda.min")
      data.frame(names=rownames(coef_vec), coef=coef_vec[,1], stringsAsFactors=FALSE)
    } else {
      NULL
    }
  }
  
  stopCluster(cl)  
  
  ## 结果
  clean_iter <- iterative_res[!sapply(iterative_res, is.null)]
  coef.df <- Reduce(function(...){ merge(..., by="names", all=TRUE) },
                    lapply(seq_along(clean_iter), function(i){
                      temp <- clean_iter[[i]]
                      colnames(temp)[2] <- paste0("coef.iter", i)
                      temp
                    }))
  rownames(coef.df) <- coef.df$names
  coef.df <- coef.df[, -match("names", colnames(coef.df))]
  
  ### 选中次数
  for(i in seq_along(clean_iter)){
    coef_vec <- clean_iter[[i]]
    selected <- coef_vec$names[coef_vec$coef != 0]
    selected <- setdiff(selected, "(Intercept)")
    feature_selected_count[selected] <- feature_selected_count[selected] + 1
  }
  
  ### VIMP
  vimp_res <- calculate_vimp(coef.df)
  
  ### 图
  if(!is.null(pdf_filename)){
    pdf(pdf_filename, width=16, height=10)
    par(mfrow=c(1,2))
    
    # 1. 系数
    coef_df_plot <- coef.df[apply(coef.df,1,function(x) any(x!=0)), ]
    coef_df_plot <- t(coef_df_plot)
    coef_df_plot <- coef_df_plot[, order(apply(coef_df_plot,2,median), decreasing=TRUE)]
    if(is.null(main_title)) main_title <- "Iterative LASSO Coefficients"
    boxplot(coef_df_plot, col="darkolivegreen3",
            main=paste(main_title,"- Coefficients"),
            ylab="Coefficient", las=2, cex.axis=0.7)
    abline(h=0,col="gray",lty=2)
    
    # 2. VIMP
    barplot(vimp_res, col="skyblue", las=2,
            main=paste(main_title,"- VIMP"), ylab="VIMP")
    
    par(mfrow=c(1,1))
    dev.off()
  }
  
  return(list(
    iterative_res = iterative_res,
    coef.df = coef.df,
    selection_count = feature_selected_count,
    vimp = vimp_res
  ))
}


predictors <- c("sex","age","ethnicity","work","marriage","education","alco",
                "smk","sport","still","dNa","mental","Fhis_hype","hype_duration_5",
                "dbp_arv","sbp_arv",
                "bmi","HR","AF","cvd","hf","stroke_w_o_TIA","DM","ckdepi",
                "ua","dyslipidemia","drug_time")

sivs_results <- list()
for(df_name in names(data_list_imp)){
  df <- data_list_imp[[df_name]]
  df_sub <- df[, c("sbp_control_bad", predictors)]
  
  pdf_filename <- paste0("SIVS_", df_name, ".pdf")
  sivs_results[[df_name]] <- run_iterative_lasso(
    df_sub=df_sub,
    response_name="sbp_control_bad",
    predictors=predictors,
    main_title=df_name,
    pdf_filename=pdf_filename
  )
  
  vimp_df <- data.frame(
    Feature = names(sivs_results[[df_name]]$vimp),
    VIMP = as.numeric(sivs_results[[df_name]]$vimp),
    Selection_Count = sivs_results[[df_name]]$selection_count[names(sivs_results[[df_name]]$vimp)]
  )
  write.csv(vimp_df, paste0("VIMP_", df_name, ".csv"), row.names=FALSE)
}



# XGBoost - gbtree ----
library(xgboost)
library(dplyr)

predictors <- c("sex","age","ethnicity","work","marriage","education","alco",
                "smk","sport","still","dNa","mental","Fhis_hype","hype_duration_5",
                "dbp_arv","sbp_arv",
                "bmi","HR","AF","cvd","hf","stroke_w_o_TIA","DM","ckdepi",
                "ua","dyslipidemia","drug_time")

xgb_params <- list(
  booster = "gbtree",
  eta = 0.05,
  max_depth = 4,
  gamma = 4,
  subsample = 0.75,
  colsample_bytree = 0.75,
  objective = "binary:logistic",
  eval_metric = "logloss"
)

# 循环data_list_imp
for(name in names(data_list_imp)){
  df <- data_list_imp[[name]]
  y_train <- df$sbp_control_bad
  X_train <- df %>% select(all_of(predictors))
  # one-hot 编码
  X_train_mat <- model.matrix(~ . -1, data = X_train)
  
  xgb_train <- xgb.DMatrix(data = X_train_mat, label = y_train)
  xgb_model <- xgb.train(
    params = xgb_params,
    data = xgb_train,
    nrounds = 1000,
    verbose = 0
  )
  
  # 特征重要性
  importance_matrix <- xgb.importance(feature_names = colnames(X_train_mat), model = xgb_model)
  importance_matrix$OriginalFeature <- sapply(importance_matrix$Feature, function(x) {
    matched <- predictors[sapply(predictors, function(p) grepl(p, x))]
    if(length(matched) > 0) matched[1] else x
  })
  
  importance_summary <- importance_matrix %>%
    group_by(OriginalFeature) %>%
    summarise(TotalGain = sum(Gain)) %>%
    arrange(desc(TotalGain))
  
  write.csv(importance_summary, file = paste0("XGB_importance_", name, ".csv"), row.names = FALSE)
  
  pdf(file = paste0("XGB_importance_", name, ".pdf"), width = 10, height = 8)
  importance_plot <- importance_matrix[order(-importance_matrix$Gain), ]
  xgb.plot.importance(importance_plot, rel_to_first = TRUE, xlab = "Relative Gain",
                      main = paste("XGBoost Feature Importance -", name))
  dev.off()
}












# RFE - RF+XBG+glm ----
library(plyr)
library(caret)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

## importance 排序 ----
load("boruta.Rdata")
str(boruta_results$ACEi)
VIMP_ACEi_sivs <- read_csv("VIMP_ACEi.csv")
VIMP_beta_sivs <- read_csv("VIMP_beta.csv")
VIMP_CCB_sivs <- read_csv("VIMP_CCB.csv")
VIMP_diuretic_sivs <- read_csv("VIMP_diuretic.csv")
str(VIMP_ACEi_sivs)
XGB_importance_CCB <- read_csv("XGB_importance_CCB.csv")
XGB_importance_ACEi <- read_csv("XGB_importance_ACEi.csv")
XGB_importance_beta <- read_csv("XGB_importance_beta.csv")
XGB_importance_diuretic <- read_csv("XGB_importance_diuretic.csv")
str(XGB_importance_ACEi)

merge_three_importance <- function(boruta_df, vimp_df, xgb_df) {
  
  boruta_clean <- boruta_df %>%
    dplyr::select(Variable, importance) %>%
    dplyr::rename(Feature = Variable, Boruta_importance = importance)
  
  vimp_clean <- vimp_df %>%
    dplyr::select(Feature, LASSO_VIMP = VIMP)
  
  xgb_clean <- xgb_df %>%
    dplyr::rename(Feature = OriginalFeature, XGB_TotalGain = TotalGain) %>%
    dplyr::select(Feature, XGB_TotalGain)
  
  merged <- boruta_clean %>%
    full_join(vimp_clean, by = "Feature") %>%
    full_join(xgb_clean, by = "Feature")
  
  return(merged)
}

merged_all_list <- list(
  ACEi     = merge_three_importance(boruta_results$ACEi, VIMP_ACEi_sivs, XGB_importance_ACEi),
  beta     = merge_three_importance(boruta_results$beta, VIMP_beta_sivs, XGB_importance_beta),
  CCB      = merge_three_importance(boruta_results$CCB, VIMP_CCB_sivs, XGB_importance_CCB),
  diuretic = merge_three_importance(boruta_results$diuretic, VIMP_diuretic_sivs, XGB_importance_diuretic)
)

merged_all <- bind_rows(
  lapply(names(merged_all_list), function(drug) {
    merged_all_list[[drug]] %>% mutate(Drug = drug)
  })
)

merged_all$XGB_TotalGain[is.na(merged_all$XGB_TotalGain)] <- 0

merged_all <- merged_all %>%
  mutate(
    Boruta_norm = (Boruta_importance - min(Boruta_importance, na.rm = TRUE)) /
      (max(Boruta_importance, na.rm = TRUE) - min(Boruta_importance, na.rm = TRUE)),
    LASSO_norm  = (LASSO_VIMP - min(LASSO_VIMP, na.rm = TRUE)) /
      (max(LASSO_VIMP, na.rm = TRUE) - min(LASSO_VIMP, na.rm = TRUE)),
    XGB_norm    = (XGB_TotalGain - min(XGB_TotalGain, na.rm = TRUE)) /
      (max(XGB_TotalGain, na.rm = TRUE) - min(XGB_TotalGain, na.rm = TRUE)),
    Comprehensive_Importance = Boruta_norm + LASSO_norm + XGB_norm
  )

write.csv(merged_all,"importance_total.csv")


merged_long <- merged_all %>%
  pivot_longer(
    cols = c(Boruta_importance, LASSO_VIMP, XGB_TotalGain),
    names_to = "Importance_Type",
    values_to = "Importance"
  ) %>%
  mutate(Importance = tidyr::replace_na(Importance, 0))

# 对 LASSO_VIMP 取 log1p(sqrt
merged_long <- merged_long %>%
  mutate(
    Importance_trans = if_else(Importance_Type == "LASSO_VIMP",
                               log1p(sqrt(pmax(Importance, 0))), Importance)
  ) %>%
  group_by(Drug, Importance_Type) %>%
  mutate(
    minv = min(Importance_trans, na.rm = TRUE),
    maxv = max(Importance_trans, na.rm = TRUE),
    Importance_norm = if_else(
      maxv > minv,
      (Importance_trans - minv) / (maxv - minv),
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


## RFE ----
library(pROC)
library(xgboost)
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
importance_total <- read.csv("importance_total.csv", stringsAsFactors = FALSE)

predictors <- c("sex","age","ethnicity","work","marriage","education","alco",
                "smk","sport","still","dNa","mental","Fhis_hype","hype_duration_5",
                "bmi","HR","AF","cvd","hf","stroke_w_o_TIA","DM","ckdepi","CKD_risk","ua","dyslipidemia")

response <- "sbp_control_bad"













