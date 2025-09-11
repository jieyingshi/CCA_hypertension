load_packages <- function(packages) {
  to_install <- packages[!(packages %in% installed.packages()[, "Package"])]
  if (length(to_install) > 0) {
    install.packages(to_install, dependencies = TRUE)
  }
  invisible(lapply(packages, library, character.only = TRUE))
}

load("~/data/ACEi_1.RData")
load("~/data/beta_1.RData")
load("~/data/CCB_1.RData")
load("~/data/diuretic_1.RData")
load("~/data/B.RData")
load("~/data_list_imp.Rdata")


# 缺失值处理 ----
load_packages(c("VIM","mice"))

print(names(B))
cols_str <- paste0('"', names(B), '"', collapse = ",");cat(cols_str)

dd <- aggr(B[,c("sex","age","dwell","ethnicity","work","marriage","education","alco",
                "smk","sport","still","dNa","mental","Fhis_hype","hype_duration_5",
                "bmi","HR","AF","cvd","hf","stroke","DM","ckdepi","CKD_risk","b_K",
                "b_Na","hcy","ua","dyslipidemia")],
           col=c("navyblue",'yellow'))


pdf("缺失值.pdf", width = 22, height = 14)  
aggr(B[,c("sex","age","dwell","ethnicity","work","marriage","education","alco",
          "smk","sport","still","dNa","mental","Fhis_hype","hype_duration_5",
          "bmi","HR","AF","cvd","hf","stroke","DM","ckdepi","CKD_risk","b_K",
          "b_Na","hcy","ua","dyslipidemia")],
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
              "bmi","HR","AF","cvd","hf","stroke","DM","ckdepi","CKD_risk","b_K",
              "b_Na","hcy","ua","dyslipidemia"), 
     cex.axis=0.6, gap=0, 
     ylab=c("Histogram  of missing data"," "))
dev.off()


##直接删除
pMiss <- function(x){round(sum(is.na(x))/length(x),3)}
cols_to_delete <- which(apply(B[,c("sex","age","dwell","ethnicity","work","marriage","education","alco",
                                   "smk","sport","still","dNa","mental","Fhis_hype","hype_duration_5",
                                   "bmi","HR","AF","cvd","hf","stroke","DM","ckdepi","CKD_risk","b_K",
                                   "b_Na","hcy","ua","dyslipidemia")], 2, pMiss) > 0.2);cols_to_delete


# baseline -----
load_packages(c("tableone"))
total <- bind_rows(ACEi_1, beta_1, CCB_1, diuretic_1);nrow(total) #138479

vars <- setdiff(names(total), c("patient_uid", "drug","dwell","dwell_c"))  
catVars <- c("sex",  "ethnicity", "work", "marriage", "education", 
             "alco", "smk", "sport", "still", "mental", "Fhis_hype", 
             "hype_duration_5", "hype_duration_10", "AF", "ccvd", "anti_glu", 
             "dm_report", "dm_2", "DM", "CKD_risk", "hua", "anti_chol", "dyslipidemia")

tab <- CreateTableOne(vars = vars, factorVars = catVars,
                      data = total, strata = "drug", addOverall = TRUE)

tab1 <- print(tab, showAllLevels = TRUE)
write.csv(tab1,"table1.csv")

# 结局描述 -----
load_packages(c("ggplot2","ggExtra"))
ggplot(total, aes(x=drug, y=sbpTTR1)) +
  geom_violin(
    width = 0.8,
    trim = TRUE,
    fill = "snow1",  
    color = "black",
    alpha = 0.6
  ) +
  geom_boxplot(
    aes(fill=drug),     
    width = 0.2,
    outlier.shape = NA,
    color = "black"
  ) +
  scale_fill_manual(values=c("lightgoldenrod1","darkseagreen3","aquamarine3","cadetblue3")) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size=12, angle=0, hjust=0.5),
    axis.text.y = element_text(size=12),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  ) +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Drug class", y = "SBP TTR")


# 特征筛选 ----
load_packages(c("Boruta","openxlsx"))

## boruta_complete.cases ----
# clean_data1 <-total[complete.cases(total[, c("sex","age","ethnicity","work","marriage","education","alco",
#                                              "smk","sport","still","dNa","mental","Fhis_hype","hype_duration_5",
#                                              "bmi","HR","AF","cvd","hf","stroke","DM","ckdepi","CKD_risk","ua","dyslipidemia")]), ]
# predictors <- c("sex","age","ethnicity","work","marriage","education","alco",
#                 "smk","sport","still","dNa","mental","Fhis_hype","hype_duration_5",
#                 "bmi","HR","AF","cvd","hf","stroke","DM","ckdepi","CKD_risk","ua","dyslipidemia")
# 
# formula <- as.formula(paste("sbpTTR1 ~", paste(predictors, collapse = " + ")))
# 
# KKK <- Boruta(formula ,
#               data = clean_data1,
#               ntree = 10)
# attStats(KKK)
# plot(KKK)

responses <- c("sbpTTR1", "SBPload1", "meanSBP")
data_list <- list(
  total = total,
  ACEi_1 = ACEi_1,
  beta_1 = beta_1,
  CCB_1 = CCB_1,
  diuretic_1 = diuretic_1
)
predictors <- c("sex","age","ethnicity","work","marriage","education","alco",
                "smk","sport","still","dNa","mental","Fhis_hype","hype_duration_5",
                "bmi","HR","AF","cvd","hf","stroke","DM","ckdepi","CKD_risk","ua","dyslipidemia")


Boruta_results <- list()
for(df_name in names(data_list)) {
  df <- data_list[[df_name]]
  
  for(resp in responses) {
    vars_needed <- c(predictors, resp)
    clean_df <- df[complete.cases(df[, vars_needed]), ]
    
    formula <- as.formula(paste(resp, "~", paste(predictors, collapse = " + ")))
    # 树 ----
    KKK <- Boruta(formula, data = clean_df, ntree = 1000)
    
    Boruta_results[[paste(df_name, resp, sep = "_")]] <- attStats(KKK)
    pdf_filename <- paste0("BorutaPlot_", df_name, "_", resp, ".pdf")
    pdf(pdf_filename, width = 16, height = 10)
    plot(KKK, main = paste(df_name, "-", resp), cex.axis = 1.5, las = 2)
    dev.off()
  }
}

wb <- createWorkbook()

for(name in names(Boruta_results)) {
  df <- Boruta_results[[name]]
  df_out <- cbind(Variable = rownames(df), as.data.frame(df))
  rownames(df_out) <- NULL
  
  addWorksheet(wb, name)
  writeData(wb, name, df_out)
}

saveWorkbook(wb, "Boruta_results.xlsx", overwrite = TRUE)


## boruta_mice ----
load_packages(c("parallel","mice","dplyr"))

mice_impute <- function(df) {
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
  complete(imp, 1)  
}
save(data_list_imp, file = "data_list_imp.Rdata")

## boruta_mice_Category_outcomes and predictors ----
add_categorical_vars <- function(df) {
  df <- df %>%
    mutate(
      age_c = ifelse(age >= 65, "old", "young"),
      
      bmi_c = case_when(
        bmi < 28 ~ "normal",
        bmi >= 28 ~ "obese",
        TRUE ~ NA_character_
      ),
      dNa_c = ifelse(dNa >= 15,"Excessive","Normal"),
      HR_c = ifelse(HR >= 80, "high","normal"),
      ckdepi_c = ifelse(ckdepi < 60, "ckd", "normal"),
      ua_c = case_when(
        sex == "Male" & ua >= 420 / 59.48 ~ 1,     
        sex == "Female" & ua >= 360 / 59.48 ~ 1, 
        TRUE ~ 0                        
      ),
      sbp_sd_c  = ifelse(sbp_sd  >= quantile(df$sbp_sd,   2/3, na.rm = TRUE),  "high", "low"),
      sbp_cv_c  = ifelse(sbp_cv  >= quantile(df$sbp_cv,   2/3, na.rm = TRUE),  "high", "low"),
      sbp_arv_c = ifelse(sbp_arv >= quantile(df$sbp_arv,  2/3, na.rm = TRUE), "high", "low"))
  
  return(df)
}
data_list_imp_c <- map(data_list_imp, add_categorical_vars)

predictors_c <- c("sex","age_c","ethnicity","work","marriage","education","alco",
                  "smk","sport","still","dNa_c","mental","Fhis_hype","hype_duration_5",
                  "bmi_c","HR_c","AF","cvd","hf","stroke","DM","ckdepi_c","CKD_risk","ua_c","dyslipidemia",
                  "sbp_sd_c","sbp_cv_c","sbp_arv_c")


## bootstrap + bruta ----
# df <-prepared_data$SBPload1[sample(1:nrow(prepared_data$SBPload1), size = 0.4 * nrow(prepared_data$SBPload1)), ]
# formula <- as.formula(paste("SBPload1_bin ~", paste(predictors_c, collapse = " + ")))
# start_time <- Sys.time()
# KKK <- Boruta(formula,
#               data = df,
#               ntree = 500,
#               maxRuns = 20,
#               doTrace = 3)
# end_time <- Sys.time()
# print(end_time - start_time)
# attStats(KKK)
# plot(KKK)

dir.create("Boruta_stepwise", showWarnings = FALSE)
dir.create("Boruta_summary", showWarnings = FALSE)

# 准备数据 ----
responses <- c("sbpTTR", "SBPload", "meanSBP")
predictors_c <- c("sex","age_c","ethnicity","work","marriage","education","alco",
                  "smk","sport","still","dNa_c","mental","Fhis_hype","hype_duration_5",
                  "bmi_c","HR_c","AF","cvd","hf","stroke","DM","ckdepi_c","CKD_risk","ua_c","dyslipidemia",
                  "sbp_sd_c","sbp_cv_c","sbp_arv_c")
df <- data_list_imp_c$total
df <- df %>% 
  rename(sbpTTR = sbpTTR1,
         SBPload = SBPload1)

prepared_data <- list()
for (resp in responses) {
  qs <- quantile(df[[resp]], probs = c(1/3, 2/3), na.rm = TRUE)
  df_temp <- df
  df_temp[[paste0(resp, "_bin")]] <- cut(df_temp[[resp]],
                                         breaks = c(-Inf, qs[1], qs[2], Inf),
                                         labels = c("Low", "Mid","High"),
                                         include.lowest = TRUE)
  df_sub <- df_temp %>% filter(.data[[paste0(resp, "_bin")]] != "Mid")
  df_sub[[paste0(resp, "_bin")]] <- droplevels(df_sub[[paste0(resp, "_bin")]])
  prepared_data[[resp]] <- df_sub
}

B <- 100
frac <- 0.5
save_every <- 5

Boruta_stability_results <- list()

for (resp in responses) {
  
  # 检查已完成的批次
  existing_files <- list.files("Boruta_stepwise", pattern = paste0("^", resp, "_bootstrap_.*\\.csv$"), full.names = TRUE)
  done_b <- c()
  if (length(existing_files) > 0) {
    done_b <- unlist(lapply(existing_files, function(f) {
      nums <- as.integer(unlist(regmatches(f, gregexpr("\\d+", f))))
      seq(nums[1], nums[2])
    }))
  }
  
  todo_b <- setdiff(seq_len(B), done_b)
  if (length(todo_b) == 0) {
    message(resp, " all bootstraps already done, skipping bootstrap.")
    next
  }
  
  df_sub <- prepared_data[[resp]]
  formula <- as.formula(paste0(paste0(resp, "_bin"), " ~ ", paste(predictors_c, collapse = " + ")))
  
  res_list <- list()
  
  for (b in todo_b) {
    set.seed(123 + b * 1000)
    id <- sample(nrow(df_sub), size = floor(frac * nrow(df_sub)), replace = FALSE)
    df_b <- df_sub[id, , drop = FALSE]
    
    K <- Boruta(formula, data = df_b, ntree = 500, maxRuns = 30, doTrace = 0, holdHistory = TRUE)
    att <- attStats(K)
    
    res_list[[b]] <- data.frame(
      var = rownames(att),
      decision = att$decision,
      normHits = att$normHits,
      medianImp = att$medianImp,
      run = b,
      stringsAsFactors = FALSE
    )
    
    message("-------------Bootstrap ", b, "/", B, " for ", resp,"---------------")
    
    # 每次按 save_every 保存
    if (b %% save_every == 0 || b == max(todo_b)) {
      df_save <- bind_rows(res_list)
      write.csv(df_save,
                file = paste0("Boruta_stepwise/", resp, "_bootstrap_", 
                              min(df_save$run), "_to_", max(df_save$run), ".csv"),
                row.names = FALSE)
      message("Saved bootstrap ", min(df_save$run), " to ", max(df_save$run), " for ", resp)
      res_list <- list()  # 清空
    }
  }
}



# 汇总所有 bootstrap
for (resp in responses) {
  files <- list.files("Boruta_stepwise", pattern = paste0("^", resp, "_bootstrap_.*\\.csv$"), full.names = TRUE)
  res_all <- bind_rows(lapply(files, read.csv))
  
  summary <- res_all %>%
    group_by(var) %>%
    summarise(
      sel_rate = mean(decision == "Confirmed"),
      mean_normHits = mean(normHits, na.rm = TRUE),
      mean_medianImp = mean(medianImp, na.rm = TRUE),
      Confirmed = sum(decision == "Confirmed"),
      Tentative = sum(decision == "Tentative"),
      Rejected = sum(decision == "Rejected"),
      .groups = "drop"
    ) %>%
    arrange(desc(sel_rate), desc(mean_normHits), desc(mean_medianImp))
  
  Boruta_stability_results[[resp]] <- summary
  
  # 写 Excel
  wb_file <- paste0("Boruta_summary/Boruta_stability_summary.xlsx")
  if (!file.exists(wb_file)) wb <- createWorkbook() else wb <- loadWorkbook(wb_file)
  if (resp %in% names(wb)) removeWorksheet(wb, resp)
  addWorksheet(wb, resp)
  writeData(wb, resp, summary)
  saveWorkbook(wb, wb_file, overwrite = TRUE)
  
  # 写 PDF
  pdf_file <- paste0("Boruta_summary/Boruta_stability_", resp, ".pdf")
  pdf(pdf_file, width = 16, height = 10)
  summary_plot <- summary[order(summary$mean_medianImp, decreasing = TRUE), ]
  colors <- colorRampPalette(c("lightblue", "navy"))(100)
  bar_colors <- colors[ceiling(summary_plot$sel_rate * 99) + 1]
  barplot(summary_plot$mean_medianImp,
          names.arg = summary_plot$var,
          las = 2,
          cex.names = 0.8,
          main = paste("Boruta stability -", resp),
          col = bar_colors,
          border = "grey")
  legend("topright", legend = c("sel_rate low", "sel_rate high"),
         fill = c("lightblue", "navy"), border = "grey", cex = 0.8)
  dev.off()
}
























# 信息量/雷达图 -----




















