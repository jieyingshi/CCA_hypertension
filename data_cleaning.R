load_packages <- function(packages) {
  to_install <- packages[!(packages %in% installed.packages()[, "Package"])]
  if (length(to_install) > 0) {
    install.packages(to_install, dependencies = TRUE)
  }
  invisible(lapply(packages, library, character.only = TRUE))
}

my_packages <- c(
  "tidyverse",
  "rms",
  "MASS",
  "mice",
  "caret",
  "uuid",
  "digest",
  "data.table",
  "lubridate",
  "rlang",
  "pracma"
)

load_packages(my_packages)



# data_clean -----
baseline <- fread("~/data/基线.csv", header = TRUE, encoding = "UTF-8")
print(dim(baseline)) #[1] 640418   216
length(unique(baseline$patient_uid)) #[1] 640418

baseline_drug <- fread("~/data/基线用药.csv", header = TRUE, encoding = "UTF-8")
print(dim(baseline_drug)) 
#[1] 2629971      13

print(names(baseline_drug))
summary(as.factor(baseline_drug$H0182))
# [1] "drug_id"     "patient_uid" "hospital_id" 
# "H0182"(药物类型)  
# 1=血管紧张素转化酶抑制剂（ACEI）
# 2=血管紧张素II受体拮抗剂（ARB）
# 3=α受体阻滞剂
# 4=β受体阻滞剂
# 5=利尿剂
# 6=沙库巴曲缬沙坦（ARNI）
# 7=钙拮抗剂
# 8=单片固定复方制剂
# 9=中药
# 10=其它
#"H0183"       "H0184"       "H0185"   "H0186"      "H0188"       "created_at"  "locked_time" "sync_at"
#"H0187"  (药物使用频率)      
length(unique(baseline_drug$patient_uid)) 
#[1] 637988
summary(as.factor(baseline_drug$H0182))

baseline_his_drug <- fread("~/data/基线病史用药.csv", header = TRUE, encoding = "UTF-8")
print(dim(baseline_his_drug)) 
#[1] 682680      8
print(names(baseline_his_drug))
#高血压用药史

follow <- fread("~/data/随访.csv", header = TRUE, encoding = "UTF-8")
print(dim(follow)) 
# [1] 4577024     105
print(names(follow))
# H0225 随访日期
# H0226 随访状态 1=存活;2=失访;3=死亡
# H0233 (电话sbp), H0235 (门诊sbp), H0234 (电话dbp), H0236 (门诊dbp)
# H0255 临床情况处理 1=继续原药物治疗; 2=调整治疗方案
sum(is.na(follow$H0233)&is.na(follow$H0235)) #[1] 163637
sum(is.na(follow$H0233) & is.na(follow$H0235)) / nrow(follow) #[1] 0.03575183
sum(is.na(follow$H0234)&is.na(follow$H0236)) #[1] 163700
sum(is.na(follow$H0234)&is.na(follow$H0236)) / nrow(follow) #[1] 0.0357656
length(unique(follow$followup_id)) 
# [1] 4577024
length(subset(follow,follow$patient_uid=="d5e32731-477a-238a-0c2f-9f993190e6f7")$followup_id)

follow_drug <- fread("~/data/随访用药.csv", header = TRUE, encoding = "UTF-8")
print(dim(follow_drug)) 
# [1] 18807899       14
print(names(follow_drug))
# H0257 药物类型
# 1=血管紧张素转化酶抑制剂（ACEI） 
# 2=血管紧张素II受体拮抗剂（ARB） 
# 3=α受体阻滞剂
# 4=β受体阻滞剂
# 5=利尿剂
# 6=沙库巴曲缬沙坦（ARNI）
# 7=钙拮抗剂
# 8=单片固定复方制剂
# 9=中药
# 10=其它
# H0262 (药物使用频率)


# 合并基线和随访 ----
## baseline ----
baseline_1 <- baseline %>%  #640418
  mutate(
    day = H0004,
    sbp = H0094,
    dbp = H0095,
    # secondary_htn ----
    secondary_htn = if_else(
      grepl("\\b[1-8]\\b", H0179), 
      "secondary_htn", 
      "primary_htn"
    )
  ) %>%
  filter(secondary_htn == "primary_htn") %>%
  dplyr::select(patient_uid, day, sbp,dbp) #629895

baseline_drug_1 <- baseline_drug %>% 
  mutate(
    drug = case_when(
      H0182 %in% c("ACEI", "ARB","沙库巴曲缬沙坦(ARNI)") ~ 1,
      H0182 == "β受体阻滞剂" ~ 2,
      H0182 == "钙拮抗剂" ~ 3,
      H0182 == "利尿剂" ~ 4,
      H0182 == "α受体阻滞剂" ~ 5,
      TRUE ~ 6 ),
    drug_time = H0187,
    drugnochange = 1
  )%>%
  dplyr::select(patient_uid, drug, drug_time,drugnochange)

### follow ----
follow_1 <- follow %>% #4577024
  mutate(
    day = H0225,
    sbp = pmax(H0233, H0235, na.rm = TRUE),
    dbp = pmax(H0234, H0236, na.rm = TRUE),
    drugnochange = H0255
  ) %>%
  dplyr::select(patient_uid, followup_id, day, drugnochange, sbp, dbp)

follow_drug_1 <- follow_drug %>% 
  mutate(
    drug = case_when(
      H0257 %in% c("ACEI", "ARB","沙库巴曲缬沙坦(ARNI)") ~ 1,
      H0257 == "β受体阻滞剂" ~ 2,
      H0257 == "钙拮抗剂" ~ 3,
      H0257 == "利尿剂" ~ 4,
      H0257 == "α受体阻滞剂" ~ 5,
      TRUE ~ 6 ),
    drug_time = H0262
  )%>%
  dplyr::select(patient_uid,followup_id, drug, drug_time)



# 所有单药用药者 -----
baseline_single <- baseline_drug_1 %>%
  filter(drug != 6) %>%
  group_by(patient_uid) %>%
  filter(n() == 1) %>%
  ungroup() #194828

follow_single <- follow_drug_1 %>%
  filter(drug != 6) %>%
  group_by(patient_uid,followup_id) %>%
  filter(n() == 1) %>%
  ungroup() #1262885

nomo_baseline <- left_join(baseline_single,baseline_1,by="patient_uid")
nomo_follow <- left_join(follow_single,follow_1,by=c("patient_uid","followup_id"))
nomo <- bind_rows(nomo_baseline,nomo_follow) #1457713
nomo <- nomo %>% 
  filter(is.na(day)==F) #1449643

setDT(nomo)
setorder(nomo, patient_uid, day)
nomo <- nomo[!is.na(sbp) & sbp >= 70 & sbp <= 250] #1440986

nomo[, rleid_drug := rleid(patient_uid, drug)]
segments <- nomo[, .(
  start_day = min(day),
  end_day = max(day),
  duration = as.numeric(difftime(max(day), min(day), units = "days")),
  n_sbp = .N  
), by = .(patient_uid, rleid_drug, drug)]
long_segments <- segments[duration >= 180 & n_sbp >= 3]

nomo_data <- nomo[long_segments, on = .(patient_uid, rleid_drug), nomatch = 0]
nrow(nomo_data) #1317004
length(unique(nomo_data$patient_uid))  #138568

save(nomo_data, file = "data/nomo_data.RData")

# sbp ttr -----
load("data/nomo_data.RData")

A <- B <- C <- D <- α <- data.frame()
for (i in unique(nomo_data$drug)) {
  subset_data <- nomo_data[nomo_data$drug == i, ]
  if (i == 1) {
    A <- rbind(A, subset_data)
  } else if (i == 2) {
    B <- rbind(B, subset_data)
  } else if (i == 3) {
    C <- rbind(C, subset_data)
  } else if (i == 4) {
    D <- rbind(D, subset_data)
  } else if (i == 5) {
    α <- rbind(α, subset_data)
  } 
}


length(unique(A$patient_uid)) #35357
length(unique(B$patient_uid)) #21679
length(unique(C$patient_uid)) #77214
length(unique(D$patient_uid)) #4229
length(unique(α$patient_uid)) #471


# 计算TTR/load/mean ----
my_packages <- c(
  "tidyverse",  
  "progressr",
  "pracma"
)
load_packages(my_packages)

data_frames <- list(A, B, C, D)
with_progress({
  p <- progressor(along = data_frames)
  
  data_frames <- map(data_frames, function(df) {
    p()  
    df %>%
      group_by(patient_uid) %>%
      summarise(
        day_num  = list(as.numeric(as.Date(day) - min(as.Date(day)))),
        sbp_vals = list(sbp),
        duration = first(duration),
        drug_time = first(drug_time),
        drug = first(drug),
        .groups = "drop"
      ) %>%
      rowwise() %>%
      mutate(
        x_seq = list(seq(min(unlist(day_num)), max(unlist(day_num)), length.out = max(duration * 2, 100))),
        SBP_interp = list(approx(unlist(day_num), unlist(sbp_vals), xout = unlist(x_seq))$y),
        
        days_sbpttr1 = sum(unlist(SBP_interp) < 130, na.rm = TRUE),
        sbpTTR1 = days_sbpttr1 / length(unlist(SBP_interp)),
        
        days_sbpttr2 = sum(unlist(SBP_interp) < 140, na.rm = TRUE),
        sbpTTR2 = days_sbpttr2 / length(unlist(SBP_interp)),
        
        meanSBP = mean(unlist(sbp_vals), na.rm = TRUE),
        
        total_sbp_area = trapz(unlist(x_seq), unlist(SBP_interp)),
        total_over_130_area = trapz(unlist(x_seq), pmax(unlist(SBP_interp) - 130, 0)),
        total_over_140_area = trapz(unlist(x_seq), pmax(unlist(SBP_interp) - 140, 0)),
        
        SBPload1 = ifelse(total_sbp_area > 0, total_over_130_area / total_sbp_area, NA_real_),
        SBPload2 = ifelse(total_sbp_area > 0, total_over_140_area / total_sbp_area, NA_real_),
        sbp_sd    = sd(unlist(sbp_vals), na.rm = TRUE),
        sbp_cv    = sd(unlist(sbp_vals), na.rm = TRUE) / mean(unlist(sbp_vals), na.rm = TRUE),
        sbp_arv   = mean(abs(diff(unlist(sbp_vals))), na.rm = TRUE)
      ) %>%
      ungroup()
  })
})


Ai <- data_frames[[1]]
Bi <- data_frames[[2]]
Ci <- data_frames[[3]]
Di <- data_frames[[4]]

ACEi <- distinct(subset(Ai,select = c("patient_uid","sbpTTR1","sbpTTR2","SBPload1","SBPload2","meanSBP","drug", "drug_time","sbp_sd","sbp_cv","sbp_arv")))
beta <- distinct(subset(Bi,select = c("patient_uid","sbpTTR1","sbpTTR2","SBPload1","SBPload2","meanSBP","drug", "drug_time","sbp_sd","sbp_cv","sbp_arv")))
CCB <- distinct(subset(Ci,select = c("patient_uid","sbpTTR1","sbpTTR2","SBPload1","SBPload2","meanSBP","drug", "drug_time","sbp_sd","sbp_cv","sbp_arv")))
diuretic <- distinct(subset(Di,select = c("patient_uid","sbpTTR1","sbpTTR2","SBPload1","SBPload2","meanSBP","drug", "drug_time","sbp_sd","sbp_cv","sbp_arv")))

diuretic <- diuretic %>%
  mutate(
    drug_time = as.character(drug_time),
    drug_time = case_when(
      drug_time %in% c("0", "20", "25") ~ NA_character_,  
      drug_time == "qd" ~ "1",                   
      TRUE ~ drug_time             
    ),
    drug_time = as.numeric(drug_time)
  )
CCB <- CCB %>%
  mutate(
    drug_time = str_trim(drug_time),           
    drug_time = case_when(
      drug_time %in% c( "5", "7","98","41","30","20", "间断","166","10","11","1146","12","13.54","109","14","14.38","0") ~ NA_character_, 
      drug_time %in% c("-", "一","qd", "1", "１","`1") ~ "1",                   
      drug_time %in% c( "二", "bid")  ~ "2",           
      TRUE ~ drug_time                                       
    ),
    drug_time = as.numeric(drug_time)         
  )

beta <- beta %>%
  mutate(
    drug_time = str_trim(drug_time),          
    drug_time = case_when(
      drug_time %in% c("-", "1`","qd") ~ "1",                   
      drug_time %in% c( "12.5", "23.75")  ~ "0.5", 
      drug_time  == "0" ~ NA_character_,  
      TRUE ~ drug_time                                       
    ),
    drug_time = as.numeric(drug_time)         
  )
ACEi <- ACEi %>%
  mutate(
    drug_time = str_trim(drug_time),          
    drug_time = case_when(
      drug_time %in% c("一","qd") ~ "1",                   
      drug_time == "bid"  ~ "2", 
      drug_time  %in% c("0","10","150","16","20","40","41","7","75","8","9","80","未服降压药") ~ NA_character_,  
      TRUE ~ drug_time                                       
    ),
    drug_time = as.numeric(drug_time)         
  )

save(ACEi, file = "data/ACEi.RData")
save(beta, file = "data/beta.RData")
save(CCB, file = "data/CCB.RData")
save(diuretic, file = "data/diuretic.RData")


# 合并基线特征 ----
load("data/ACEi.RData")
load("data/beta.RData")
load("data/CCB.RData")
load("data/diuretic.RData")
uids <- c(ACEi$patient_uid, CCB$patient_uid, beta$patient_uid,diuretic$patient_uid) #138479

B <- baseline %>% 
  filter(patient_uid %in% uids) %>% #138103
  mutate(
    # 基本人口学 ----
    sex   = recode(H0006, `1` = "Male", `2` = "Female"),
    age      = H0009,
    dwell    = case_when(
      H0017 %in% c("安徽", "安徽阜阳", "安徽省") ~ "安徽",
      H0017 %in% c("北京", "北京市") ~ "北京",
      H0017 %in% c("福建", "福建福州") ~ "福建",
      H0017 %in% c("甘肃", "甘肃省") ~ "甘肃",
      H0017 %in% c("广东") ~ "广东",
      H0017 %in% c("广西") ~ "广西",
      H0017 %in% c("贵州") ~ "贵州",
      H0017 %in% c("哈尔滨") ~ "黑龙江",
      H0017 %in% c("海南") ~ "海南",
      H0017 %in% c("河北", "河北省") ~ "河北",
      H0017 %in% c("河南", "河南省") ~ "河南",
      H0017 %in% c("湖北", "湖北省") ~ "湖北",
      H0017 %in% c("湖南", "湖南省") ~ "湖南",
      H0017 %in% c("吉林", "吉林省") ~ "吉林",
      H0017 %in% c("江苏", "江苏省") ~ "江苏",
      H0017 %in% c("江西", "江西省") ~ "江西",
      H0017 %in% c("辽宁", "辽宁省") ~ "辽宁",
      H0017 %in% c("内蒙古", "内蒙古自治区") ~ "内蒙古",
      H0017 %in% c("宁夏") ~ "宁夏",
      H0017 %in% c("青海", "青海省") ~ "青海",
      H0017 %in% c("山东", "山东省") ~ "山东",
      H0017 %in% c("山西", "山西省") ~ "山西",
      H0017 %in% c("陕西", "陕西省") ~ "陕西",
      H0017 %in% c("上海") ~ "上海",
      H0017 %in% c("四川", "四川省") ~ "四川",
      H0017 %in% c("台湾") ~ "台湾",
      H0017 %in% c("天津", "天津市") ~ "天津",
      H0017 %in% c("西藏") ~ "西藏",
      H0017 %in% c("香港") ~ "香港",
      H0017 %in% c("新疆", "新疆维吾尔自治区") ~ "新疆",
      TRUE ~ NA_character_
    ),
    dwell_c = case_when(
      dwell %in% c("北京", "天津", "河北", "山西", "内蒙古") ~ "华北",
      dwell %in% c("辽宁", "吉林", "黑龙江") ~ "东北",
      dwell %in% c("上海", "江苏", "浙江", "安徽", "福建", "江西", "山东") ~ "华东",
      dwell %in% c("河南", "湖北", "湖南") ~ "华中",
      dwell %in% c("广东", "广西", "海南") ~ "华南",
      dwell %in% c("重庆", "四川", "贵州", "云南", "西藏") ~ "西南",
      dwell %in% c("陕西", "甘肃", "青海", "宁夏", "新疆") ~ "西北",
      TRUE ~ NA_character_
    ),
    ethnicity = case_when(
      H0011 == "汉族" ~ "汉族",
      # 西南 + 西北
      H0011 %in% c("彝族","苗族","布依族","傣族","拉佤族","佤族","傈僳族","景颇族",
                   "羌族","普米族","怒族","纳西族","瑶族","白族",
                   "维吾尔族","哈萨克族","裕固族","回族","撒拉族","土族","东乡族") ~ "西南/西北民族",
      # 其他少数民族（北方、藏区、其他、外国血统）
      H0011 %in% c("蒙古族","达斡尔族","鄂温克族","满族","朝鲜族",
                   "藏族","门巴族","珞巴族",
                   "仡佬族","布朗族","保安族","阿昌族","基诺族","锡伯族",
                   "外国血统") ~ "其他少数民族",
      TRUE ~ NA_character_
    ),
    work = case_when(
      H0013 %in% c(1,2) ~ "High SES",
      H0013 %in% c(3,4,5) ~ "Medium/Low SES",
      H0013 %in% c(6,7,8,9) ~ "Non-working",
      H0013 == 10 ~ "Other",
      TRUE ~ NA_character_
    ),
    marriage = case_when(
      H0016 == 2 ~ "Married",
      H0016 %in% c(1,3,4) ~ "Not married",
      TRUE ~ NA_character_
    ),
    education = case_when(
      H0015 %in% 1:3 ~ "Below High School",
      H0015 %in% 4:6 ~ "High School or Above",
      TRUE ~ NA_character_
    ),
    
    # 饮酒与吸烟 ----
    alco = H0025,
    smk  = case_when(
      H0033 == 1 ~ 1,
      H0033 == 0 & H0034 == 2 ~ 2,
      H0033 == 0 & H0034 == 1 ~ 0,
      is.na(H0034) ~ as.numeric(H0033)
    ),
    
    # 运动 ----
    sport = H0038,
    still = H0037,
    dNa    = H0043,
    mental = case_when(
      H0044 %in% 1:3 ~"Anxiety or Depression",
      H0044 == 4 ~ "Good Mental Health",
      TRUE ~ NA_character_
    ),
    
    # 身体指标 ----
    bmi = if_else(H0105 >= 10 & H0105 <= 60, H0105, NA_real_),
    HR_1 = if_else(H0107 >= 30 & H0107 <= 220, H0107, NA_integer_),
    HR_ecg = if_else(H0141 >= 30 & H0141 <= 220, H0141, NA_integer_),
    HR = round(rowMeans(cbind(HR_1, HR_ecg), na.rm = TRUE)),
    
    # 合并症 ----
    Fhis_hype   = case_when(
      H0045==2 ~ 1,
      H0045==1 ~ 0,
      TRUE ~ NA_real_
    ),
    hype_duration = H0050,
    hype_duration_5 = case_when(
      H0050 >=5 ~ 1,
      TRUE ~ 0
    ),
    hype_duration_10 = case_when(
      H0050 >=10 ~ 1,
      TRUE ~ 0
    ),
    AF = case_when(
      (H0142 == 4) | (H0194 ==1) | (H0072 ==1) ~ 1,
      is.na(H0142) & is.na(H0194) & is.na(H0072)  ~ NA_real_,
      TRUE ~ 0
    ),
    
    ## 糖尿病及血糖 ----
    anti_glu  = H0207,
    dm_report = case_when(
      H0063 == 3 ~ "preDM",
      H0062 == 2 | H0063 %in% c(1, 2) | (H0207 == 1) ~ "DM",
      H0062 == 1 ~ "no",
      is.na(H0062) & is.na(H0063) &is.na(H0207) ~ NA_character_  
    ),
    fpg = case_when(
      H0126 == 1 ~ H0125 * 100,      
      H0126 == 2 ~ H0125 * 18,       
      TRUE ~ NA_real_
    ),
    fpg = if_else(fpg >= 40 & fpg <= 500, fpg, NA_real_),   
    HbA1c = if_else(H0128 >= 3 & H0128 <=30, H0128, NA_real_),
    dm_2 = case_when(
      (anti_glu == 1) | (HbA1c >= 6.5) | (fpg >= 126) ~ "DM",              
      ((HbA1c >= 5.7 & HbA1c < 6.5) | (fpg >= 100 & fpg <= 125)) ~ "preDM", 
      (is.na(anti_glu) & is.na(HbA1c) & is.na(fpg)) ~ NA_character_,        
      TRUE ~ "no" 
    ), 
    DM = case_when(
      is.na(dm_2) & is.na(dm_report) ~ NA_character_,   
      dm_report == "DM" | dm_2 =="DM" ~ "DM",
      dm_report == "preDM" | dm_2 == "preDM" ~ "preDM",
      dm_report == "no" & dm_2 =="no" ~ "no"
    ),
    
    ## 肾脏指标 ----
    Scr = if_else(H0110 >= 10 & H0110 <= 1500, H0110, NA_real_),
    uacr = case_when(
      H0131 == 1 ~ H0130,                 
      H0131 == 2 ~ H0130 * 8.84,        
      H0131 == 3 ~ H0130 * 8841,         
      TRUE ~ NA_real_
    ),
    uacr = if_else(uacr >= 0 & uacr <= 3000, uacr, NA_real_),
    
    ckdepi = case_when(
      sex == "Female" & Scr / 88.4 <= 0.7 ~ 142 * (Scr / 88.4 / 0.7)^-0.241 * 0.9938^age,
      sex == "Female" & Scr / 88.4 > 0.7 ~ 142 * (Scr / 88.4 / 0.7)^-1.200 * 0.9938^age,
      sex == "Male" & Scr / 88.4 <= 0.9 ~ 142 * (Scr / 88.4 / 0.9)^-0.302 * 0.9938^age,
      sex == "Male" & Scr / 88.4 > 0.9 ~ 142 * (Scr / 88.4 / 0.9)^-1.200 * 0.9938^age,
      TRUE ~ NA_real_
    ),
    
    gfr_category = as.character(case_when(
      ckdepi >= 90 | is.na(ckdepi) ~ "G1",
      ckdepi >= 60 & ckdepi < 90 ~ "G2",
      ckdepi >= 45 & ckdepi < 60 ~ "G3a",
      ckdepi >= 30 & ckdepi < 45 ~ "G3b",
      ckdepi >= 15 & ckdepi < 30 ~ "G4",
      ckdepi < 15 ~ "G5"
    )),
    
    alb_category = as.character(case_when(
      uacr < 30 | is.na(uacr) ~ "A1",
      uacr >= 30 & uacr <= 299 ~ "A2",
      uacr >= 300 ~ "A3"
    )),
    
    KDIGO = case_when(
      gfr_category %in% c("G1","G2") & alb_category == "A1" ~ 0L,
      (gfr_category == "G1" & alb_category == "A2") |
        (gfr_category == "G2" & alb_category == "A2") |
        (gfr_category == "G3a" & alb_category == "A1") ~ 1L,
      (gfr_category %in% c("G1","G2") & alb_category == "A3") |
        (gfr_category == "G3b" & alb_category == "A1") |
        (gfr_category == "G3a" & alb_category == "A2") ~ 2L,
      (gfr_category == "G3a" & alb_category == "A3") |
        (gfr_category == "G3b" & alb_category %in% c("A2","A3")) |
        (gfr_category %in% c("G4","G5")) ~ 3L,
      TRUE ~ NA_integer_
    ),
    
    CKD_risk = factor(case_when(
      KDIGO == 0L ~ "low risk",
      KDIGO == 1L ~ "moderately increased",
      KDIGO == 2L ~ "high risk",
      KDIGO == 3L ~ "very high risk",
      TRUE ~ NA_character_
    ), levels = c("low risk", "moderately increased", "high risk", "very high risk")),
    
    # 血清学 ----
    b_K  = if_else(H0108 >= 2 & H0108 <= 8, H0108, NA_real_),
    b_Na = if_else(H0109 >= 100 & H0109 <= 200, H0109, NA_real_),
    hcy = if_else(H0115 >= 2 & H0115 <= 100, H0115, NA_real_),
    ua =  (if_else(H0112 >= 60 & H0112 <= 1200, H0112, NA_real_)) / 59.48, 
    hua = case_when(
      sex == "Male" & ua >= 420 / 59.48 ~ 1,     
      sex == "Female" & ua >= 360 / 59.48 ~ 1,     
      H0065 == 2 ~ 1,                  
      is.na(ua) & (is.na(H0065) | H0065 ==3) ~ NA_real_, 
      TRUE ~ 0                        
    ),
    
    ## 血脂 ----
    tc   = case_when(
      H0118 == 1 ~ H0117,               
      H0118 == 2 ~ H0117 * 38.67,       
      TRUE ~ NA_real_
    ),
    ldlc = case_when(
      H0122 == 1 ~ H0121,               
      H0122 == 2 ~ H0121 * 38.67,       
      TRUE ~ NA_real_
    ),
    hdlc = case_when(
      H0124 == 1 ~ H0123,               
      H0124 == 2 ~ H0123 * 38.67,       
      TRUE ~ NA_real_
    ),
    tg   = case_when(
      H0120 == 1 ~ H0119,               
      H0120 == 2 ~ H0119 * 88.6,       
      TRUE ~ NA_real_
    ),
    tc = if_else(!is.na(tc) & tc >= 50 & tc <= 500, tc, NA_real_),
    ldlc = if_else(!is.na(ldlc) & ldlc >= 10 & ldlc <= 400, ldlc, NA_real_),
    hdlc = if_else(!is.na(hdlc) & hdlc >= 5 & hdlc <= 150, hdlc, NA_real_),
    tg = if_else(!is.na(tg) & tg >= 10 & tg <= 1000, tg, NA_real_),
    
    anti_chol = H0196,
    
    dyslipidemia = case_when(
      is.na(tc) & is.na(ldlc)  ~ NA_real_,
      tc >= 240 | ldlc >= 160   ~ 1,
      TRUE ~ 0
    ),
    
    # cvd ----
    cvd = case_when(
      H0068 == 1 | (H0169 == 3 & H0170 %in% c(2, 3)) ~ 1,
      is.na(H0068) & (is.na(H0169) | is.na(H0170) )~ NA_real_,
      TRUE ~ 0
    ),
    # HF ----
    hf = as.character(H0074),
    # stroke ----
    stroke = as.character(H0066),
    stroke_w_o_TIA = as.character(case_when(
      (H0066 ==1 & H0067 != 1) | H0067 %in% c(2, 3) ~ 1,
      is.na(H0066) & is.na(H0067)~ NA_real_,
      TRUE ~ 0
    ))) %>%  
  # 选择最终变量
  dplyr::select(patient_uid, sex, age, dwell, dwell_c, ethnicity, work, marriage, education,
                alco, smk, sport, still, dNa, mental,
                Fhis_hype, hype_duration_5,hype_duration_10,
                bmi, HR_1, HR_ecg,HR,
                AF,cvd,hf,stroke,
                anti_glu, dm_report, fpg, HbA1c, dm_2,DM,
                Scr, uacr, ckdepi, CKD_risk,
                b_K, b_Na, hcy, ua, hua,
                tc, ldlc, hdlc, tg, anti_chol, dyslipidemia)
nrow(B)
save(B, file = "data/B.RData")


load( "data/B.RData")

ACEi_1 <- left_join(ACEi,B,by="patient_uid")
beta_1 <- left_join(beta,B,by="patient_uid")
CCB_1 <- left_join(CCB,B,by="patient_uid")
diuretic_1 <- left_join(diuretic,B,by="patient_uid")


save(ACEi_1, file = "data/ACEi_1.RData")
save(beta_1, file = "data/beta_1.RData")
save(CCB_1, file = "data/CCB_1.RData")
save(diuretic_1, file = "data/diuretic_1.RData")
