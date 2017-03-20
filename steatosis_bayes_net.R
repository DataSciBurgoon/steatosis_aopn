library(gRain)

#Building the Bayes Net
yn <- c("yes", "no")
nrf2 <- cptable(~ nrf2, values=c(0.5, 0.5), levels=yn)
fxr <- cptable(~ fxr | nrf2, values=c(0.95, 0.05, 0.05, 0.95), levels=yn)
shp <- cptable(~ shp | fxr, values=c(0.95, 0.05, 0.05, 0.95), levels=yn)
lxr <- cptable(~ lxr | shp, values=c(0.05, 0.95, 0.95, 0.05), levels=yn)
ppar_alpha <- cptable(~ ppara | fxr+shp+lxr, 
                      values=c(0.50, 0.95, 0.50, 0.50, 0.50, 0.50, 0.05, 0.50,
                               0.50, 0.05, 0.50, 0.50, 0.50, 0.50, 0.95, 0.50),
                      levels=yn)
hsd17b4 <- cptable(~ hsd17b4 | ppara, values=c(0.95, 0.05, 0.05, 0.95), levels=yn)
fatty_acid_beta_oxidation <- cptable(~ fatty_acid_beta_oxidation | hsd17b4,
                                     values = c(1, 0, 0, 1), levels=yn)
steatosis <-  cptable(~steatosis | fatty_acid_beta_oxidation, 
                      values=c(0, 1, 1, 0), levels=yn)
lrh1 <- cptable(~lrh1 | shp,
                values=c(0.05, 0.95, 0.95, 0.05), levels=yn)
pi3k <- cptable(~pi3k, values=c(0.50, 0.50), levels=yn)
mtorc2 <- cptable(~mtorc2 | ir, values=c(1, 0, 0, 1), levels=yn)
ir <- cptable(~ir, values=c(0.50, 0.50), levels=yn)
akt <- cptable(~akt | pi3k+mtorc2, values=c(0.95, 0.95, 0.95, 0.05,
                                            0.05, 0.05, 0.05, 0.95),
               levels=yn)
lfabp <- cptable(~lfabp | akt+pi3k, values=c(0.95, 0.95, 0.95, 0.05,
                                             0.05, 0.05, 0.05, 0.95),
                 levels=yn)
pparg <- cptable(~pparg | lfabp, values=c(0.95, 0.05, 0.05, 0.95), levels=yn)
fas <- cptable(~fas | lrh1+lxr+pparg, values=c(0.95, 0.75, 0.75, 0.50, 0.75, 0.50, 0.50, 0.05,
                                         0.05, 0.25, 0.25, 0.50, 0.25, 0.50, 0.50, 0.95),
               levels=yn)
mtorc1 <- cptable(~mtorc1 | akt, values=c(0.95, 0.05, 0.05, 0.95), levels=yn)
apkc <- cptable(~apkc | pi3k, values=c(1, 0, 0, 1), levels=yn)
srebp1 <- cptable(~srebp1 | mtorc1 + apkc, values=c(0.95, 0.95, 0.95, 0.05,
                                                    0.05, 0.05, 0.05, 0.95),
                  levels=yn)
scd1 <- cptable(~scd1 | srebp1, values=c(1, 0, 0, 1), levels=yn)
lipogenesis <- cptable(~lipogenesis | scd1 + fas, values=c(0.95, 0.95, 0.66, 0,
                                                           0.05, 0.05, 0.33, 1),
                       levels=yn)
cytosolic_fatty_acids <- cptable(~cytosolic_fatty_acids | lipogenesis + lfabp + fatty_acid_beta_oxidation, 
                                 values=c(0.01, 1, 0.01, 0.95, 0.01, 0.66, 0.01, 0.50,
                                          0.99, 0, 0.99, 0.05, 0.99, 0.34, 0.99, 0.50),
                                 levels=yn)






plist <- compileCPT(list(nrf2, fxr, shp, lxr, ppar_alpha, hsd17b4, 
                         fatty_acid_beta_oxidation, steatosis, lrh1, pi3k, 
                         mtorc2, ir, akt, lfabp, pparg, fas, mtorc1, apkc, srebp1, scd1, lipogenesis,
                         cytosolic_fatty_acids))
bn_test <- grain(plist)
plot(bn_test)
querygrain(setEvidence(bn_test, evidence=list(hsd17b4="no")) )



a1 <- cptable(~accuracy_atg_era_trans_up, values=c(posterior_atg_era_trans_up_assay_pos_er_pos + 
                                                     posterior_atg_era_trans_up_assay_neg_er_neg, 
                                                   1 - (posterior_atg_era_trans_up_assay_pos_er_pos + 
                                                          posterior_atg_era_trans_up_assay_neg_er_neg)), levels=yn)
e1 <- cptable(~atg_era_trans_up|accuracy_atg_era_trans_up:chem_activates_er, 
              values=c(posterior_atg_era_trans_up_assay_pos_er_pos,
                       posterior_atg_era_trans_up_assay_neg_er_pos,
                       0.5,
                       0.5,
                       posterior_atg_era_trans_up_assay_pos_er_neg,
                       posterior_atg_era_trans_up_assay_neg_er_neg,
                       0.5,
                       0.5),
              levels=yn)
a2 <- cptable(~accuracy_nvs_nr_her, values=c(posterior_nvs_nr_her_assay_pos_er_pos +
                                               posterior_nvs_nr_her_assay_neg_er_neg,
                                             1 - (posterior_nvs_nr_her_assay_pos_er_pos +
                                                    posterior_nvs_nr_her_assay_neg_er_neg)), levels=yn)
e2 <- cptable(~nvs_nr_her|accuracy_nvs_nr_her:chem_activates_er,
              values = c(posterior_nvs_nr_her_assay_pos_er_pos,
                         posterior_nvs_nr_her_assay_neg_er_pos,
                         0.5,
                         0.5,
                         posterior_nvs_nr_her_assay_pos_er_neg,
                         posterior_nvs_nr_her_assay_neg_er_neg,
                         0.5,
                         0.5),
              levels=yn)
a3 <- cptable(~accuracy_ot_er_eraera_0480, values=c(posterior_OT_ER_ERaERa_0480_assay_pos_er_pos +
                                                      posterior_OT_ER_ERaERa_0480_assay_neg_er_neg,
                                                    1 - (posterior_OT_ER_ERaERa_0480_assay_pos_er_pos +
                                                           posterior_OT_ER_ERaERa_0480_assay_neg_er_neg)), levels=yn)
e3 <- cptable(~ot_er_eraera_0480|accuracy_ot_er_eraera_0480:chem_activates_er,
              values = c(posterior_OT_ER_ERaERa_0480_assay_pos_er_pos,
                         posterior_OT_ER_ERaERa_0480_assay_neg_er_pos,
                         0.5,
                         0.5,
                         posterior_OT_ER_ERaERa_0480_assay_pos_er_neg,
                         posterior_OT_ER_ERaERa_0480_assay_neg_er_neg,
                         0.5,
                         0.5),
              levels=yn)
a4 <- cptable(~accuracy_ot_er_eraera_1440, values=c(posterior_OT_ER_ERaERa_1440_assay_pos_er_pos +
                                                      posterior_OT_ER_ERaERa_1440_assay_neg_er_neg,
                                                    1 - (posterior_OT_ER_ERaERa_1440_assay_pos_er_pos +
                                                           posterior_OT_ER_ERaERa_1440_assay_neg_er_neg)), levels=yn)
d4 <- cptable(~ot_er_eraera_dep, values=c(0.5, 0.5), levels=yn)
#e4 <- cptable(~ot_er_eraera_1440|accuracy_ot_er_eraera_1440:ot_er_eraera_dep:ot_er_eraera_0480:chem_activates_er,
#              )

e4 <- cptable(~ot_er_eraera_1440|accuracy_ot_er_eraera_1440:chem_activates_er,
              values = c(posterior_OT_ER_ERaERa_1440_assay_pos_er_pos,
                         posterior_OT_ER_ERaERa_1440_assay_neg_er_pos,
                         0.5,
                         0.5,
                         posterior_OT_ER_ERaERa_1440_assay_pos_er_neg,
                         posterior_OT_ER_ERaERa_1440_assay_neg_er_neg,
                         0.5,
                         0.5),
              levels=yn)
a5 <- cptable(~accuracy_ot_er_eraerb_0480, values=c(posterior_OT_ER_ERaERb_0480_assay_pos_er_pos +
                                                      posterior_OT_ER_ERaERb_0480_assay_neg_er_neg,
                                                    1 - (posterior_OT_ER_ERaERb_0480_assay_pos_er_pos +
                                                           posterior_OT_ER_ERaERb_0480_assay_neg_er_neg)), levels=yn)
e5 <- cptable(~ot_er_eraerb_0480|accuracy_ot_er_eraerb_0480:chem_activates_er,
              values = c(posterior_OT_ER_ERaERb_0480_assay_pos_er_pos,
                         posterior_OT_ER_ERaERb_0480_assay_neg_er_pos,
                         0.5,
                         0.5,
                         posterior_OT_ER_ERaERb_0480_assay_pos_er_neg,
                         posterior_OT_ER_ERaERb_0480_assay_neg_er_neg,
                         0.5,
                         0.5),
              levels=yn)
a6 <- cptable(~accuracy_ot_er_eraerb_1440, values=c(posterior_OT_ER_ERaERb_1440_assay_pos_er_pos +
                                                      posterior_OT_ER_ERaERb_1440_assay_neg_er_neg,
                                                    1 - (posterior_OT_ER_ERaERb_1440_assay_pos_er_pos +
                                                           posterior_OT_ER_ERaERb_1440_assay_neg_er_neg)), levels=yn)
e6 <- cptable(~ot_er_eraerb_1440|accuracy_ot_er_eraerb_1440:chem_activates_er,
              values = c(posterior_OT_ER_ERaERb_1440_assay_pos_er_pos,
                         posterior_OT_ER_ERaERb_1440_assay_neg_er_pos,
                         0.5,
                         0.5,
                         posterior_OT_ER_ERaERb_1440_assay_pos_er_neg,
                         posterior_OT_ER_ERaERb_1440_assay_neg_er_neg,
                         0.5,
                         0.5),
              levels=yn)
a7 <- cptable(~accuracy_ot_er_erberb_0480, values=c(posterior_OT_ER_ERbERb_0480_assay_pos_er_pos +
                                                      posterior_OT_ER_ERbERb_0480_assay_neg_er_neg,
                                                    1-(posterior_OT_ER_ERbERb_0480_assay_pos_er_pos +
                                                         posterior_OT_ER_ERbERb_0480_assay_neg_er_neg)), levels=yn)
e7 <- cptable(~ot_er_erberb_0480|accuracy_ot_er_erberb_0480:chem_activates_er,
              values = c(posterior_OT_ER_ERbERb_0480_assay_pos_er_pos,
                         posterior_OT_ER_ERbERb_0480_assay_neg_er_pos,
                         0.5,
                         0.5,
                         posterior_OT_ER_ERbERb_0480_assay_pos_er_neg,
                         posterior_OT_ER_ERbERb_0480_assay_neg_er_neg,
                         0.5,
                         0.5),
              levels=yn)
a8 <- cptable(~accuracy_ot_er_erberb_1440, values=c(posterior_OT_ER_ERbERb_1440_assay_pos_er_pos +
                                                      posterior_OT_ER_ERbERb_1440_assay_neg_er_neg,
                                                    1-(posterior_OT_ER_ERbERb_1440_assay_pos_er_pos +
                                                         posterior_OT_ER_ERbERb_1440_assay_neg_er_neg)), levels=yn)
e8 <- cptable(~ot_er_erberb_1440|accuracy_ot_er_erberb_1440:chem_activates_er,
              values = c(posterior_OT_ER_ERbERb_1440_assay_pos_er_pos,
                         posterior_OT_ER_ERbERb_1440_assay_neg_er_pos,
                         0.5,
                         0.5,
                         posterior_OT_ER_ERbERb_1440_assay_pos_er_neg,
                         posterior_OT_ER_ERbERb_1440_assay_neg_er_neg,
                         0.5,
                         0.5),
              levels=yn)
a9 <- cptable(~accuracy_tox21_era_bla_agonist_ratio, values=c(posterior_TOX21_ERa_BLA_Agonist_ratio_assay_pos_er_pos +
                                                                posterior_TOX21_ERa_BLA_Agonist_ratio_assay_neg_er_neg,
                                                              1-(posterior_TOX21_ERa_BLA_Agonist_ratio_assay_pos_er_pos +
                                                                   posterior_TOX21_ERa_BLA_Agonist_ratio_assay_neg_er_neg)),
              levels=yn)
e9 <- cptable(~tox21_era_bla_agonist_ratio|accuracy_tox21_era_bla_agonist_ratio:chem_activates_er,
              values=c(posterior_TOX21_ERa_BLA_Agonist_ratio_assay_pos_er_pos,
                       posterior_TOX21_ERa_BLA_Agonist_ratio_assay_neg_er_pos,
                       0.5,
                       0.5,
                       posterior_TOX21_ERa_BLA_Agonist_ratio_assay_pos_er_neg,
                       posterior_TOX21_ERa_BLA_Agonist_ratio_assay_neg_er_neg,
                       0.5,
                       0.5),
              levels=yn)
a10 <- cptable(~accuracy_atg_era_trans_dn, values=c(posterior_ATG_ERa_TRANS_dn_assay_pos_er_pos +
                                                      posterior_ATG_ERa_TRANS_dn_assay_neg_er_neg,
                                                    1-(posterior_ATG_ERa_TRANS_dn_assay_pos_er_pos +
                                                         posterior_ATG_ERa_TRANS_dn_assay_neg_er_neg)),
               levels=yn)
e10 <- cptable(~atg_era_trans_dn|accuracy_atg_era_trans_dn:chem_activates_er,
               values = c(posterior_ATG_ERa_TRANS_dn_assay_pos_er_pos,
                          posterior_ATG_ERa_TRANS_dn_assay_neg_er_pos,
                          0.5,
                          0.5,
                          posterior_ATG_ERa_TRANS_dn_assay_pos_er_neg,
                          posterior_ATG_ERa_TRANS_dn_assay_neg_er_neg,
                          0.5,
                          0.5),
               levels=yn)
a11 <- cptable(~accuracy_atg_ere_cis_dn, values=c(posterior_ATG_ERE_CIS_dn_assay_pos_er_pos +
                                                    posterior_ATG_ERE_CIS_dn_assay_neg_er_neg,
                                                  1-(posterior_ATG_ERE_CIS_dn_assay_pos_er_pos +
                                                       posterior_ATG_ERE_CIS_dn_assay_neg_er_neg)),
               levels=yn)
e11 <- cptable(~atg_ere_cis_dn|accuracy_atg_ere_cis_dn:chem_activates_er,
               values = c(posterior_ATG_ERE_CIS_dn_assay_pos_er_pos,
                          posterior_ATG_ERE_CIS_dn_assay_neg_er_pos,
                          0.5,
                          0.5,
                          posterior_ATG_ERE_CIS_dn_assay_pos_er_neg,
                          posterior_ATG_ERE_CIS_dn_assay_neg_er_neg,
                          0.5,
                          0.5),
               levels=yn)
h2 <- cptable(~chem_estrogenic|chem_activates_er,
              values = c(0.9, 0.1, 0.01, 0.99),
              levels=yn)

plist <- compileCPT(list(h1, a1, e1, a2, e2, a3, e3, a4, e4, a5, e5, a6, e6, a7, e7,
                         a8, e8, a9, e9, a10, e10, a11, e11, h2))
estrogenic_bn <- grain(plist)

querygrain(setEvidence(estrogenic_bn, evidence=list(chem_activates_er="yes")) )

estrogenic_bn_run <- function(ev_list, estrogenic_bn){
  evidence_list <- list(atg_era_trans_up=ev_list$atg_era_trans_up, nvs_nr_her=ev_list$nvs_nr_her, 
                        ot_er_eraera_0480=ev_list$ot_er_eraera_0480, ot_er_eraera_1440=ev_list$ot_er_eraera_1440,
                        ot_er_eraerb_0480=ev_list$ot_er_eraerb_0480, ot_er_eraerb_1440=ev_list$ot_er_eraerb_1440, 
                        ot_er_erberb_0480=ev_list$ot_er_erberb_0480, ot_er_erberb_1440=ev_list$ot_er_erberb_1440, 
                        tox21_era_bla_agonist_ratio=ev_list$tox21_era_bla_agonist_ratio, 
                        atg_era_trans_dn=ev_list$atg_era_trans_dn, atg_ere_cis_dn=ev_list$atg_ere_cis_dn)
  grain_results <- querygrain(setEvidence(estrogenic_bn, evidence=evidence_list))
  return(grain_results)
}

convert_nums_to_yn <- function(x){
  if(length(x)==1 && x == 1){
    return("yes")
  }
  else{
    return("no")
  }
}

unique_chems_er_pos <- unique(tc_gt_er_pos$casn)
unique_chems_er_neg <- unique(tc_gt_er_neg$casn)

tc_gt_er_pos_list <- list()
tc_gt_er_neg_list <- list()

for(i in 1:length(unique_chems_er_pos)){
  tc_gt_er_pos_chem_sub <- subset(tc_gt_er_pos, casn == unique_chems_er_pos[i])
  chem_assay_list <- list()
  chem_assay_list$atg_era_trans_up <- convert_nums_to_yn(subset(tc_gt_er_pos_chem_sub, aenm == "ATG_ERa_TRANS_up")$hitc)
  chem_assay_list$nvs_nr_her <- convert_nums_to_yn(subset(tc_gt_er_pos_chem_sub, aenm == "NVS_NR_hER")$hitc)
  chem_assay_list$ot_er_eraera_0480 <- convert_nums_to_yn(subset(tc_gt_er_pos_chem_sub, aenm == "OT_ER_ERaERa_0480")$hitc)
  chem_assay_list$ot_er_eraera_1440 <- convert_nums_to_yn(subset(tc_gt_er_pos_chem_sub, aenm == "OT_ER_ERaERa_1440")$hitc)
  chem_assay_list$ot_er_eraerb_0480 <- convert_nums_to_yn(subset(tc_gt_er_pos_chem_sub, aenm == "OT_ER_ERaERb_0480")$hitc)
  chem_assay_list$ot_er_eraerb_1440 <- convert_nums_to_yn(subset(tc_gt_er_pos_chem_sub, aenm == "OT_ER_ERaERb_1440")$hitc)
  chem_assay_list$ot_er_erberb_0480 <- convert_nums_to_yn(subset(tc_gt_er_pos_chem_sub, aenm == "OT_ER_ERbERb_0480")$hitc)
  chem_assay_list$ot_er_erberb_1440 <- convert_nums_to_yn(subset(tc_gt_er_pos_chem_sub, aenm == "OT_ER_ERbERb_1440")$hitc)
  chem_assay_list$tox21_era_bla_agonist_ratio <- convert_nums_to_yn(subset(tc_gt_er_pos_chem_sub, aenm == "TOX21_ERa_BLA_Agonist_ratio")$hitc)
  chem_assay_list$atg_era_trans_dn <- convert_nums_to_yn(subset(tc_gt_er_pos_chem_sub, aenm == "ATG_ERa_TRANS_dn")$hitc)
  chem_assay_list$atg_ere_cis_dn <- convert_nums_to_yn(subset(tc_gt_er_pos_chem_sub, aenm == "ATG_ERE_CIS_dn")$hitc)
  print(chem_assay_list)
  tc_gt_er_pos_list[[unique_chems_er_pos[i]]] <- chem_assay_list
}

toxcast_er_pos_grain_results <- lapply(tc_gt_er_pos_list, estrogenic_bn_run, estrogenic_bn)
toxcast_er_pos_grain_results_probs <- NULL

for(i in 1:length(toxcast_er_pos_grain_results)){
  toxcast_er_pos_grain_results_probs[i] <- toxcast_er_pos_grain_results[[i]]$chem_estrogenic[[1]]
}

for(i in 1:length(unique_chems_er_neg)){
  tc_gt_er_neg_chem_sub <- subset(tc_gt_er_pos, casn == unique_chems_er_pos[i])
  chem_assay_list <- list()
  chem_assay_list$atg_era_trans_up <- convert_nums_to_yn(subset(tc_gt_er_neg_chem_sub, aenm == "ATG_ERa_TRANS_up")$hitc)
  chem_assay_list$nvs_nr_her <- convert_nums_to_yn(subset(tc_gt_er_neg_chem_sub, aenm == "NVS_NR_hER")$hitc)
  chem_assay_list$ot_er_eraera_0480 <- convert_nums_to_yn(subset(tc_gt_er_neg_chem_sub, aenm == "OT_ER_ERaERa_0480")$hitc)
  chem_assay_list$ot_er_eraera_1440 <- convert_nums_to_yn(subset(tc_gt_er_neg_chem_sub, aenm == "OT_ER_ERaERa_1440")$hitc)
  chem_assay_list$ot_er_eraerb_0480 <- convert_nums_to_yn(subset(tc_gt_er_neg_chem_sub, aenm == "OT_ER_ERaERb_0480")$hitc)
  chem_assay_list$ot_er_eraerb_1440 <- convert_nums_to_yn(subset(tc_gt_er_neg_chem_sub, aenm == "OT_ER_ERaERb_1440")$hitc)
  chem_assay_list$ot_er_erberb_0480 <- convert_nums_to_yn(subset(tc_gt_er_neg_chem_sub, aenm == "OT_ER_ERbERb_0480")$hitc)
  chem_assay_list$ot_er_erberb_1440 <- convert_nums_to_yn(subset(tc_gt_er_neg_chem_sub, aenm == "OT_ER_ERbERb_1440")$hitc)
  chem_assay_list$tox21_era_bla_agonist_ratio <- convert_nums_to_yn(subset(tc_gt_er_neg_chem_sub, aenm == "TOX21_ERa_BLA_Agonist_ratio")$hitc)
  chem_assay_list$atg_era_trans_dn <- convert_nums_to_yn(subset(tc_gt_er_neg_chem_sub, aenm == "ATG_ERa_TRANS_dn")$hitc)
  chem_assay_list$atg_ere_cis_dn <- convert_nums_to_yn(subset(tc_gt_er_neg_chem_sub, aenm == "ATG_ERE_CIS_dn")$hitc)
  print(chem_assay_list)
  tc_gt_er_neg_list[[unique_chems_er_neg[i]]] <- chem_assay_list
}

toxcast_er_neg_grain_results <- lapply(tc_gt_er_neg_list, estrogenic_bn_run, estrogenic_bn)
toxcast_er_neg_grain_results_probs <- NULL

for(i in 1:length(toxcast_er_neg_grain_results)){
  toxcast_er_neg_grain_results_probs[i] <- toxcast_er_neg_grain_results[[i]]$chem_estrogenic[[1]]
}

true_positives <- length(which(toxcast_er_pos_grain_results_probs > 0.50)) / length(toxcast_er_pos_grain_results)
true_negatives <- length(which(toxcast_er_neg_grain_results_probs < 0.50)) / length(toxcast_er_neg_grain_results)
false_positives <- length(which(toxcast_er_neg_grain_results_probs > 0.50)) / length(toxcast_er_neg_grain_results)
false_negatives <- length(which(toxcast_er_pos_grain_results_probs < 0.50)) / length(toxcast_er_pos_grain_results)

toxcast_preds <- c(toxcast_er_pos_grain_results_probs, toxcast_er_neg_grain_results_probs)
toxcast_labs <- c(rep(1, length(toxcast_er_pos_grain_results_probs)), 
                  rep(0, length(toxcast_er_neg_grain_results_probs)))
er_tc_preds <- prediction(toxcast_preds, toxcast_labs)
er_tc_perf <- performance(er_tc_preds, "tpr", "fpr")
plot(er_tc_perf, colorize=TRUE)

plot(performance(er_tc_preds, "sens", "fpr"), colorize=TRUE)
abline(0,1)

# #Clearly this isn't doing very well with all these assays. We need to do better.
# #So, I'm going to remove some of the lower accuracy assays and see if we get
# #better results.
# querygrain(estrogenic_bn)
# 
# #We can see that the accuracy on the following 2 assays is not much better than
# #guessing or 50%:
# # 1) ATG_ERa_TRANS_dn
# # 2) ATG_ERE_CIS_dn
# #
# #In addition, some of the assays, are accurate, but not very helpful. The results
# #that they give tend to be close to 50%:
# # 1) OT_ER_ERaERa_1440
# # 2) OT_ER_ERaERb_0480
# # 3) OT_ER_ERaERb_1440
# # 4) OT_ER_ERbERb_0480
# # 5) OT_ER_ERbERb_1440
# # 6) TOX21_ERa_BLA_Agonist_ratio
# 
# 
# #Building the Streamlined Bayes Net
# yn <- c("yes", "no")
# h1 <-  cptable(~chem_activates_er, values=c(round(prior_er_pos,2), 1-round(prior_er_pos,2)),levels=yn)
# a1 <- cptable(~accuracy_atg_era_trans_up, values=c(posterior_atg_era_trans_up_assay_pos_er_pos + 
#                                                      posterior_atg_era_trans_up_assay_neg_er_neg, 
#                                                    1 - (posterior_atg_era_trans_up_assay_pos_er_pos + 
#                                                           posterior_atg_era_trans_up_assay_neg_er_neg)), levels=yn)
# e1 <- cptable(~atg_era_trans_up|accuracy_atg_era_trans_up:chem_activates_er, 
#               values=c(posterior_atg_era_trans_up_assay_pos_er_pos,
#                        posterior_atg_era_trans_up_assay_neg_er_pos,
#                        0.5,
#                        0.5,
#                        posterior_atg_era_trans_up_assay_pos_er_neg,
#                        posterior_atg_era_trans_up_assay_neg_er_neg,
#                        0.5,
#                        0.5),
#               levels=yn)
# a2 <- cptable(~accuracy_nvs_nr_her, values=c(posterior_nvs_nr_her_assay_pos_er_pos +
#                                                posterior_nvs_nr_her_assay_neg_er_neg,
#                                              1 - (posterior_nvs_nr_her_assay_pos_er_pos +
#                                                     posterior_nvs_nr_her_assay_neg_er_neg)), levels=yn)
# e2 <- cptable(~nvs_nr_her|accuracy_nvs_nr_her:chem_activates_er,
#               values = c(posterior_nvs_nr_her_assay_pos_er_pos,
#                          posterior_nvs_nr_her_assay_neg_er_pos,
#                          0.5,
#                          0.5,
#                          posterior_nvs_nr_her_assay_pos_er_neg,
#                          posterior_nvs_nr_her_assay_neg_er_neg,
#                          0.5,
#                          0.5),
#               levels=yn)
# a3 <- cptable(~accuracy_ot_er_eraera_0480, values=c(posterior_OT_ER_ERaERa_0480_assay_pos_er_pos +
#                                                       posterior_OT_ER_ERaERa_0480_assay_neg_er_neg,
#                                                     1 - (posterior_OT_ER_ERaERa_0480_assay_pos_er_pos +
#                                                            posterior_OT_ER_ERaERa_0480_assay_neg_er_neg)), levels=yn)
# e3 <- cptable(~ot_er_eraera_0480|accuracy_ot_er_eraera_0480:chem_activates_er,
#               values = c(posterior_OT_ER_ERaERa_0480_assay_pos_er_pos,
#                          posterior_OT_ER_ERaERa_0480_assay_neg_er_pos,
#                          0.5,
#                          0.5,
#                          posterior_OT_ER_ERaERa_0480_assay_pos_er_neg,
#                          posterior_OT_ER_ERaERa_0480_assay_neg_er_neg,
#                          0.5,
#                          0.5),
#               levels=yn)
# 
# h2 <- cptable(~chem_estrogenic|chem_activates_er,
#               values = c(0.9, 0.1, 0.01, 0.99),
#               levels=yn)
# 
# plist <- compileCPT(list(h1, a1, e1, a2, e2, a3, e3, h2))
# estrogenic_bn <- grain(plist)
# 
# querygrain(setEvidence(estrogenic_bn, evidence=list(chem_activates_er="yes")) )
# 
# estrogenic_bn_run <- function(ev_list, estrogenic_bn){
#   evidence_list <- list(atg_era_trans_up=ev_list$atg_era_trans_up, nvs_nr_her=ev_list$nvs_nr_her, 
#                         ot_er_eraera_0480=ev_list$ot_er_eraera_0480)
#   grain_results <- querygrain(setEvidence(estrogenic_bn, evidence=evidence_list))
#   return(grain_results)
# }
# 
# convert_nums_to_yn <- function(x){
#   if(length(x)==1 && x == 1){
#     return("yes")
#   }
#   else{
#     return("no")
#   }
# }
# 
# unique_chems_er_pos <- unique(tc_gt_er_pos$casn)
# unique_chems_er_neg <- unique(tc_gt_er_neg$casn)
# 
# tc_gt_er_pos_list <- list()
# tc_gt_er_neg_list <- list()
# 
# for(i in 1:length(unique_chems_er_pos)){
#   tc_gt_er_pos_chem_sub <- subset(tc_gt_er_pos, casn == unique_chems_er_pos[i])
#   chem_assay_list <- list()
#   chem_assay_list$atg_era_trans_up <- convert_nums_to_yn(subset(tc_gt_er_pos_chem_sub, aenm == "ATG_ERa_TRANS_up")$hitc)
#   chem_assay_list$nvs_nr_her <- convert_nums_to_yn(subset(tc_gt_er_pos_chem_sub, aenm == "NVS_NR_hER")$hitc)
#   chem_assay_list$ot_er_eraera_0480 <- convert_nums_to_yn(subset(tc_gt_er_pos_chem_sub, aenm == "OT_ER_ERaERa_0480")$hitc)
#   print(chem_assay_list)
#   tc_gt_er_pos_list[[unique_chems_er_pos[i]]] <- chem_assay_list
# }
# 
# toxcast_er_pos_grain_results <- lapply(tc_gt_er_pos_list, estrogenic_bn_run, estrogenic_bn)
# toxcast_er_pos_grain_results_probs <- NULL
# 
# for(i in 1:length(toxcast_er_pos_grain_results)){
#   toxcast_er_pos_grain_results_probs[i] <- toxcast_er_pos_grain_results[[i]]$chem_estrogenic[[1]]
# }
# 
# for(i in 1:length(unique_chems_er_neg)){
#   tc_gt_er_neg_chem_sub <- subset(tc_gt_er_pos, casn == unique_chems_er_pos[i])
#   chem_assay_list <- list()
#   chem_assay_list$atg_era_trans_up <- convert_nums_to_yn(subset(tc_gt_er_neg_chem_sub, aenm == "ATG_ERa_TRANS_up")$hitc)
#   chem_assay_list$nvs_nr_her <- convert_nums_to_yn(subset(tc_gt_er_neg_chem_sub, aenm == "NVS_NR_hER")$hitc)
#   chem_assay_list$ot_er_eraera_0480 <- convert_nums_to_yn(subset(tc_gt_er_neg_chem_sub, aenm == "OT_ER_ERaERa_0480")$hitc)
#   print(chem_assay_list)
#   tc_gt_er_neg_list[[unique_chems_er_neg[i]]] <- chem_assay_list
# }
# 
# toxcast_er_neg_grain_results <- lapply(tc_gt_er_neg_list, estrogenic_bn_run, estrogenic_bn)
# toxcast_er_neg_grain_results_probs <- NULL
# 
# for(i in 1:length(toxcast_er_neg_grain_results)){
#   toxcast_er_neg_grain_results_probs[i] <- toxcast_er_neg_grain_results[[i]]$chem_estrogenic[[1]]
# }
# 
# true_positives <- length(which(toxcast_er_pos_grain_results_probs > 0.50)) / length(toxcast_er_pos_grain_results)
# true_negatives <- length(which(toxcast_er_neg_grain_results_probs < 0.50)) / length(toxcast_er_neg_grain_results)
# false_positives <- length(which(toxcast_er_neg_grain_results_probs > 0.50)) / length(toxcast_er_neg_grain_results)
# false_negatives <- length(which(toxcast_er_pos_grain_results_probs < 0.50)) / length(toxcast_er_pos_grain_results)
# 
# toxcast_preds <- c(toxcast_er_pos_grain_results_probs, toxcast_er_neg_grain_results_probs)
# toxcast_labs <- c(rep(1, length(toxcast_er_pos_grain_results_probs)), 
#                   rep(0, length(toxcast_er_neg_grain_results_probs)))
# er_tc_preds <- prediction(toxcast_preds, toxcast_labs)
# er_tc_perf <- performance(er_tc_preds, "tpr", "fpr")
# plot(er_tc_perf, colorize=TRUE)
# 
# plot(performance(er_tc_preds, "sens", "fpr"), colorize=TRUE)
# abline(0,1)

#######
# Scenario 3X: Assay(+), But posterior TP rate changes. What Should Player 2 Do? Go vs No-Go
#######
#Posterior Probabilities (across all conditions)
##P(ER+ | ATG_ERa_TRANS_up +)
amalgum_pos_er_pos_rate <- .73
amalgum_pos_er_neg_rate <- .65
amalgum_neg_er_neg_rate <- .35
amalgum_neg_er_pos_rate <- .27

posterior_amalgum_pos_er_pos <- 
  (amalgum_pos_er_pos_rate * prior_er_pos) / ((amalgum_pos_er_pos_rate * prior_er_pos) +
                                                (amalgum_pos_er_neg_rate * prior_er_neg) +
                                                (amalgum_neg_er_neg_rate * prior_er_neg) +
                                                (amalgum_neg_er_pos_rate * prior_er_pos))
##P(ER+ | ATG_ERa_TRANS_up -)
posterior_amalgum_neg_er_pos <-
  (atg_era_trans_up_assay_neg_er_pos_rate * prior_er_pos) / ((atg_era_trans_up_assay_pos_er_pos_rate * prior_er_pos) +
                                                               (atg_era_trans_up_assay_pos_er_neg_rate * prior_er_neg) +
                                                               (atg_era_trans_up_assay_neg_er_neg_rate * prior_er_neg) +
                                                               (atg_era_trans_up_assay_neg_er_pos_rate * prior_er_pos))
##P(ER- | ATG_ERa_TRANS_up +)
posterior_amalgum_pos_er_neg <-
  (atg_era_trans_up_assay_pos_er_neg_rate * prior_er_neg) / ((atg_era_trans_up_assay_pos_er_pos_rate * prior_er_pos) +
                                                               (atg_era_trans_up_assay_pos_er_neg_rate * prior_er_neg) +
                                                               (atg_era_trans_up_assay_neg_er_neg_rate * prior_er_neg) +
                                                               (atg_era_trans_up_assay_neg_er_pos_rate * prior_er_pos))
##P(ER- | ATG_ERa_TRANS_up -)
posterior_amalgum_neg_er_neg <-
  (atg_era_trans_up_assay_neg_er_neg_rate * prior_er_neg) / ((atg_era_trans_up_assay_pos_er_pos_rate * prior_er_pos) +
                                                               (atg_era_trans_up_assay_pos_er_neg_rate * prior_er_neg) +
                                                               (atg_era_trans_up_assay_neg_er_neg_rate * prior_er_neg) +
                                                               (atg_era_trans_up_assay_neg_er_pos_rate * prior_er_pos))


posterior_amalgum_pos_er_pos <- seq(from=0.01, to=0.50, by = 0.01)
posterior_amalgum_neg_er_pos <- 0.5 - posterior_amalgum_pos_er_pos
results_mat <- data.frame(posterior_true_pos_rate = numeric(), expected_utility.go = numeric(), expected_utility.nogo = numeric())

for(i in 1:length(posterior_amalgum_pos_er_pos)){
  print(posterior_amalgum_pos_er_pos[i])
  eu_go_decision_amalgum_pos_temp <- ((player2_payoff_estrogenic_hit_go * 
                                         posterior_amalgum_pos_er_pos[i]) +
                                        (player2_payoff_not_estrogenic_hit_go *
                                           posterior_amalgum_pos_er_neg))
  
  eu_nogo_decision_amalgum_pos_temp <- ((player2_payoff_estrogenic_hit_nogo *
                                           posterior_amalgum_pos_er_pos[i]) +
                                          player2_payoff_not_estrogenic_hit_nogo *
                                          posterior_amalgum_pos_er_neg)
  
  print(eu_go_decision_amalgum_pos_temp)
  print(eu_nogo_decision_amalgum_pos_temp)
  
  results_mat[i,] <- c(posterior_amalgum_pos_er_pos[i], eu_go_decision_amalgum_pos_temp,
                       eu_nogo_decision_amalgum_pos_temp)
}

results_mat
results_mat <- reshape(results_mat, direction="long",varying=c(2,3))

#This shows the equilibrium point with respect to P(Assay+ | ER+)
png("p_amalgum+_er+_equilibrium_plot.png")
ggplot(data=results_mat, aes(x=posterior_true_pos_rate, y=expected_utility, shape=time)) +
  geom_point() +
  xlab("P(Assay+ | ER+)")
dev.off()


#The results are still pretty awful. What I did was use the ToxCast hit term.
#Clearly that didn't work.
#
#I think this is a case where using EPA's
#AUC approach may work well. Basically, EPA calculates the AUC using the 
#trapezoidal rule for each dose-response curve, and then compares that to the
#dose-response curve for 17-beta estradiol in that same assay. That AUC ratio
#would give me an indication of where to set cut-offs to say that the assay
#is giving me an estrogenic chemical.
#
#So, I would calculate the AUC ratios, look at the Sensitivity vs False Positive Rate
#ROC curve on a per-assay basis, and determine a cut-off. 

n <- 51
x <- seq(0, pi/2, len = n)
y <- sin(x)
trapz(x, y)  

