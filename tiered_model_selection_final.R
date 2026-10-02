#######################################################################
######### TIERED MODEL SELECTION — MOOSE iSSF ANALYSIS ################
######### Yukon-Charley Rivers National Preserve (2020–2023) ##########
#######################################################################
#
# Tier 1 — Additive covariate selection
#   1a: Topography (ruggedness vs elevation, linear vs quadratic)
#   1b: River distance (include or exclude)
#   1c: Habitat covariates (subsets of salix, wetland, spruce, deciduous)
#
# Tier 2 — CalfStatus × habitat interactions (which habitats does calf
#          status modify?)
#
# Tier 3 — Wolf × CalfStatus × habitat three-way interactions (which
#          habitat × calf interactions are further modulated by wolf
#          density?)
#
# NOTE: Tier 1a/1b decisions feed into 1c. Tier 1c winner feeds into
# Tier 2. Tier 2 winner feeds into Tier 3. 
#######################################################################

pkgs = c('tidyverse', 'glmmTMB', 'AICcmodavg')
sapply(pkgs, require, character = TRUE)
conflicted::conflict_prefer("filter", "dplyr")

here::here() %>% setwd()

mse_ssf = readRDS("data/mse_ssf_full2.rds") %>%
  mutate(Step = paste0(Step, "_", AnimalId))

#### seasonal comparison ####
mse_ssf_sum_cf = mse_ssf %>% dplyr::filter(month(Time_2) %in% c(5,6,7,8)) %>% na.omit() %>% filter(CalfStatus == 1)
mse_ssf_sum_nc = mse_ssf %>% dplyr::filter(month(Time_2) %in% c(5,6,7,8)) %>% na.omit() %>% filter(CalfStatus == 0)
mse_ssf_win_cf = mse_ssf %>% dplyr::filter(month(Time_2) %in% c(11,12,1,2)) %>% na.omit() %>% filter(CalfStatus == 1)
mse_ssf_win_nc = mse_ssf %>% dplyr::filter(month(Time_2) %in% c(11,12,1,2)) %>% na.omit() %>% filter(CalfStatus == 0)

mse_ssf_sum_cs = rbind(mse_ssf_sum_cf, mse_ssf_sum_nc)
mse_ssf_win_cs = rbind(mse_ssf_win_cf, mse_ssf_win_nc)


#######################################################################
############## TIER 1a: TOPOGRAPHY ####################################
#######################################################################
# Hold habitat (all 4) + wolf + river constant. Vary topo metric.
# 4 models per season.
#######################################################################

#### SUMMER ####

t1a_1_sum <- readRDS("model_outputs/t1a_1_sum.rds")
pars <- getME(t1a_1_sum, "theta")
betas <- fixef(t1a_1_sum)$cond

# t1a_2_sum: Ruggedness linear
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + salshr_sc + wtlnd_sc +
                      picgla_sc + dectre_sc + riv_dst_sc + rugg_sc +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+salshr_sc|AnimalId) + (0+wtlnd_sc|AnimalId) +
                      (0+picgla_sc|AnimalId) + (0+dectre_sc|AnimalId) + (0+riv_dst_sc|AnimalId) +
                      (0+rugg_sc|AnimalId),
                    family=poisson, data=mse_ssf_sum_cs,
                    control=glmmTMBControl(parallel=4, optimizer=optim, optArgs=list(method="BFGS")),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(0, 10))),
                    map=list(theta=factor(c(NA, 1:10))))
t1a_1_sum = glmmTMB:::fitTMB(mod_struc)


# t1a_2_sum: Ruggedness quadratic
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + salshr_sc + wtlnd_sc +
                      picgla_sc + dectre_sc + riv_dst_sc + rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+salshr_sc|AnimalId) + (0+wtlnd_sc|AnimalId) +
                      (0+picgla_sc|AnimalId) + (0+dectre_sc|AnimalId) + (0+riv_dst_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_sum_cs,
                    control=glmmTMBControl(parallel=4, optimizer=optim, optArgs=list(method="BFGS")),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(log(0.5), 11))),
                    map=list(theta=factor(c(NA, 1:11))))
t1a_2_sum = glmmTMB:::fitTMB(mod_struc)
saveRDS(t1a_2_sum, "model_outputs/t1a_2_sum.rds"); rm(t1a_2_sum)

# t1a_3_sum: Elevation linear
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + salshr_sc + wtlnd_sc +
                      picgla_sc + dectre_sc + riv_dst_sc + elev_sc +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+salshr_sc|AnimalId) + (0+wtlnd_sc|AnimalId) +
                      (0+picgla_sc|AnimalId) + (0+dectre_sc|AnimalId) + (0+riv_dst_sc|AnimalId) +
                      (0+elev_sc|AnimalId),
                    family=poisson, data=mse_ssf_sum_cs,
                    control=glmmTMBControl(parallel=4),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(log(0.5), 10))),
                    map=list(theta=factor(c(NA, 1:10))))
t1a_3_sum = glmmTMB:::fitTMB(mod_struc)
saveRDS(t1a_3_sum, "model_outputs/t1a_3_sum.rds"); rm(t1a_3_sum)

# t1a_4_sum: Elevation quadratic
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + salshr_sc + wtlnd_sc +
                      picgla_sc + dectre_sc + riv_dst_sc + elev_sc + I(elev_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+salshr_sc|AnimalId) + (0+wtlnd_sc|AnimalId) +
                      (0+picgla_sc|AnimalId) + (0+dectre_sc|AnimalId) + (0+riv_dst_sc|AnimalId) +
                      (0+elev_sc|AnimalId) + (0+I(elev_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_sum_cs,
                    control=glmmTMBControl(parallel=4),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(log(0.5), 11))),
                    map=list(theta=factor(c(NA, 1:11))))
t1a_4_sum = glmmTMB:::fitTMB(mod_struc)
saveRDS(t1a_4_sum, "model_outputs/t1a_4_sum.rds"); rm(t1a_4_sum)


#### WINTER ####

# t1a_1_win: Ruggedness linear
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + salshr_sc + wtlnd_sc +
                      picgla_sc + dectre_sc + riv_dst_sc + rugg_sc +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+salshr_sc|AnimalId) + (0+wtlnd_sc|AnimalId) +
                      (0+picgla_sc|AnimalId) + (0+dectre_sc|AnimalId) + (0+riv_dst_sc|AnimalId) +
                      (0+rugg_sc|AnimalId),
                    family=poisson, data=mse_ssf_win_cs,
                    control=glmmTMBControl(parallel=4),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(log(0.5), 10))),
                    map=list(theta=factor(c(NA, 1:10))))
t1a_1_win = glmmTMB:::fitTMB(mod_struc)
saveRDS(t1a_1_win, "model_outputs/t1a_1_win.rds"); rm(t1a_1_win)

# t1a_2_win: Ruggedness quadratic
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + salshr_sc + wtlnd_sc +
                      picgla_sc + dectre_sc + riv_dst_sc + rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+salshr_sc|AnimalId) + (0+wtlnd_sc|AnimalId) +
                      (0+picgla_sc|AnimalId) + (0+dectre_sc|AnimalId) + (0+riv_dst_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_win_cs,
                    control=glmmTMBControl(parallel=4),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(log(0.5), 11))),
                    map=list(theta=factor(c(NA, 1:11))))
t1a_2_win = glmmTMB:::fitTMB(mod_struc)
saveRDS(t1a_2_win, "model_outputs/t1a_2_win.rds"); rm(t1a_2_win)

# t1a_3_win: Elevation linear
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + salshr_sc + wtlnd_sc +
                      picgla_sc + dectre_sc + riv_dst_sc + elev_sc +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+salshr_sc|AnimalId) + (0+wtlnd_sc|AnimalId) +
                      (0+picgla_sc|AnimalId) + (0+dectre_sc|AnimalId) + (0+riv_dst_sc|AnimalId) +
                      (0+elev_sc|AnimalId),
                    family=poisson, data=mse_ssf_win_cs,
                    control=glmmTMBControl(parallel=4),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(log(0.5), 10))),
                    map=list(theta=factor(c(NA, 1:10))))
t1a_3_win = glmmTMB:::fitTMB(mod_struc)
saveRDS(t1a_3_win, "model_outputs/t1a_3_win.rds"); rm(t1a_3_win)

# t1a_4_win: Elevation quadratic
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + salshr_sc + wtlnd_sc +
                      picgla_sc + dectre_sc + riv_dst_sc + elev_sc + I(elev_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+salshr_sc|AnimalId) + (0+wtlnd_sc|AnimalId) +
                      (0+picgla_sc|AnimalId) + (0+dectre_sc|AnimalId) + (0+riv_dst_sc|AnimalId) +
                      (0+elev_sc|AnimalId) + (0+I(elev_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_win_cs,
                    control=glmmTMBControl(parallel=4),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(log(0.5), 11))),
                    map=list(theta=factor(c(NA, 1:11))))
t1a_4_win = glmmTMB:::fitTMB(mod_struc)
saveRDS(t1a_4_win, "model_outputs/t1a_4_win.rds"); rm(t1a_4_win)


#### TIER 1a COMPARISON ####
t1a_1_sum=readRDS("model_outputs/t1a_1_sum.rds"); t1a_2_sum=readRDS("model_outputs/t1a_2_sum.rds")
t1a_3_sum=readRDS("model_outputs/t1a_3_sum.rds"); t1a_4_sum=readRDS("model_outputs/t1a_4_sum.rds")
aictab(lst(t1a_1_sum, t1a_2_sum, t1a_3_sum, t1a_4_sum))
rm(t1a_1_sum, t1a_2_sum, t1a_3_sum, t1a_4_sum)

t1a_1_win=readRDS("model_outputs/t1a_1_win.rds"); t1a_2_win=readRDS("model_outputs/t1a_2_win.rds")
t1a_3_win=readRDS("model_outputs/t1a_3_win.rds"); t1a_4_win=readRDS("model_outputs/t1a_4_win.rds")
aictab(lst(t1a_1_win, t1a_2_win, t1a_3_win, t1a_4_win))
rm(t1a_1_win, t1a_2_win, t1a_3_win, t1a_4_win)


#######################################################################
############## TIER 1b: RIVER DISTANCE ################################
#######################################################################
# Hold habitat (all 4) + wolf + best topo constant. Compare with/without
# river distance. 2 models per season.
# Assumes: rugg_sc + I(rugg_sc^2) from Tier 1a.
#######################################################################

#### SUMMER ####

# t1b_1_sum: With river distance
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + salshr_sc + wtlnd_sc +
                      picgla_sc + dectre_sc + riv_dst_sc + rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+salshr_sc|AnimalId) + (0+wtlnd_sc|AnimalId) +
                      (0+picgla_sc|AnimalId) + (0+dectre_sc|AnimalId) + (0+riv_dst_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_sum_cs,
                    control=glmmTMBControl(parallel=4),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(0, 11))),
                    map=list(theta=factor(c(NA, 1:11))))
t1b_1_sum = glmmTMB:::fitTMB(mod_struc)
saveRDS(t1b_1_sum, "model_outputs/t1b_1_sum.rds"); rm(t1b_1_sum)

# t1b_2_sum: Without river distance
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + salshr_sc + wtlnd_sc +
                      picgla_sc + dectre_sc + rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+salshr_sc|AnimalId) + (0+wtlnd_sc|AnimalId) +
                      (0+picgla_sc|AnimalId) + (0+dectre_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_sum_cs,
                    control=glmmTMBControl(parallel=4),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(0, 10))),
                    map=list(theta=factor(c(NA, 1:10))))
t1b_2_sum = glmmTMB:::fitTMB(mod_struc)
saveRDS(t1b_2_sum, "model_outputs/t1b_2_sum.rds"); rm(t1b_2_sum)

#### WINTER ####

# t1b_1_win: With river distance
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + salshr_sc + wtlnd_sc +
                      picgla_sc + dectre_sc + riv_dst_sc + rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+salshr_sc|AnimalId) + (0+wtlnd_sc|AnimalId) +
                      (0+picgla_sc|AnimalId) + (0+dectre_sc|AnimalId) + (0+riv_dst_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_win_cs,
                    control=glmmTMBControl(parallel=4, optimizer=optim, optArgs=list(method="BFGS")),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(log(0.5), 11))),
                    map=list(theta=factor(c(NA, 1:11))))
t1b_1_win = glmmTMB:::fitTMB(mod_struc)
saveRDS(t1b_1_win, "model_outputs/t1b_1_win.rds"); rm(t1b_1_win)

# t1b_2_win: Without river distance
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + salshr_sc + wtlnd_sc +
                      picgla_sc + dectre_sc + rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+salshr_sc|AnimalId) + (0+wtlnd_sc|AnimalId) +
                      (0+picgla_sc|AnimalId) + (0+dectre_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_win_cs,
                    control=glmmTMBControl(parallel=4, optimizer=optim, optArgs=list(method="BFGS")),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(log(0.5), 10))),
                    map=list(theta=factor(c(NA, 1:10))))
t1b_2_win = glmmTMB:::fitTMB(mod_struc)
saveRDS(t1b_2_win, "model_outputs/t1b_2_win.rds"); rm(t1b_2_win)


#### TIER 1b COMPARISON ####
t1b_1_sum=readRDS("model_outputs/t1b_1_sum.rds"); t1b_2_sum=readRDS("model_outputs/t1b_2_sum.rds")
aictab(lst(t1b_1_sum, t1b_2_sum))
rm(t1b_1_sum, t1b_2_sum)
# 
t1b_1_win=readRDS("model_outputs/t1b_1_win.rds"); t1b_2_win=readRDS("model_outputs/t1b_2_win.rds")
aictab(lst(t1b_1_win, t1b_2_win))
rm(t1b_1_win, t1b_2_win)


#######################################################################
############## TIER 1c: HABITAT COVARIATES ############################
#######################################################################
# All non-empty subsets of {salix, wetland, spruce, deciduous}.
# Base: wlf_kde_sea_sc + riv_dst_sc + rugg_sc + I(rugg_sc^2)
# 15 models per season.
# Assumes: rugg_quad + river from Tier 1a/1b decisions.
#######################################################################

#### SUMMER ####

# t1c_01_sum: salix only
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + riv_dst_sc + salshr_sc +
                      rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+riv_dst_sc|AnimalId) + (0+salshr_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_sum_cs,
                    control=glmmTMBControl(parallel=4),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(log(0.5), 8))),
                    map=list(theta=factor(c(NA, 1:8))))
t1c_01_sum = glmmTMB:::fitTMB(mod_struc)
saveRDS(t1c_01_sum, "model_outputs/t1c_01_sum.rds"); rm(t1c_01_sum)

# t1c_02_sum: wetland only
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + riv_dst_sc + wtlnd_sc +
                      rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+riv_dst_sc|AnimalId) + (0+wtlnd_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_sum_cs,
                    control=glmmTMBControl(parallel=4, optimizer=optim, optArgs=list(method="BFGS")),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(log(0.5), 8))),
                    map=list(theta=factor(c(NA, 1:8))))
t1c_02_sum = glmmTMB:::fitTMB(mod_struc)
saveRDS(t1c_02_sum, "model_outputs/t1c_02_sum.rds"); rm(t1c_02_sum)

# t1c_03_sum: spruce only
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + riv_dst_sc + picgla_sc +
                      rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+riv_dst_sc|AnimalId) + (0+picgla_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_sum_cs,
                    control=glmmTMBControl(parallel=4),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(log(0.5), 8))),
                    map=list(theta=factor(c(NA, 1:8))))
t1c_03_sum = glmmTMB:::fitTMB(mod_struc)
saveRDS(t1c_03_sum, "model_outputs/t1c_03_sum.rds"); rm(t1c_03_sum)

# t1c_04_sum: deciduous only
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + riv_dst_sc + dectre_sc +
                      rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+riv_dst_sc|AnimalId) + (0+dectre_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_sum_cs,
                    control=glmmTMBControl(parallel=4),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(log(0.5), 8))),
                    map=list(theta=factor(c(NA, 1:8))))
t1c_04_sum = glmmTMB:::fitTMB(mod_struc)
saveRDS(t1c_04_sum, "model_outputs/t1c_04_sum.rds"); rm(t1c_04_sum)

# t1c_05_sum: salix + wetland
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + riv_dst_sc + salshr_sc + wtlnd_sc +
                      rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+riv_dst_sc|AnimalId) + (0+salshr_sc|AnimalId) +
                      (0+wtlnd_sc|AnimalId) + (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_sum_cs,
                    control=glmmTMBControl(parallel=4),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(0, 9))),
                    map=list(theta=factor(c(NA, 1:9))))
t1c_05_sum = glmmTMB:::fitTMB(mod_struc)
saveRDS(t1c_05_sum, "model_outputs/t1c_05_sum.rds"); rm(t1c_05_sum)

# t1c_06_sum: salix + spruce
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + riv_dst_sc + salshr_sc + picgla_sc +
                      rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+riv_dst_sc|AnimalId) + (0+salshr_sc|AnimalId) +
                      (0+picgla_sc|AnimalId) + (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_sum_cs,
                    control=glmmTMBControl(parallel=4, optimizer=optim, optArgs=list(method="BFGS")),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(log(0.5), 9))),
                    map=list(theta=factor(c(NA, 1:9))))
t1c_06_sum = glmmTMB:::fitTMB(mod_struc)
saveRDS(t1c_06_sum, "model_outputs/t1c_06_sum.rds"); rm(t1c_06_sum)

# t1c_07_sum: salix + deciduous
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + riv_dst_sc + salshr_sc + dectre_sc +
                      rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+riv_dst_sc|AnimalId) + (0+salshr_sc|AnimalId) +
                      (0+dectre_sc|AnimalId) + (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_sum_cs,
                    control=glmmTMBControl(parallel=4),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(0, 9))),
                    map=list(theta=factor(c(NA, 1:9))))
t1c_07_sum = glmmTMB:::fitTMB(mod_struc)
saveRDS(t1c_07_sum, "model_outputs/t1c_07_sum.rds"); rm(t1c_07_sum)

# t1c_08_sum: wetland + spruce
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + riv_dst_sc + wtlnd_sc + picgla_sc +
                      rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+riv_dst_sc|AnimalId) + (0+wtlnd_sc|AnimalId) +
                      (0+picgla_sc|AnimalId) + (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_sum_cs,
                    control=glmmTMBControl(parallel=4, optimizer=optim, optArgs=list(method="BFGS")),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(log(0.5), 9))),
                    map=list(theta=factor(c(NA, 1:9))))
t1c_08_sum = glmmTMB:::fitTMB(mod_struc)
saveRDS(t1c_08_sum, "model_outputs/t1c_08_sum.rds"); rm(t1c_08_sum)


t1c_09_sum <- readRDS("model_outputs/t1c_09_sum.rds")
pars <- getME(t1c_09_sum, "theta")
betas <- fixef(t1c_09_sum)$cond

# t1c_09_sum: wetland + deciduous
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + riv_dst_sc + wtlnd_sc + dectre_sc +
                      rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+riv_dst_sc|AnimalId) + (0+wtlnd_sc|AnimalId) +
                      (0+dectre_sc|AnimalId) + (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_sum_cs,
                    control=glmmTMBControl(parallel=4),
                    doFit=FALSE, start=list(beta = betas, theta = pars),
                    map=list(theta=factor(c(NA, 1:9))))
t1c_09_sum = glmmTMB:::fitTMB(mod_struc)
saveRDS(t1c_09_sum, "model_outputs/t1c_09_sum.rds"); rm(t1c_09_sum)

# t1c_10_sum: spruce + deciduous
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + riv_dst_sc + picgla_sc + dectre_sc +
                      rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+riv_dst_sc|AnimalId) + (0+picgla_sc|AnimalId) +
                      (0+dectre_sc|AnimalId) + (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_sum_cs,
                    control=glmmTMBControl(parallel=4, optimizer=optim, optArgs=list(method="BFGS")),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(log(0.5), 9))),
                    map=list(theta=factor(c(NA, 1:9))))
t1c_10_sum = glmmTMB:::fitTMB(mod_struc)
saveRDS(t1c_10_sum, "model_outputs/t1c_10_sum.rds"); rm(t1c_10_sum)

# t1c_11_sum: salix + wetland + spruce
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + riv_dst_sc + salshr_sc + wtlnd_sc + picgla_sc +
                      rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+riv_dst_sc|AnimalId) + (0+salshr_sc|AnimalId) +
                      (0+wtlnd_sc|AnimalId) + (0+picgla_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_sum_cs,
                    control=glmmTMBControl(parallel=4, optimizer=optim, optArgs=list(method="BFGS")),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(log(0.5), 10))),
                    map=list(theta=factor(c(NA, 1:10))))
t1c_11_sum = glmmTMB:::fitTMB(mod_struc)
saveRDS(t1c_11_sum, "model_outputs/t1c_11_sum.rds"); rm(t1c_11_sum)

# t1c_12_sum: salix + wetland + deciduous
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + riv_dst_sc + salshr_sc + wtlnd_sc + dectre_sc +
                      rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+riv_dst_sc|AnimalId) + (0+salshr_sc|AnimalId) +
                      (0+wtlnd_sc|AnimalId) + (0+dectre_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_sum_cs,
                    control=glmmTMBControl(parallel=4, optimizer=optim, optArgs=list(method="BFGS")),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(log(0.5), 10))),
                    map=list(theta=factor(c(NA, 1:10))))
t1c_12_sum = glmmTMB:::fitTMB(mod_struc)
saveRDS(t1c_12_sum, "model_outputs/t1c_12_sum.rds"); rm(t1c_12_sum)

# t1c_13_sum: salix + spruce + deciduous
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + riv_dst_sc + salshr_sc + picgla_sc + dectre_sc +
                      rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+riv_dst_sc|AnimalId) + (0+salshr_sc|AnimalId) +
                      (0+picgla_sc|AnimalId) + (0+dectre_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_sum_cs,
                    control=glmmTMBControl(parallel=4),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(0, 10))),
                    map=list(theta=factor(c(NA, 1:10))))
t1c_13_sum = glmmTMB:::fitTMB(mod_struc)
saveRDS(t1c_13_sum, "model_outputs/t1c_13_sum.rds"); rm(t1c_13_sum)

# t1c_14_sum: wetland + spruce + deciduous
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + riv_dst_sc + wtlnd_sc + picgla_sc + dectre_sc +
                      rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+riv_dst_sc|AnimalId) + (0+wtlnd_sc|AnimalId) +
                      (0+picgla_sc|AnimalId) + (0+dectre_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_sum_cs,
                    control=glmmTMBControl(parallel=4, optimizer=optim, optArgs=list(method="BFGS")),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(log(0.5), 10))),
                    map=list(theta=factor(c(NA, 1:10))))
t1c_14_sum = glmmTMB:::fitTMB(mod_struc)
saveRDS(t1c_14_sum, "model_outputs/t1c_14_sum.rds"); rm(t1c_14_sum)

# t1c_15_sum: salix + wetland + spruce + deciduous
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + riv_dst_sc + salshr_sc + wtlnd_sc + picgla_sc + dectre_sc +
                      rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+riv_dst_sc|AnimalId) + (0+salshr_sc|AnimalId) +
                      (0+wtlnd_sc|AnimalId) + (0+picgla_sc|AnimalId) + (0+dectre_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_sum_cs,
                    control=glmmTMBControl(parallel=4),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(0, 11))),
                    map=list(theta=factor(c(NA, 1:11))))
t1c_15_sum = glmmTMB:::fitTMB(mod_struc)
saveRDS(t1c_15_sum, "model_outputs/t1c_15_sum.rds"); rm(t1c_15_sum)


#### WINTER — same 15 models, data = mse_ssf_win_cs ####

# t1c_01_win: salix only
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + riv_dst_sc + salshr_sc +
                      rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+riv_dst_sc|AnimalId) + (0+salshr_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_win_cs,
                    control=glmmTMBControl(parallel=4, optimizer=optim, optArgs=list(method="BFGS")),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(log(0.5), 8))),
                    map=list(theta=factor(c(NA, 1:8))))
t1c_01_win = glmmTMB:::fitTMB(mod_struc)
saveRDS(t1c_01_win, "model_outputs/t1c_01_win.rds"); rm(t1c_01_win)

# t1c_02_win: wetland only
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + riv_dst_sc + wtlnd_sc +
                      rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+riv_dst_sc|AnimalId) + (0+wtlnd_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_win_cs,
                    control=glmmTMBControl(parallel=4, optimizer=optim, optArgs=list(method="BFGS")),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(log(0.5), 8))),
                    map=list(theta=factor(c(NA, 1:8))))
t1c_02_win = glmmTMB:::fitTMB(mod_struc)
saveRDS(t1c_02_win, "model_outputs/t1c_02_win.rds"); rm(t1c_02_win)

# t1c_03_win: spruce only
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + riv_dst_sc + picgla_sc +
                      rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+riv_dst_sc|AnimalId) + (0+picgla_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_win_cs,
                    control=glmmTMBControl(parallel=4, optimizer=optim, optArgs=list(method="BFGS")),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(log(0.5), 8))),
                    map=list(theta=factor(c(NA, 1:8))))
t1c_03_win = glmmTMB:::fitTMB(mod_struc)
saveRDS(t1c_03_win, "model_outputs/t1c_03_win.rds"); rm(t1c_03_win)

# t1c_04_win: deciduous only
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + riv_dst_sc + dectre_sc +
                      rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+riv_dst_sc|AnimalId) + (0+dectre_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_win_cs,
                    control=glmmTMBControl(parallel=4, optimizer=optim, optArgs=list(method="BFGS")),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(log(0.5), 8))),
                    map=list(theta=factor(c(NA, 1:8))))
t1c_04_win = glmmTMB:::fitTMB(mod_struc)
saveRDS(t1c_04_win, "model_outputs/t1c_04_win.rds"); rm(t1c_04_win)

# t1c_05_win: salix + wetland
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + riv_dst_sc + salshr_sc + wtlnd_sc +
                      rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+riv_dst_sc|AnimalId) + (0+salshr_sc|AnimalId) +
                      (0+wtlnd_sc|AnimalId) + (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_win_cs,
                    control=glmmTMBControl(parallel=4, optimizer=optim, optArgs=list(method="BFGS")),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(log(0.5), 9))),
                    map=list(theta=factor(c(NA, 1:9))))
t1c_05_win = glmmTMB:::fitTMB(mod_struc)
saveRDS(t1c_05_win, "model_outputs/t1c_05_win.rds"); rm(t1c_05_win)

# t1c_06_win: salix + spruce
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + riv_dst_sc + salshr_sc + picgla_sc +
                      rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+riv_dst_sc|AnimalId) + (0+salshr_sc|AnimalId) +
                      (0+picgla_sc|AnimalId) + (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_win_cs,
                    control=glmmTMBControl(parallel=4, optimizer=optim, optArgs=list(method="BFGS")),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(log(0.5), 9))),
                    map=list(theta=factor(c(NA, 1:9))))
t1c_06_win = glmmTMB:::fitTMB(mod_struc)
saveRDS(t1c_06_win, "model_outputs/t1c_06_win.rds"); rm(t1c_06_win)

# t1c_07_win: salix + deciduous
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + riv_dst_sc + salshr_sc + dectre_sc +
                      rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+riv_dst_sc|AnimalId) + (0+salshr_sc|AnimalId) +
                      (0+dectre_sc|AnimalId) + (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_win_cs,
                    control=glmmTMBControl(parallel=4, optimizer=optim, optArgs=list(method="BFGS")),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(log(0.5), 9))),
                    map=list(theta=factor(c(NA, 1:9))))
t1c_07_win = glmmTMB:::fitTMB(mod_struc)
saveRDS(t1c_07_win, "model_outputs/t1c_07_win.rds"); rm(t1c_07_win)

# t1c_08_win: wetland + spruce
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + riv_dst_sc + wtlnd_sc + picgla_sc +
                      rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+riv_dst_sc|AnimalId) + (0+wtlnd_sc|AnimalId) +
                      (0+picgla_sc|AnimalId) + (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_win_cs,
                    control=glmmTMBControl(parallel=4, optimizer=optim, optArgs=list(method="BFGS")),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(log(0.5), 9))),
                    map=list(theta=factor(c(NA, 1:9))))
t1c_08_win = glmmTMB:::fitTMB(mod_struc)
saveRDS(t1c_08_win, "model_outputs/t1c_08_win.rds"); rm(t1c_08_win)

# t1c_09_win: wetland + deciduous
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + riv_dst_sc + wtlnd_sc + dectre_sc +
                      rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+riv_dst_sc|AnimalId) + (0+wtlnd_sc|AnimalId) +
                      (0+dectre_sc|AnimalId) + (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_win_cs,
                    control=glmmTMBControl(parallel=4, optimizer=optim, optArgs=list(method="BFGS")),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(log(0.5), 9))),
                    map=list(theta=factor(c(NA, 1:9))))
t1c_09_win = glmmTMB:::fitTMB(mod_struc)
saveRDS(t1c_09_win, "model_outputs/t1c_09_win.rds"); rm(t1c_09_win)

# t1c_10_win: spruce + deciduous
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + riv_dst_sc + picgla_sc + dectre_sc +
                      rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+riv_dst_sc|AnimalId) + (0+picgla_sc|AnimalId) +
                      (0+dectre_sc|AnimalId) + (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_win_cs,
                    control=glmmTMBControl(parallel=4, optimizer=optim, optArgs=list(method="BFGS")),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(log(0.5), 9))),
                    map=list(theta=factor(c(NA, 1:9))))
t1c_10_win = glmmTMB:::fitTMB(mod_struc)
saveRDS(t1c_10_win, "model_outputs/t1c_10_win.rds"); rm(t1c_10_win)

# t1c_11_win: salix + wetland + spruce
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + riv_dst_sc + salshr_sc + wtlnd_sc + picgla_sc +
                      rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+riv_dst_sc|AnimalId) + (0+salshr_sc|AnimalId) +
                      (0+wtlnd_sc|AnimalId) + (0+picgla_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_win_cs,
                    control=glmmTMBControl(parallel=4, optimizer=optim, optArgs=list(method="BFGS")),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(log(0.5), 10))),
                    map=list(theta=factor(c(NA, 1:10))))
t1c_11_win = glmmTMB:::fitTMB(mod_struc)
saveRDS(t1c_11_win, "model_outputs/t1c_11_win.rds"); rm(t1c_11_win)

# t1c_12_win: salix + wetland + deciduous
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + riv_dst_sc + salshr_sc + wtlnd_sc + dectre_sc +
                      rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+riv_dst_sc|AnimalId) + (0+salshr_sc|AnimalId) +
                      (0+wtlnd_sc|AnimalId) + (0+dectre_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_win_cs,
                    control=glmmTMBControl(parallel=4, optimizer=optim, optArgs=list(method="BFGS")),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(log(0.5), 10))),
                    map=list(theta=factor(c(NA, 1:10))))
t1c_12_win = glmmTMB:::fitTMB(mod_struc)
saveRDS(t1c_12_win, "model_outputs/t1c_12_win.rds"); rm(t1c_12_win)

# t1c_13_win: salix + spruce + deciduous
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + riv_dst_sc + salshr_sc + picgla_sc + dectre_sc +
                      rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+riv_dst_sc|AnimalId) + (0+salshr_sc|AnimalId) +
                      (0+picgla_sc|AnimalId) + (0+dectre_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_win_cs,
                    control=glmmTMBControl(parallel=4, optimizer=optim, optArgs=list(method="BFGS")),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(log(0.5), 10))),
                    map=list(theta=factor(c(NA, 1:10))))
t1c_13_win = glmmTMB:::fitTMB(mod_struc)
saveRDS(t1c_13_win, "model_outputs/t1c_13_win.rds"); rm(t1c_13_win)

# t1c_14_win: wetland + spruce + deciduous
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + riv_dst_sc + wtlnd_sc + picgla_sc + dectre_sc +
                      rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+riv_dst_sc|AnimalId) + (0+wtlnd_sc|AnimalId) +
                      (0+picgla_sc|AnimalId) + (0+dectre_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_win_cs,
                    control=glmmTMBControl(parallel=4, optimizer=optim, optArgs=list(method="BFGS")),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(log(0.5), 10))),
                    map=list(theta=factor(c(NA, 1:10))))
t1c_14_win = glmmTMB:::fitTMB(mod_struc)
saveRDS(t1c_14_win, "model_outputs/t1c_14_win.rds"); rm(t1c_14_win)

# t1c_15_win: salix + wetland + spruce + deciduous
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + riv_dst_sc + salshr_sc + wtlnd_sc + picgla_sc + dectre_sc +
                      rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+riv_dst_sc|AnimalId) + (0+salshr_sc|AnimalId) +
                      (0+wtlnd_sc|AnimalId) + (0+picgla_sc|AnimalId) + (0+dectre_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_win_cs,
                    control=glmmTMBControl(parallel=4, optimizer=optim, optArgs=list(method="BFGS")),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(log(0.5), 11))),
                    map=list(theta=factor(c(NA, 1:11))))
t1c_15_win = glmmTMB:::fitTMB(mod_struc)
saveRDS(t1c_15_win, "model_outputs/t1c_15_win.rds"); rm(t1c_15_win)


#### TIER 1c COMPARISON ####
t1c_sum_nms = paste0("t1c_", sprintf("%02d", 1:15), "_sum")
for(nm in t1c_sum_nms) assign(nm, readRDS(paste0("model_outputs/", nm, ".rds")))
aictab(mget(t1c_sum_nms))
rm(list = t1c_sum_nms)

t1c_win_nms = paste0("t1c_", sprintf("%02d", 1:15), "_win")
for(nm in t1c_win_nms) assign(nm, readRDS(paste0("model_outputs/", nm, ".rds")))
aictab(mget(t1c_win_nms))
rm(list = t1c_win_nms)


# t1c_15_sum: salix + wetland + spruce + deciduous
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_full_sc + riv_dst_sc + salshr_sc +
                      wtlnd_sc + picgla_sc + dectre_sc + rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_full_sc|AnimalId) + (0+riv_dst_sc|AnimalId) + (0+salshr_sc|AnimalId) +
                      (0+wtlnd_sc|AnimalId) + (0+picgla_sc|AnimalId) + (0+dectre_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_sum_cs,
                    control=glmmTMBControl(parallel=4),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(0, 11))),
                    map=list(theta=factor(c(NA, 1:11))))
t1c_15_full = glmmTMB:::fitTMB(mod_struc)
saveRDS(t1c_15_full, "model_outputs/t1c_15_full.rds"); rm(t1c_15_full)

# t1c_15_sum: salix + wetland + spruce + deciduous
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wolf_kde_pack_sc + riv_dst_sc + salshr_sc +
                      wtlnd_sc + picgla_sc + dectre_sc + rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wolf_kde_pack_sc|AnimalId) + (0+riv_dst_sc|AnimalId) + (0+salshr_sc|AnimalId) +
                      (0+wtlnd_sc|AnimalId) + (0+picgla_sc|AnimalId) + (0+dectre_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_sum_cs,
                    control=glmmTMBControl(parallel=4),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(0, 11))),
                    map=list(theta=factor(c(NA, 1:11))))
t1c_15_pack = glmmTMB:::fitTMB(mod_struc)
saveRDS(t1c_15_pack, "model_outputs/t1c_15_pack.rds"); rm(t1c_15_pack)

# t1c_15_sum: salix + wetland + spruce + deciduous
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_wt_sc + riv_dst_sc + salshr_sc +
                      wtlnd_sc + picgla_sc + dectre_sc + rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_wt_sc|AnimalId) + (0+riv_dst_sc|AnimalId) + (0+salshr_sc|AnimalId) +
                      (0+wtlnd_sc|AnimalId) + (0+picgla_sc|AnimalId) + (0+dectre_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_sum_cs,
                    control=glmmTMBControl(parallel=4),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(0, 11))),
                    map=list(theta=factor(c(NA, 1:11))))
t1c_15_wtpk = glmmTMB:::fitTMB(mod_struc)
saveRDS(t1c_15_wtpk, "model_outputs/t1c_15_wtpk.rds"); rm(t1c_15_wtpk)


t1c_15_wtpk = readRDS("model_outputs/t1c_15_wtpk.rds")
t1c_15_pack = readRDS("model_outputs/t1c_15_pack.rds")
t1c_15_full = readRDS("model_outputs/t1c_15_full.rds")
t1c_15_sum = readRDS("model_outputs/t1c_15_sum.rds")


aictab(lst(t1c_15_wtpk, t1c_15_pack, t1c_15_full, t1c_15_sum))


#######################################################################
############## TIER 2: CALF STATUS × HABITAT ##########################
#######################################################################
# Base = Tier 1 winner (additive). Test all subsets of CalfStatus ×
# {salix, wetland, spruce} interactions. CalfStatus main effect
# included whenever any interaction is present.
#######################################################################

#### SUMMER ####
t2_0_sum <- readRDS("model_outputs/t2_0_sum.rds")
pars <- getME(t2_0_sum, "theta")
betas <- fixef(t2_0_sum)$cond

# t2_0_sum: CalfStatus × all four habitats
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + salshr_sc + wtlnd_sc + picgla_sc + dectre_sc +
                      riv_dst_sc + wlf_kde_sea_sc:picgla_sc +
                      rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+salshr_sc|AnimalId) + (0+wtlnd_sc|AnimalId) +
                      (0+picgla_sc|AnimalId) + (0+dectre_sc|AnimalId) + (0+riv_dst_sc|AnimalId) +
                      (0+wlf_kde_sea_sc:picgla_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_sum_cs,
                    control=glmmTMBControl(parallel=4,
                                           optCtrl = list(iter.max = 2000, eval.max = 4000)),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(0, 12))),
                    map=list(theta=factor(c(NA, 1:12))))
t2_0_sum = glmmTMB:::fitTMB(mod_struc)
saveRDS(t2_0_sum, "model_outputs/t2_0_sum.rds"); rm(t2_0_sum)


# t2_1_sum <- readRDS("model_outputs/t2_1_sum.rds")
# pars <- getME(t2_1_sum, "theta")
# betas <- fixef(t2_1_sum)$cond

# t2_0_sum: CalfStatus × all four habitats
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + salshr_sc + wtlnd_sc + picgla_sc + dectre_sc +
                      riv_dst_sc + wlf_kde_sea_sc:wtlnd_sc + rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+salshr_sc|AnimalId) + (0+wtlnd_sc|AnimalId) +
                      (0+picgla_sc|AnimalId) + (0+dectre_sc|AnimalId) + (0+riv_dst_sc|AnimalId) +
                      (0+wlf_kde_sea_sc:wtlnd_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_sum_cs,
                    control=glmmTMBControl(parallel=4,
                                           optCtrl = list(iter.max = 2000, eval.max = 4000)),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(0, 12))),
                    map=list(theta=factor(c(NA, 1:12))))
t2_1_sum = glmmTMB:::fitTMB(mod_struc)
saveRDS(t2_1_sum, "model_outputs/t2_1_sum.rds"); rm(t2_1_sum)

# 
t2_2_sum <- readRDS("model_outputs/t2_2_sum.rds")
pars <- getME(t2_2_sum, "theta")
betas <- fixef(t2_2_sum)$cond

# t2_0_sum: CalfStatus × all four habitats
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + salshr_sc + wtlnd_sc + picgla_sc + dectre_sc +
                      riv_dst_sc +
                      wlf_kde_sea_sc:dectre_sc + rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+salshr_sc|AnimalId) + (0+wtlnd_sc|AnimalId) +
                      (0+picgla_sc|AnimalId) + (0+dectre_sc|AnimalId) + (0+riv_dst_sc|AnimalId) +
                      (0+wlf_kde_sea_sc:dectre_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_sum_cs,
                    control=glmmTMBControl(parallel=4,
                                           optCtrl = list(iter.max = 2000, eval.max = 4000)),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(0, 12))),
                    map=list(theta=factor(c(NA, 1:12))))
t2_2_sum = glmmTMB:::fitTMB(mod_struc)
saveRDS(t2_2_sum, "model_outputs/t2_2_sum.rds"); rm(t2_2_sum)


# t2_0_sum: CalfStatus × all four habitats
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + salshr_sc + wtlnd_sc + picgla_sc + dectre_sc +
                      riv_dst_sc + wlf_kde_sea_sc:wtlnd_sc +
                      wlf_kde_sea_sc:dectre_sc + rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+salshr_sc|AnimalId) + (0+wtlnd_sc|AnimalId) +
                      (0+picgla_sc|AnimalId) + (0+dectre_sc|AnimalId) + (0+riv_dst_sc|AnimalId) +
                      (0+wlf_kde_sea_sc:wtlnd_sc|AnimalId) +
                      (0+wlf_kde_sea_sc:dectre_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_sum_cs,
                    control=glmmTMBControl(parallel=4,
                                           optCtrl = list(iter.max = 2000, eval.max = 4000)),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(0, 13))),
                    map=list(theta=factor(c(NA, 1:13))))
t2_3_sum = glmmTMB:::fitTMB(mod_struc)
saveRDS(t2_3_sum, "model_outputs/t2_3_sum.rds"); rm(t2_3_sum)


# t2_0_sum: CalfStatus × all four habitats
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + salshr_sc + wtlnd_sc + picgla_sc + dectre_sc +
                      riv_dst_sc +
                      wlf_kde_sea_sc:salshr_sc + rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+salshr_sc|AnimalId) + (0+wtlnd_sc|AnimalId) +
                      (0+picgla_sc|AnimalId) + (0+dectre_sc|AnimalId) + (0+riv_dst_sc|AnimalId) +
                      (0+wlf_kde_sea_sc:salshr_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_sum_cs,
                    control=glmmTMBControl(parallel=4,
                                           optCtrl = list(iter.max = 2000, eval.max = 4000)),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(0, 12))),
                    map=list(theta=factor(c(NA, 1:12))))
t2_4_sum = glmmTMB:::fitTMB(mod_struc)
saveRDS(t2_4_sum, "model_outputs/t2_4_sum.rds"); rm(t2_4_sum)


# t2_0_sum: CalfStatus × all four habitats
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + salshr_sc + wtlnd_sc + picgla_sc + dectre_sc +
                      riv_dst_sc + wlf_kde_sea_sc:salshr_sc +
                      wlf_kde_sea_sc:wtlnd_sc + rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+salshr_sc|AnimalId) + (0+wtlnd_sc|AnimalId) +
                      (0+picgla_sc|AnimalId) + (0+dectre_sc|AnimalId) + (0+riv_dst_sc|AnimalId) +
                      (0+wlf_kde_sea_sc:salshr_sc|AnimalId) + (0+wlf_kde_sea_sc:wtlnd_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_sum_cs,
                    control=glmmTMBControl(parallel=4,
                                           optCtrl = list(iter.max = 2000, eval.max = 4000)),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(0, 13))),
                    map=list(theta=factor(c(NA, 1:13))))
t2_5_sum = glmmTMB:::fitTMB(mod_struc)
saveRDS(t2_5_sum, "model_outputs/t2_5_sum.rds"); rm(t2_5_sum)



# t2_0_sum: CalfStatus × all four habitats
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + salshr_sc + wtlnd_sc + picgla_sc + dectre_sc +
                      riv_dst_sc + wlf_kde_sea_sc:salshr_sc + wlf_kde_sea_sc:dectre_sc +
                      wlf_kde_sea_sc:wtlnd_sc + rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+salshr_sc|AnimalId) + (0+wtlnd_sc|AnimalId) +
                      (0+picgla_sc|AnimalId) + (0+dectre_sc|AnimalId) + (0+riv_dst_sc|AnimalId) +
                      (0+wlf_kde_sea_sc:salshr_sc|AnimalId) + (0+wlf_kde_sea_sc:wtlnd_sc|AnimalId) +
                      (0+wlf_kde_sea_sc:dectre_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_sum_cs,
                    control=glmmTMBControl(parallel=4,
                                           optCtrl = list(iter.max = 2000, eval.max = 4000)),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(0, 14))),
                    map=list(theta=factor(c(NA, 1:14))))
t2_6_sum = glmmTMB:::fitTMB(mod_struc)
saveRDS(t2_6_sum, "model_outputs/t2_6_sum.rds"); rm(t2_6_sum)



# t2_0_sum: CalfStatus × all four habitats
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + salshr_sc + wtlnd_sc + picgla_sc + dectre_sc +
                      riv_dst_sc + wlf_kde_sea_sc:salshr_sc + wlf_kde_sea_sc:dectre_sc + wlf_kde_sea_sc:picgla_sc +
                      wlf_kde_sea_sc:wtlnd_sc + rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+salshr_sc|AnimalId) + (0+wtlnd_sc|AnimalId) +
                      (0+picgla_sc|AnimalId) + (0+dectre_sc|AnimalId) + (0+riv_dst_sc|AnimalId) +
                      (0+wlf_kde_sea_sc:salshr_sc|AnimalId) + (0+wlf_kde_sea_sc:wtlnd_sc|AnimalId) +
                      (0+wlf_kde_sea_sc:dectre_sc|AnimalId) + (0+wlf_kde_sea_sc:picgla_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_sum_cs,
                    control=glmmTMBControl(parallel=4,
                                           optCtrl = list(iter.max = 2000, eval.max = 4000)),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(0, 15))),
                    map=list(theta=factor(c(NA, 1:15))))
t2_7_sum = glmmTMB:::fitTMB(mod_struc)
saveRDS(t2_7_sum, "model_outputs/t2_7_sum.rds"); rm(t2_7_sum)


#### WINTER ####
# t2_0_win: CalfStatus × all four habitats
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + salshr_sc + wtlnd_sc + picgla_sc + dectre_sc +
                      riv_dst_sc + wlf_kde_sea_sc:picgla_sc +
                      rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+salshr_sc|AnimalId) + (0+wtlnd_sc|AnimalId) +
                      (0+picgla_sc|AnimalId) + (0+dectre_sc|AnimalId) + (0+riv_dst_sc|AnimalId) +
                      (0+wlf_kde_sea_sc:picgla_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_win_cs,
                    control=glmmTMBControl(parallel=4,
                                           optCtrl = list(iter.max = 2000, eval.max = 4000)),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(0, 12))),
                    map=list(theta=factor(c(NA, 1:12))))
t2_0_win = glmmTMB:::fitTMB(mod_struc)
saveRDS(t2_0_win, "model_outputs/t2_0_win.rds"); rm(t2_0_win)


# t2_0_win: CalfStatus × all four habitats
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + salshr_sc + wtlnd_sc + picgla_sc + dectre_sc +
                      riv_dst_sc + wlf_kde_sea_sc:wtlnd_sc + rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+salshr_sc|AnimalId) + (0+wtlnd_sc|AnimalId) +
                      (0+picgla_sc|AnimalId) + (0+dectre_sc|AnimalId) + (0+riv_dst_sc|AnimalId) +
                      (0+wlf_kde_sea_sc:wtlnd_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_win_cs,
                    control=glmmTMBControl(parallel=4,
                                           optCtrl = list(iter.max = 2000, eval.max = 4000)),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(0, 12))),
                    map=list(theta=factor(c(NA, 1:12))))
t2_1_win = glmmTMB:::fitTMB(mod_struc)
saveRDS(t2_1_win, "model_outputs/t2_1_win.rds"); rm(t2_1_win)


# t2_0_win: CalfStatus × all four habitats
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + salshr_sc + wtlnd_sc + picgla_sc + dectre_sc +
                      riv_dst_sc +
                      wlf_kde_sea_sc:dectre_sc + rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+salshr_sc|AnimalId) + (0+wtlnd_sc|AnimalId) +
                      (0+picgla_sc|AnimalId) + (0+dectre_sc|AnimalId) + (0+riv_dst_sc|AnimalId) +
                      (0+wlf_kde_sea_sc:dectre_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_win_cs,
                    control=glmmTMBControl(parallel=4,
                                           optCtrl = list(iter.max = 2000, eval.max = 4000)),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(0, 12))),
                    map=list(theta=factor(c(NA, 1:12))))
t2_2_win = glmmTMB:::fitTMB(mod_struc)
saveRDS(t2_2_win, "model_outputs/t2_2_win.rds"); rm(t2_2_win)


# t2_0_win: CalfStatus × all four habitats
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + salshr_sc + wtlnd_sc + picgla_sc + dectre_sc +
                      riv_dst_sc + wlf_kde_sea_sc:wtlnd_sc +
                      wlf_kde_sea_sc:dectre_sc + rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+salshr_sc|AnimalId) + (0+wtlnd_sc|AnimalId) +
                      (0+picgla_sc|AnimalId) + (0+dectre_sc|AnimalId) + (0+riv_dst_sc|AnimalId) +
                      (0+wlf_kde_sea_sc:wtlnd_sc|AnimalId) +
                      (0+wlf_kde_sea_sc:dectre_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_win_cs,
                    control=glmmTMBControl(parallel=4,
                                           optCtrl = list(iter.max = 2000, eval.max = 4000)),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(0, 13))),
                    map=list(theta=factor(c(NA, 1:13))))
t2_3_win = glmmTMB:::fitTMB(mod_struc)
saveRDS(t2_3_win, "model_outputs/t2_3_win.rds"); rm(t2_3_win)


# t2_0_win: CalfStatus × all four habitats
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + salshr_sc + wtlnd_sc + picgla_sc + dectre_sc +
                      riv_dst_sc +
                      wlf_kde_sea_sc:salshr_sc + rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+salshr_sc|AnimalId) + (0+wtlnd_sc|AnimalId) +
                      (0+picgla_sc|AnimalId) + (0+dectre_sc|AnimalId) + (0+riv_dst_sc|AnimalId) +
                      (0+wlf_kde_sea_sc:salshr_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_win_cs,
                    control=glmmTMBControl(parallel=4,
                                           optCtrl = list(iter.max = 2000, eval.max = 4000)),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(0, 12))),
                    map=list(theta=factor(c(NA, 1:12))))
t2_4_win = glmmTMB:::fitTMB(mod_struc)
saveRDS(t2_4_win, "model_outputs/t2_4_win.rds"); rm(t2_4_win)


# t2_0_win: CalfStatus × all four habitats
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + salshr_sc + wtlnd_sc + picgla_sc + dectre_sc +
                      riv_dst_sc + wlf_kde_sea_sc:salshr_sc +
                      wlf_kde_sea_sc:wtlnd_sc + rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+salshr_sc|AnimalId) + (0+wtlnd_sc|AnimalId) +
                      (0+picgla_sc|AnimalId) + (0+dectre_sc|AnimalId) + (0+riv_dst_sc|AnimalId) +
                      (0+wlf_kde_sea_sc:salshr_sc|AnimalId) + (0+wlf_kde_sea_sc:wtlnd_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_win_cs,
                    control=glmmTMBControl(parallel=4,
                                           optCtrl = list(iter.max = 2000, eval.max = 4000)),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(0, 13))),
                    map=list(theta=factor(c(NA, 1:13))))
t2_5_win = glmmTMB:::fitTMB(mod_struc)
saveRDS(t2_5_win, "model_outputs/t2_5_win.rds"); rm(t2_5_win)



# t2_0_win: CalfStatus × all four habitats
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + salshr_sc + wtlnd_sc + picgla_sc + dectre_sc +
                      riv_dst_sc + wlf_kde_sea_sc:salshr_sc + wlf_kde_sea_sc:dectre_sc +
                      wlf_kde_sea_sc:wtlnd_sc + rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+salshr_sc|AnimalId) + (0+wtlnd_sc|AnimalId) +
                      (0+picgla_sc|AnimalId) + (0+dectre_sc|AnimalId) + (0+riv_dst_sc|AnimalId) +
                      (0+wlf_kde_sea_sc:salshr_sc|AnimalId) + (0+wlf_kde_sea_sc:wtlnd_sc|AnimalId) +
                      (0+wlf_kde_sea_sc:dectre_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_win_cs,
                    control=glmmTMBControl(parallel=4,
                                           optCtrl = list(iter.max = 2000, eval.max = 4000)),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(0, 14))),
                    map=list(theta=factor(c(NA, 1:14))))
t2_6_win = glmmTMB:::fitTMB(mod_struc)
saveRDS(t2_6_win, "model_outputs/t2_6_win.rds"); rm(t2_6_win)


# t2_0_win: CalfStatus × all four habitats
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + salshr_sc + wtlnd_sc + picgla_sc + dectre_sc +
                      riv_dst_sc + wlf_kde_sea_sc:salshr_sc + wlf_kde_sea_sc:dectre_sc + wlf_kde_sea_sc:picgla_sc +
                      wlf_kde_sea_sc:wtlnd_sc + rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+salshr_sc|AnimalId) + (0+wtlnd_sc|AnimalId) +
                      (0+picgla_sc|AnimalId) + (0+dectre_sc|AnimalId) + (0+riv_dst_sc|AnimalId) +
                      (0+wlf_kde_sea_sc:salshr_sc|AnimalId) + (0+wlf_kde_sea_sc:wtlnd_sc|AnimalId) +
                      (0+wlf_kde_sea_sc:dectre_sc|AnimalId) + (0+wlf_kde_sea_sc:picgla_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_win_cs,
                    control=glmmTMBControl(parallel=4,
                                           optCtrl = list(iter.max = 2000, eval.max = 4000)),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(0, 15))),
                    map=list(theta=factor(c(NA, 1:15))))
t2_7_win = glmmTMB:::fitTMB(mod_struc)
saveRDS(t2_7_win, "model_outputs/t2_7_win.rds"); rm(t2_7_win)


t2_sum_nms = paste0("t2_", 0:7, "_sum")
for(nm in t2_sum_nms) assign(nm, readRDS(paste0("model_outputs/", nm, ".rds")))
aictab(mget(t2_sum_nms))
rm(list = t2_sum_nms)


t2_win_nms = paste0("t2_", 0:7, "_win")
for(nm in t2_win_nms) assign(nm, readRDS(paste0("model_outputs/", nm, ".rds")))
aictab(mget(t2_win_nms))
rm(list = t2_win_nms)


#######################################################################
####### TIER 3: CALF STATUS × WOLF DENSITY x HABITAT ##################
#######################################################################
# Base = Tier 2 winner (two way interactions). Test all subsets of 
# CalfStatus × WolfDensity interactions. CalfStatus and WolfDensity main
# effect included whenever any interaction is present.
#######################################################################

######## tier 3 ##########
# # t3_7_sum: Global — three-way on dectre and wetlnd
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + salshr_sc + wtlnd_sc + picgla_sc +
                      riv_dst_sc + dectre_sc + CalfStatus + CalfStatus:wtlnd_sc + CalfStatus:dectre_sc +
                      wtlnd_sc:wlf_kde_sea_sc + CalfStatus:wlf_kde_sea_sc + dectre_sc:wlf_kde_sea_sc +
                      CalfStatus:wtlnd_sc:wlf_kde_sea_sc + CalfStatus:dectre_sc:wlf_kde_sea_sc +
                      rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+salshr_sc|AnimalId) + (0+wtlnd_sc|AnimalId) + (0+dectre_sc|AnimalId) +
                      (0+picgla_sc|AnimalId) + (0+riv_dst_sc|AnimalId) + (0+CalfStatus|AnimalId) +
                      (0+CalfStatus:wtlnd_sc|AnimalId) + (0+CalfStatus:dectre_sc|AnimalId) +
                      (0+wtlnd_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+CalfStatus:wlf_kde_sea_sc|AnimalId) + (0+dectre_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+CalfStatus:wtlnd_sc:wlf_kde_sea_sc|AnimalId) + (0+CalfStatus:dectre_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_sum_cs,
                    control=glmmTMBControl(parallel=4,
                                           optCtrl = list(iter.max = 2000, eval.max = 4000)),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(0, 19))),
                    map=list(theta=factor(c(NA, 1:19))))
t3_0_sum = glmmTMB:::fitTMB(mod_struc)
saveRDS(t3_0_sum, "model_outputs/t3_0_sum.rds"); rm(t3_0_sum)


# # t3_7_sum: Global — three-way on all three habitats
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + salshr_sc + wtlnd_sc + picgla_sc +
                      riv_dst_sc + dectre_sc + CalfStatus + CalfStatus:dectre_sc +
                      wtlnd_sc:wlf_kde_sea_sc + CalfStatus:wlf_kde_sea_sc + dectre_sc:wlf_kde_sea_sc +
                      CalfStatus:dectre_sc:wlf_kde_sea_sc +
                      rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+salshr_sc|AnimalId) + (0+wtlnd_sc|AnimalId) + (0+dectre_sc|AnimalId) +
                      (0+picgla_sc|AnimalId) + (0+riv_dst_sc|AnimalId) + (0+CalfStatus|AnimalId) +
                      (0+CalfStatus:dectre_sc|AnimalId) +
                      (0+wtlnd_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+CalfStatus:wlf_kde_sea_sc|AnimalId) + (0+dectre_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+CalfStatus:dectre_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_sum_cs,
                    control=glmmTMBControl(parallel=4,
                                           optCtrl = list(iter.max = 2000, eval.max = 4000)),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(0, 17))),
                    map=list(theta=factor(c(NA, 1:17))))
t3_1_sum = glmmTMB:::fitTMB(mod_struc)
saveRDS(t3_1_sum, "model_outputs/t3_1_sum.rds"); rm(t3_1_sum)


# # t3_7_sum: Global — three-way on all three habitats
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + salshr_sc + wtlnd_sc + picgla_sc +
                      riv_dst_sc + dectre_sc + CalfStatus + CalfStatus:wtlnd_sc +
                      wtlnd_sc:wlf_kde_sea_sc + CalfStatus:wlf_kde_sea_sc + dectre_sc:wlf_kde_sea_sc +
                      CalfStatus:wtlnd_sc:wlf_kde_sea_sc +
                      rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+salshr_sc|AnimalId) + (0+wtlnd_sc|AnimalId) + (0+dectre_sc|AnimalId) +
                      (0+picgla_sc|AnimalId) + (0+riv_dst_sc|AnimalId) + (0+CalfStatus|AnimalId) +
                      (0+CalfStatus:wtlnd_sc|AnimalId) +
                      (0+wtlnd_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+CalfStatus:wlf_kde_sea_sc|AnimalId) + (0+dectre_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+CalfStatus:wtlnd_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_sum_cs,
                    control=glmmTMBControl(parallel=4,
                                           optCtrl = list(iter.max = 2000, eval.max = 4000)),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(0, 17))),
                    map=list(theta=factor(c(NA, 1:17))))
t3_2_sum = glmmTMB:::fitTMB(mod_struc)
saveRDS(t3_2_sum, "model_outputs/t3_2_sum.rds"); rm(t3_2_sum)


# # t3_7_sum: Global — three-way on all three habitats
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + salshr_sc + wtlnd_sc + picgla_sc +
                      riv_dst_sc + dectre_sc + CalfStatus + CalfStatus:dectre_sc +
                      CalfStatus:wlf_kde_sea_sc + dectre_sc:wlf_kde_sea_sc +
                      CalfStatus:dectre_sc:wlf_kde_sea_sc +
                      rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+salshr_sc|AnimalId) + (0+wtlnd_sc|AnimalId) + (0+dectre_sc|AnimalId) +
                      (0+picgla_sc|AnimalId) + (0+riv_dst_sc|AnimalId) + (0+CalfStatus|AnimalId) +
                      (0+CalfStatus:dectre_sc|AnimalId) +
                      (0+CalfStatus:wlf_kde_sea_sc|AnimalId) + (0+dectre_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+CalfStatus:dectre_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_sum_cs,
                    control=glmmTMBControl(parallel=4,
                                           optCtrl = list(iter.max = 2000, eval.max = 4000)),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(0, 16))),
                    map=list(theta=factor(c(NA, 1:16))))
t3_3_sum = glmmTMB:::fitTMB(mod_struc)
saveRDS(t3_3_sum, "model_outputs/t3_3_sum.rds"); rm(t3_3_sum)


# # t3_7_sum: Global — three-way on all three habitats
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + salshr_sc + wtlnd_sc + picgla_sc +
                      riv_dst_sc + dectre_sc + CalfStatus + CalfStatus:wtlnd_sc +
                      wtlnd_sc:wlf_kde_sea_sc + CalfStatus:wlf_kde_sea_sc +
                      CalfStatus:wtlnd_sc:wlf_kde_sea_sc +
                      rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+salshr_sc|AnimalId) + (0+wtlnd_sc|AnimalId) + (0+dectre_sc|AnimalId) +
                      (0+picgla_sc|AnimalId) + (0+riv_dst_sc|AnimalId) + (0+CalfStatus|AnimalId) +
                      (0+CalfStatus:wtlnd_sc|AnimalId) +
                      (0+wtlnd_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+CalfStatus:wlf_kde_sea_sc|AnimalId) +
                      (0+CalfStatus:wtlnd_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_sum_cs,
                    control=glmmTMBControl(parallel=4,
                                           optCtrl = list(iter.max = 2000, eval.max = 4000)),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(0, 16))),
                    map=list(theta=factor(c(NA, 1:16))))
t3_4_sum = glmmTMB:::fitTMB(mod_struc)
saveRDS(t3_4_sum, "model_outputs/t3_4_sum.rds"); rm(t3_4_sum)




# # t3_7_sum: Global — three-way on all three habitats
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + salshr_sc + wtlnd_sc + picgla_sc +
                      riv_dst_sc + dectre_sc + CalfStatus +
                      CalfStatus:salshr_sc + CalfStatus:wtlnd_sc + CalfStatus:picgla_sc + CalfStatus:dectre_sc +
                      salshr_sc:wlf_kde_sea_sc + wtlnd_sc:wlf_kde_sea_sc + picgla_sc:wlf_kde_sea_sc +
                      CalfStatus:wlf_kde_sea_sc + dectre_sc:wlf_kde_sea_sc +
                      CalfStatus:salshr_sc:wlf_kde_sea_sc + CalfStatus:wtlnd_sc:wlf_kde_sea_sc + CalfStatus:picgla_sc:wlf_kde_sea_sc +
                      CalfStatus:dectre_sc:wlf_kde_sea_sc +
                      rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+salshr_sc|AnimalId) + (0+wtlnd_sc|AnimalId) + (0+dectre_sc|AnimalId) +
                      (0+picgla_sc|AnimalId) + (0+riv_dst_sc|AnimalId) + (0+CalfStatus|AnimalId) +
                      (0+CalfStatus:salshr_sc|AnimalId) + (0+CalfStatus:wtlnd_sc|AnimalId) + (0+CalfStatus:picgla_sc|AnimalId) + (0+CalfStatus:dectre_sc|AnimalId) +
                      (0+salshr_sc:wlf_kde_sea_sc|AnimalId) + (0+wtlnd_sc:wlf_kde_sea_sc|AnimalId) + (0+picgla_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+CalfStatus:wlf_kde_sea_sc|AnimalId) + (0+dectre_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+CalfStatus:salshr_sc:wlf_kde_sea_sc|AnimalId) + (0+CalfStatus:wtlnd_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+CalfStatus:picgla_sc:wlf_kde_sea_sc|AnimalId) + (0+CalfStatus:dectre_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_sum_cs,
                    control=glmmTMBControl(parallel=4,
                                           optCtrl = list(iter.max = 2000, eval.max = 4000)),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(0, 25))),
                    map=list(theta=factor(c(NA, 1:25))))
t3_5_sum = glmmTMB:::fitTMB(mod_struc)
saveRDS(t3_5_sum, "model_outputs/t3_5_sum.rds"); rm(t3_5_sum)



# # t3_7_sum: Global — three-way on all three habitats
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + salshr_sc + wtlnd_sc + picgla_sc +
                      riv_dst_sc + dectre_sc + CalfStatus +
                      CalfStatus:salshr_sc + CalfStatus:wtlnd_sc + CalfStatus:dectre_sc +
                      salshr_sc:wlf_kde_sea_sc + wtlnd_sc:wlf_kde_sea_sc + picgla_sc:wlf_kde_sea_sc +
                      CalfStatus:wlf_kde_sea_sc + dectre_sc:wlf_kde_sea_sc +
                      CalfStatus:salshr_sc:wlf_kde_sea_sc + CalfStatus:wtlnd_sc:wlf_kde_sea_sc +
                      CalfStatus:dectre_sc:wlf_kde_sea_sc +
                      rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+salshr_sc|AnimalId) + (0+wtlnd_sc|AnimalId) + (0+dectre_sc|AnimalId) +
                      (0+picgla_sc|AnimalId) + (0+riv_dst_sc|AnimalId) + (0+CalfStatus|AnimalId) +
                      (0+CalfStatus:salshr_sc|AnimalId) + (0+CalfStatus:wtlnd_sc|AnimalId) + (0+CalfStatus:dectre_sc|AnimalId) +
                      (0+salshr_sc:wlf_kde_sea_sc|AnimalId) + (0+wtlnd_sc:wlf_kde_sea_sc|AnimalId) + (0+picgla_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+CalfStatus:wlf_kde_sea_sc|AnimalId) + (0+dectre_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+CalfStatus:salshr_sc:wlf_kde_sea_sc|AnimalId) + (0+CalfStatus:wtlnd_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+CalfStatus:dectre_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_sum_cs,
                    control=glmmTMBControl(parallel=4,
                                           optCtrl = list(iter.max = 2000, eval.max = 4000)),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(0, 23))),
                    map=list(theta=factor(c(NA, 1:23))))
t3_6_sum = glmmTMB:::fitTMB(mod_struc)
saveRDS(t3_6_sum, "model_outputs/t3_6_sum.rds"); rm(t3_6_sum)


# # t3_7_sum: Global — three-way on all three habitats
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + salshr_sc + wtlnd_sc + picgla_sc +
                      riv_dst_sc + dectre_sc + CalfStatus +
                      CalfStatus:salshr_sc + CalfStatus:wtlnd_sc + CalfStatus:dectre_sc +
                      salshr_sc:wlf_kde_sea_sc + wtlnd_sc:wlf_kde_sea_sc + picgla_sc:wlf_kde_sea_sc +
                      CalfStatus:wlf_kde_sea_sc + dectre_sc:wlf_kde_sea_sc +
                      CalfStatus:wtlnd_sc:wlf_kde_sea_sc +
                      CalfStatus:dectre_sc:wlf_kde_sea_sc +
                      rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+salshr_sc|AnimalId) + (0+wtlnd_sc|AnimalId) + (0+dectre_sc|AnimalId) +
                      (0+picgla_sc|AnimalId) + (0+riv_dst_sc|AnimalId) + (0+CalfStatus|AnimalId) +
                      (0+CalfStatus:salshr_sc|AnimalId) + (0+CalfStatus:wtlnd_sc|AnimalId) + (0+CalfStatus:dectre_sc|AnimalId) +
                      (0+salshr_sc:wlf_kde_sea_sc|AnimalId) + (0+wtlnd_sc:wlf_kde_sea_sc|AnimalId) + (0+picgla_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+CalfStatus:wlf_kde_sea_sc|AnimalId) + (0+dectre_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+CalfStatus:wtlnd_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+CalfStatus:dectre_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_sum_cs,
                    control=glmmTMBControl(parallel=4,
                                           optCtrl = list(iter.max = 2000, eval.max = 4000)),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(0, 22))),
                    map=list(theta=factor(c(NA, 1:22))))
t3_7_sum = glmmTMB:::fitTMB(mod_struc)
saveRDS(t3_7_sum, "model_outputs/t3_7_sum.rds"); rm(t3_7_sum)


# # t3_7_sum: Global — three-way on all three habitats
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + salshr_sc + wtlnd_sc + picgla_sc +
                      riv_dst_sc + dectre_sc + CalfStatus +
                      CalfStatus:salshr_sc + CalfStatus:wtlnd_sc + CalfStatus:dectre_sc +
                      wtlnd_sc:wlf_kde_sea_sc +
                      CalfStatus:wlf_kde_sea_sc + dectre_sc:wlf_kde_sea_sc +
                      CalfStatus:wtlnd_sc:wlf_kde_sea_sc +
                      CalfStatus:dectre_sc:wlf_kde_sea_sc +
                      rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+salshr_sc|AnimalId) + (0+wtlnd_sc|AnimalId) + (0+dectre_sc|AnimalId) +
                      (0+picgla_sc|AnimalId) + (0+riv_dst_sc|AnimalId) + (0+CalfStatus|AnimalId) +
                      (0+CalfStatus:salshr_sc|AnimalId) + (0+CalfStatus:wtlnd_sc|AnimalId) + (0+CalfStatus:dectre_sc|AnimalId) +
                      (0+wtlnd_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+CalfStatus:wlf_kde_sea_sc|AnimalId) + (0+dectre_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+CalfStatus:wtlnd_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+CalfStatus:dectre_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_sum_cs,
                    control=glmmTMBControl(parallel=4,
                                           optCtrl = list(iter.max = 2000, eval.max = 4000)),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(0, 20))),
                    map=list(theta=factor(c(NA, 1:20))))
t3_8_sum = glmmTMB:::fitTMB(mod_struc)
saveRDS(t3_8_sum, "model_outputs/t3_8_sum.rds"); rm(t3_8_sum)





######## tier 3 ##########
# # t3_7_win: Global — three-way on dectre and wetlnd
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + salshr_sc + wtlnd_sc + picgla_sc +
                      riv_dst_sc + dectre_sc + CalfStatus + CalfStatus:wtlnd_sc + CalfStatus:dectre_sc +
                      wtlnd_sc:wlf_kde_sea_sc + CalfStatus:wlf_kde_sea_sc + dectre_sc:wlf_kde_sea_sc +
                      CalfStatus:wtlnd_sc:wlf_kde_sea_sc + CalfStatus:dectre_sc:wlf_kde_sea_sc +
                      rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+salshr_sc|AnimalId) + (0+wtlnd_sc|AnimalId) + (0+dectre_sc|AnimalId) +
                      (0+picgla_sc|AnimalId) + (0+riv_dst_sc|AnimalId) + (0+CalfStatus|AnimalId) +
                      (0+CalfStatus:wtlnd_sc|AnimalId) + (0+CalfStatus:dectre_sc|AnimalId) +
                      (0+wtlnd_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+CalfStatus:wlf_kde_sea_sc|AnimalId) + (0+dectre_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+CalfStatus:wtlnd_sc:wlf_kde_sea_sc|AnimalId) + (0+CalfStatus:dectre_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_win_cs,
                    control=glmmTMBControl(parallel=4,
                                           optCtrl = list(iter.max = 2000, eval.max = 4000)),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(0, 19))),
                    map=list(theta=factor(c(NA, 1:19))))
t3_0_win = glmmTMB:::fitTMB(mod_struc)
saveRDS(t3_0_win, "model_outputs/t3_0_win.rds"); rm(t3_0_win)


# # t3_7_win: Global — three-way on all three habitats
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + salshr_sc + wtlnd_sc + picgla_sc +
                      riv_dst_sc + dectre_sc + CalfStatus + CalfStatus:dectre_sc +
                      wtlnd_sc:wlf_kde_sea_sc + CalfStatus:wlf_kde_sea_sc + dectre_sc:wlf_kde_sea_sc +
                      CalfStatus:dectre_sc:wlf_kde_sea_sc +
                      rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+salshr_sc|AnimalId) + (0+wtlnd_sc|AnimalId) + (0+dectre_sc|AnimalId) +
                      (0+picgla_sc|AnimalId) + (0+riv_dst_sc|AnimalId) + (0+CalfStatus|AnimalId) +
                      (0+CalfStatus:dectre_sc|AnimalId) +
                      (0+wtlnd_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+CalfStatus:wlf_kde_sea_sc|AnimalId) + (0+dectre_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+CalfStatus:dectre_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_win_cs,
                    control=glmmTMBControl(parallel=4,
                                           optCtrl = list(iter.max = 2000, eval.max = 4000)),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(0, 17))),
                    map=list(theta=factor(c(NA, 1:17))))
t3_1_win = glmmTMB:::fitTMB(mod_struc)
saveRDS(t3_1_win, "model_outputs/t3_1_win.rds"); rm(t3_1_win)


# # t3_7_win: Global — three-way on all three habitats
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + salshr_sc + wtlnd_sc + picgla_sc +
                      riv_dst_sc + dectre_sc + CalfStatus + CalfStatus:wtlnd_sc +
                      wtlnd_sc:wlf_kde_sea_sc + CalfStatus:wlf_kde_sea_sc + dectre_sc:wlf_kde_sea_sc +
                      CalfStatus:wtlnd_sc:wlf_kde_sea_sc +
                      rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+salshr_sc|AnimalId) + (0+wtlnd_sc|AnimalId) + (0+dectre_sc|AnimalId) +
                      (0+picgla_sc|AnimalId) + (0+riv_dst_sc|AnimalId) + (0+CalfStatus|AnimalId) +
                      (0+CalfStatus:wtlnd_sc|AnimalId) +
                      (0+wtlnd_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+CalfStatus:wlf_kde_sea_sc|AnimalId) + (0+dectre_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+CalfStatus:wtlnd_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_win_cs,
                    control=glmmTMBControl(parallel=4,
                                           optCtrl = list(iter.max = 2000, eval.max = 4000)),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(0, 17))),
                    map=list(theta=factor(c(NA, 1:17))))
t3_2_win = glmmTMB:::fitTMB(mod_struc)
saveRDS(t3_2_win, "model_outputs/t3_2_win.rds"); rm(t3_2_win)


# # t3_7_win: Global — three-way on all three habitats
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + salshr_sc + wtlnd_sc + picgla_sc +
                      riv_dst_sc + dectre_sc + CalfStatus + CalfStatus:dectre_sc +
                      CalfStatus:wlf_kde_sea_sc + dectre_sc:wlf_kde_sea_sc +
                      CalfStatus:dectre_sc:wlf_kde_sea_sc +
                      rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+salshr_sc|AnimalId) + (0+wtlnd_sc|AnimalId) + (0+dectre_sc|AnimalId) +
                      (0+picgla_sc|AnimalId) + (0+riv_dst_sc|AnimalId) + (0+CalfStatus|AnimalId) +
                      (0+CalfStatus:dectre_sc|AnimalId) +
                      (0+CalfStatus:wlf_kde_sea_sc|AnimalId) + (0+dectre_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+CalfStatus:dectre_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_win_cs,
                    control=glmmTMBControl(parallel=4,
                                           optCtrl = list(iter.max = 2000, eval.max = 4000)),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(0, 16))),
                    map=list(theta=factor(c(NA, 1:16))))
t3_3_win = glmmTMB:::fitTMB(mod_struc)
saveRDS(t3_3_win, "model_outputs/t3_3_win.rds"); rm(t3_3_win)


# # t3_7_win: Global — three-way on all three habitats
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + salshr_sc + wtlnd_sc + picgla_sc +
                      riv_dst_sc + dectre_sc + CalfStatus + CalfStatus:wtlnd_sc +
                      wtlnd_sc:wlf_kde_sea_sc + CalfStatus:wlf_kde_sea_sc +
                      CalfStatus:wtlnd_sc:wlf_kde_sea_sc +
                      rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+salshr_sc|AnimalId) + (0+wtlnd_sc|AnimalId) + (0+dectre_sc|AnimalId) +
                      (0+picgla_sc|AnimalId) + (0+riv_dst_sc|AnimalId) + (0+CalfStatus|AnimalId) +
                      (0+CalfStatus:wtlnd_sc|AnimalId) +
                      (0+wtlnd_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+CalfStatus:wlf_kde_sea_sc|AnimalId) +
                      (0+CalfStatus:wtlnd_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_win_cs,
                    control=glmmTMBControl(parallel=4,
                                           optCtrl = list(iter.max = 2000, eval.max = 4000)),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(0, 16))),
                    map=list(theta=factor(c(NA, 1:16))))
t3_4_win = glmmTMB:::fitTMB(mod_struc)
saveRDS(t3_4_win, "model_outputs/t3_4_win.rds"); rm(t3_4_win)




# # t3_7_win: Global — three-way on all three habitats
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + salshr_sc + wtlnd_sc + picgla_sc +
                      riv_dst_sc + dectre_sc + CalfStatus +
                      CalfStatus:salshr_sc + CalfStatus:wtlnd_sc + CalfStatus:picgla_sc + CalfStatus:dectre_sc +
                      salshr_sc:wlf_kde_sea_sc + wtlnd_sc:wlf_kde_sea_sc + picgla_sc:wlf_kde_sea_sc +
                      CalfStatus:wlf_kde_sea_sc + dectre_sc:wlf_kde_sea_sc +
                      CalfStatus:salshr_sc:wlf_kde_sea_sc + CalfStatus:wtlnd_sc:wlf_kde_sea_sc + CalfStatus:picgla_sc:wlf_kde_sea_sc +
                      CalfStatus:dectre_sc:wlf_kde_sea_sc +
                      rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+salshr_sc|AnimalId) + (0+wtlnd_sc|AnimalId) + (0+dectre_sc|AnimalId) +
                      (0+picgla_sc|AnimalId) + (0+riv_dst_sc|AnimalId) + (0+CalfStatus|AnimalId) +
                      (0+CalfStatus:salshr_sc|AnimalId) + (0+CalfStatus:wtlnd_sc|AnimalId) + (0+CalfStatus:picgla_sc|AnimalId) + (0+CalfStatus:dectre_sc|AnimalId) +
                      (0+salshr_sc:wlf_kde_sea_sc|AnimalId) + (0+wtlnd_sc:wlf_kde_sea_sc|AnimalId) + (0+picgla_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+CalfStatus:wlf_kde_sea_sc|AnimalId) + (0+dectre_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+CalfStatus:salshr_sc:wlf_kde_sea_sc|AnimalId) + (0+CalfStatus:wtlnd_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+CalfStatus:picgla_sc:wlf_kde_sea_sc|AnimalId) + (0+CalfStatus:dectre_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_win_cs,
                    control=glmmTMBControl(parallel=4,
                                           optCtrl = list(iter.max = 2000, eval.max = 4000)),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(0, 25))),
                    map=list(theta=factor(c(NA, 1:25))))
t3_5_win = glmmTMB:::fitTMB(mod_struc)
saveRDS(t3_5_win, "model_outputs/t3_5_win.rds"); rm(t3_5_win)



# # t3_7_win: Global — three-way on all three habitats
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + salshr_sc + wtlnd_sc + picgla_sc +
                      riv_dst_sc + dectre_sc + CalfStatus +
                      CalfStatus:salshr_sc + CalfStatus:wtlnd_sc + CalfStatus:dectre_sc +
                      salshr_sc:wlf_kde_sea_sc + wtlnd_sc:wlf_kde_sea_sc + picgla_sc:wlf_kde_sea_sc +
                      CalfStatus:wlf_kde_sea_sc + dectre_sc:wlf_kde_sea_sc +
                      CalfStatus:salshr_sc:wlf_kde_sea_sc + CalfStatus:wtlnd_sc:wlf_kde_sea_sc +
                      CalfStatus:dectre_sc:wlf_kde_sea_sc +
                      rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+salshr_sc|AnimalId) + (0+wtlnd_sc|AnimalId) + (0+dectre_sc|AnimalId) +
                      (0+picgla_sc|AnimalId) + (0+riv_dst_sc|AnimalId) + (0+CalfStatus|AnimalId) +
                      (0+CalfStatus:salshr_sc|AnimalId) + (0+CalfStatus:wtlnd_sc|AnimalId) + (0+CalfStatus:dectre_sc|AnimalId) +
                      (0+salshr_sc:wlf_kde_sea_sc|AnimalId) + (0+wtlnd_sc:wlf_kde_sea_sc|AnimalId) + (0+picgla_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+CalfStatus:wlf_kde_sea_sc|AnimalId) + (0+dectre_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+CalfStatus:salshr_sc:wlf_kde_sea_sc|AnimalId) + (0+CalfStatus:wtlnd_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+CalfStatus:dectre_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_win_cs,
                    control=glmmTMBControl(parallel=4,
                                           optCtrl = list(iter.max = 2000, eval.max = 4000)),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(0, 23))),
                    map=list(theta=factor(c(NA, 1:23))))
t3_6_win = glmmTMB:::fitTMB(mod_struc)
saveRDS(t3_6_win, "model_outputs/t3_6_win.rds"); rm(t3_6_win)


# # t3_7_win: Global — three-way on all three habitats
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + salshr_sc + wtlnd_sc + picgla_sc +
                      riv_dst_sc + dectre_sc + CalfStatus +
                      CalfStatus:salshr_sc + CalfStatus:wtlnd_sc + CalfStatus:dectre_sc +
                      salshr_sc:wlf_kde_sea_sc + wtlnd_sc:wlf_kde_sea_sc + picgla_sc:wlf_kde_sea_sc +
                      CalfStatus:wlf_kde_sea_sc + dectre_sc:wlf_kde_sea_sc +
                      CalfStatus:wtlnd_sc:wlf_kde_sea_sc +
                      CalfStatus:dectre_sc:wlf_kde_sea_sc +
                      rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+salshr_sc|AnimalId) + (0+wtlnd_sc|AnimalId) + (0+dectre_sc|AnimalId) +
                      (0+picgla_sc|AnimalId) + (0+riv_dst_sc|AnimalId) + (0+CalfStatus|AnimalId) +
                      (0+CalfStatus:salshr_sc|AnimalId) + (0+CalfStatus:wtlnd_sc|AnimalId) + (0+CalfStatus:dectre_sc|AnimalId) +
                      (0+salshr_sc:wlf_kde_sea_sc|AnimalId) + (0+wtlnd_sc:wlf_kde_sea_sc|AnimalId) + (0+picgla_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+CalfStatus:wlf_kde_sea_sc|AnimalId) + (0+dectre_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+CalfStatus:wtlnd_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+CalfStatus:dectre_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_win_cs,
                    control=glmmTMBControl(parallel=4,
                                           optCtrl = list(iter.max = 2000, eval.max = 4000)),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(0, 22))),
                    map=list(theta=factor(c(NA, 1:22))))
t3_7_win = glmmTMB:::fitTMB(mod_struc)
saveRDS(t3_7_win, "model_outputs/t3_7_win.rds"); rm(t3_7_win)


# # t3_7_win: Global — three-way on all three habitats
mod_struc = glmmTMB(Type ~ log_sl + sl_sc + cos_ta + wlf_kde_sea_sc + salshr_sc + wtlnd_sc + picgla_sc +
                      riv_dst_sc + dectre_sc + CalfStatus +
                      CalfStatus:salshr_sc + CalfStatus:wtlnd_sc + CalfStatus:dectre_sc +
                      wtlnd_sc:wlf_kde_sea_sc +
                      CalfStatus:wlf_kde_sea_sc + dectre_sc:wlf_kde_sea_sc +
                      CalfStatus:wtlnd_sc:wlf_kde_sea_sc +
                      CalfStatus:dectre_sc:wlf_kde_sea_sc +
                      rugg_sc + I(rugg_sc^2) +
                      (1|Step) + (0+log_sl|AnimalId) + (0+sl_sc|AnimalId) + (0+cos_ta|AnimalId) +
                      (0+wlf_kde_sea_sc|AnimalId) + (0+salshr_sc|AnimalId) + (0+wtlnd_sc|AnimalId) + (0+dectre_sc|AnimalId) +
                      (0+picgla_sc|AnimalId) + (0+riv_dst_sc|AnimalId) + (0+CalfStatus|AnimalId) +
                      (0+CalfStatus:salshr_sc|AnimalId) + (0+CalfStatus:wtlnd_sc|AnimalId) + (0+CalfStatus:dectre_sc|AnimalId) +
                      (0+wtlnd_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+CalfStatus:wlf_kde_sea_sc|AnimalId) + (0+dectre_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+CalfStatus:wtlnd_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+CalfStatus:dectre_sc:wlf_kde_sea_sc|AnimalId) +
                      (0+rugg_sc|AnimalId) + (0+I(rugg_sc^2)|AnimalId),
                    family=poisson, data=mse_ssf_win_cs,
                    control=glmmTMBControl(parallel=4,
                                           optCtrl = list(iter.max = 2000, eval.max = 4000)),
                    doFit=FALSE, start=list(theta=c(log(1e3), rep(0, 20))),
                    map=list(theta=factor(c(NA, 1:20))))
t3_8_win = glmmTMB:::fitTMB(mod_struc)
saveRDS(t3_8_win, "model_outputs/t3_8_win.rds"); rm(t3_8_win)


##
t3_sum_nms = paste0("t3_", 0:8, "_sum")
for(nm in t3_sum_nms) assign(nm, readRDS(paste0("model_outputs/", nm, ".rds")))
aictab(mget(t3_sum_nms))
rm(list = t3_sum_nms)


t3_win_nms = paste0("t3_", 0:8, "_win")
for(nm in t3_win_nms) assign(nm, readRDS(paste0("model_outputs/", nm, ".rds")))
aictab(mget(t3_win_nms))
rm(list = t3_win_nms)


