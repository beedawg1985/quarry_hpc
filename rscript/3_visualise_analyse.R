setwd('/home/barneyharris/projects/quarry')
source(paste0(getwd(),'/rscript/general_functions.R'))

# load('data/intpols.RDS')

# results <- list()
# results_offset_nonoise <- analyseDat(filePattern = 'march_offset_nonoise_nodiffs')
# results_offset_50noise <- analyseDat(filePattern = 'intdat_march_offset_wnoise')
# results$`diffDat_feb28_nodiffs` <- analyseDat(filePattern = 'diffDat_feb28_nodiffs')
# results_nooffset_50noise <- analyseDat(filePattern = 'diffDat_feb28_nodiffs')


# with bicubic
a <- results$prepData_unmod_locs_NYSE_2013_2009_1to50_smpper0_all$bestRas.l$testErr.r$`1`$ras$`6`$intras

results <- list()

sessionTag <- 'prepData_unmod_locs_NYSE_2013_2009_1to50_smpper0_all'
results[[sessionTag]] <- 
  analyseDat(filePattern = sessionTag,
             plot2pdf = F,
             bestPlots = F,
             descPlots = F,
             oldSave = T)
saveRDS(results,file=paste0('data/results_',sessionTag,'.RDS'))

ex <- exploreFurther(results[[1]])

# saveRDS(results,file='data/results_prepData_bestOffSets_buffpol10_smper0_all.RDS')
sessionTag <- 'prepData_bestOffSets_buffpol10_smper0_all'
results[[sessionTag]] <- 
  analyseDat(filePattern = sessionTag,
             plot2pdf = F,
             bestPlots = F,
             descPlots = F,
             oldSave = T)

sessionTag <- 'prepData_bestOffSets_buffpol10_smper50_all'
results[[sessionTag]] <- 
  analyseDat(filePattern = sessionTag,
             plot2pdf = F,
             bestPlots = F,
             descPlots = F,
             oldSave = T)

 
# resultsTimes <- [['individual result RDS file']]$orig.maps$intTimes[3:7] %>% 
#   map(.f = function(m) map_df(m, ~data.frame(val = .x$val, u = .x$unit_chr))) %>% 
#   bind_rows(.id='int_method') %>% 
#   dplyr::group_by(int_method) %>% 
#   dplyr::summarise(total_time = lubridate::make_difftime(sum(val),
#                                                          units='mins'))


results$`prepData_bestOffSets_buffpol10_smper50_all` <- 
  analyseDat(filePattern = 'prepData_bestOffSets_buffpol10_smper50_all',
             plot2pdf = F,
             bestPlots = F,
             descPlots = F,
             oldSave = T)


# comparison tables
bestRuns.inc.r <- results[[1]]$bestRuns %>% 
  dplyr::filter(error_var %in% 
                  c('compareDiff.inc.r',
                    'testErr.inc.r')) %>% 
  dplyr::select(-c(bid,qv,quantile)) %>%
  tidyr::pivot_wider(names_from = error_var,
                     values_from = c(run_no,error_value)) %>% 
  mutate(intpol_fid = as.character(intpol_fid))

bestCompareDiff.inc.r <- bestRuns.inc.r %>% 
  group_by(intpol_fid) %>% 
  slice_min(order_by = error_value_compareDiff.inc.r)

 extra
results$`prepData_bestOffSets_buffpol10_smper50_all` <- 
  analyseDat(filePattern = "prepData_bestOffSets_buffpol10_smper50_all",
             plot2pdf = T)


saveRDS(l, file='data/results_prepData_bestOffSets_buffpol10_smper0.RDS')
a <- results$`intdat_march5_offset_w50noise_all`$bestRuns %>% 
  filter(int_method == 'Random Forest SP' && 
           error_var == 'testErr.inc.r')

results$`march19_offset_w50noise_rfsp` <- 
  analyseDat(filePattern = 'march19_offset_w50noise_rfsp',
             plot2pdf = F)

aa <- results$`march19_offset_w50noise_rfsp`$bestRuns %>% 
  filter(int_method == 'Random Forest SP' && 
           error_var == 'testErr.inc.r')

results$`intdat_march6_nooffset_w50noise_all` <- 
  analyseDat(filePattern = 'intdat_march6_nooffset_w50noise_all',
             plot2pdf = F)

results$`intdat_march6_offset_nonoise_all` <- 
  analyseDat(filePattern = 'intdat_march6_offset_nonoise_all',
             plot2pdf = F)

results$`intdat_march5_offset_w50noise_all`$plots$crossPlots$testErr.inc.r$`GRASS Regularized Splines Tension`

results$`intdat_march5_offset_w50noise_all`$plots$bestRasPlots$testErr.inc.r$`5`$diff_b

f1 <- results$`intdat_march5_offset_w50noise_all`$plots$crossPlots$testErr.inc.r$`GRASS Regularized Splines Tension`
f2 <- results$`intdat_march5_offset_w50noise_all`$plots$crossPlots$compareDiff.inc.r$`GRASS Regularized Splines Tension`
f3 <- results$`intdat_march6_nooffset_w50noise_all`$plots$crossPlots$compareDiff.inc.r$`GRASS Bicubic Spline`

reportPlots <- list(f1 = f1,
                    f2 = f2,
                    f3 = f3)

saveRDS(reportPlots, file='manuscript/reportplots.RDS')

# save(results,file='/media/mal/working_files/quarry/results_march5-6.RDS')




