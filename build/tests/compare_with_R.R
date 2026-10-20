# Runs the original R model on the same inputs and prints the largest difference from the JS port.
#   Rscript tests/compare_with_R.R path/to/topmodel_r/R/topmodel_core.R     (run from build/, after compare_with_R.js)
args <- commandArgs(TRUE); source(args[1]); core_use_cpp <- function() FALSE
J <- jsonlite::fromJSON('tests/out/in.json', simplifyVector=TRUE)
JS <- jsonlite::fromJSON('tests/out/js.json')
F <- J$F; TI <- J$TI
for (k in names(J$cases)) {
  c <- J$cases[[k]]; p <- c$p
  if (!is.null(p$luc)) { if (identical(p$luc$compRecDays,'Inf')) p$luc$compRecDays <- Inf
    for (s in c('pre','start','final')) p$luc[[s]] <- as.numeric(ifelse(is.na(p$luc[[s]]),NA,p$luc[[s]])) }
  o <- enhanced_topmodel_core(F, TI, p, if (isTRUE(c$R$on)) c$R else NULL)
  d <- sapply(c('Qsim','ET','SWE','Sbar','satFrac','Qie','Qimp','Qssf','storage','Qcomp','Qdev','PET','Ea'), function(f) max(abs(o[[f]] - JS[[k]][[f]])))
  cat(k, ' meanQ', round(mean(o$Qsim),4), ' maxdiff:', paste(names(d), signif(d,2), collapse=' '), '\n')
}
