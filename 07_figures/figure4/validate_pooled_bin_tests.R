args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1)
root <- dirname(args[[1]])
d <- read.delim(args[[1]],check.names=FALSE)
stopifnot(nrow(d)==28)
for (i in which(!is.na(d$p_one_sided_less))) {
  z <- matrix(c(d$high_N_obs[i],d$low_N_obs[i],d$high_S_obs[i],d$low_S_obs[i]),nrow=2)
  p <- fisher.test(z,alternative="less")$p.value
  stopifnot(isTRUE(all.equal(p,d$p_one_sided_less[i],tolerance=1e-12)))
}
for (sp in unique(d$species)) {
  idx <- which(d$species==sp)
  p <- d$p_one_sided_less[idx]
  valid <- which(!is.na(p))
  if (length(valid)) {
    q <- p.adjust(p[valid],method="BH",n=4)
    stopifnot(isTRUE(all.equal(q,d$BH_q_within_species_four_planned[idx][valid],tolerance=1e-12)))
  }
}
valid <- which(!is.na(d$p_one_sided_less))
stopifnot(isTRUE(all.equal(p.adjust(d$p_one_sided_less[valid],method="BH",n=28),
                           d$BH_q_all_species_28_planned[valid],tolerance=1e-12)))
cat(sprintf("PASS: %d one-sided Fisher P values; seven 4-comparison BH families and global 28-comparison BH.\n",length(valid)))

rows <- list()
for (comparison in unique(d$comparison)) {
  sub <- d[d$comparison==comparison & !is.na(d$p_one_sided_less),,drop=FALSE]
  n <- nrow(sub)
  high_total <- sum(sub$high_n_sites)
  low_total <- sum(sub$low_n_sites)
  low_dir <- sum(sub$effect_direction=="high_bin_lower_N_over_S")
  high_dir <- sum(sub$effect_direction=="high_bin_higher_N_over_S")
  equal_dir <- sum(sub$effect_direction %in% c("equal_sample_N_over_S","undefined_no_N_S_contrast"))
  rec <- data.frame(comparison=comparison,n_included_species=n,
      included_species=paste(sub$species,collapse=";"),high_bin_sites=high_total,low_bin_sites=low_total,
      n_species_high_lower=low_dir,n_species_high_higher=high_dir,
      n_species_equal_or_no_N_S_contrast=equal_dir,
      test="exact stratified Mantel-Haenszel; H1 common OR<1",
      common_odds_ratio_conditional_MLE=NA_real_,p_one_sided_exact=NA_real_,
      status="not_estimated",stringsAsFactors=FALSE)
  if (n>=2) {
    arr <- array(0L,dim=c(2,2,n))
    for (i in seq_len(n)) {
      arr[,,i] <- matrix(c(sub$high_N_obs[i],sub$low_N_obs[i],
                           sub$high_S_obs[i],sub$low_S_obs[i]),nrow=2)
    }
    fit <- tryCatch(mantelhaen.test(arr,alternative="less",exact=TRUE),error=function(e)e)
    if (inherits(fit,"error")) {
      rec$status <- paste0("exact_test_error: ",conditionMessage(fit))
    } else {
      rec$common_odds_ratio_conditional_MLE <- unname(fit$estimate)
      rec$p_one_sided_exact <- fit$p.value
      rec$status <- if (high_total<10) "exact_test_computable_but_high_bin_sites_under_10" else "exact_test_computable"
    }
  } else {
    rec$status <- "fewer_than_two_nonempty_species_strata"
  }
  rows[[length(rows)+1]] <- rec
}
out <- do.call(rbind,rows)
out$BH_q_four_stratified_comparisons <- p.adjust(out$p_one_sided_exact,method="BH")
write.table(out,file.path(root,"cross_species_exact_stratified_tests.tsv"),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
print(out)
