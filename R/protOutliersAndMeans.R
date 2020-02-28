#'================================================================
#' scoresMod to identify outliers in a vector
#' Outliers as exceeding the normal probability of proba
#'
#' @param             x: a vector of numbers
#' @param  reject.vec.i: pre-specified rejection locations
#' @param         proba: pre-specified normal probability
#' @return             : TRUE (outliers) or FALSE (keeping)
#' @note               : when sd(x)=0, scores() generates NA, mark
#'                       NA as FALSE
#'================================================================
scoresMod <- function(x, reject.vec = NULL, proba = 0.99) {
  x[reject.vec] <- NA
  non.na.indices <- (1:length(x))[!is.na(x)]
  ind <- scores(x[!is.na(x)], prob = proba)
  ind[is.na(ind)] <- F
  outlier.ind <- rep(F, length(x))
  outlier.ind[non.na.indices] <- ind

  return(outlier.ind*1)
}



#'================================================================
#' boxplotMod to identify outliers in a vector
#' Outliers that are 3 times the interquartile range from either
#' quartile; note that the default for boxplots is 1.5 times the
#' interquartile range
#'
#' @param             x: a vector of numbers
#' @param  reject.vec.i: pre-specified rejection locations
#' @param         range: pre-specified times of IQR range
#' @return             : TRUE (outliers) or FALSE (keeping)
#'================================================================
boxplotMod <- function(x, reject.vec = NULL, range = 3) {
    x[reject.vec] <- NA
    outlier.x <- boxplot(x, plot = F, range = range)$out
    outlier.ind <- (x %in% outlier.x)

    return(outlier.ind*1)
}



#'================================================================
#' log2eps to log transformation with epsilon adjustment
#'
#' @param             x: a vector of numbers
#' @param           eps: little value to avoid log2(0)
#' @return             : log transformed x
#'================================================================
log2eps <- function(x, .eps = 0.029885209) log2(x + .eps)



#'================================================================
#' OutlierIdentify to identify outliers in a matrix
#' using both scoresMod and boxplotMod
#'
#' @param    input.Data: a matrix of numbers
#' @param        refCol: number of columns before Mass Spectrometry
#'                       data columns
#' @param          nCol: how many columns in MS data
#' @param     save.data: save output data externally?
#' @param     byPeptide: outliers by peptide(T) or gene(F)?
#' @return             : Input.Data plus outlier indicators
#' @examples
#' protsCombineCnew <- Outlier.Identify(protClass=tmtMS2orig, refCol=4, nCol=9, Save.data=F)

#'================================================================
OutlierIdentify <- function(protClass, refCol, nCol, save.data = FALSE, byPeptide = FALSE) {
  #get warning message status, ON: 0; OFF: -1
  oldw <- getOption("warn")
  #suppress warning message
  options(warn = -1)

  #transformation using scoresMod
  output_score <- protClass %>%
    arrange(geneId, geneSeqId) %>%
    mutate_at(vars((refCol + 1):(refCol + nCol)), log2eps) %>%
    nest(-geneId) %>%
    mutate(
      outlier = map(data, ~ apply(.x[,(refCol + 1):(refCol + nCol)], 2, scoresMod)),
      tidied = map(outlier, tidy)) %>%
    unnest(tidied) %>%
    select(-data,
           -outlier,
           -names,
           -x)

  name <- colnames(protClass)[(refCol + 1):(refCol + nCol)]

  for(i in 1:(length(name))) {
    channel <- quo(!!as.name(name[i]))

    output_score <- output_score %>%
      rename(!!paste0(quo_name(channel), ".score") := !!as.name(quo_name(channel)))
  }

  output_score <- output_score %>%
    filter(complete.cases(.)) #gene with only one spectrum is removed

  #transformation using boxplotMod
  output_box <- protClass %>%
    arrange(geneId, geneSeqId) %>%
    mutate_at(vars((refCol + 1):(refCol + nCol)), log2eps) %>%
    nest(-geneId) %>%
    mutate(
      outlier = map(data, ~ apply(.x[,(refCol + 1):(refCol + nCol)], 2, boxplotMod)),
      tidied = map(outlier, tidy)) %>%
    unnest(tidied) %>%
    select(-data,
           -outlier,
           -names,
           -x)

  for(i in 1:(length(name))) {
    channel <- quo(!!as.name(name[i]))

    output_box <- output_box %>%
      rename(!!paste0(quo_name(channel), ".box") := !!as.name(quo_name(channel)))
  }

  output_box <- output_box %>%
    filter(complete.cases(.))

  #id for gene with more than one spectrum
  id <- protClass %>%
    arrange(geneId, geneSeqId) %>%
    filter(geneId %in% c(unique(output_score$geneId))) %>%
    select(Filename..scan)

  id2 <- protClass %>%
    arrange(geneId, geneSeqId) %>%
    filter(geneId %in% c(unique(output_box$geneId))) %>%
    select(Filename..scan)

  if(!all.equal(id, id2)) {
    cat("Observations identified using scoresMod and boxplotMod are DIFFERENT. \n")
  } else {
    output_score <- output_score %>%
      bind_cols(id) %>%
      select(-geneId)

    output_box <- output_box %>%
      bind_cols(id) %>%
      select(-geneId)

    output <- protClass %>%
      left_join(output_score, by = "Filename..scan") %>%
      left_join(output_box, by = "Filename..scan") %>%
      mutate_at(vars((refCol+nCol+2+1):(refCol+nCol+2+nCol*2)), ~replace_na(., 0)) %>%
      mutate(outlier.score = rowSums(select(., contains("score"))),
             outlier.box = rowSums(select(., contains("box"))))
  }

  if(save.data){
    if(!byPeptide) {
      write.csv(output, file = paste("OutlierList", dataUse, ".txt", sep = ""), row.names = F)
    } else {
      setwd("..")
      setwd("..")
      setwd("dataModPeptide")
      setwd(dataUse)
      write.csv(output, file = paste("OutlierListPeptide", dataUse, ".txt", sep = ""), row.names = F)
    }
  }

  #resume waring message status
  options(warn = oldw)

  return(output)
}



#'================================================================
#' geneProfileEstimate to calculate mean and SE for each channel in
#' each gene using random effect model or arithmetic mean
#' random effect model can avoid dominance of a sequence with
#' too many spectra
#'
#' @param  protsCombineCnew: a matrix of numbers
#' @param           outlier: outlier-removing method: 1)boxplot,
#'                           2)score, 3)none
#' @param            refCol: number of columns before Mass
#'                           Spectrometry data columns
#' @param              nCol: how many columns in MS data
#' @param               eps: epsilon to avoid log(0)
#' @param           dataUse: name of data set
#' @param         save.data: save output data externally?
#' @return                 : calculated mean and SE, raw means
#'                           transformed back (inverse of log2eps)
#' @examples
#' gProfSummaryTMTms2 <- geneProfileEstimate(protsCombineCnew = protsCombineCnew, outlier="boxplot",
#' refCol=4, nCol=9, eps=0.029885209, dataUse="tmtMS2", save.data=F)

#' head(gProfSummaryTMTms2)
#'================================================================
geneProfileEstimate <- function(protsCombineCnew, outlier, refCol, nCol, eps, dataUse, save.data = FALSE) {
  #get warning message status, ON: 0; OFF: -1
  oldw <- getOption("warn")
  #suppress warning message
  options(warn = -1)

  #remove outliers identified by the pre-specified method and log transform raw data
  finalList <- protsCombineCnew %>%
    arrange(geneId, geneSeqId) %>%
    filter(if     (outlier == "boxplot") outlier.box == 0
           else if(outlier == "score")   outlier.score == 0
           else if(outlier == "none")    gene != "") %>%
    mutate_at(vars((refCol + 1):(refCol + nCol)), log2eps)

  #scale log transformed data to better fit random effect model
  scaled_data <- finalList %>%
    group_by(geneId) %>%
    mutate_at(vars((refCol + 1):(refCol + nCol)), scale) %>%
    ungroup() %>%
    rename_at((refCol + 1):(refCol + nCol), list( ~ paste0("scaled.", .))) %>%
    select(Filename..scan,
           starts_with("scaled."))

  #combine scaled data with raw data for model fitting
  #filter good data for random effect model, lmerGood means 1) 3 or more spectra per sequence and 2) 2 or more sequences and 3) there's variation among geneId (SD!=0)
  model_data <- finalList %>%
    rename_at((refCol + 1):(refCol + nCol), list( ~ paste0("origin.", .))) %>%
    left_join(scaled_data, by = "Filename..scan") %>%
    group_by(geneId) %>%
    mutate(Nspectra = n(),
           Nseq = n_distinct(geneSeqId)) %>%
    ungroup() %>%
    mutate(lmerGood = ifelse(((Nspectra/Nseq) >= 3 & Nseq >= 2 & rowSums(is.na(select(., starts_with("scaled.")))) == 0), T, F)) %>%
    select(geneId,
           geneSeqId,
           Nspectra,
           Nseq,
           lmerGood,
           gene,
           peptide,
           modifications,
           Filename..scan,
           starts_with("scaled."),
           starts_with("origin.")) %>%
    filter(lmerGood == T)

  #fit random effect model for each channel of data
  output <- finalList %>%
    group_by(geneId) %>%
    summarise(geneName = unique(gene),
              Nspectra = n(),
              Nseq = n_distinct(geneSeqId)) %>%
    ungroup()

  name <- colnames(finalList)[(refCol + 1):(refCol + nCol)]

  for(i in 1:(length(name))) {
    channel <- quo(!!as.name(name[i]))

    model_fit <- model_data %>%
      nest(-geneId) %>%
      mutate(fit      = tryCatch(map(data, ~ lmer(!!as.name(paste0("scaled.", quo_name(channel))) ~ 1 + (1 | geneSeqId), data = .x)),
                                 error = function(e) {cat("ERROR :", conditionMessage(e), "\n")}),
             tidied   = tryCatch(map(fit, tidy),
                                 error = function(e) {cat("ERROR :", conditionMessage(e), "\n")}),
             singular = tryCatch(map(fit, isSingular),
                                 error = function(e) {cat("ERROR :", conditionMessage(e), "\n")})) %>%
      unnest(tidied) %>%
      filter(term == "(Intercept)") %>%
      select(-data,
             -fit,
             -term,
             -statistic,
             -group)

    model_mu_sigma <- model_data %>%
      group_by(geneId) %>%
      summarise(mu    = mean(!!as.name(paste0("origin.", quo_name(channel)))),
                sigma = sd(  !!as.name(paste0("origin.", quo_name(channel))))) %>%
      ungroup()

    model_fit <- model_fit %>%
      left_join(model_mu_sigma, by = "geneId") %>%
      mutate(!!paste0(quo_name(channel), ".coef")     := estimate*sigma + mu,
             !!paste0(quo_name(channel), ".se")       := std.error*sigma,
             !!paste0(quo_name(channel), ".singular") := singular) %>%
      select(-estimate,
             -std.error,
             -singular,
             -mu,
             -sigma)

    output <- full_join(output, model_fit, by = "geneId")

    all_mu_sigma <- finalList %>%
      group_by(geneId) %>%
      summarise(!!paste0(quo_name(channel), ".all.mu") := mean(!!channel),
                !!paste0(quo_name(channel), ".all.se") := sd(!!channel)/n()) %>%
      ungroup()

    output <- full_join(output, all_mu_sigma, by = "geneId")

    output <- output %>%
      mutate(!!paste0("log.", quo_name(channel), ".MU") := ifelse(!!as.name(paste0(quo_name(channel), ".singular")) == "FALSE",
                                                                  !!as.name(paste0(quo_name(channel), ".coef")),
                                                                  !!as.name(paste0(quo_name(channel), ".all.mu"))),
             !!paste0("log.", quo_name(channel), ".SE") := ifelse(!!as.name(paste0(quo_name(channel), ".singular")) == "FALSE",
                                                                  !!as.name(paste0(quo_name(channel), ".se")),
                                                                  !!as.name(paste0(quo_name(channel), ".all.se"))),
             !!paste0("raw.", quo_name(channel), ".MU") := 2^!!as.name(paste0("log.", quo_name(channel), ".MU")) - eps,
             !!paste0("raw.", quo_name(channel), ".SE") := 2^!!as.name(paste0("log.", quo_name(channel), ".MU"))*log(2)*!!as.name(paste0("log.", quo_name(channel), ".SE"))) %>%
      select(-!!paste0(quo_name(channel), ".coef"),
             -!!paste0(quo_name(channel), ".all.mu"),
             -!!paste0(quo_name(channel), ".se"),
             -!!paste0(quo_name(channel), ".all.se"),
             -!!paste0(quo_name(channel), ".singular"))
  }

  output <- output %>%
    select(geneId,
           geneName,
           Nspectra,
           Nseq,
           matches("(raw.*MU)"),
           matches("(raw.*SE)"),
           matches("(log.*MU)"),
           matches("(log.*SE)"))

  #adjust raw mu and raw se, make raw mu sum up to one, raw se is divided by the total sum
  raw.sum <- rowSums(output[, (4+1):(4+nCol)])
  output[,(4+nCol+1):(4+nCol*2)] <- output[,(4+nCol+1):(4+nCol*2)]/raw.sum
  output[,(4+1):(4+nCol)] <- output[,(4+1):(4+nCol)]/raw.sum

  if(save.data) write.csv(output, file = paste("geneProfileSummary", dataUse, ".txt", sep = ""), row.names = F)

  #resume waring message status
  options(warn = oldw)

  return(output)
}






























