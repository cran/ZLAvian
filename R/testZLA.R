#' @name testZLA
#' @title Assess evidence for Zipf's law of abbreviation
#' @description Assesses evidence for Zipf's Law of Abbreviation in a population where samples from the population repertoire can be assigned to individuals.
#'
#' @param data a dataframe containing columns "note" (factor/character; identifies the note/phrase type of each token), "duration" (numeric; describes the duration of each token), and "ID" (factor; identifies the individual that produced each token). Other columns in the dataframe are ignored.
#' @param minimum the minimum number of times a note type must appear in the data set to be included in the analysis. Must be a positive integer.
#' @param null the number of permutations used to estimate the null distribution. Must be a positive integer 99 or greater.
#' @param est takes values "mixed," "mean," or "median." If est = "mixed" then the expected logged duration for each note type in the population is computed as the intercept from an intercept-only mixed effects model, fit using the lmer() function of lme4 (Bates et al. 2015), that includes a random effect of individual ID. This computes a weighted mean across birds, and accords greater weights to birds that produce the note type more frequently. If est = "mean" then the expected logged duration for each note type in the population is computed as the mean of the means for the individual birds, with each individual bird weighted equally. This is faster than the "mixed" method. If est = "median" then the expected logged duration for each note type within birds is taken to be the median logged duration of the note type when produced by that bird, and the expected logged duration for each note type in the population is taken to be the median of the medians for the birds that produced that note type. The expected durations for note types are used in the permutation algorithm.
#' @param cores divides (parallelizes) computation of the null distribution among cores. Cores must be an integer between 1 and the number of cores available on the users machine, inclusive.
#' @param transform takes values "log" or "none." Indicates how duration data should be transformed prior to analysis. Gilman and colleagues (2023) argue that log transformation may often be appropriate for duration data, but some other measures might be better analysed as raw (untransformed) values.
#'
#' @return a matrix that reports Kendall's tau and the p-value associated with Kendall's tau computed at both the population and individual levels.
#'
#' @author CD Durrant and R. Tucker Gilman (2024)
#'
#' @import lme4
#' @import performance
#' @import parallel
#' @import doParallel
#' @export

testZLA<-function (data, minimum = 1, null = 999, est = "mean", cores = 2, transform = "log")
{

  #initialise objects that will hold outputs
  ZLAIndWeighted = list(tau = NULL, pValue = NULL)
  ZLAIndUnweighted = list(tau = NULL, pValue = NULL)
  ZLAPopulation = list(tau = NULL, pValue = NULL)
  outputList <- list(overview = NULL, shannon = NULL, stats = NULL,
                     unweighted = NULL, plotObject = NULL, thresholds = NULL)

  ###conduct some preliminary checks
  if (minimum != round(minimum) | minimum < 1) {
    warning(paste("invalid input for \"minimum\" (", minimum,
                  "). Will be updated to closest valid input: ", max(1,
                                                                     round(minimum, 0))))
    minimum <- max(1, round(minimum, 0))
  }

  if (null != round(null) | null < 99) {
    warning(paste("invalid input for \"null\" (", null, "). Will be updated to closest valid input: ",
                  max(99, round(null, 0))))
    null <- max(99, round(null, 0))
  }

  if (cores != round(cores) | cores < 1 | cores > detectCores()) {
    warning(paste("invalid input for \"cores\" (", cores,
                  "). Will be updated to closest valid input: ", min(max(1,
                                                                         round(cores, 0)), detectCores())))
    cores <- min(max(1, round(cores, 0)), detectCores())
  }

  if (est != "mean" & est != "mixed" & est != "median") {
    validOpts <- c("mean", "mixed", "median")
    closestMatch <- agrep(est, validOpts, value = TRUE)
    if (length(closestMatch) > 0) {
      warning(paste0("invalid argument for \"est\".\n",
                     "Closest valid option is \"", closestMatch, "\".\n",
                     "Setting \"est\" to \"", closestMatch))
      est <- closestMatch
    }else {
      warning("invalid argument for \"est\" with no close match. \"est\" set to default \"mean\"")
      est <- "mean"
    }
  }

  if (!(transform %in% c("log","none"))) {
    warning("invalid input for transform. Test conducted using raw durations.")
    transform <- "none"
  }

  expectedColumns <- c("note", "duration", "ID")
  extraColumns <- setdiff(names(data), expectedColumns)
  if (!("data.frame" %in% class(data))) {
    stop("input is not a dataframe")
  }
  else if (!all(expectedColumns %in% names(data))) {
    missingColumns <- expectedColumns[!expectedColumns %in%
                                        names(data)]
    stop(paste("the following columns are missing or potentially misspelt:",
               paste(missingColumns, collapse = ", ")))
  }else if (!("character" %in% class(data$note)) | !("numeric" %in%
                                                    class(data$duration)) | !(class(data$ID) %in% c("character",
                                                                                                    "factor"))) {
    stop("invalid data type in one or more columns. \"note\" and \"ID\" must be the data type \"character\". \"duration\" must be the data type \"double\"")
  }

  if (transform == "log" & min(data$duration) < 0){
    stop("duration data includes negative numbers.\nLog transform not meaningful.\nCorrect data errors or analyse raw duration data.", call. = FALSE)
  }
  ###

  #log transform duration (unless transform == "none")
  if (transform == "log"){
    data$duration=log(data$duration)
  }

  #get the mean of median for durations of note types within birds
  if (est=="median"){
    fxkDuration <- aggregate(data$duration, by = list(data$ID,
                                                           data$note), median)
    meanDuration <- aggregate(fxkDuration$x, by = list(fxkDuration$Group.2),
                              median)
  }else{
    fxkDuration <- aggregate(data$duration, by = list(data$ID,
                                                           data$note), mean)
    meanDuration <- aggregate(fxkDuration$x, by = list(fxkDuration$Group.2),
                            mean)
  }
  fxkCount <- aggregate(rep(1, nrow(data)), by = list(data$ID,
                                                      data$note), sum)

  #build the data set vocalDB. If est="mixed," overwrite the unweighted mean with the intercept from a mixed effects model
  oldWarning <- getOption("warn")
  options(warn=-1)
  if (est %in% c("mean","median")) {
    vocalDB <- data.frame(paste(fxkCount$Group.1, fxkCount$Group.2,
                                sep = "."), fxkCount$x, fxkDuration$x, fxkCount$Group.1,
                          fxkCount$Group.2, meanDuration$x[match(fxkCount$Group.2,
                                                                 meanDuration$Group.1)])
    colnames(vocalDB) <- c("IDNote", "count", "logDuration",
                           "ID", "note", "meanDuration")
    vocalDB$indDev <- vocalDB$logDuration - vocalDB$meanDuration
  } else {
    model2 <- meanDuration
    suppressMessages({
      for (i in unique(data$note)) {
        focData <- data[which(data$note == i), ]
        if (length(unique(focData$ID)) > 1 & length(unique(focData$ID)) <
            nrow(focData)) {

          tryCatch({
          linearModel = lmer(duration ~ 1 + (1 | ID), data = focData, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e+05)))},
          error = function(msg){
            options(warn=oldWarning)
            stop("The lmer function cannot fit mixed effects models\nto the duration data. Try transforming\nthe data, or use the 'mean' rather than 'mixed' method.")})

          if (is.null(linearModel@optinfo$conv$lme4$code)) {
            lmeConv <- 0
          }
          else {
            lmeConv <- linearModel@optinfo$conv$lme4$code
          }
          if (lmeConv <= 0 & linearModel@optinfo$conv$opt <=
              0 & check_singularity(linearModel) == FALSE) {
            model2$x[which(model2$Group.1 == i)] <- summary(linearModel)[[10]][1]
          }
        }
      }
    })
    vocalDB <- data.frame(paste(fxkCount$Group.1, fxkCount$Group.2,
                                sep = "."), fxkCount$x, fxkDuration$x, fxkCount$Group.1,
                          fxkCount$Group.2, model2$x[match(fxkCount$Group.2,
                                                           model2$Group.1)])
    colnames(vocalDB) <- c("IDNote", "count", "logDuration",
                           "ID", "note", "meanDuration")
    vocalDB$indDev <- vocalDB$logDuration - vocalDB$meanDuration
  }
  options(warn=oldWarning)


  if (minimum > 1) {
    classCount <- aggregate(vocalDB$count, by = list(vocalDB$note),
                            sum)
    classRemove <- classCount$Group.1[which(classCount$x <
                                              minimum)]
    if (length(classRemove) > 0) {
      vocalDB <- vocalDB[-which(vocalDB$note %in% classRemove),
      ]
    }
    if (dim(vocalDB)[1] == 0) {
      stop("\"minimum\" too high; all notes removed")
    }
  }

  #begin analyses
  fullPopNotes <- aggregate(vocalDB$count, by = list(vocalDB$note),
                            sum)
  fullPopNotes$logDuration <- vocalDB$meanDuration[match(fullPopNotes$Group.1,
                                                         vocalDB$note)]
  fpKendall <- cor.test(fullPopNotes$x, fullPopNotes$logDuration,
                        method = "kendall", alternative = "less", exact = FALSE)
  totalIndividuals <- length(unique(data$ID))
  totalNotes <- nrow(data)
  totalNoteClasses <- length(unique(data$note))
  notesPerIndividual <- (totalNotes/totalIndividuals)
  ZLAPopulation$tau <- fpKendall$estimate
  ZLAPopulation$pValue <- fpKendall$p.value
  classesPerIndividual <- nrow(fxkCount)/totalIndividuals
  ndPop <- aggregate(rep(1, nrow(data)), by = list(data$note),
                     sum)
  ndProp <- ndPop$x/nrow(data)
  ShannonPopulation <- -sum(ndProp * log(ndProp))
  ndInd <- aggregate(rep(1, nrow(data)), by = list(data$note,
                                                   data$ID), sum)
  ncInd <- aggregate(rep(1, nrow(data)), by = list(data$ID),
                     sum)
  indProp <- ndInd$x/ncInd$x[match(ndInd$Group.2, ncInd$Group.1)]
  ShannonIndividual <- -sum(indProp * log(indProp))/totalIndividuals
  for (i in unique(vocalDB$ID)) {
    if (min(vocalDB$count[which(vocalDB$ID == i)]) == max(vocalDB$count[which(vocalDB$ID ==
                                                                              i)])) {
      vocalDB <- vocalDB[-which(vocalDB$ID == i), ]
    }
  }
  recIDs <- unique(vocalDB$ID)
  recTaus <- NULL
  tauVar <- NULL
  for (i in 1:length(recIDs)) {
    thisRecData <- vocalDB[which(vocalDB$ID == recIDs[i]),
    ]
    tau <- cor.test(thisRecData$count, thisRecData$logDuration,
                    method = "kendall", exact = FALSE)
    recTaus <- rbind(recTaus, tau$estimate)
    tauN <- nrow(thisRecData)
    tauVar <- rbind(tauVar, 2 * (2 * tauN + 5)/(9 * tauN *
                                                  (tauN - 1)))
  }
  if (is.null(recTaus)) {
    stop("All birds sing all notes in their repertoires with equal frequencies. Concordance cannot be computed.")
  }
  ZLAIndUnweighted$tau <- mean(recTaus)
  ZLAIndWeighted$tau <- sum(recTaus/tauVar)/sum(1/tauVar)
  nullTaus <- NULL
  nullTausW <- NULL
  noteList <- unique(vocalDB$note)
  computeNullTau <- function(j) {
    noteMap <- data.frame(noteList, sample(vocalDB$meanDuration[match(noteList,
                                                                      vocalDB$note)]))
    vocalDB$nullDuration <- noteMap[match(vocalDB$note, noteMap[,
                                                                1]), 2] + sign(rnorm(1)) * vocalDB$indDev
    sampleTau <- NULL
    for (i in 1:length(recIDs)) {
      thisRecData <- vocalDB[which(vocalDB$ID == recIDs[i]),
      ]
      tau <- cor.test(thisRecData$count, thisRecData$nullDuration,
                      method = "kendall", exact = FALSE)
      sampleTau <- rbind(sampleTau, tau$estimate)
    }
    return(list(nullTau = mean(sampleTau), nullTauW = sum(sampleTau/tauVar)/sum(1/tauVar)))
  }
  cl <- makeCluster(cores)
  clusterExport(cl = cl, varlist = c("vocalDB", "noteList",
                                     "null", "recIDs", "tauVar"), envir = environment())
  results <- parLapply(cl, 1:null, computeNullTau)
  stopCluster(cl)
  nullTaus <- sapply(results, function(res) res$nullTau)
  nullTausW <- sapply(results, function(res) res$nullTauW)
  ZLAIndUnweighted$pValue <- (1 + sum(nullTaus < ZLAIndUnweighted$tau))/(1 +
                                                                           length(nullTaus))
  ZLAIndWeighted$pValue <- (1 + sum(nullTausW < ZLAIndWeighted$tau))/(1 +
                                                                        length(nullTausW))
  suppressWarnings(outputList$stats <- matrix(c(ZLAIndWeighted$tau,
                                                ZLAPopulation$tau, ZLAIndWeighted$pValue, ZLAPopulation$pValue),
                                              nrow = 2))
  rownames(outputList$stats) <- c("individual level", "population level")
  colnames(outputList$stats) <- c("tau", "p")

  ###find the threshold tau at which ZLA would be detected (alpha = 0.05)
  if (null>98){

    SortedNullTaus=sort(c(nullTaus,ZLAIndUnweighted$tau))
    SortedNullTausW=sort(c(nullTausW,ZLAIndWeighted$tau))

    minNS=SortedNullTaus[ceiling(0.05*(null+1))]
    maxNS=SortedNullTaus[floor(0.95*(null+1))]
    minNSW=SortedNullTausW[ceiling(0.05*(null+1))]
    maxNSW=SortedNullTausW[floor(0.95*(null+1))]

    minRank=rank(unique(SortedNullTaus))[which(unique(SortedNullTaus)==minNS)]-1
    maxRank=rank(unique(SortedNullTaus))[which(unique(SortedNullTaus)==maxNS)]+1
    minRankW=rank(unique(SortedNullTausW))[which(unique(SortedNullTausW)==minNSW)]-1
    maxRankW=rank(unique(SortedNullTausW))[which(unique(SortedNullTausW)==maxNSW)]+1

    minThresh=NA
    if (minRank %in% rank(unique(SortedNullTaus))){
      minThresh=unique(SortedNullTaus)[minRank]
    }
    maxThresh=NA
    if (maxRank %in% rank(unique(SortedNullTaus))){
      maxThresh=unique(SortedNullTaus)[maxRank]
    }
    minThreshW=NA
    if (minRankW %in% rank(unique(SortedNullTausW))){
      minThreshW=unique(SortedNullTausW)[minRankW]
    }
    maxThreshW=NA
    if (maxRankW %in% rank(unique(SortedNullTausW))){
      maxThreshW=unique(SortedNullTausW)[maxRankW]
    }


    outputList$thresholds=matrix(c(minThresh,minThreshW,
                                  maxThresh,maxThreshW),nrow=2)

    rownames(outputList$thresholds)=c("unweighted","weighted")
    colnames(outputList$thresholds)=c("lower","upper")

  }
  ###

  suppressWarnings(outputList$unweighted <- matrix(c(ZLAIndUnweighted$tau,
                                                     ZLAIndUnweighted$pValue)))
  rownames(outputList$unweighted) <- c("tau", "p")
  suppressWarnings(outputList$overview <- matrix(c(totalNotes,
                                                   totalIndividuals, totalNoteClasses, notesPerIndividual,
                                                   classesPerIndividual)))
  rownames(outputList$overview) <- c("total notes", "total individuals",
                                     "total note classes", "notes per individual", "classes per individual")
  suppressWarnings(outputList$shannon <- matrix(c(ShannonPopulation,
                                                  ShannonIndividual)))
  rownames(outputList$shannon) <- c("population index", "individual index")
  outputList$plotObject <- list(figureData = vocalDB, allData = data, transform = transform)
  print(outputList$stats)
  return(invisible(outputList))
}
