#' @name plotZLA
#' @title Web plots to illustrate Zipf's law of abbreviation
#' @description Produces a webplot. It requires data calculated from the testZLA function.
#'
#' @param inputObject output from the function testZLA
#'
#' @param title xyz
#'
#' @param ylab xyz
#'
#' @param x.scale xyz
#'
#' @param y.base xyz
#'
#' @return Produces a webplot illustrating concordance between note duration and frequency of use within individuals. Requires output from the testZLA function.
#' @details
#' In the figure produced by plotZLA, each point represents a note or phrase type in the population repertoire. Note types are joined by a line if both note types are produced by the same individual. The weight of the line is proportional to the number of individuals that produce both note types. The color of the lines indicates whether there is a positive (blue) or negative (red) concordance between the duration and frequency of use of the note types. Negative concordances are consistent with Zipf's law of abbreviation. Shades between blue and red indicate that the concordance is positive in some individuals and negative in others. For example, this can happen if some individuals use the note types more frequently than others, such that the rank order of frequency of use varies among individuals. Grey crosses centered on each point show the longest and shortest durations of the note type (vertical) and the highest and lowest frequencies of use (horizontal) in the population.
#'
#' @references Davis, M. K. and Chen, G. (2007) Graphic Kendall's tau. Computational Statistics & Data Analysis, 51(5), 2373-2378. doi: 10.2307/2346786
#'
#' @author CD Durrant and R. Tucker Gilman (2024)
#'
#' @export

plotZLA<-function (inputObject, title="", ylab="duration (s)", x.scale = "log", y.base = 2)
{

  #catch and correct unrecognised arguments
  if (!(x.scale %in% c("log","linear"))) {
    warning("x.scale must be either 'log' or 'linear'.\nPlot produced using log scaling on the x-axis.")
    x.scale <- "log"
  }

  #convert logged durations to base 2 logs (b/c this is more intuitive to readers)
  if (inputObject$plotObject$transform == "log"){

    if (!(y.base %in% c(2,10))) {
      warning("y.base must be either 2 or 10.\nplotZLA will attempt the best choice for the data.")
      if (max(inputObject$plotObject$figureData$logDuration)-min(inputObject$plotObject$figureData$logDuration)>log(100)){
        y.base = 10
      } else {
        y.base = 2
      }
    }

    inputObject$plotObject$figureData$logDuration=inputObject$plotObject$figureData$logDuration/log(y.base)
    inputObject$plotObject$figureData$meanDuration=inputObject$plotObject$figureData$meanDuration/log(y.base)
    inputObject$plotObject$figureData$indDev=inputObject$plotObject$figureData$indDev/log(y.base)

  }


  dataW <- inputObject$plotObject$figureData
  dataFull <- inputObject$plotObject$allData
  totalNotesPerBird <- aggregate(dataW$count, by = list(dataW$ID),
                                 sum)
  dataW$noteProp <- dataW$count/totalNotesPerBird$x[match(dataW$ID,
                                                          totalNotesPerBird$Group.1)]
  noteClassData <- aggregate(dataW[, c(6, 8)], by = list(dataW$note),
                             mean)
  noteClassMax <- aggregate(dataW[, c(3, 8)], by = list(dataW$note),
                            max)
  noteClassMin <- aggregate(dataW[, c(3, 8)], by = list(dataW$note),
                            min)
  #proportions will be stores in
  if (x.scale=="log"){
    noteClassData$plotProp <- log10(noteClassData$noteProp)
    noteClassMax$plotProp <- log10(noteClassMax$noteProp)
    noteClassMin$plotProp <- log10(noteClassMin$noteProp)
  } else {
    noteClassData$plotProp <- noteClassData$noteProp
    noteClassMax$plotProp <- noteClassMax$noteProp
    noteClassMin$plotProp <- noteClassMin$noteProp
  }

  stochSup <- NULL
  noteList <- unique(dataW$note)
  for (b in unique(dataW$ID)) {
    thisBirdNotes <- dataFull[which(dataFull$ID == b), ]
    for (i in 1:(length(noteList) - 1)) {
      for (j in (i + 1):length(noteList)) {
        if (sum(noteList[c(i, j)] %in% unique(thisBirdNotes$note)) ==
            2 & sum(thisBirdNotes$note == noteList[i]) !=
            sum(thisBirdNotes$note == noteList[j])) {
          suppressWarnings(W <- wilcox.test(thisBirdNotes$duration[which(thisBirdNotes$note ==
                                                                           noteList[i])], thisBirdNotes$duration[which(thisBirdNotes$note ==
                                                                                                                         noteList[j])]))
          ss <- W$statistic/(sum(thisBirdNotes$note ==
                                   noteList[i]) * sum(thisBirdNotes$note ==
                                                        noteList[j]))
          if (sum(thisBirdNotes$note == noteList[i]) <
              sum(thisBirdNotes$note == noteList[j])) {
            ss <- 1 - ss
          }
          bt <- binom.test(sum(thisBirdNotes$note ==
                                 noteList[i]), sum(thisBirdNotes$note == noteList[i]) +
                             sum(thisBirdNotes$note == noteList[j]))
          stochSup <- rbind(stochSup, c(b, noteList[i],
                                        noteList[j], ss, 1 - bt$p.value))
        }
      }
    }
  }
  stochSup = as.data.frame(stochSup)
  colnames(stochSup) = c("bird.id", "note1", "note2", "ss",
                         "w")
  stochSup[, 4] = as.numeric(stochSup[, 4])
  stochSup[, 5] = as.numeric(stochSup[, 5])
  stochSup = stochSup[which(stochSup$w > 0), ]
  npW = aggregate(stochSup$w, by = list(stochSup$note1, stochSup$note2),
                  sum)
  npSS = aggregate(stochSup$w * stochSup$ss, by = list(stochSup$note1,
                                                       stochSup$note2), sum)
  npSS$x = npSS$x/npW$x

  #get x limits for figure
  xFbB = 0.05 * (max(noteClassMax$plotProp) - min(noteClassMin$plotProp))
  xUbFb = max(noteClassMax$plotProp) + xFbB
  if (x.scale=="log"){
    xLbFb = min(noteClassMin$plotProp) - xFbB
  } else {
    xLbFb = 0
  }

  #ensure x axis includes at least 2 ticks
  if (x.scale=="log" & (ceiling(xUbFb)-floor(xLbFb))<3){
    xUbFb=ceiling(xUbFb)
    xLbFb=floor(xLbFb)
  }

  #get y limits for figure
  yFbB = 0.05 * (max(noteClassMax$logDuration) - min(noteClassMin$logDuration))
  yLbFb = min(noteClassMin$logDuration) - yFbB
  yUbFb = max(noteClassMax$logDuration) + yFbB
  #ensure y axis includes at least 3 ticks
  if ((ceiling(yUbFb)-floor(yLbFb))<4 & inputObject$plotObject$transform=="log"){
    yUbFb=ceiling(yUbFb)
    yLbFb=floor(yLbFb)
  }

  par(mar=c(5.1, 5.1, 4.1, 1.1),las=1)
  plot(0, 0, col = 0, xlab = "", ylab = "", xlim = c(xLbFb,
                                                     xUbFb), ylim = c(yLbFb, yUbFb), main = title, cex.axis = 1.4,axes=F)
  box()
  if (x.scale=="log"){
    all.x.labs=c(parse(text='10^-10'),parse(text='10^-9'),parse(text='10^-8'),parse(text='10^-7'),parse(text='10^-6'),parse(text='10^-5'),parse(text='10^-4'),parse(text='10^-3'),parse(text='10^-2'),parse(text='10^-1'),parse(text='10^0'))
    axis(1,at=(-10:0),labels=all.x.labs,cex.axis=1)
  } else {
    axis(1)
  }

  yat=seq(ceiling(yLbFb),floor(yUbFb))
  if (length(yat)==2 & y.base==2){
    yat=log2(2^(yat[1])+(0:5)/5)
  }
  if (length(yat)==2 & y.base==10){
    yat=log10(10^(yat[1])+(0:5)/5)
  }

  if (inputObject$plotObject$transform=="log" & y.base==2){
    axis(2,at=yat,labels=round(y.base^yat,3),cex.axis=1)
  } else if (inputObject$plotObject$transform=="log" & y.base==10){
    all.y.labs=c(parse(text='10^-10'),parse(text='10^-9'),parse(text='10^-8'),parse(text='10^-7'),parse(text='10^-6'),parse(text='10^-5'),parse(text='10^-4'),parse(text='10^-3'),parse(text='10^-2'),parse(text='10^-1'),parse(text='10^0'),parse(text='10^1'),parse(text='10^2'),parse(text='10^3'),parse(text='10^4'),parse(text='10^5'),parse(text='10^6'),parse(text='10^7'),parse(text='10^8'),parse(text='10^9'),parse(text='10^10'))
    axis(2,at=(-10:10),labels=all.y.labs,cex.axis=1)
  } else {
    axis(2,cex.axis=1)
  }
  title(ylab = ylab, cex.lab = 1.2, line = 3.5)
  title(xlab = "relative frequency", cex.lab = 1.2,
        line = 2.5)
  lws = 0.25
  for (i in 1:nrow(npSS)) {
    n1 = noteClassData[which(noteClassData$Group.1 == npSS[i,
                                                           1]), ]
    n2 = noteClassData[which(noteClassData$Group.1 == npSS[i,
                                                           2]), ]
    segments(n1$plotProp, n1$meanDuration, n2$plotProp, n2$meanDuration,
             col = rgb(1 - npSS$x[i], 0, npSS$x[i]), lwd = lws *
               (1 + 5 * npW$x[i]/max(npW$x))/2)
  }
  points(noteClassData$plotProp, noteClassData$meanDuration,
         pch = 19, col = 0, cex = 0.5)
  points(noteClassData$plotProp, noteClassData$meanDuration,
         pch = 1, cex = 0.5)
  for (i in 1:nrow(noteClassData)) {
    segments(noteClassMin$plotProp[i], noteClassData$meanDuration[i],
             noteClassMax$plotProp[i], noteClassData$meanDuration[i],
             col = "darkgrey", lwd = 2 * lws)
    segments(noteClassData$plotProp[i], noteClassMin$logDuration[i],
             noteClassData$plotProp[i], noteClassMax$logDuration[i],
             col = "darkgrey", lwd = 2 * lws)
  }
}
