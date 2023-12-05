#' @name plotZLA
#' @title Web plots to illustrate Zipf's law of abbreviation
#' @description Produces a webplot. It requires data calculated from the testZLA function.
#'
#' @param inputObject output from the function testZLA
#'
#' @return Produces a webplot illustrating concordance between note duration and frequency of use within individuals. Requires output from the testZLA function.
#' @details
#' In the figure produced by plotZLA, each point represents a note or phrase type in the population repertoire. Note types are joined by a line if both note types are produced by the same individual. The weight of the line is proportional to the number of individuals that produce both note types. The color of the lines indicates whether there is a positive (blue) or negative (red) concordance between the duration and frequency of use of the note types. Negative concordances are consistent with Zipf's law of abbreviation. Shades between blue and red indicate that the concordance is positive in some individuals and negative in others. For example, this can happen if some individuals use the note types more frequently than others, such that the rank order of frequency of use varies among individuals. Grey crosses centered on each point show the longest and shortest durations of the note type (vertical) and the highest and lowest frequencies of use (horizontal) in the population.
#'
#' @references Davis, M. K. and Chen, G. (2007) Graphic Kendall's tau. Computational Statistics & Data Analysis, 51(5), 2373-2378. doi: 10.2307/2346786
#'
#' @author CD Durrant and R. Tucker Gilman (2023)
#'
#' @export

plotZLA <- function(inputObject){

  #Import data from the vocalZLA object
  dataW <- inputObject$plotObject$figureData
  dataFull <- inputObject$plotObject$allData

  #Add relative frequency of use within birds
  totalNotesPerBird <- aggregate(dataW$count,
                                by = list(dataW$ID),
                                sum)
  dataW$noteProp <- dataW$count/totalNotesPerBird$x[match(dataW$ID,
                                                         totalNotesPerBird$Group.1)]

  #Get mean duration and proportion by note class with birds weighted equally
  noteClassData <- aggregate(dataW[,c(6,8)],
                            by = list(dataW$note),
                            mean)
  noteClassData$logProp <- log(noteClassData$noteProp)

  #Add max and min
  noteClassMax <- aggregate(dataW[,c(3,8)],
                           by = list(dataW$note),
                           max)
  noteClassMax$logProp <- log(noteClassMax$noteProp)
  noteClassMin <- aggregate(dataW[,c(3,8)],
                           by = list(dataW$note),
                           min)
  noteClassMin$logProp <- log(noteClassMin$noteProp)

  #For every note class combination produced by every bird, let's compute the
  #stochastic superiority of the more frequently used note class
  stochSup <- NULL
  noteList <- unique(dataW$note)
  for(b in unique(dataW$ID)){
    thisBirdNotes <- dataFull[which(dataFull$ID == b),]
    for(i in 1:(length(noteList)-1)){
      for(j in (i+1):length(noteList)){
        #Get ss only if bird sings both notes and unequal number of times
        if(sum(noteList[c(i,j)] %in% unique(thisBirdNotes$note)) == 2 & sum(thisBirdNotes$note == noteList[i]) != sum(thisBirdNotes$note == noteList[j])){
          suppressWarnings(W <- wilcox.test(thisBirdNotes$duration[which(thisBirdNotes$note == noteList[i])], thisBirdNotes$duration[which(thisBirdNotes$note == noteList[j])]))

          #Get measure of stochastic superiority
          ss <- W$statistic/(sum(thisBirdNotes$note == noteList[i])*sum(thisBirdNotes$note == noteList[j]))

          #Correct measure for direction
          if(sum(thisBirdNotes$note == noteList[i]) < sum(thisBirdNotes$note == noteList[j])){
            ss <- 1 - ss
          }
          #Get weight for proportional difference
          bt <- binom.test(sum(thisBirdNotes$note == noteList[i]), sum(thisBirdNotes$note == noteList[i]) + sum(thisBirdNotes$note == noteList[j]))
          stochSup <- rbind(stochSup, c(b, noteList[i], noteList[j], ss, 1-bt$p.value))
        }
      }
    }
  }
  stochSup = as.data.frame(stochSup)
  colnames(stochSup) = c("bird.id",
                         "note1",
                         "note2",
                         "ss",
                         "w")
  stochSup[,4] = as.numeric(stochSup[,4])
  stochSup[,5] = as.numeric(stochSup[,5])
  stochSup = stochSup[which(stochSup$w > 0),] #remove note pairs where we cannot tell which is larger

  #Find note pair weight and mean stochastic superiority
  npW = aggregate(stochSup$w,
                  by = list(stochSup$note1,
                            stochSup$note2),
                  sum)
  npSS = aggregate(stochSup$w*stochSup$ss,
                   by = list(stochSup$note1,
                             stochSup$note2),
                   sum)
  npSS$x = npSS$x/npW$x

  #make plot area
  xFbB = .05*(max(noteClassMax$logProp) - min(noteClassMin$logProp))
  xLbFb = min(noteClassMin$logProp) - xFbB
  xUbFb = max(noteClassMax$logProp) + xFbB
  yFbB = .05*(max(noteClassMax$logDuration) - min(noteClassMin$logDuration))
  yLbFb = min(noteClassMin$logDuration) - yFbB
  yUbFb = max(noteClassMax$logDuration) + yFbB
  plot(0,
       0,
       col = 0,
       xlab = "",
       ylab = "",
       xlim = c(xLbFb,
                xUbFb),
       ylim = c(yLbFb,
                yUbFb),
       main = "",
       cex.axis = 1.4)
  title(ylab = "mean logged duration",
        cex.lab = 1.4,
        line = 2.75)
  title(xlab = "logged mean relative frequency",
        cex.lab = 1.4,
        line = 2.75)

  #Draw in the lines
  lws = 0.25 #allows us to scale the line width. reducing line width in "busy" figures allows us to see more lines
  for(i in 1:nrow(npSS)){
    n1 = noteClassData[which(noteClassData$Group.1 == npSS[i,1]),]
    n2 = noteClassData[which(noteClassData$Group.1 == npSS[i,2]),]
    segments(n1$logProp,
             n1$meanDuration,
             n2$logProp,
             n2$meanDuration,
             col = rgb(1-npSS$x[i],
                       0,
                       npSS$x[i]),
             lwd = lws*(1+5*npW$x[i]/max(npW$x))/2)
  }
  #Add point and crosses
  points(noteClassData$logProp,
         noteClassData$meanDuration,
         pch = 19,
         col = 0,
         cex = 0.5)
  points(noteClassData$logProp,
         noteClassData$meanDuration,
         pch = 1,
         cex = 0.5)
  for(i in 1:nrow(noteClassData)){
    segments(noteClassMin$logProp[i],
             noteClassData$meanDuration[i],
             noteClassMax$logProp[i],
             noteClassData$meanDuration[i],
             col = "darkgrey",
             lwd = 2*lws)
    segments(noteClassData$logProp[i],
             noteClassMin$logDuration[i],
             noteClassData$logProp[i],
             noteClassMax$logDuration[i],
             col = "darkgrey",
             lwd = 2*lws)
  }
}
