## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  warning = FALSE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(ZLAvian)

## ----warning = FALSE----------------------------------------------------------
data(testdata, package = "ZLAvian")

## ----include = FALSE----------------------------------------------------------
data = Java.sparrow.notes

## ----warning = FALSE----------------------------------------------------------
test.ZLA.output = testZLA(data, minimum = 1, null = 999, est = "mixed", cores = 2)

## ----warning = FALSE----------------------------------------------------------
store <- testZLA(data, minimum = 1, null = 999, est = "mixed", cores = 2)

plotZLA(store, ylab = "duration (ms)", x.scale = "linear")

