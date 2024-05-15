\name{Java.sparrow.notes}
\alias{Java.sparrow.notes}
\docType{data}
\title{Java sparrow note duration and frequency.}
\description{
This dataset reports the durations and frequencies of use for 22,970 notes in 676 undirected songs produced by 73 Java sparrows (\emph{Padda oryzivora}). Songs were recorded by Masayo Soma and colleagues and notes were assigned to types by Rebecca Lewis (Lewis et al. 2021; Lewis et al. 2023).}
\usage{data("Java.sparrow.notes")}
\format{
  A list with 22,970 observations of the following 3 variables.
  \describe{
    \item{\code{duration}}{duration, in ms, of each note in the dataset}
    \item{\code{note}}{note type to which each note was assigned}
    \item{\code{ID}}{identifier of the individual bird that produced each note}
  }
}
\source{<https://figshare.manchester.ac.uk/articles/dataset/Code_and_Data_from_Lewis_et_al_2023_Animal_Behavior/22550149>}
\references{
Lewis, R.N., Soma, M., de Kort, S.R., Gilman, R.T. (2021) Like father like son: cultural and genetic contributions to song inheritance in an Estrildid finch. \emph{Frontiers in Psychology} 12, 654198. (doi.org/10.3389/fpsyg.2021.654198)

Lewis, R.N., Kwong, A., Soma, M., de Kort, S.R., Gilman, R.T. (2023) Inheritance of temporal song features in Java sparrows. \emph{Animal Behaviour} 206, 61-74. (doi.org/10.1016/j.anbehav.2023.09.012)}
\keyword{datasets}
\examples{
  head(Java.sparrow.notes)
}
