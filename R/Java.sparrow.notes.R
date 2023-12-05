#' Java sparrow note duration and frequency.
#'
#' This dataset reports the durations and frequencies of use for 22,970 notes in
#' 676 undirected songs produced by 73 Java sparrows (Padda oryzivora).
#' Songs were recorded by Masayo Soma and colleagues and notes were assigned to
#' classes by Rebecca Lewis (Lewis et al. 2021; Lewis et al. 2023).
#'
#' @docType data
#' @usage data(Java.sparrow.notes)
#' @format A list of 22970 rows and 3 variables{
#'  \describe{
#'    \item{duration}{Amount of seconds the note was sung for}
#'    \item{note}{Identifier of note type}
#'    \item{ID}{Identifier of individual bird}
#'  }
#'}
#' @source \href{https://figshare.manchester.ac.uk/articles/dataset/Code_and_Data_from_Lewis_et_al_2023_Animal_Behavior/22550149}
#' @keywords dataset
#' @references Lewis RN, Soma M, de Kort SR, Gilman RT (2021) Like father like son: Cultural and genetic contributions to song inheritance in an estrildid finch. Frontiers in Psychology 12, 20-30.
#' @references Lewis RN, Kwong A, Soma M, de Kort SR, Gilman RT (2023) Inheritance of temporal song features in Java sparrows. Animal Behaviour 206, 61-74
