\name{NEWS}
\alias{NEWS}

\title{NEWS for Package \pkg{pogit}}


\section{Changes in version 1.3.0}{
  \itemize{
    \item Fixed empty labels in Rd files of 'logitBvs', 'negbinBvs', 'pogitBvs', 'poissonBvs'.
  }
}

\section{Changes in version 1.2.0}{
  \itemize{
    \item New (invisible) function 'compute_mixture' since 'compute.mixture' from
    package BayesLogit (archived) is no longer supported. 
    \item DESCRIPTION file: Imports changed (BayesLogit removed)
  }
}

\section{Changes in version 1.1.0}{
  \itemize{
    \item New function 'negbinBvs' to deal with overdispersion of count data (and all affected functions). 
    \item Minor changes in 'poissonBvs' (uses now 'select_poissonOD' for observation-specific random intercept)
    \item Minor changes in 'pogitBvs' (uses now 'select_poissonOD' for observation-specific random intercept)
    \item Updated CITATION file to refer to newly published article in Statistical Modelling 16(1), 24 - 46, doi: 10.1177/1471082x15588398.
  }
}

\section{Changes in version 1.0.1}{
  \itemize{
    \item Fixed missing imports from non-base default packages.
    \item Updated CITATION file.
  }
}


