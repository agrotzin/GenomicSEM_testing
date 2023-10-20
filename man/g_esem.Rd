\name{g_esem}
\alias{g_esem}
\title{Genomic Exploratory Structural Equation Modeling}
\description{
Function to estimating exploratory structural equation modeling (ESEM) using genomic data as input
}
\usage{
g_esem(covstruc, model = "", rotation = "geomin", rotation.args=list(NULL),std.lv=TRUE, toler=NULL \dots)

}
\arguments{
   \item{covstruc}{Output from the multivariable LDSC function of Genomic SEM}
   \item{model}{ESEM model to be specified using lavaan notation}  
  \item{rotation}{What type of rotation to use for the EFA. Default is geomin}
  \item{rotation.args}{Optional arguments used to run the EFA}
   \item{std.lv}{Optional argument to denote whether all latent variables are standardized using unit variance identification (default = FALSE)}
   \item{toler}{Optional argument to set lower tolerance for matrix inversion used to produce sadnwich corrected standard errors. (default = FALSE)}
}

\value{
  The function estimates an ESEM model, along with model fit indices, using output from GenomicSEM LDSC.
}

\examples{

}
