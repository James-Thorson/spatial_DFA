\name{InjectParams}
\alias{InjectParams}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to input previous run into tagged-list of input parameters for TMB
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
InjectParams(previous, skeleton)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{previous}{
%%     ~~Describe \code{previous} here~~
}
  \item{skeleton}{
%%     ~~Describe \code{skeleton} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (previous, skeleton) 
{
    new = skeleton
    for (i in 1:4) {
        parname = c("delta_i", "gamma_k", "log_sigma_p", "zinfl_pz")[i]
        new[[parname]] = previous[[parname]]
    }
    for (i in 1:5) {
        parname = c("alpha_j", "phi_j", "loglambda_j", "rho_j", 
            "L_val")[i]
        new[[parname]][1:length(previous[[parname]])] = previous[[parname]]
    }
    new[["logkappa_jz"]][1:nrow(previous[["logkappa_jz"]]), ][] = as.vector(previous[["logkappa_jz"]])
    new[["logkappa_jz"]][nrow(new[["logkappa_jz"]]), ] = previous[["logkappa_jz"]][nrow(previous[["logkappa_jz"]]), 
        ]
    new[["Omega_input"]][, 1:ncol(previous[["Omega_input"]])][1:length(as.vector(previous[["Omega_input"]]))] = as.vector(previous[["Omega_input"]])
    new[["Epsilon_input"]][, 1:dim(previous[["Epsilon_input"]])[2], 
        ][1:length(as.vector(previous[["Epsilon_input"]]))] = as.vector(previous[["Epsilon_input"]])
    return(new)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
