\name{Rotate_Fn}
\alias{Rotate_Fn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to implement varimax or PCA-inspired rotation, and perform test on whether transformation matrix is permissible
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Rotate_Fn(L_pj, Psi, RotationMethod = "Varimax", testcutoff = 1e-10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{L_pj}{
%%     ~~Describe \code{L_pj} here~~
}
  \item{Psi}{
%%     ~~Describe \code{Psi} here~~
}
  \item{RotationMethod}{
%%     ~~Describe \code{RotationMethod} here~~
}
  \item{testcutoff}{
%%     ~~Describe \code{testcutoff} here~~
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
function (L_pj, Psi, RotationMethod = "Varimax", testcutoff = 1e-10) 
{
    if (RotationMethod == "Varimax") {
        Hinv = varimax(L_pj, normalize = FALSE)
        L_pj_rot = L_pj \%*\% Hinv$rotmat
        Psi_rot = array(NA, dim = dim(Psi))
        for (n in 1:Nknots) Psi_rot[n, , ] = solve(Hinv$rotmat) \%*\% 
            Psi[n, , ]
    }
    if (RotationMethod == "PCA") {
        Eigen = eigen(L_pj \%*\% t(L_pj))
        L_pj_rot = (Eigen$vectors \%*\% diag(sqrt(trunc_machineprec(Eigen$values))))[, 
            1:Nfactors, drop = FALSE]
        rownames(L_pj_rot) = rownames(L_pj)
        Hinv = list(rotmat = corpcor::pseudoinverse(L_pj_rot) \%*\% 
            L_pj)
        Psi_rot = array(NA, dim = dim(Psi))
        for (n in 1:Nknots) Psi_rot[n, , ] = Hinv$rotmat \%*\% 
            Psi[n, , ]
    }
    approx_equal = function(m1, m2, d = 1e-10) (2 * abs(m1 - 
        m2)/mean(m1 + m2)) < d
    if (!all(approx_equal(L_pj \%*\% t(L_pj), L_pj_rot \%*\% t(L_pj_rot), 
        d = testcutoff))) 
        stop("Covariance matrix is changed by rotation")
    for (i in 1:dim(Psi)[[1]]) {
        for (j in 1:dim(Psi)[[3]]) {
            if (!all(approx_equal(L_pj \%*\% Psi[i, , j], L_pj_rot \%*\% 
                Psi_rot[i, , j], d = testcutoff))) 
                stop(paste0("Linear predictor is wrong for site ", 
                  i, " and time ", j))
        }
    }
    if (!all(approx_equal(Hinv$rotmat \%*\% t(Hinv$rotmat), diag(Nfactors), 
        d = testcutoff))) 
        stop("Rotation matrix is not a rotation")
    Return = list(L_pj_rot = L_pj_rot, Psi_rot = Psi_rot, Hinv = Hinv)
    return(Return)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
