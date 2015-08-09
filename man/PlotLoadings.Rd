\name{PlotLoadings}
\alias{PlotLoadings}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to plot loadings in nonstandard but useful format
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PlotLoadings(L_pj, whichfactor, Nspecies, addtitle = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{L_pj}{
%%     ~~Describe \code{L_pj} here~~
}
  \item{whichfactor}{
%%     ~~Describe \code{whichfactor} here~~
}
  \item{Nspecies}{
%%     ~~Describe \code{Nspecies} here~~
}
  \item{addtitle}{
%%     ~~Describe \code{addtitle} here~~
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
function (L_pj, whichfactor, Nspecies, addtitle = TRUE) 
{
    plot(1, type = "n", xlim = c(1, nrow(L_pj) + 0.5), ylim = range(L_pj) + 
        diff(range(L_pj)) * c(0, 0.1), xlab = "", ylab = "", 
        xaxt = "n")
    if (addtitle == TRUE) 
        mtext(text = paste("Factor", whichfactor), side = 3, 
            line = 0.1, adj = 0)
    abline(h = 0)
    for (p in 1:Nspecies) {
        lines(y = c(0, L_pj[p, whichfactor]), x = rep(p, 2), 
            lwd = 5)
        text(x = (p + 0.2), y = 0 + (0.02 * max(L_pj[p, whichfactor])), 
            label = rownames(L_pj)[p], srt = 90, pos = 4, cex = 1.2)
    }
    legend("top", legend = paste0("Proportion of explained variance= ", 
        round(100 * sum(L_pj[, whichfactor]^2)/sum(L_pj^2), 1), 
        "\%"), bty = "n")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
