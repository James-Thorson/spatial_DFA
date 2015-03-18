\name{plotLake}
\alias{plotLake}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot lake in Alaska with color-coded sites
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotLake(mapObj, col, site_set, cex = 2, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mapObj}{
%%     ~~Describe \code{mapObj} here~~
}
  \item{col}{
%%     ~~Describe \code{col} here~~
}
  \item{site_set}{
%%     ~~Describe \code{site_set} here~~
}
  \item{cex}{
%%     ~~Describe \code{cex} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (mapObj, col, site_set, cex = 2, ...) 
{
    require(maptools)
    plot(1, type = "n", xlim = attr(mapObj, "bbox")["x", ], ylim = attr(mapObj, 
        "bbox")["y", ], ...)
    plot(mapObj, add = TRUE)
    stMercCoord <- read.table(paste0(DataFile, "bsSiteMercCoord.txt"), 
        head = TRUE)
    LatLong2UTM <- match(site_set, table = as.character(stMercCoord[, 
        "site"]))
    points(stMercCoord$v[LatLong2UTM], stMercCoord$e[LatLong2UTM], 
        pch = 20, cex = cex, col = col)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
