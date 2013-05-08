\name{createBins}
\alias{createBins}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
createBins(binsize, genome = "hg19")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{binsize}{
%%     ~~Describe \code{binsize} here~~
}
  \item{genome}{
%%     ~~Describe \code{genome} here~~
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
function (binsize, genome = "hg19") 
{
    genome.build <- as.integer(gsub("[^0-9]", "", genome))
    if (genome.build \%in\% c(19, 37)) {
        genome <- c(249250621, 243199373, 198022430, 191154276, 
            180915260, 171115067, 159138663, 146364022, 141213431, 
            135534747, 135006516, 133851895, 115169878, 107349540, 
            102531392, 90354753, 81195210, 78077248, 59128983, 
            63025520, 48129895, 51304566, 155270560, 59373566)
        names(genome) <- c(1:22, "X", "Y")
    }
    else if (genome.build \%in\% c(18, 36)) {
        genome <- c(247249719, 242951149, 199501827, 191273063, 
            180857866, 170899992, 158821424, 146274826, 140273252, 
            135374737, 134452384, 132349534, 114142980, 106368585, 
            100338915, 88827254, 78774742, 76117153, 63811651, 
            62435964, 46944323, 49691432, 154913754, 57772954)
        names(genome) <- c(1:22, "X", "Y")
    }
    else {
        stop("Unknown genome: ", genome)
    }
    start <- end <- integer()
    for (chromosome in names(genome)) {
        chromosome.size <- genome[chromosome]
        chromosome.starts <- seq(from = 1, to = chromosome.size, 
            by = binsize * 1000)
        chromosome.ends <- chromosome.starts + binsize * 1000 - 
            1
        chromosome.ends[length(chromosome.ends)] <- chromosome.size
        start <- c(start, chromosome.starts)
        end <- c(end, chromosome.ends)
    }
    bins <- data.frame(chromosome = rep(names(genome), ceiling(genome/(binsize * 
        1000))), start, end, stringsAsFactors = FALSE)
    rownames(bins) <- paste(bins$chromosome, ":", bins$start, 
        "-", bins$end, sep = "")
    bins
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
