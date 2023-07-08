\name{ars}
\alias{ars}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%% 
Adaptive Rejection Sampling 
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
A function to compute the density of the distribution of interest
computed using Adaptive Rejection Sampling.
}
\usage{
main(old_density_function, starting_points = NA, bounds = c(-Inf, Inf), n = 20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{old_density_function}{
  a density function of which to compute sample the density of with 
  the Adaptive Rejection Sampling.
%%     ~~Describe \code{x} here~~
}
\item{starting_points}{
  which points shall the sampling start at. The default is set as NA to 
  indicate if no points are inputted. Then starting points are evaluated
  by the function. 
}
\item{bounds}{
  which bounds the sampling shall use. The default is -Inf and Inf to indicate
  if no bounds are inputted. Then the bounds are evaluated by the function. 
}
\item{n}{
  how many points shall sample. The default is set at 20 points. 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
%% do we want to write any thing here??
}
\value{

`ars` returns a numerical vector of length n with the accepted sampled points.
These points are samples from the density function of interest. 

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\source{
The Adaptive Rejection Sampling method implemented is based on  

W. R. Gilks and P. Wild 1992, 'Adaptive RejectionSamplingforGibbs Sampling', 
_Journal of the Royal Statistical Society,_ Series C (Applied Statistics), 
Vol. 41, No. 2 (1992), pp. 337-348.
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
Arm Wonghirundacha, Daniel Sunko and Inga Huld Ármann. 
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
function (dnorm) 
{
  }
}

main(dnorm)

% add here code for plotting the points 

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
