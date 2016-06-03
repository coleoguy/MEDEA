\name{modelMedea}
\alias{modelMedea}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulates a MEDEA system similar to that observed in Tribolium castaneum
}
\description{
Returns the results of forward time individual based population genetic simulations.
}
\usage{
modelMedea(max, N, sims, loss, gain, bottleneck, bottleneck.freq, mon.both, fate, pop.init, s, mutation, verbose}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{max}{numeric - maximum number of generations to simulate}
  \item{N}{population size to simulate}
  \item{sims}{the number of simulations to run}
  \item{loss}{the per generation per locus probability of losing a MEDEA allele}
  \item{gain}{the per generation per locus probability of gaining a MEDEA allele}
  \item{bottleneck}{the proportion of the population that is lost in a bottleneck event}
  \item{bottleneck.freq}{the frequency of bottleneck events i.e. 7 would mean every 7th genertion a bottleneck occurs}
  \item{mon.both}{if TRUE then the function will return a list the elements of which are matrices reporting the frequency of both MEDEA alleles in every generation of each simulation.}
  \item{fate}{if TRUE then the function will return a matrix with the columns representing the iterations and the rows representing the fate of MEDEA factor A and B.}
  \item{pop.init}{a numeric vector of length 16 describing the initial population.  The values should be in the order: ABAB, ABAb, ABaB, ABab, AbAB, AbAb, AbaB, Abab, aBAB, aBAb, aBaB, aBab, abAB, abAb, abaB, abab.}
  \item{s}{this is the selection coefficient of the Medea allele fitness is modeled to be additive: AABB fitness = 1-4*s, AABb fitness = 1-3*s, AAbb fitness = 1-2*s; Aabb fitness = 1-s; aabb fitness = 1}
  \item{mutation}{if TRUE then gain and loss mutations of MEDEA factors are modeled.. currently too slow to use.}
  \item{verbose}{print many print statements for troubleshooting}
}

}
\details{
If this ever is to be published lots of info here about how we imagine the model to work
}
\value{
returns either a matrix or list as described above}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Heath Blackmon}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}