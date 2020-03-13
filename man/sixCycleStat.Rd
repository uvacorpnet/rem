\name{sixCycleStat}
\alias{sixCycleStat}
\alias{sixCycle}
\title{Calculate six cycle statistics}
%
\description{Calculate the endogenous network statistic \code{sixCycle} that 
measures the tendency for events to close four cycles in two-mode event sequences.}
\usage{
sixCycleStat(data, time, sender, target, halflife, 
    weight = NULL, 
    eventtypevar = NULL, 
    eventtypevalue = 'standard', 
    eventfiltervar = NULL, 
    eventfilterAB = NULL, eventfilterAJ = NULL, 
    eventfilterIB = NULL, eventfilterIJ = NULL,
    eventvar = NULL,
    variablename = 'fourCycle',
    returnData = FALSE,
    dataPastEvents = NULL,
    showprogressbar = FALSE, 
    inParallel = FALSE, cluster = NULL
)
}
\arguments{
\item{data}{ A data frame containing all the variables.}

\item{time}{ Numeric variable that represents the event sequence. The variable
has to be sorted in ascending order.}

\item{sender}{ A string (or factor or numeric) variable that represents the sender of the event.}

\item{target}{ A string (or factor or numeric) variable that represents the target of the event.}

\item{halflife}{ A numeric value that is used in the decay function.
The vector of past events is weighted by an exponential decay function using the specified halflife. The halflife parameter determins after how long a period the event weight should be halved. E.g. if \code{halflife = 5}, the weight of an event that occured 5 units in the past is halved. Smaller halflife values give more importance to more recent events, while larger halflife values should be  used if time does not affect the sequence of events that much.}

\item{weight}{ An optional numeric variable that represents the weight of each event. If \code{weight = NULL} each event is given an event weight of \code{1}.
}

\item{eventtypevar}{ An optional variable that represents the type of the event. Use \code{eventtypevalue} to specify how the \code{eventtypevar} should be used to filter past events.
}

\item{eventtypevalue}{ An optional value (or set of values) used to specify how paste events should be filtered depending on their type. \code{'standard'} is implemented.}

\item{eventfiltervar}{ An optinoal variable that allows filtering of past events using an event attribute; not implemented.}

\item{eventfilterAB}{ An optional value used to specify how 
paste events should be filtered depending on their attribute; not implemented.}

\item{eventfilterAJ}{ see \code{eventfilterAB}.}

\item{eventfilterIB}{see \code{eventfilterAB}.}

\item{eventfilterIJ}{see \code{eventfilterAB}.}

\item{eventvar}{ An optional dummy variable with 0 values for null-events and 1 values for true events. If the \code{data} is in the form of counting process data, use the \code{eventvar}-option to specify which variable contains the 0/1-dummy for event occurrence. If this variable is not specified, all events in the past will be considered for the calulation of the four cycle statistic, regardless if they occurred or not (= are null-events). Misspecification could result in grievous errors in the calculation of the network statistic.}

\item{variablename}{ An optional value (or values) with the name the four cycle statistic variable should be given. To be used if \code{returnData = TRUE}.}

\item{returnData}{ \code{TRUE/FALSE}. Set to \code{FALSE} by default. The new  variable(s) are bound directly to the \code{data.frame} provided and the data frame is returned in full.}

\item{dataPastEvents}{ An optional \code{data.frame} with the following variables:
	column 1 = time variable, 
	column 2 = sender variable, 
	column 3 = target on other variable (or all "1"), 
	column 4 = weight variable (or all "1"), 
	column 5 = event type variable (or all "1"), 
	column 6 = event filter variable (or all "1"). Make sure that the data frame does not contain null events. Filter it out for true events only.}

\item{showprogressbar}{\code{TRUE/FALSE}. To be implemented.}

\item{inParallel}{ \code{TRUE/FALSE}. An optional boolean to specify if the loop should be run in parallel.}

\item{cluster}{ An optional numeric or character value that defines the cluster. By specifying a single number, the cluster option uses the provided number of nodes to parallellize. By specifying a cluster using the \code{makeCluster}-command in the \code{doParallel}-package, the loop can be run on multiple nodes/cores. E.g., \code{cluster = makeCluster(12, type="FORK")}.}

}
\details{
The \code{sixCycleStat()}-function calculates an endogenous statistic that measures whether events have a tendency to form six cycles.

The effect is further described in the following paper:

D. Valeeva, F.W. Takes and E.M. Heemskerk, The duality of firms and directors in board interlock networks: A relational event modeling approach, Social Networks 62: 68-79, Elsevier, 2020.

}
% \value{
% 
% }
% \references{
% 
% }
% \note{
% 
% }
\author{
Diliara Valeeva, Frank Takes and Eelke Heemskerk of the University of Amsterdam's CORPNET group  \email{corpnet@uva.nl}
}
\seealso{
\link{rem-package}
}
\examples{See fourCycleStat() examples.
}
%\keyword{key}




