\name{RIOM}
\alias{RIOM}
\title{Function for calculating warning accuracy and sensitivity weights in randomization inference with outcome misclassification}
\usage{
RIOM(treat.ind, outcome.ind, index, alpha = 0.05, null.hypothesis, type.random, timelimit = 1000, gap = 0.00001)
}
\description{
Function for calculating warning accuracy and sensitivity weights in randomization inference with outcome misclassification, which works for
both Fisher's sharp null and Neyman's weak null.
}
\arguments{
 \item{treat.ind}{A N-length vector of the treatment indicators of all N units in the study, 1 if treated and 0 if not. The n-th entry is the treatment indicator of unit n.}
\item{outcome.ind}{A N-length vector of the measured binary outcomes. The n-th entry is the measured outcome of unit n.}
\item{index}{A N-length vector: the n-th entry is the index of the stratum of unit n. Suppose that there are I strata in total, then the value of each entry is an integer between 1 and I.}
\item{alpha}{The alpha level of the two-sided test. The default is 0.05.}
\item{null.hypothesis}{The null hypothesis of interest. Two options: "sharp" = Fisher's sharp null; "weak" = Neyman's weak null.}
\item{type.random}{The type of the randomization design. Two options: "1" = a type 1 randomization design; "2" = a type 2 randomization design. See Heng and Shaw (2022) for details.}
\item{timelimit}{The limit for the runtime in seconds of the function. The default is 1000 seconds.}
\item{gap}{The tolerable gap of the upper bound and lower bound of the optimial solution. The default is 0.00001 (i.e., we want an exact warning accuracy).}
}
\value{
\item{p-value}{Two-sided p-value based on measured outcomes.}
\item{Difference-in-means Estimate}{The point estimate of the average treatment effect based on measured outcomes (only for weak null).}
\item{Confidence Interval}{The confidence interval of the average treatment effect based on measured outcomes (only for weak null).}
\item{Warning Accuracy}{The warning accuracy given the observed data and alpha level.}
\item{Minimal Alteration Number}{The minimal alteration number given the observed data and alpha level.}
\item{Sensitivity Weights}{The four sensitivity weights given the observed data and alpha level.}
\item{Runtime}{The total computation time in seconds.}
}
