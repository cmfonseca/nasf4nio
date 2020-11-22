In addition folder you can find documentation for segment operations implemented in symmtsp.c(initSegment,
randomMoveTowards, findSuitableBreakPoint and randomMoveTowards).

Also there are some defined problem configurations taken from:
https://people.sc.fsu.edu/~jburkardt/datasets/tsp/tsp.html


To start a program compile the project and run on of symmtsp-ils, symmtsp-sga or symmtsp-shd
with program Arguments.
For symmtsp-ils and symmtsp-shd example of arguments are: addition\att48.txt 100     where first argument is path to the problem configuration file and second argument
                                                                                     is number of iterations
For  symmtsp-sga example of argumens are:    addition\att48.txt 10 100  here second argument is population size and third argument is number of iterations