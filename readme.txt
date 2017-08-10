---- LATENT CLASS CHOICE MODEL ----

Written by: Neeraj Saxena
emailid: n.saxena@unsw.edu.au


The application of this model can be found in:
1. In the paper "Modelling the Route Choice Behaviour under Stop-&-go Traffic for Different Car Driver Segments" by Neeraj Saxena, Taha H. Rashidi, Vinayak V. Dixit & S. Travis Waller --> under review
2. Chapter 4 of the thesis "Modelling the effect of the number of stop-&-gos on the route choice behaviour of car drivers" by Neeraj Saxena


main.m --> main function which initiates the procedure... This file must be run everytime
Change the NOC field, which represents the number of latent classes

doitr.m --> Calls the fminunc function to carry out non-linear unconstrained minimisation routine
Change the upper and lower bound matrices every time you change NOC. Also activate the starting value matrix (b0) accordingly.

loglik.m --> evaluates log-likelihood.. It in turn calls llgrad2.m

llgrad2.m --> The function which evaluates the likelihood function and the analytical gradient.

makedraws.m --> Takes MLHS draws for the random coefficient



****** Written by Neeraj Saxena; University of NSW; Date: 12/06/2016 ******