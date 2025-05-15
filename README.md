# SCR-Count-Prior-Data-Augmentation
Count Prior Data Augmentation for SCR models.

This repository demonstrates 3 different approaches to jointly update N and z when using what I am calling
"count prior data augmentation". The distinction of this data augmentation approach is that instead of
inducing a distribution on the z's where N is a derived variable, we work with the distribution of N directly.
However, we still need z indicator variables to represent individual effects.
Let's say N ~ Poisson(lambda) and we assume N individuals can be
allocated to M (data augmentation limit) z indices at random. Then the prior for z[1:M] is 1/(M choose N). We can initialize N and z such
that N <- sum(z) and update N and z jointly and maintaining this sum constraint.

I provide three approaches described below. To date, I have been using approach 1 where the asymmetric proposal probabilities
are canceled by the prior ratio of z[1:M], but approach 2 is somewhat more efficient. These zSamplers are in the NimbleFunctions file.


Approach 1: Update 1 z at a time, but multiple updates per iteration.
add update, select one of M-N z=0 inds, 
subtract update, select from N z=1 inds, autoreject if select a captured individual.
This is what I have used historically, proposal probs cancel with z prior, 1/choose(M,N)

Approach 2: Update 1 z at a time, but multiple updates per iteration.
add update, select one of M-N z=0 inds,
subtract update, select from N-n.det z=1 inds, so detected individuals are not selected
since they are always rejected. Asymmetric proposal probs need to be accounted for here.
Approach 2 produces a greater effective sample size per unit time than Approach 1. 
I like this one best.

Approach 3: Same as Approach 2, but we update multiple z's at a time and this number is tuned.
We can also do multiple multi z updates per iteration
add update, select one of M-N z=0 inds, 
subtract update, select from N-n.det z=1 inds, so detected individuals are not selected
since they are always rejected. Asymmetric proposal probs need to be accounted for here.
z.ups <- round(M*0.25) # how many N/z proposals per iteration? Not sure what is optimal, setting to 25% of M here.
reversible jump approaches typically only update N once per iteration. This is not optimal,
particularly with large individual heterogeneity in detection probability as is typical in SCR.
Optimal number to update at once is greater as cumulative capture probability decreases.