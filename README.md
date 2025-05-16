# SCR-Count-Prior-Data-Augmentation
Count Prior Data Augmentation for SCR models.

This repository demonstrates 3 different approaches to jointly update N and z when using what I am calling
"count prior data augmentation". The distinction of this data augmentation approach is that instead of
specifying independent distributions on the individual z's that induce a distribution on N, we work with the distribution of N directly.
However, we still need z indicator variables to represent individual effects.
Let's say N ~ Poisson(lambda) and we assume N individuals can be
allocated to M (data augmentation limit) z indices at random. Then the prior for z[1:M] is 1/(M choose N). We can initialize N and z such
that N <- sum(z) and update N and z jointly and maintaining this sum constraint.

I provide three approaches described below. To date, I have been using approach 1 where the asymmetric proposal probabilities
are canceled by the prior ratio of z[1:M], but approach 2 is somewhat more efficient. These zSamplers are in the NimbleFunctions file.

In all 3 approaches, we choose to add or subtract individuals with probability 0.5. Each approach differs in how we choose individuals
to add/subtract.

Approach 1: Update 1 z at a time, but use multiple updates per iteration, say 25% of M.
add update: select one of M-N z=0 inds, 
subtract update: select from N z=1 inds, autoreject if select a captured individual that cannot be turned off.
This is what I have used historically, proposal probs cancel with z prior, 1/choose(M,N)

Approach 2: Update 1 z at a time, but multiple updates per iteration.
add update: select one of M-N z=0 inds,
subtract update: select from N-n.det z=1 inds, so detected individuals are not selected
since they are always rejected. Asymmetric proposal probs need to be accounted for here.
Approach 2 produces a greater effective sample size per unit time than Approach 1 since we never select
individuals that cannot be turned off. I like this one best.

Approach 3: Same as Approach 2, but we update multiple z's at a time and this number is tuned.
We can also do multiple multi z updates per iteration

This is similar to reversible jump approaches that typically only update N once per iteration. This is not optimal,
particularly with large individual heterogeneity in detection probability as is typical in SCR.
Optimal number to update at once is greater as cumulative capture probability decreases.