# SCR-N-Prior-Data-Augmentation
N-Prior Data Augmentation for SCR models.

This repository contains methods for fitting SCR models via data augmentation where the prior is placed on N instead of the inclusion
indicators \bm{z}. Schofield and Barker (2014) distinguish the complete data likelihoods where N or z (they use w) are parameters by CDL-N and CDL-w.
Typical data augmentation (Royle et al. 2007) is an example of CDL-z. S&B also introduced a data augmentation approach for CDL-N.
I stumbled across a similar approach, but only recently noticed how it was related to the approach of S&B.

https://link.springer.com/article/10.1007/s10651-013-0262-3

Both my approach and that of S&B assume that N ~ Poisson(lambda) and require a distribution for [\bm{z}|N] to explain how the N individuals
are allocated to the indices of \bm{z}. S&B assume that the individuals truly in the population occupy indices 1,...,N and the "false" individuals
occupy indices N+1,...,M. Further, due to the mechanics of their N/z update as specified in the BUGS code, they need to account for the ordering of the true and false individuals.
This leads to P(\bm{z}|N) = 1/((M choose N)(N! - n!) (M - N)!), accounting for the number of ways to choose N true individuals from M total individuals,
the number of ways to order the uncaptured, true individuals and the number of ways to order the false individuals.

My approach does not assume any ordering of the individuals. I assume individuals are randomly allocated to the z indices and P(\bm{z}|N) = 1/((M choose N)).
However, this requires a custom update that turns individuals on/off at random with respect to their indices. I provide three slightly different ways to do this.
To date, I have been using approach 1 where the asymmetric proposal probabilities for selecting the individuals to turn on or off are canceled 
by the prior combinatorial terms, but approach 2 is somewhat more efficient. These zSamplers are in the NimbleFunctions file.

What are the advantages of using this approach? First, this approach is required (I believe) for latent ID SCR models where individuals 
allocated samples do not stay at the top of the capture history. These individuals could not be turned off with the S&B approach. Second, it yields a 
the effective sample size per iteration is larger, and as implemented, this approach runs about twice as fast. The greater effective sample size is likely related to the
fact that individuals must be turned on/off consecutively. For example, if the next individual to turn on has a capture probability that is too high to be consistent with the model parameters
(e.g., it is very close to one or more traps), you will not turn on this individual or any subsequent individuals on this iteration. In my approach, we randomly select individuals
to turn on/off individually, so 1 individual does not effect the acceptance of other individuals' indicators. 
Then the faster run time is probably fully explained by the fact that I am using an efficient Metropolis-Hastings approach instead
of the nimble default sampler for N, which is a slice sampler.

What are the advantages of the S&B approach? It works without a custom update. The way the step function is used to fill in the z values is very slick!

What are the advantages of using a prior on N instead of z? Mainly, an exact Poisson specification without relying on M being very large relative to N, which in practice never happens.
It also introduces a way to do Jolly-Seber estimation much faster and use Poisson distributions for recruits.


Some Augustine approach details:

In all 3 approaches, we choose to add or subtract individuals with probability 0.5. Each approach differs in how we choose individuals
to add/subtract.

Approach 1: Update 1 z at a time, but use multiple updates per iteration, say 25% of M.
add update: select one of M-N z=0 inds, 
subtract update: select from N z=1 inds, autoreject if select a captured individual that cannot be turned off.
This is what I have used historically, proposal probs cancel with z prior combinatorial terms.

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

Some Schofield and Barker approach details:

The S&B approach can be fit using "testscript SCR Poisson SB.R" which uses the model file "NimbleModel SCR SB.R".
You can compare our approaches on the same data sets with "testscript SCR Poisson SB.R"

S&B also show how to use a Poisson prior using CDL-z by specifying individual distributions on each z indicators that induce a Poisson
distribution on N. I don't have that in this repository.
