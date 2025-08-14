# SCR-N-Prior-Data-Augmentation
N Prior Data Augmentation for SCR models.

This repository contains methods for fitting SCR models via data augmentation where the prior is placed on N instead of the inclusion
indicators, z. Schofield and Barker (2014) distinguish the complete data likelihoods where N or z (they use w) are parameters by CDL-N and CDL-z.
Typical data augmentation (Royle et al. 2007) is an example of CDL-z. Here, I will call this "Bernoulli-z data augmentation" (ignoring the inconvenience that
a Binomial-N prior is induced). S&B introduced a data augmentation approach for CDL-N and an alternative one for CDL-z that specifies independent z priors that 
induce a Poisson prior on N. I stumbled across an approach similar to their CDL-N approach, but only recently noticed how it was related to the approach of S&B.
I explore this below and provide nimble code to fit both approaches.

https://link.springer.com/article/10.1007/s10651-013-0262-3

Both my approach and that of S&B assume that N ~ Poisson(lambda) and require a distribution for [\bm{z}|N] to explain how the N individuals
are allocated to the indices of \bm{z}. S&B assume that the individuals truly in the population occupy indices 1,...,N and the "false" individuals
occupy indices N+1,...,M. Further, due to the mechanics of their N/z update as specified in the BUGS code, 
they need to account for the ordering of the true and false individuals.
This leads to P(\bm{z}|N) = 1/((M choose N)(N - n)!(M - N)!), accounting for the number of ways to choose N true individuals from M total individuals,
the number of ways to order the uncaptured, true individuals and the number of ways to order the false individuals.

My approach does not assume an ordering of the individuals. I assume individuals are randomly allocated to the z indices and P(\bm{z}|N) = 1/((M choose N)).
However, this requires a custom update that turns individuals on/off at random with respect to their indices. I provide three slightly different ways to do this.
To date, I have been using approach 1 where the asymmetric proposal probabilities for selecting the individuals to turn on or off are canceled 
by the prior combinatorial terms, but approach 2 is somewhat more efficient. These zSamplers are in the NimbleFunctions file.

Perhaps a useful distinction in terminology is ordered-z vs. random-z N prior data augmentation.

What are the advantages of using random-z N prior data augmentation? First, it can be used for latent ID SCR models where individuals 
allocated samples do not stay at the top of the capture history and the number detected is not fixed, which is not compatible with the ordered-z prior. Second, it yields an 
effective sample size per iteration that is larger, and as implemented, this approach runs about twice as fast. The greater effective sample size is likely related to the
fact that individuals must be turned on/off consecutively in the ordered-z approach. 
For example, if the next individual to turn on has a capture probability that is too high to be consistent with the model parameters
(e.g., it is very close to one or more traps), you will very likely not turn on this individual or any subsequent individuals on this iteration. In effect, if individual N+1
is a poor individual to propose to turn on, it can block individuals to propose to turn on that are much more likely to be accepted. You cannot skip over individual N+1,
you have to wait until the next iteration(s) where its new activity center location (and the new values of other parameters) may make it more likely to be turned on.

Using random-z N prior data augmentation, we randomly select individuals
to turn on/off individually, so one individual does not affect the acceptance of other individuals' indicators. This is also a feature of Bernoulli-z data augmentation, and 
S&B's w-CDL data augmentation, which I believe is desirable. Then the faster run time is probably fully explained by the fact that I am using an efficient Metropolis-Hastings approach instead
of the nimble default sampler for N, which is a slice sampler.

What are the advantages of ordered-z data augmentation? It works without a custom update which requires more MCMC knowledge and effort to modify for new model structures.
The way the step function is used to fill in the z values is very slick!

What are the advantages of using a prior on N instead of z? Mainly, an exact Poisson specification without relying on M being very large relative to N, which in practice never happens.
It also introduces a way to do Jolly-Seber estimation much faster and use Poisson distributions for recruits.


Some random-z approach details:

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
We can also do multiple multi-z updates per iteration.

Some ordered-z approach approach details:

The ordered-z approach can be fit using "testscript SCR Poisson SB.R" which uses the model file "NimbleModel SCR SB.R".
You can compare our approaches on the same data sets with "testscript SCR Poisson SB.R"