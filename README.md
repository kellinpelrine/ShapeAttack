# ShapeAttack

Author: Kellin Pelrine

This is an algorithm to produce adversarial examples satisfying a shape constraint, e.g. an adversarial monotonic function.

The core is a genetic algorithm. It is a gradient-free method and only requires querying a black box. It benefits substantially from parallelization, which is trivial to implement.

One application: testing robustness of macroeconomic models. Shape constraints allow one to search, for example, for adversarial perturbations of a utility function that keep the utility increasing in consumption.
