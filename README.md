# ShapeAttack
This is an algorithm to produce adversarial examples satisfying a shape constraint, e.g. an adversarial monotonic function.

Author: Kellin Pelrine. Code for RBC model modified from code provided by Fabrizio Zilibotti.

<br/><br/>
The core is a genetic algorithm. It is a gradient-free method and only requires querying a black box. It benefits substantially from parallelization, which is trivial to implement.

![image](https://github.com/kellinpelrine/ShapeAttack/blob/master/Algorithm.JPG)

<br/><br/>
One application: testing robustness of macroeconomic models. Shape constraints allow one to search, for example, for adversarial perturbations of a utility function that keep the utility increasing in consumption. 

This is illustrated by the following figure. Perturbing the utility function in a simple RBC model from blue (log utility) to red, with a monotone increasing constraint, changes the correlation between output and consumption (as determined by the model) from about 0.75 to 0.01. 

![image](https://github.com/kellinpelrine/ShapeAttack/blob/master/Utility_Perturbation.JPG)

More information is contained in the report.
