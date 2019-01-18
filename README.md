# RBF-chaos

Context:

1) Steady-state performance constraints for dynamical models based on RBF networks. Engineering Applications of Artificial Intelligence, vol. 20, issue 7.
2) Chaos from a buck switching regulator operating in discontinuous mode. International Journal of Circuit Theory and Applications, vol. 22. 

Bifurcacao: Builds the bifurcation diagram of an RBF model of the sine map.

Buckmap: Simulates a model of the Poincaré map of a buck switching regulator.

ident_buck_fp: Builds an RBF model of the Poincaré map of a buck switching regulator, constrained to have the same fixed point as the original system. 

ident_sine: Builds constrained RBF models of the sine map with cubic nonlinearities. Models have the same trivial fixed point of the original system and reproduce the symmetries of the remaining fixed points.

ident_sine1: Same, but models now reproduce all symmetrical fixed points.

ident_sine2: Same, but now varying basis function widths.

lyap_rbf: Calculates the variation of the largest Lyapunov exponent of an RBF model of the sine map, given different basis function widths.

modelo_buck: Model of the Poincaré map of a buck switching regulator.

modelo_seno: Model of a sine map with cubic nonlinearities.
