* This is a declare file which should be included (with `#include declare.inc`) in your FORM program before using analytic topology files or inserted Feynman rules.

* s1m1 etc are massive propagators 1/(M1^2-p1^2) when using Minkowski momenta
auto symbol s1,...,s9;
* Ancillary symbols in analytic topology files
auto symbol n0,...,n9, D0,...,D9;
* This flag is mutlipled to integrals with an unity integrand. In DimReg they are set to zero.
symbol FLAGSCALELESS;
* Internal momenta
auto vector p0,...,p9, P0,...,P9;
* External momenta
auto vector q1,...,q9;

* add_index_1 replacements
auto index al1,...,al9;
* add_index_2 replacements
auto index be1,...,be9;

* lorentz_index replacements (external legs)
auto index mu1,...,mu9;
* lorentz_index replacements (internal legs)
auto index nu1,...,nu9;

* spinor_index replacements (external legs)
auto index i1,...,i9;
* spinor_index replacements (internal legs)
auto index j1,...,j9;

* fundamental_index replacements (external legs)
auto index fe1,...,fe9;
* fundamental_index replacements (internal legs)
auto index fi1,...,fi9;

* colour_index replacements (external legs)
cf a;
* colour_index replacements (external legs)
cf b;
