#-
Off Statistics;
On HighFirst;

#include- form-declarations.h

#ifndef `LOOPS'
	#message Error, please specify the number of loops as LOOPS.
	#terminate
#endif
#ifndef `DIA'
	#message Error, please specify the diagram number as DIA.
	#terminate
#endif

* Load the diagram from the tapir file:
Local d`LOOPS'l`DIA'amp =
	#include- ../ew_tapir/dia/munuenu-`LOOPS'l.dia # d`LOOPS'l`DIA'

#ifndef `BRIDGEMOMENTA'
	#define BRIDGEMOMENTA ""
#endif

* Load the mapped topology, and make the necessary momentum replacements
#include ../ew_tapir/topo/mapping-`LOOPS'l.h # d`LOOPS'l`DIA'
`MOMREPLACEMENT'
.sort

Argument;
`BRIDGEMOMENTA'
EndArgument;

Argument;
#include ../ew_tapir/topo/`INT1' # NUMERATORMOMENTA
EndArgument;
.sort

Identify Mom(?a) = Vec(?a);

*Feynman rules
* auxGamma are gamma matrices, their argument is a Lorentz index and the spinor1 and spinor2.
* Internal spinor indices are called iufo and j. External are called i.

Identify auxPL(i1?,i2?) = cFT(g7,XX,i1,i2)/2;
Identify auxPR(i1?,i2?) = cFT(g6,XX,i1,i2)/2;
Identify auxSlash(?a,i1?,i2?) = cFT(?a,XX,i1,i2);
Identify auxGamma(?a,i1?,i2?) = cFT(?a,XX,i1,i2);
Identify cFT(- p1?, ?a) = - cFT(p1, ?a);

* Start to build spinor line at some arbitrary point:
#do i = 1,4
	Identify, once cFT(?a) = FT`i'(?a);
	Repeat;
		Identify FT`i'(?a,i1?,i2?) * cFT(?b,i2?,i3?) = FT`i'(?a,XX,?b,XX,i1,i3);
		Identify cFT(?b,i3?,i1?) * FT`i'(?a,i1?,i2?) = FT`i'(?b,XX,?a,XX,i3,i2);
		Identify FT`i'(?a,i1?,i1?) = FT`i'(?a);
	EndRepeat;

	Repeat Identify FT`i'(?a,XX,?b) = FT`i'(?a) * FT`i'(?b);
	Identify FT`i' = 1;
#enddo

* The external fermions now result in an FT`i'(i<k>,i<l>), where
* i<k> and i<l> are spinor indices. There may be a g6 or g7:
* take it out, and then keep the external spinor indices:
#do i = 1,4
	#do j = 1,10
	#do k = 1,10
		Identify FT`i'(?a,i`j',i`k') = FT`i'(?a) * FT`i'(i`j',i`k');
		Identify FT`i'(i`j',i`k') = spinIndExt(`i',i`j',i`k');
	#enddo
	#enddo
	Identify FT`i' = 1;
#enddo

*propagators substitutions
* Dtran, Dlong are gauge boson propagators, their arguments are 2 Lorenx indices, the runnin momentum, the gauge and mass.
* Dgoldst is the goldstone propagator, that should have a mass
Identify Dph(ind1?,ind2?,mom?,gaug?) = (d_(ind1, ind2) - (1-gaug) * Vec(ind1,mom) * Vec(ind2,mom) * Den(mom,0,0))* Den(mom,0,0);
Identify Dgoldst(mom?,gaug?,mass?) = Den(mom, gaug, mass);
Identify Dlong(ind1?,ind2?,mom?,gaug?,mass?) = (gaug*Vec(ind1,mom)*Vec(ind2,mom)*Den(mom,gaug,mass))* Den(mom,gaug,mass);
Identify Dtran(ind1?,ind2?,mom?,gaug?,mass?) = (d_(ind1, ind2)-Vec(ind1,mom)*Vec(ind2,mom)*Den(mom,gaug,mass))* Den(mom,gaug,mass);
Identify DH(mom?, mass?) = Den(mom, 0, mass);
Identify Dghost(mom?, gaug?) = Den(mom,0,0);
Identify Dghost(mom?, gaug?, mass?) = Den(mom,gaug,mass);

#do i = 1,`NUMTRACES'
	Identify FT`i'(g6) = g6_(`i');
	Identify FT`i'(g7) = g7_(`i');
	Identify FT`i'(nu?) = g_(`i',nu);
	Identify FT`i'(mom?,mass?) = g_(`i',mom) + gi_(`i')*mass;
#enddo

*split the momenta in the numerator
SplitArg Vec;
Repeat;
	Identify Vec(ind?,mom?,?a) = Vecr(ind,mom) + Vec(ind,?a);
	Identify Vec(ind?) = 0;
EndRepeat;
FactArg Vecr;
Identify Vecr(ind?,x?number_,mom?) = x * Vecr(ind,mom);
Identify Vecr(?a) = Vec(?a);
.sort

*change the label in the loop momenta so that only internal momenta enter tensor reduction
* p1 to p4 are truly internal momenta because we have previously applied the bridges
#do i = 1,4
	Identify Vec(ind?,p`i') = Vecr(ind,p`i');
#enddo

*TENSOR REDUCTION : highest ranks first
* contract all the scalar product pi^2, etc
Identify Vecr(ind?,momen1?)*Vecr(ind?,momen2?) = momen1.momen2;


*tensor reduction for rank 4 - no symmetries, most generic (commented for now as it might not be needed)
*Identify Vecr(ind1?,momen1?)*Vecr(ind2?,momen2?)*Vecr(ind3?,momen3?)*Vecr(ind4?,momen4?)=1/(d* (-2 + d + d^2))*(d_(ind1, ind4)* d_(ind2, ind3) *((1 + d)* momen1.momen4 * momen2.momen3 - momen1.momen3 * momen2.momen4 - momen1.momen2 * momen3.momen4) +   d_(ind1, ind3) * d_(ind2, ind4)* (-momen1.momen4 momen2.momen3 +  (1 + d)* momen1.momen3 * momen2.momen4 - momen1.momen2 * momen3.momen4) + d_(ind1, ind2) * d_(ind3, ind4)* (-momen1.momen4 * momen2.momen4 - momen1.momen3 * momen2.momen4 + (1 + d) *momen1.momen2 * momen3.momen4));


*tensor reduction for rank 4 with symmetries
Identify Vecr(ind1?,momen?)*Vecr(ind2?,momen?)*Vecr(ind3?,momen1?)*Vecr(ind4?,momen1?) = 1/(d *(d^2 + d - 2))*(- momen.momen * momen1.momen1 + d* (momen.momen1)^2)*(d_(ind1,ind3)*d_(ind2,ind4) + d_(ind1,ind4)*d_(ind2,ind3));
Identify Vecr(ind1?,momen?)*Vecr(ind3?,momen1?)*Vecr(ind2?,momen?)*Vecr(ind4?,momen1?) = 1/(d *(d^2 + d - 2))*(- momen.momen * momen1.momen1 + d* (momen.momen1)^2)*(d_(ind1,ind3)*d_(ind2,ind4) + d_(ind1,ind4)*d_(ind2,ind3));
Identify Vecr(ind1?,momen?)*Vecr(ind2?,momen1?)*Vecr(ind3?,momen1?)*Vecr(ind4?,momen1?) =1/(d^2+2*d)*(momen.momen1)*(momen1.momen1)*(d_(ind1,ind2)*d_(ind3,ind4)+d_(ind1,ind4)*d_(ind3,ind2)+d_(ind1,ind3)*d_(ind2,ind4));


*tensor reduction for rank 3
Identify Vecr(ind1?,momen1?)*Vecr(ind2?,momen2?)*Vecr(ind3?,momen3?) = 0;

*tensor reduction for rank 2
Identify Vecr(ind1?,momen1?)*Vecr(ind2?,momen2?) = d_(ind1,ind2) * 1/d * momen1.momen2;

*tensor reduction for rank 1
Identify Vecr(ind1?,momen1?) = 0;
.sort

*split args in denominator 
SplitArg Den;
*set ext mom to zero in the denominators
*Multiply replace_(q1,0);
Identify Den(?a,q1,?b)=Den(?a,?b);
Identify Den(?a,q2,?b)=Den(?a,?b);
Identify Den(?a,q3,?b)=Den(?a,?b);
Identify Den(?a,q4,?b)=Den(?a,?b);
Identify Den(?a,-q1,?b)=Den(?a,?b);
Identify Den(?a,-q2,?b)=Den(?a,?b);
Identify Den(?a,-q3,?b)=Den(?a,?b);
Identify Den(?a,-q4,?b)=Den(?a,?b);

*adjust signs in denominators (-p)^2 = p^2
Identify Den(-p3,gaug?,mass?)=Den(p3,gaug,mass);
Identify Den(-p4,gaug?,mass?)=Den(p4,gaug,mass);

* Simplify the d-dependent coefficients which come from the tensor reduction:
*PolyRatFun prf;
*Denominators dentmp;
Identify dentmp(x?) = prf(1,x);
Identify d^x? = prf(d^x,1);
Identify 1/x?number_ = prf(1,x);
Identify Den(gaug?, mass?) = prf(1,gaug * mass^2);
*join all prf functions into one only -  NOT WORKING
Identify prf(a1?, a2?)*prf(a3?, a4?) = prf(a1* a3, a2* a4);

*set ext mom to zero in the vectors
#do i = 1,4
	Identify Vec(ind?,q`i') = 0;
#enddo

*project in operator basis
Identify g_(1,6_,ind?)*g_(2,6_,ind?)*spinIndExt(1,i1,i4)*spinIndExt(2,i2,i3) = Op;

*simplify denominators for partial fractioning
Identify Den(mom?, 0, 0) = Den(mom, 0);
Identify Den(mom?, 0, mass?) = Den(mom, mass);
Identify Den(mom?, gaug?, mass?) = Den(mom, gaug*mass);

* rewrite the scalar products in terms of denominators to have cancellations between numerator and denominator
#do i = 3,4
	Identify p`i'.p`i'*Den(p`i',0) = 1;
	Identify p`i'.p`i'^x? = 1/(Den(p`i',0))^x;
#enddo

*join all dens in list for partial fractioning 
Identify Den(mom?, mass?)^x? = Int(mom, mass, x);
Identify Den(mom?, mass?) = Int(mom, mass, 1);

Repeat;
Identify Int(?a)*Int(?b) = Int(?a, ?b);
EndRepeat;

*if the same den appear, sum the powers (regardless of the order they appear)
Identify Int(?a, mom?, mass?, x?, mom?, mass?, x1?, ?b) = Int(?a, mom, mass, x+x1, ?b);
Identify Int(?a, mom?, mass?, x?, ?c, mom?, mass?, x1?, ?b) = Int(?a, mom, mass, x+x1,?c, ?b);

*Partial fraction decomposition
Repeat;
Identify Int(mom?, mass?, x?pos_, mom?, mass1?, x1?pos_) = prf(1,mass^2-mass^2)*(Int(mom, mass, x, mom, mass1, x1-1)-Int(mom, mass, x-1, mom, mass1, x1));
EndRepeat;

Bracket Int;
*Print[];
Print +s;
.end


* Write the result to a file:
*.sort
*#write <results/d`LOOPS'l`DIA'.h> "%E", d`LOOPS'l`DIA'amp

*Bracket d2l1,...,d2l270;
*Print +s;
*.end
