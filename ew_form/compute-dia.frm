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

*On code;

Identify Mom(?a) = Vec(?a);

*Feynman rules
* auxGamma are gamma matrices, their argument is a Lorentz index ans the spinor1 and spinor2.

*Identify flavourTag(?f,i1?,i2?) = delta_(i1,i2);

Identify auxPL(i1?,i2?) = cFT(g7,XX,i1,i2)/2;
Identify auxPR(i1?,i2?) = cFT(g6,XX,i1,i2)/2;
Identify auxSlash(?a,i1?,i2?) = cFT(?a,XX,i1,i2);
Identify auxGamma(?a,i1?,i2?) = cFT(?a,XX,i1,i2);

Identify cFT(- p1?, ?a) = - cFT(p1, ?a);

#do i = 1,2

  Identify once cFT(?a) = FT`i'(?a);
  Repeat;
    Identify FT`i'(?a,i1?,i2?) * cFT(?b,i2?,i3?) = FT`i'(?a,XX,?b,XX,i1,i3);
    Identify cFT(?b,i3?,i1?) * FT`i'(?a,i1?,i2?) = FT`i'(?b,XX,?a,XX,i3,i2);
    Identify FT`i'(?a,i1?,i2?) * cFT(?b,i3?,i1?) = FT`i'(?a,XX,?b,XX,i2,i3);
    Identify cFT(?b,i3?,i1?) * FT`i'(?a,i2?,i3?) = FT`i'(?b,XX,?a,XX,i1,i2);
    Identify FT`i'(?a,i1?,i1?) = FT`i'(?a);
  EndRepeat;
  Repeat;
    Identify FT`i'(XX,?a) = FT`i'(?a);
    Identify FT`i'(?a,XX) = FT`i'(?a);
    Identify FT`i'(?a,XX,?b) = FT`i'(?a) * FT`i'(?b);
  EndRepeat;

  .sort

#enddo

*propagators substitutions
* Dtran, Dlong are fermion propagators, their arguments are 2 Lorenx indices, the runnin momentum, the gauge and mass.
*Dgoldst is the goldstone propagator, that should have a mass
Identify Dph(ind1?,ind2?,?mom,gaug?) = (d_(ind1, ind2) - (1-gaug) * Vec(ind1,?mom) * Vec(ind2,?mom) * Den(?mom,0,0))* Den(?mom,0,0);
Identify Dgoldst(?mom,gaug?,mass?) = Den(?mom, gaug, mass);
Identify Dlong(ind1?,ind2?,?mom,gaug?,mass?) = (gaug*Vec(ind1,?mom)*Vec(ind2,?mom)*Den(?mom,gaug,mass))* Den(?mom,gaug,mass);
Identify Dtran(ind1?,ind2?,?mom,gaug?,mass?) = (d_(ind1, ind2)-Vec(ind1,?mom)*Vec(ind2,?mom)*Den(?mom,gaug,mass))* Den(?mom,gaug,mass);


#do i = 1,`NUMTRACES'
  Identify FT`i'(g7) = g7_(`i');
#enddo

#do i = 1,`NUMTRACES'
	Identify FT`i'(nu?) = g_(`i',nu);
#enddo

#do i = 1,`NUMTRACES'
	Identify FT`i'(i1?,i2?) = 1;
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

*Bracket Den;
*Print[];
*.end


*split args in denominator 
SplitArg Den;

*set external momenta to zero
* this can be safely done as all denominators are either massive or
* have internal momenta going in the loop and therefore there are no singularities (there will be spurious IR divergences though...)
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

Bracket Den;
Print[];
.end

* rewrite the scalar products in terms of denominators
Identify mom.mom1 = mom.mom


* Compute the traces:
*Tracen,1;
*Tracen,2;
*.sort
*Print +s;
*.end

*project in operator basis
Identify g_(1,6_,ind?)*g_(2,6_,ind?) = Op;
.sort

*Print +s;
*.end



* Clean up the notation:
*PolyRatFun prf;
* Kinematics
*Identify q1.q1 = 0;
*Identify q2.q2 = 0;
*Identify q1.q2 = prf(s,2);
*Identify M1^2 = prf(mts,1);
*Identify s^n1? = prf(s^n1,1);
*Identify d^n1? = prf(d^n1,1);
.sort

** Insert the IBP reduction rules from Kira:
*#include ../kira-`LOOPS'l/results/`INT1'/kira_`INT1'.inc
*Identify num(s?) = prf(s,1);
*Identify den(s?) = prf(1,s);
*Identify mts^n1? = prf(mts^n1,1);
*Identify s^n1? = prf(s^n1,1);
*Identify d^n1? = prf(d^n1,1);



* Write the result to a file:
*.sort
*#write <results/d`LOOPS'l`DIA'.h> "%E", d`LOOPS'l`DIA'amp

*Bracket d2l1,...,d2l21;
*Print +s;
*.end
