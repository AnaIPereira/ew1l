* Declarations:

* Work in d = 4-2*ep dimensions:
Symbol d,ep;
* Indices declared following this will have dimension d:
Dimension d;


* Symbols used in tapir-generated topology files and defined by a
* file in the tapir source:
* #include- Software/src/tapir/declare.inc
#include- declare.inc

* Functions to hold the propagator indices. "AutoDeclare" means that
* any variable which starts with these strings, are automatically
* defined to be a CFunction:
AutoDeclare CFunction d1l, d2l;


AutoDeclare Function cFT, FT;
Symbol g7, XX, g6

* Kinematic variables and masses
Symbol M1,M2,M3,M4,mass;
Vector mom, mom1;
Cfunction flavourTag;
Symbol fve,fvm,fe,fmu,fta,fvt,fu,fd,fs,fb,fc,ft;

CFunction Mom;

AutoDeclare Vector momen;

*Coupling constants and gauge parameters
Symbol gaug, xia, xiw, xiz;

AutoDeclare Index iufo;
AutoDeclare Symbol ufoGC;

*indices
Symbol mu,spin;
Indices ind,ind1,ind2,ind3,ind4,i,j,nu;

*gamma matrices
Function auxGamma, auxPL, auxPR,auxSlash;

*momentum
CFunction Vec, Vecr;
Vector p1,...,p4;
Vector q1,...,q4;
Symbol p,q;
*CFunction to hold rational polynomials:
CFunction prf;

* Propagator denominators:
CFunction Den,Deng,Dlong,Dph,Dgoldst,Dtran,DH,Dghost;

* gauge parameter
*#define GAUGEa "xia"
*Symbol xia;
*#define GAUGEw "xiw"
*Symbol xiw;
*#define GAUGEz "xiz"
*Symbol xiz;

*Symbol operator basis
Symbol Op, Ev3, Ev5;


* Symbols used in IBP reduction rules from Kira:
CFunction num, den;


* For computing colour factors
Symbol NA,NC,ca,cf;
CFunction GM,prop,a,b;
CFunction colF,colT,colD;
CFunction colFac;


* Maximal number of fermion lines which appear in the diagrams
#define NUMTRACES "2"
* Maximal number of propagators appearing in the diagrams
#define NUMPROPS "4"

