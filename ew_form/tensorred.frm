*TENSOR REDUCTION
* contract all the scalar product pi^2, etc
Identify Vec(ind?,momen?)*Vec(ind?,momen?) = momen.momen;
Identify Vec(ind?,momen?)^2 = momen.momen;

*tensor reduction for rank 4 with symmetries
Identify Vec(ind1?,momen?)*Vec(ind2?,momen?)*Vec(ind3?,momen1?)*Vec(ind4?,momen1?) =
1/(d *(d^2 + d - 2))*(- momen.momen * momen1.momen1 + d* (momen.momen1)^2)*(d_(ind1,ind3)*d_(ind2,ind4) + d_(ind1,ind4)*d_(ind2,ind3));
Identify Vec(ind1?,momen?)*Vec(ind3?,momen1?)*Vec(ind2?,momen?)*Vec(ind4?,momen1?) =
1/(d *(d^2 + d - 2))*(- momen.momen * momen1.momen1 + d* (momen.momen1)^2)*(d_(ind1,ind3)*d_(ind2,ind4) + d_(ind1,ind4)*d_(ind2,ind3));
Identify Vec(ind1?,momen?)*Vec(ind2?,momen1?)*Vec(ind3?,momen1?)*Vec(ind4?,momen1?) =
1/(d^2+2*d)*(momen.momen1)*(momen1.momen1)*(d_(ind1,ind2)*d_(ind3,ind4)+d_(ind1,ind4)*d_(ind3,ind2)+d_(ind1,ind3)*d_(ind2,ind4));


*tensor reduction for rank 4 - no symmetries, most generic 
*Identify Vec(ind1?,momen1?)*Vec(ind2?,momen2?)*Vec(ind3?,momen3?)*Vec(ind4?,momen4?)=
1/(d* (-2 + d + d^2))*(d_(ind1, ind4)* d_(ind2, ind3) *((1 + d)* momen1.momen4 * momen2.momen3 - 
   momen1.momen3 * momen2.momen4 - momen1.momen2 * momen3.momen4) + 
  d_(ind1, ind3) * d_(ind2, ind4)* (-momen1.momen4 momen2.momen3 + 
  (1 + d)* momen1.momen3 * momen2.momen4 - momen1.momen2 * momen3.momen4) + 
  d_(ind1, ind2) * d_(ind3, ind4)* (-momen1.momen4 * momen2.momen4 - 
     momen1.momen3 * momen2.momen4 + (1 + dim) *momen1.momen2 * momen3.momen4))




 