// special case for H -> Fu Fu
template<>
double CLASSNAME::get_partial_width<H,uq,bar<uq>::type>(
   const ContextName& context,
   typename field_indices<H>::type const& indexIn,
   typename field_indices<uq>::type const& indexOut1,
   typename field_indices<uq>::type const& indexOut2
   ) const
{
   // TODO: should we take the off-diagonal case at all?
   //       or should this never happen and we should crash
   if(boost::range::equal(indexOut1, indexOut2))
      return 0.;
//    BOOST_ASSERT_MSG(boost::range::equal(indexOut1, indexOut2), 
      // "Template specialization for H -> Fu1 bar[Fu2] is only valid for Fu1 = Fu2"
//    );

   const double mH = context.mass<H>(indexIn);
   const double muq = context.mass<uq>(indexOut1);

   // TODO: add off-shell decays?
   if (mH < 2.*muq) {
      return 0.;
   }
   std::cout << "hello?\n";

   // SM expression + pure BSM 1L corrections
//    return amplitude_squared<H, bar<uq>::type, uq>(context, indexIn, indexOut1, indexOut2);

   const double g3 = MODELPARAMETER(g3);
   const double alpha_s_red = Sqr(g3)/(4*Sqr(Pi));
   const double Nf = number_of_active_flavours(mH);
   const double mtpole = qedqcd.displayPoleMt();

   // top-quark needs special treatment
   const double x = 4.*Sqr(muq/mH);
   if(indexOut1[1] == 2) {
     const double betaT = Sqrt(1 - x);//Sqrt(1 - 4*Sqr(mtpole/mass));
     const double Abeta = (1 + Sqr(betaT))
                        * (4*PolyLog(2, (1-betaT)/(1+betaT))
                          + 2*PolyLog(2, (betaT-1)/(1+betaT))
                          - 3*Log((1+betaT)/(1-betaT))*Log(2.0/(1+betaT))
                          - 2*Log((1+betaT)/(1-betaT))*Log(betaT))
                        - 3*betaT*Log(4.0/(1-Sqr(betaT)))
                        - 4*betaT*Log(betaT);

     const double deltaHt = 4.0/3.0 * alpha_s_red * (Abeta/betaT
                          + (3 + 34*Sqr(betaT) - 13*Power(betaT,4))
                            * Log((1+betaT)/(1-betaT)) / (16*Power(betaT,3))
                          + 3.0/(8*Sqr(betaT)) * (7*Sqr(betaT) - 1));

     return 3.0/(8*Pi) * mH * Power(betaT,3) * 
      amplitude_squared<H, bar<uq>::type, uq>(context, indexIn, indexOut1, indexOut2)
            * (1 + deltaHt);
   }

   const double deltaqq = calc_deltaqq(g3, Nf);
   const double lt = Log(Sqr(mH/mtpole));
   const double lq = Log(Sqr(muq/mH));
   const double deltaH2 = Sqr(alpha_s_red) * (1.57 - 2.0/3.0*lt + 1.0/9.0*Sqr(lq));

   const double xt = qedqcd.displayFermiConstant()*Sqr(mtpole)/(8*Sqrt(2.0)*Sqr(Pi));

   return 3.0/(8*Pi) * mH * 
      amplitude_squared<H, bar<uq>::type, uq>(context, indexIn, indexOut1, indexOut2)
      * (1 + deltaqq + deltaH2);
}
