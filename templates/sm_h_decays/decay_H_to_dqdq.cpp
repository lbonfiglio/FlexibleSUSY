// template specialization for the H -> Fd Fd case


template<>
double CLASSNAME::get_partial_width<H,bar<dq>::type,dq>(
   const ContextName& context,
   typename field_indices<H>::type const& indexIn,
   typename field_indices<dq>::type const& indexOut1,
   typename field_indices<dq>::type const& indexOut2
   ) const
{
   // TODO: should we take the off-diagonal case at all?
   //       or should this never happen and we should crash
   if(!boost::range::equal(indexOut1, indexOut2))
      return 0.;
//    BOOST_ASSERT_MSG(boost::range::equal(indexOut1, indexOut2), 
      // "Template specialization for H -> Fd1 bar[Fd2] is only valid for Fd1 = Fd2"
//    );

   const double mHOS = context.physical_mass<H>(indexIn);
   const double mdqDR = context.mass<dq>(indexOut1);
   const double mdqOS = context.physical_mass<dq>(indexOut1);
   BOOST_ASSERT_MSG(is_zero(mdqDR) || is_zero(mdqOS),
                    "Quarks should not be massless");
      const auto x = 4.*mdqOS*mdqOS/(mHOS*mHOS);

   // TODO: add off-shell decays?
   if (mHOS < 2.*mdqDR) {
      return 0.;
   }

   const double alpha_s = get_alphas(context);
   const double alpha_s_red = alpha_s/(Pi);
   const double Nf = number_of_active_flavours(mHOS);
   const double mtpole = qedqcd.displayPoleMt();

   const double deltaqqOS = 
      4./3 * alpha_s_red * calc_DeltaH(sqrt(1. - x)) +
      calc_deltaqq(alpha_s_red, Nf);

   const double deltaqqDR = calc_deltaqq(alpha_s_red, Nf) + 4./3*alpha_s_red* calc_DeltaH(sqrt(1. - 4.*mdqDR*mdqDR/(mHOS*mHOS)));

   // chiral breaking correctios
   // TODO: probably shouldn't be applied in case of CP-breaking 
   double deltaH2 = 0.;
   if(!info::CPViolationInHiggsSector) {
      const double lt = Log(Sqr(mHOS/mtpole));
      const double lq = Log(Sqr(mdqDR/mHOS));
      deltaH2 = Sqr(alpha_s_red) * (1.57 - 2.0/3.0*lt + 1.0/9.0*Sqr(lq));
   }

   const double flux = 1./(2.*mHOS);
   const double color_factor = 3;
   const double phase_spaceDR = 1./(8.*Pi) * beta(mHOS, mdqDR, mdqDR);
   const double phase_spaceOS = 1./(8.*Pi) * beta(mHOS, mdqOS, mdqOS);
   
   // DRbar or MSbar coupling
   const double amp2DR = amplitude_squared<H, bar<dq>::type, dq>(context, indexIn, indexOut1, indexOut2);
   const double amp2OS = amplitude_squared<H, bar<dq>::type, dq>(context, indexIn, indexOut1, indexOut2)* pow(mdqOS / mdqDR, 2);

   return flux * color_factor *
          (
             // low mass limit
             (1 - x) * phase_spaceDR * amp2DR * (1. + deltaqqDR + deltaH2) +
             x * phase_spaceOS * amp2OS  *
                (1. + deltaqqOS + deltaH2));
}
