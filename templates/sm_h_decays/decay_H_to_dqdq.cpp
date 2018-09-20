// template specialization for the H -> Fd Fd case

// TODO: we need to distinguish between scalar and scalar-pseudoscalar-mixture Higgses
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

   // TODO: add off-shell decays?
   if (mHOS < 2.*mdqDR) {
      return 0.;
   }

   const double alpha_s = get_alphas(context);
   const double alpha_s_red = alpha_s/(Pi);
   const double Nf = number_of_active_flavours(mHOS);
   const double mtpole = qedqcd.displayPoleMt();

   // higher order corrections in the chiral limit
   const double deltaqq = calc_deltaqq(alpha_s_red, Nf);

   // chiral breaking correctios
   // TODO: probably shouldn't be applied in case of CP-breaking 
   double deltaH2 = 0.;
   if(!info::CPViolationInHiggsSector) {
      const double lt = Log(Sqr(mHOS/mtpole));
      const double lq = Log(Sqr(mdqDR/mHOS));
      deltaH2 = Sqr(alpha_s_red) * (1.57 - 2.0/3.0*lt + 1.0/9.0*Sqr(lq));
   }

   const double flux = 1./(2.*mHOS);
   const double phase_space = 1./(8.*Pi) * beta(mHOS, mdqDR, mdqDR);
   const double color_factor = 3;

   return flux * phase_space * color_factor *
      amplitude_squared<H, bar<dq>::type, dq>(context, indexIn, indexOut1, indexOut2) *
      (1. + deltaqq);
}
