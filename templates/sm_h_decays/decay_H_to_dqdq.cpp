
// special case for H -> Fd Fd
template<>
double CLASSNAME::get_partial_width<H,dq,bar<dq>::type>(
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

   const double mH = context.mass<H>(indexIn);
   const double mdq = context.mass<dq>(indexOut1);

   // TODO: add off-shell decays?
   if (mH < 2.*mdq) {
      return 0.;
   }

   const double g3 = MODELPARAMETER(g3);
   const double Nf = number_of_active_flavours(mH);
   const double alpha_s_red = Sqr(g3)/(4*Sqr(Pi));
   const double mtpole = qedqcd.displayPoleMt();

   const double deltaqq = calc_deltaqq(alpha_s_red, Nf);

   const double lt = Log(Sqr(mH/mtpole));
   const double lq = Log(Sqr(mdq/mH));
   const double deltaH2 = Sqr(alpha_s_red) * (1.57 - 2.0/3.0*lt + 1.0/9.0*Sqr(lq));

   // SM expression + pure BSM 1L corrections

   const double phase_space = 1./(8.*Pi) * beta(mH, mdq, mdq);
   const double flux = 1./(2.*mH);
   const double color_factor = 3;


   return flux * phase_space * color_factor *
      amplitude_squared<H, bar<dq>::type, dq>(context, indexIn, indexOut1, indexOut2) *
      (1. + deltaqq);
}