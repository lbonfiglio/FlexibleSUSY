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
   if (!boost::range::equal(indexOut1, indexOut2))
      return 0.;
//    BOOST_ASSERT_MSG(boost::range::equal(indexOut1, indexOut2), 
      // "Template specialization for H -> Fd1 bar[Fd2] is only valid for Fd1 = Fd2"
//    );

   const double mHOS = context.physical_mass<H>(indexIn);
   const double mdqOS = context.physical_mass<dq>(indexOut1);
   const auto xOS = std::pow(mdqOS / mHOS, 2);
   const auto betaOS = std::sqrt(1. - 4. * xOS);

   const double mdqDR = context.mass<dq>(indexOut1);
   const auto xDR = std::pow(mdqDR / mHOS, 2);
   const auto betaDR = std::sqrt(1. - 4. * xDR);

   BOOST_ASSERT_MSG(!is_zero(mdqDR) && !is_zero(mdqOS),
                    "Quarks should not be massless");

   // TODO: add off-shell decays?
   if (mHOS < 2. * mdqDR) {
      return 0.;
   }

   // get HBBbar vertex
   // we don't use amplitude_squared here because we need both this vertex
   // both with running and pole masses
   const auto indices = concatenate(indexIn, indexOut1, indexOut2);
   const auto HBBbarVertexDR =
      Vertex<H, bar<dq>::type, dq>::evaluate(indices, context);
   BOOST_ASSERT_MSG(
      is_zero(HBBbarVertexDR.left() - HBBbarVertexDR.right()),
      "Left and right coupling of CP-even Higgs to fermions should be equal");

   const double alpha_s = get_alphas(context);
   const double alpha_s_red = alpha_s / Pi;
   const double Nf = number_of_active_flavours(mHOS);
   const double mtpole = qedqcd.displayPoleMt();

   const double flux = 1. / (2. * mHOS);
   const double color_factor = 3;

   // on-shell calculation
   const double phase_spaceOS = 1. / (8. * Pi) * beta(mHOS, mdqOS, mdqOS);
   const double deltaqqOS = 4. / 3. * alpha_s_red * calc_DeltaH(betaOS) +
                            calc_deltaqq(alpha_s_red, Nf);

   // chiral breaking correctios
   // TODO: probably shouldn't be applied in case of CP-breaking
   double deltaH2 = 0.;
   if (!info::CPViolationInHiggsSector) {
      const double lt = std::log(std::pow(mHOS / mtpole, 2));
      const double lq = std::log(std::pow(mdqDR / mHOS, 2));
      deltaH2 =
         std::pow(alpha_s_red, 2) * (1.57 - 2.0 / 3.0 * lt + 1.0 / 9.0 * std::pow(lq, 2));
   }
   // this is the same (as far as analytic expression is concerned)
   // that would come from amplitude_squared but written in terms
   // of OS parameters
   const auto amp2OS = std::pow(mHOS, 2) * std::pow(betaOS, 2) * 2. *
                       std::norm(HBBbarVertexDR.left()) *
                       std::pow(mdqOS / mdqDR, 2);
   const auto resOS =
      flux * color_factor * phase_spaceOS * amp2OS * (1. + deltaqqOS);

   // MSbar/DRbar calculation
   const double phase_spaceDR = 1. / (8. * Pi) * beta(mHOS, mdqDR, mdqDR);
   const double deltaqqDR = 2. * (1. - 10. * xDR) / (1 - 4. * xDR) *
                               (4. / 3. - std::log(xDR)) * alpha_s_red +
                            4. / 3. * alpha_s_red * calc_DeltaH(betaDR) +
                            calc_deltaqq(alpha_s_red, Nf);
   // this is the same (as far as analytic expression is concerned)
   // that would come from amplitude_squared but written in terms
   // of DRbar/MSbar parameters
   const auto amp2DR = std::pow(mHOS, 2) * std::pow(betaDR, 2) * 2. *
                       std::norm(HBBbarVertexDR.left());
   const auto resDR =
      flux * color_factor * phase_spaceDR * amp2DR * (1. + deltaqqDR + deltaH2);

   // interpolate between on-shell formula in x ~ 1 regime
   // and MSbar/DRbar formula for x<<1
   return
      // low x limit
      (1 - 4. * xOS) * resDR +
      // high x limit
      4 * xOS * resOS;
}
