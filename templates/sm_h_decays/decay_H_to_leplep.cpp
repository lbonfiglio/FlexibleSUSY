// CP-even Higgs to charged leptons

template<>
double CLASSNAME::get_partial_width<H,bar<lep>::type,lep>(
      const cxx_qft::SM_evaluation_context& context,
      typename cxx_qft::field_indices<H>::type const& indexIn,
      typename cxx_qft::field_indices<bar<lep>::type>::type const& indexOut1,
      typename cxx_qft::field_indices<lep>::type const& indexOut2) const
   {
      using effective_couplings::beta;

      // onshell masses
      const double mIn = context.physical_mass<H>(indexIn);
      const double mOut1 = context.physical_mass<bar<lep>::type>(indexOut1);
      const double mOut1DR = context.mass<bar<lep>::type>(indexOut1);
      const double mOut2 = context.physical_mass<lep>(indexOut2);

      // phase space without symmetry factor
      const double ps = 1. / (8. * Pi) * beta(mIn, mOut1, mOut2);

      // phase space symmetry factor
      const double ps_symmetry = 1.;

      // matrix element squared
   const auto xOS = std::pow(mOut1/mIn, 2);
   const auto betaOS = sqrt(1. - 4.*xOS);
      const auto indices = concatenate(indexIn, indexOut1, indexOut2);
      const auto HLLbarVertexDR = Vertex<H, bar<lep>::type, lep>::evaluate(indices, context);
      const auto amp2OS = std::pow(mIn, 2) * std::pow(betaOS, 2) *
                2.*std::norm(HLLbarVertexDR.left()) * std::pow(mOut1 / mOut1DR, 2);

      // flux * phase space factor * symmetry factor * matrix element^2
      return 0.5 * ps * ps_symmetry * amp2OS/ mIn;
   }