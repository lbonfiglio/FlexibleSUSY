// special case for H -> W+ W-
// TODO: implement higher order corrections
template <>
double CLASSNAME::get_partial_width<H, W, conj<W>::type>(
   const ContextName& context, typename field_indices<H>::type const& indexIn,
   typename field_indices<conj<W>::type>::type const& indexOut1,
   typename field_indices<W>::type const& indexOut2) const
{

   const double mH = context.mass<H>(indexIn);
   const double mHOS = context.physical_mass<H>(indexIn);
   const double mW = context.mass<W>(indexOut1);
   const double mWOS = context.physical_mass<W>(indexOut1);
   const double x = Sqr(mW / mH);
   double res;
   // three-body-decays if below threshold
   if (4 * x > 1.0) {
      // const auto vd = MODELPARAMETER(vd);
      // const auto vu = MODELPARAMETER(vu);
      // TODO: specify the vev correctly
      const auto vev = 246.0; // sqrt(Sqr(vd) + Sqr(vu));

      // deltaV is 1 for W-bosons
      const double RT = 3 * (1 - 8 * x + 20 * Sqr(x)) / sqrt(4 * x - 1) *
                           acos(0.5 * (3 * x - 1) / pow(x, 3.0 / 2.0)) -
                        0.5 * (1 - x) / x * (2 - 13 * x + 47 * Sqr(x)) -
                        3.0 / 2.0 * (1 - 6 * x + 4 * Sqr(x)) * Log(x);

      res = 3.0 / (128 * pow(Pi, 3)) * mH / Sqr(vev) * RT;
   } else {

      const double flux = 1. / (2 * mHOS);
      // phase space without symmetry factor
      const double ps = 1. / (8. * Pi) * beta(mHOS, mWOS, mWOS);

      // matrix element squared
      const auto mat_elem = effective_coupling<H, W, conj<W>::type>(
         context, indexIn, indexOut1, indexOut2);
      const auto mat_elem_sq = mat_elem.square();

      // flux * phase space factor * matrix element squared
      return flux * ps * mat_elem_sq;
   }
   const auto indices = concatenate(indexOut2, indexOut1, indexIn);
   const auto ghWW =
      Vertex<conj<W>::type, W, H>::evaluate(indices, context).value();
   return res * std::norm(ghWW);
}
