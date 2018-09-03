
// special case for H -> W+ W-
// TODO: implement higher order corrections
template <>
double CLASSNAME::get_partial_width<H, W, conj<W>::type>(
   const ContextName& context,
   typename field_indices<H>::type const& indexIn,
   typename field_indices<conj<W>::type>::type const& indexOut1,
   typename field_indices<W>::type const& indexOut2) const
{

   const double mH = context.mass<H>(indexIn);
   const double mW = context.mass<W>(indexOut1);
   const double x = Sqr(mW / mH);
   double res;
   // three-body-decays if below threshold
   if (4 * x > 1.0) {
      // const auto vd = MODELPARAMETER(vd);
      // const auto vu = MODELPARAMETER(vu);
      // TODO: specify the vev correctly
      const auto vev = 246.0; // sqrt(Sqr(vd) + Sqr(vu));
      const double sw2 =
         Sqr(Sin(model.ThetaW())); // 1.0 - Sqr(PHYSICAL(MVWp)/PHYSICAL(MVZ));

      // deltaV is 1 for W-bosons
      const double RT = 3 * (1 - 8 * x + 20 * Sqr(x)) / sqrt(4 * x - 1) *
                           acos(0.5 * (3 * x - 1) / pow(x, 3.0 / 2.0)) -
                        0.5 * (1 - x) / x * (2 - 13 * x + 47 * Sqr(x)) -
                        3.0 / 2.0 * (1 - 6 * x + 4 * Sqr(x)) * Log(x);

      res = 3.0 / (128 * pow(Pi, 3)) * mH / Sqr(vev) * RT;
   } else {
      res = 2.0 / (128 * Pi * mH * Sqr(x)) * sqrt(1 - 4 * x) *
            (1 - 4 * x + 12 * Sqr(x));
   }
   const auto indices = concatenate(indexOut2, indexOut1, indexIn);
   const auto ghWW =
      Vertex<conj<W>::type, W, H>::evaluate(indices, context).value();
   return res * std::norm(ghWW);
}