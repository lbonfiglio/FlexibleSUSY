// special case for Higgs -> Z Z
// TODO: implement higher order corrections

template <>
double CLASSNAME::get_partial_width<H,Z,Z>(
   const ContextName& context,
   typename field_indices<H>::type const& indexIn,
   typename field_indices<Z>::type const& indexOut1,
   typename field_indices<Z>::type const& indexOut2
   ) const
{

   const double mHOS = context.physical_mass<H>(indexIn);
   const double mZ = context.mass<Z>(indexOut1);
   const double mZOS = context.physical_mass<Z>(indexOut1);
   const double x = Sqr(mZOS/mHOS);
   double res;
   // three-body-decays if below threshold
      const auto vev = 246.0; //sqrt(Sqr(vd) + Sqr(vu));
   if(4.*x > 1.0) {
      // const auto vd = MODELPARAMETER(vd);
      // const auto vu = MODELPARAMETER(vu);
      // TODO: specify the vev correctly
      const double sw2 = Sqr(Sin(model.ThetaW()));//1.0 - Sqr(PHYSICAL(MVWp)/PHYSICAL(MVZ));
      const double deltaV = 7.0/12.0 - 10.0/9.0 * sw2 + 40.0/27.0 * Sqr(sw2);

      res = 3.0/(128*pow(Pi,3)) * mHOS/Sqr(vev) * deltaV * RT(x);
   } else {

      const double flux = 1. / (2 * mHOS);
      // phase space without symmetry factor
      const double ps = 1. / (8. * Pi) * beta(mHOS, mZOS, mZOS);

      // phase space symmetry factor
      const double ps_symmetry = 1. / 2.;

      // matrix element squared
      const auto mat_elem = calculate_amplitude<H, Z, Z>(
         context, indexIn, indexOut1, indexOut2);
      const auto mat_elem_sq = mat_elem.square();

      // flux * phase space factor * matrix element squared
      return flux * ps * ps_symmetry * mat_elem_sq;
   }
   const auto indices = concatenate(indexOut1, indexOut2, indexIn);
   const auto ghZZ =
      Vertex<Z, Z, H>::evaluate(indices, context).value() * std::pow(mZOS/mZ, 2);
   return res * std::norm(ghZZ);
}
