// special case for H -> photon Z
template <>
double CLASSNAME::get_partial_width<H, A, Z>(
   const ContextName& context,
   typename field_indices<H>::type const& indexIn,
   typename field_indices<A>::type const& indexOut1,
   typename field_indices<Z>::type const& indexOut2) const
{

   const double mass = context.mass<H>(indexIn);
   const double MVZ = context.mass<Z>(indexOut2);
   if(MVZ/mass >= 1.0)
      return 0.0;
   return 0 * 1.0/(16.0*Pi) * pow(mass * (1.0 - Sqr(MVZ/mass)), 3);
      //   * effective_coupling<H, A,Z>(context, indexIn, indexOut1, indexOut2).square();
}

template<>
double CLASSNAME::get_partial_width<H,Z,A>(
   const ContextName& context,
   typename field_indices<H>::type const& indexIn,
   typename field_indices<Z>::type const& indexOut1,
   typename field_indices<A>::type const& indexOut2
   ) const
{
   return get_partial_width<H,A,Z>(context, indexIn, indexOut2, indexOut1);
}