
template<>
double CLASSNAME::get_partial_width<H,AH,AH>(
   const ContextName& context,
   typename field_indices<AH>::type const& indexIn,
   typename field_indices<AH>::type const& indexOut1,
   typename field_indices<AH>::type const& indexOut2
   ) const
{
   return 0.;
}
