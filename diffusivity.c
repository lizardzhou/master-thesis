/***********************************************************************
   UDF for defining the diffusivity of the particle mixture
 ***********************************************************************/
 #include "udf.h"

DEFINE_DIFFUSIVITY(diffusivity,c,t,i)
{
	real temperature = C_T(c, t);
	real diff = 2.252 / pow(10, 5) * pow((temperature / 273), 1.81);
	return diff;
}