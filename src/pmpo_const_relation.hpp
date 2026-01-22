#ifndef POLYMPO_CONSTITUTIVE_RELATION
#define POLYMPO_CONSTITUTIVE_RELATION

#include <stdlib.h>
#include <iostream>

namespace polyMPO{

#define PUNY  1.0e-11
static constexpr double eccentricity = 2.0;
static constexpr double eccentricitySquared = eccentricity * eccentricity;
static constexpr double dampingTimescaleParameter = 0.36;

KOKKOS_INLINE_FUNCTION
void constitutive_evp(const Vec3d& strain, Vec3d& stress, const double& icePressure, double& replacementPressure, 
                      const double& areaMP, const double& dtElastic, const double& dampingTimescale){
  auto strainDivergence = strain[0] + strain[1];
  auto strainTension = strain[0] - strain[1];
  auto strainShearing = 2*strain[2];
  
  auto stress1 = stress[0] + stress[1];
  auto stress2 = stress[0] - stress[1];
  
  auto Delta = sqrt(strainDivergence*strainDivergence + (strainTension*strainTension + strainShearing*strainShearing)/eccentricitySquared);

  auto pressureCoefficient = icePressure / Kokkos::max(Delta, PUNY);
  replacementPressure = pressureCoefficient * Delta;
  pressureCoefficient = (pressureCoefficient * dtElastic) / (2.0 * dampingTimescale);
 
  auto denominator = 1.0 + (0.5 * dtElastic) / dampingTimescale;

  stress1  = (stress1    +  pressureCoefficient                        * (strainDivergence - Delta))  / denominator;
  stress2  = (stress2    + (pressureCoefficient / eccentricitySquared) *  strainTension             ) / denominator;
  stress[2] = (stress[2] + (pressureCoefficient / eccentricitySquared) *  strainShearing * 0.5) / denominator;
 
  stress[0] = 0.5 * (stress1 + stress2);
  stress[1] = 0.5 * (stress1 - stress2);
}

}
#endif


