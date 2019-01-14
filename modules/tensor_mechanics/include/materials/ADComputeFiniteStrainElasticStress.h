//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef ADCOMPUTEFINITESTRAINELASTICSTRESS_H
#define ADCOMPUTEFINITESTRAINELASTICSTRESS_H

#include "ADComputeStressBase.h"
#include "GuaranteeConsumer.h"

template <ComputeStage>
class ADComputeFiniteStrainElasticStress;

declareADValidParams(ADComputeFiniteStrainElasticStress);

/**
 * ADComputeFiniteStrainElasticStress computes the stress following elasticity
 * theory for finite strains
 */
template <ComputeStage compute_stage>
class ADComputeFiniteStrainElasticStress : public ADComputeStressBase<compute_stage>,
                                           public GuaranteeConsumer
{
public:
  ADComputeFiniteStrainElasticStress(const InputParameters & parameters);

  void initialSetup() override;

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpStress() override;

  const ADMaterialProperty(RankTwoTensor) & _strain_increment;
  const ADMaterialProperty(RankTwoTensor) & _rotation_increment;
  const MaterialProperty<RankTwoTensor> & _stress_old;

  /**
   * The old elastic strain is used to calculate the old stress in the case
   * of variable elasticity tensors
   */
  const MaterialProperty<RankTwoTensor> & _elastic_strain_old;

  usingComputeStressBaseMembers;
};

#endif // ADCOMPUTEFINITESTRAINELASTICSTRESS_H
