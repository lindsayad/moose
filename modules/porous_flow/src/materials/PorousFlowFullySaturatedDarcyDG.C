//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "PorousFlowFullySaturatedDarcyDG.h"
#include <iostream>
using namespace std;
registerMooseObject("PorousFlowApp", PorousFlowFullySaturatedDarcyDG);

InputParameters
PorousFlowFullySaturatedDarcyDG::validParams()
{
  InputParameters params = InterfaceKernel::validParams();
  params.addParam<Real>("stabilized_para", 10.0, "DG stabilized parameter");
  params.addParam<bool>("multiply_by_density",
                        true,
                        "If true, then this Kernel is the fluid mass "
                        "flux.  If false, then this Kernel is the "
                        "fluid volume flux (which is common in "
                        "poro-mechanics)");
  params.addClassDescription("DG (IIPG) interface kernel for single phase Darcy flow");
  return params;
}

PorousFlowFullySaturatedDarcyDG::PorousFlowFullySaturatedDarcyDG(
    const InputParameters & parameters)
   : JvarMapKernelInterface<InterfaceKernel>(parameters),
     _stabilized_para(getParam<Real>("stabilized_para")),
     _multiply_by_density(getParam<bool>("multiply_by_density")),
     _permeability_L(getMaterialProperty<RealTensorValue>("PorousFlow_permeability_qp")),
     _permeability_R(getNeighborMaterialProperty<RealTensorValue>("PorousFlow_permeability_qp")),
     _density_L(getMaterialProperty<std::vector<Real>>("PorousFlow_fluid_phase_density_qp")),
     _density_R(getNeighborMaterialProperty<std::vector<Real>>("PorousFlow_fluid_phase_density_qp")),
     _viscosity_L(getMaterialProperty<std::vector<Real>>("PorousFlow_viscosity_qp")),
     _viscosity_R(getNeighborMaterialProperty<std::vector<Real>>("PorousFlow_viscosity_qp")),
     _pp_L(getMaterialProperty<std::vector<Real>>("PorousFlow_porepressure_qp")),
     _pp_R(getNeighborMaterialProperty<std::vector<Real>>("PorousFlow_porepressure_qp")),
     _grad_p_L(getMaterialProperty<std::vector<RealGradient>>("PorousFlow_grad_porepressure_qp")),
     _grad_p_R(getNeighborMaterialProperty<std::vector<RealGradient>>("PorousFlow_grad_porepressure_qp"))
{
}

Real
PorousFlowFullySaturatedDarcyDG::computeQpResidual(Moose::DGResidualType type)
{
 //c
 //c this computes the residual contributed from DG interface only
 //c based on the IIPG method;
 //c
  Real r(0.0);
  Real a_flux = computeAverageFlux();
  Real pp_jump = computePorePressureJump();
  Real flow_conduct = computeAverageConduct();
  //c
  //c Adding flux on surface based on IIPG formulation
  //c
  r -= a_flux;
  //c
  //c Adding stabilized term contribution
  //c
  r+= _stabilized_para*flow_conduct/sqrt(_current_side_volume)*pp_jump;
  switch (type)
  {
   //c
   //c [test_right-test_left]*flux;
   //c test_left==test==test_primary; test_right=test_neighbor==test_secondary.
   //c
   //c for left element
   //c
   case Moose::Element:
        r *= -_test[_i][_qp];
        break;
   //c
   //c for right element
   //c
   case Moose::Neighbor:
        r *= _test_neighbor[_i][_qp];
        break;
  }
      cerr << "_qp = " << _qp << ".\n";
      cerr << "residual = " << r << ".\n";
return r;
}

Real
PorousFlowFullySaturatedDarcyDG::computeQpJacobian(Moose::DGJacobianType type)
{
 Real jac(0.0);
   //c
   //c IIPG only;
   //c
   //c  This rouine adds Jacobian contributed from:
   //c   (1) major face integeral:
   //c    -\int 1/2{flux_left+flux_righ}_normals\cdot [\dleta p] ds;
   //c   (2) stabilized term:
   //c     \lambda\int [p] \cdot[\delta p] ds
   //c
   Real flow_conduct =  computeAverageConduct();
   Real DG_para=_stabilized_para*flow_conduct/sqrt(_current_side_volume);
   RealVectorValue delta_p_n;
   RealVectorValue grad_phi;
   switch (type)
   {
    //c
    //c  for K^(LL)
    //c
    case Moose::ElementElement:
        //c
        //c compute the vetcor delta_p_n by the test function \delta p and normal vector:
        //c
        delta_p_n = _test[_i][_qp]*_normals[_qp];
        grad_phi = mobility_L()*_permeability_L[_qp]*_grad_phi[_j][_qp];
        //c
        //c add Jacobian contributed from major face integral
        //c
        jac -= delta_p_n*grad_phi*0.5;
        //c
        //c add Jacobian contributed from DG stabilized term
        //c
        jac += DG_para*_test[_i][_qp]*_phi[_j][_qp];
     break;
     //c
     //c  for K^(LR)
     //c
    case Moose::ElementNeighbor:
         delta_p_n = _test[_i][_qp]*_normals[_qp];
         grad_phi = mobility_R()*_permeability_R[_qp]*_grad_phi_neighbor[_j][_qp];
         //c
         //c add Jacobian contributed from major face integral
         //c
         jac -= delta_p_n*grad_phi*0.5;
         jac -= DG_para*_test[_i][_qp]*_phi_neighbor[_j][_qp];
         break;
    //c
    //c  for K^(RL)
    //c
    case Moose::NeighborElement:
        delta_p_n = _test_neighbor[_i][_qp]*_normals[_qp];
        grad_phi = mobility_L()*_permeability_L[_qp]*_grad_phi[_j][_qp];
        //c
        //c add Jacobian contributed from major face integral
        //c
        jac += delta_p_n*grad_phi*0.5;
        jac -= DG_para*_test_neighbor[_i][_qp]*_phi[_j][_qp];
        break;
    //c
    //c  for K^(RR)
    //c
    case Moose::NeighborNeighbor:
        delta_p_n= _test_neighbor[_i][_qp]*_normals[_qp];
        grad_phi = mobility_R()*_permeability_R[_qp]*_grad_phi_neighbor[_j][_qp];
        //c
        //c add Jacobian contributed from major face integral
        //c
        jac += delta_p_n*grad_phi*0.5;
        jac += DG_para*_test_neighbor[_i][_qp]*_phi_neighbor[_j][_qp];
        break;
   }
   return jac;
}


Real
PorousFlowFullySaturatedDarcyDG::computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar)
{
 return 0.0;
}

Real
PorousFlowFullySaturatedDarcyDG::computeAverageFlux()
{
 //c
 //c compute the average of flux across an interface
 //c
 const unsigned ph = 0;
 Real flux;
 flux =  -(mobility_L()*_permeability_L[_qp]*_grad_p_L[_qp][ph]+
           mobility_R()*_permeability_R[_qp]*_grad_p_R[_qp][ph])
           *_normals[_qp]*0.5;
 return flux;
}

Real
PorousFlowFullySaturatedDarcyDG::computePorePressureJump()
{
 //c
 //c compute the jump of pore pressure across an interface
 //c
 const unsigned ph = 0;
 Real pp_jump;
 pp_jump = _pp_R[_qp][ph]-_pp_L[_qp][ph];
 return pp_jump;
}

Real
PorousFlowFullySaturatedDarcyDG::computeAverageConduct()
{
 //c
 //c compute the average of conductivityacross an interface
 //c
 Real conduct;
 conduct = 0.5*(mobility_L()*_permeability_L[_qp].tr()+
         mobility_R()*_permeability_R[_qp].tr())/3.0;
 return conduct;
}

Real
PorousFlowFullySaturatedDarcyDG::mobility_L() const
{
  const unsigned ph = 0;
  Real mob = 1.0 / _viscosity_L[_qp][ph];
  if (_multiply_by_density)
    mob *= _density_L[_qp][ph];
  return mob;
}

Real
PorousFlowFullySaturatedDarcyDG::mobility_R() const
{
  const unsigned ph = 0;
  Real mob = 1.0 / _viscosity_R[_qp][ph];
  if (_multiply_by_density)
    mob *= _density_R[_qp][ph];
  return mob;
}

