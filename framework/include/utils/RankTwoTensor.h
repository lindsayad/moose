//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef RANKTWOTENSOR_H
#define RANKTWOTENSOR_H

#include "Moose.h"
#include "RankFourTensor.h"

// Any requisite includes here
#include "libmesh/libmesh.h"
#include "libmesh/vector_value.h"
#include "libmesh/tensor_value.h"

#include <petscsys.h>
#include <vector>

// Forward declarations
class RankTwoTensor;
template <typename T>
class MooseArray;
typedef MooseArray<Real> VariableValue;
class ColumnMajorMatrix;

template <typename T>
void mooseSetToZero(T & v);

/**
 * Helper function template specialization to set an object to zero.
 * Needed by DerivativeMaterialInterface
 */
template <>
void mooseSetToZero<RankTwoTensor>(RankTwoTensor & v);

/**
 * RankTwoTensor is designed to handle the Stress or Strain Tensor for a fully anisotropic material.
 * It is designed to allow for maximum clarity of the mathematics and ease of use.
 * Original class authors: A. M. Jokisaari, O. Heinonen, M. R. Tonks
 *
 * RankTwoTensor holds the 9 separate Sigma_ij or Epsilon_ij entries.
 * The entries are accessed by index, with i, j equal to 1, 2, or 3, or
 * internally i, j = 0, 1, 2.
 */
class RankTwoTensor : public RealTensorValue
{
public:
  // Select initialization
  enum InitMethod
  {
    initNone,
    initIdentity
  };

  /// Default constructor; fills to zero
  RankTwoTensor();

  /// Select specific initialization pattern
  RankTwoTensor(const InitMethod);

  /**
   * To fill up the 9 entries in the 2nd-order tensor, fillFromInputVector
   * is called with one of the following fill_methods.
   * See the fill*FromInputVector functions for more details
   */
  enum FillMethod
  {
    autodetect = 0,
    isotropic1 = 1,
    diagonal3 = 3,
    symmetric6 = 6,
    general = 9
  };

  /**
   * Constructor that takes in 3 vectors and uses them to create rows
   * _coords[0][i] = row1(i), _coords[1][i] = row2(i), _coords[2][i] = row3(i)
   */
  RankTwoTensor(const TypeVector<Real> & row1,
                const TypeVector<Real> & row2,
                const TypeVector<Real> & row3);

  /// named constructor for initializing from row vectors
  static RankTwoTensor initializeFromRows(const TypeVector<Real> & row0,
                                          const TypeVector<Real> & row1,
                                          const TypeVector<Real> & row2);

  /// named constructor for initializing from column vectors
  static RankTwoTensor initializeFromColumns(const TypeVector<Real> & col0,
                                             const TypeVector<Real> & col1,
                                             const TypeVector<Real> & col2);

  /// Constructor that proxies the fillFromInputVector method
  RankTwoTensor(const std::vector<Real> & input) { this->fillFromInputVector(input); };

  /// Initialization list replacement constructors, 6 arguments
  RankTwoTensor(Real S11, Real S22, Real S33, Real S23, Real S13, Real S12);

  /// Initialization list replacement constructors, 9 arguments
  RankTwoTensor(
      Real S11, Real S21, Real S31, Real S12, Real S22, Real S32, Real S13, Real S23, Real S33);

  /// Copy constructor from RealTensorValue
  RankTwoTensor(const RealTensorValue & a) : RealTensorValue(a) {}

  /// Copy constructor from TypeTensor<Real>
  RankTwoTensor(const TypeTensor<Real> & a) : RealTensorValue(a) {}

  // Named constructors
  static RankTwoTensor Identity() { return RankTwoTensor(initIdentity); }

  /// Static method for use in validParams for getting the "fill_method"
  static MooseEnum fillMethodEnum();

  /**
   * fillFromInputVector takes 6 or 9 inputs to fill in the Rank-2 tensor.
   * If 6 inputs, then symmetry is assumed S_ij = S_ji, and
   *   _coords[0][0] = input[0]
   *   _coords[1][1] = input[1]
   *   _coords[2][2] = input[2]
   *   _coords[1][2] = input[3]
   *   _coords[0][2] = input[4]
   *   _coords[0][1] = input[5]
   * If 9 inputs then input order is [0][0], [1][0], [2][0], [0][1], [1][1], ..., [2][2]
   */
  void fillFromInputVector(const std::vector<Real> & input, FillMethod fill_method = autodetect);

  /**
   * fillFromScalarVariable takes FIRST/THIRD/SIXTH order scalar variable to fill in the Rank-2
   * tensor.
   */
  void fillFromScalarVariable(const VariableValue & scalar_variable);

public:
  /// returns _coords[i][c], ie, column c, with c = 0, 1, 2
  TypeVector<Real> column(const unsigned int c) const;

  /**
   * Returns a rotated version of the tensor data given a rank two tensor rotation tensor
   * _coords[i][j] = R_ij * R_jl * _coords[k][l]
   * @param R rotation matrix as a RealTensorValue
   */
  RankTwoTensor rotated(const RankTwoTensor & R) const;

  /**
   * rotates the tensor data given a rank two tensor rotation tensor
   * _coords[i][j] = R_ij * R_jl * _coords[k][l]
   * @param R rotation matrix as a RankTwoTensor
   */
  void rotate(const RankTwoTensor & R);

  /**
   * rotates the tensor data anticlockwise around the z-axis
   * @param a angle in radians
   */
  RankTwoTensor rotateXyPlane(Real a);

  /**
   * Returns a matrix that is the transpose of the matrix this
   * was called on.
   */
  RankTwoTensor transpose() const;

  /// sets _coords to a, and returns _coords
  RankTwoTensor & operator=(const RankTwoTensor & a);

  /// adds a to _coords
  RankTwoTensor & operator+=(const RankTwoTensor & a);

  /// returns _coords + a
  RankTwoTensor operator+(const RankTwoTensor & a) const;

  /// sets _coords -= a and returns vals
  RankTwoTensor & operator-=(const RankTwoTensor & a);

  /// returns _coords - a
  RankTwoTensor operator-(const RankTwoTensor & a) const;

  /// returns -_coords
  RankTwoTensor operator-() const;

  /// performs _coords *= a
  RankTwoTensor & operator*=(const Real a);

  /// returns _coords*a
  RankTwoTensor operator*(const Real a) const;

  /// performs _coords /= a
  RankTwoTensor & operator/=(const Real a);

  /// returns _coords/a
  RankTwoTensor operator/(const Real a) const;

  /// Defines multiplication with a vector to get a vector
  TypeVector<Real> operator*(const TypeVector<Real> & a) const;

  /// Defines multiplication with a TypeTensor<Real>
  RankTwoTensor operator*(const TypeTensor<Real> & a) const;

  /// Defines multiplication with a TypeTensor<Real>
  RankTwoTensor & operator*=(const TypeTensor<Real> & a);

  /// Defines logical equality with another RankTwoTensor
  bool operator==(const RankTwoTensor & a) const;

  /// Sets _coords to the values in a ColumnMajorMatrix (must be 3x3)
  RankTwoTensor & operator=(const ColumnMajorMatrix & a);

  /// returns _coords_ij * a_ij (sum on i, j)
  Real doubleContraction(const RankTwoTensor & a) const;

  /// returns C_ijkl = a_ij * b_kl
  RankFourTensor outerProduct(const RankTwoTensor & a) const;

  /// returns C_ijkl = a_ik * b_jl
  RankFourTensor mixedProductIkJl(const RankTwoTensor & a) const;

  /// returns C_ijkl = a_jk * b_il
  RankFourTensor mixedProductJkIl(const RankTwoTensor & a) const;

  /// returns C_ijkl = a_il * b_jk
  RankFourTensor mixedProductIlJk(const RankTwoTensor & a) const;

  /// return positive projection tensor of eigen-decomposition
  RankFourTensor positiveProjectionEigenDecomposition(std::vector<Real> & eigval,
                                                      RankTwoTensor & eigvec) const;

  /// returns A_ij - de_ij*tr(A)/3, where A are the _coords
  RankTwoTensor deviatoric() const;

  /// returns the trace of the tensor, ie _coords[i][i] (sum i = 0, 1, 2)
  Real trace() const;

  /// retuns the inverse of the tensor
  RankTwoTensor inverse() const;

  /**
   * Denote the _coords[i][j] by A_ij, then this returns
   * d(trace)/dA_ij
   */
  RankTwoTensor dtrace() const;

  /**
   * Denote the _coords[i][j] by A_ij, then
   * S_ij = A_ij - de_ij*tr(A)/3
   * Then this returns (S_ij + S_ji)*(S_ij + S_ji)/8
   * Note the explicit symmeterisation
   */
  Real generalSecondInvariant() const;

  /**
   * Calculates the second invariant (I2) of a tensor
   */
  Real secondInvariant() const;

  /**
   * Denote the _coords[i][j] by A_ij, then this returns
   * d(secondInvariant)/dA_ij
   */
  RankTwoTensor dsecondInvariant() const;

  /**
   * Denote the _coords[i][j] by A_ij, then this returns
   * d^2(secondInvariant)/dA_ij/dA_kl
   */
  RankFourTensor d2secondInvariant() const;

  /**
   * Sin(3*Lode_angle)
   * If secondInvariant() <= r0 then return r0_value
   * This is to gaurd against precision-loss errors.
   * Note that sin(3*Lode_angle) is not defined for secondInvariant() = 0
   */
  Real sin3Lode(const Real r0, const Real r0_value) const;

  /**
   * d(sin3Lode)/dA_ij
   * If secondInvariant() <= r0 then return zero
   * This is to gaurd against precision-loss errors.
   * Note that sin(3*Lode_angle) is not defined for secondInvariant() = 0
   */
  RankTwoTensor dsin3Lode(const Real r0) const;

  /**
   * d^2(sin3Lode)/dA_ij/dA_kl
   * If secondInvariant() <= r0 then return zero
   * This is to gaurd against precision-loss errors.
   * Note that sin(3*Lode_angle) is not defined for secondInvariant() = 0
   */
  RankFourTensor d2sin3Lode(const Real r0) const;

  /**
   * Denote the _coords[i][j] by A_ij, then
   * S_ij = A_ij - de_ij*tr(A)/3
   * Then this returns det(S + S.transpose())/2
   * Note the explicit symmeterisation
   */
  Real thirdInvariant() const;

  /**
   * Denote the _coords[i][j] by A_ij, then
   * this returns d(thirdInvariant()/dA_ij
   */
  RankTwoTensor dthirdInvariant() const;

  /**
   * Denote the _coords[i][j] by A_ij, then this returns
   * d^2(thirdInvariant)/dA_ij/dA_kl
   */
  RankFourTensor d2thirdInvariant() const;

  /**
   * Denote the _coords[i][j] by A_ij, then this returns
   * d(det)/dA_ij
   */
  RankTwoTensor ddet() const;

  /// Print the rank two tensor
  void print(std::ostream & stm = Moose::out) const;

  /// Add identity times a to _coords
  void addIa(const Real a);

  /// Sqrt(_coords[i][j]*_coords[i][j])
  Real L2norm() const;

  /**
   * sets _coords[0][0], _coords[0][1], _coords[1][0], _coords[1][1] to input,
   * and the remainder to zero
   */
  void surfaceFillFromInputVector(const std::vector<Real> & input);

  /**
   * computes eigenvalues, assuming tens is symmetric, and places them
   * in ascending order in eigvals
   */
  void symmetricEigenvalues(std::vector<Real> & eigvals) const;

  /**
   * computes eigenvalues and eigenvectors, assuming tens is symmetric, and places them
   * in ascending order in eigvals.  eigvecs is a matrix with the first column
   * being the first eigenvector, the second column being the second, etc.
   */
  void symmetricEigenvaluesEigenvectors(std::vector<Real> & eigvals, RankTwoTensor & eigvecs) const;

  /**
   * computes eigenvalues, and their symmetric derivatives wrt vals,
   * assuming tens is symmetric
   * @param eigvals are the eigenvalues of the matrix, in ascending order
   * @param deigvals Here digvals[i](j,k) = (1/2)*(d(eigvals[i])/dA_jk + d(eigvals[i]/dA_kj))
   * Note the explicit symmeterisation here.
   * For equal eigenvalues, these derivatives are not gauranteed to
   * be the ones you expect, since the derivatives in this case are
   * often defined by continuation from the un-equal case, and that is
   * too sophisticated for this routine.
   */
  void dsymmetricEigenvalues(std::vector<Real> & eigvals,
                             std::vector<RankTwoTensor> & deigvals) const;

  /**
   * Computes second derivatives of Eigenvalues of a rank two tensor
   * @param deriv store second derivative of the current tensor in here
   */
  void d2symmetricEigenvalues(std::vector<RankFourTensor> & deriv) const;

  /**
   * Uses the petscblaslapack.h LAPACKsyev_ routine to find, for symmetric _coords:
   *  (1) the eigenvalues (if calculation_type == "N")
   *  (2) the eigenvalues and eigenvectors (if calculation_type == "V")
   * @param calculation_type If "N" then calculation eigenvalues only
   * @param eigvals Eigenvalues are placed in this array, in ascending order
   * @param a Eigenvectors are placed in this array if calculation_type == "V".
   * See code in dsymmetricEigenvalues for extracting eigenvectors from the a output.
   */
  void syev(const char * calculation_type,
            std::vector<PetscScalar> & eigvals,
            std::vector<PetscScalar> & a) const;

  /**
   * Uses the petscblaslapack.h LAPACKsyev_ routine to perform RU decomposition and obtain the
   * rotation tensor.
   */
  void getRUDecompositionRotation(RankTwoTensor & rot) const;

  /**
   * This function initializes random seed based on a user-defined number.
   */
  static void initRandom(unsigned int);

  /**
   * This function generates a random unsymmetric rank two tensor.
   * The first real scales the random number.
   * The second real offsets the uniform random number
   */
  static RankTwoTensor genRandomTensor(Real, Real);

  /**
   * This function generates a random symmetric rank two tensor.
   * The first real scales the random number.
   * The second real offsets the uniform random number
   */
  static RankTwoTensor genRandomSymmTensor(Real, Real);

  /// RankTwoTensor from outer product of vectors
  void vectorOuterProduct(const TypeVector<Real> &, const TypeVector<Real> &);

  /// Return real tensor of a rank two tensor
  void fillRealTensor(RealTensorValue &);

  ///Assigns value to the columns of a specified row
  void fillRow(unsigned int, const TypeVector<Real> &);

  ///Assigns value to the rows of a specified column
  void fillColumn(unsigned int, const TypeVector<Real> &);

  /// returns this_ij * b_ijkl
  RankTwoTensor initialContraction(const RankFourTensor & b) const;

private:
  static constexpr unsigned int N = LIBMESH_DIM;
  static constexpr unsigned int N2 = N * N;

  friend class RankFourTensor;
  friend class RankThreeTensor;
};

#endif // RANKTWOTENSOR_H
