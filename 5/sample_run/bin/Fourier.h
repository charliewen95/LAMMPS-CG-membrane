#if !defined( FOURIER_INCLUDED )
#define FOURIER_INCLUDED

#include <math.h>

#include <complex>
#include <vector>

#if defined( _OPENMP )
	#include "omp.h"
#endif

//
// Todo:
//
// - OpenMP enabled synthesis routines (currently unused, so no rush!)
//

namespace Fourier
{

typedef std::complex<double> Complex;
typedef std::vector<Complex> Coeffs1D;
typedef std::vector< std::vector<Complex> > Coeffs2D;


//
// 1D Fourier analysis equation:
//
// F(k) = 1/N \Sum_n f[n] . exp( -i.arg )
// F(k) = 1/N \Sum_n f[n] . ( cos(arg), i.sin(arg) )
//
// Where:
//
// - arg = 2pi( k.x[n]/Lx )
// - f[n], x[n] : N samples of function f at positions x on range [0,L): f[n] = f(x[n])
// - k is an integer wave number
//
// Notes:
//
// - The prefactor, (1/N) is arbitrary. Fourier descriptions from different sources will have different prefactors.
// - We rearrange f[n].exp(-i.arg) to exp(-i.arg).f[n] to enforce correct multiplication of complex number by scalar.
// - The loop over the wave numebrs k is INCLUSIVE of max_k: 0 <= k <= max_k
//
int Analysis(
	const std::vector<double>& x,
	const std::vector<double>& f,
	double Lx,
	int max_k,
	Coeffs1D& F )
{
	double tpol = 2.0*M_PI/Lx;
	size_t N = x.size();

	if( (N<1) || (N!=f.size()) ) return -1;

	double prefactor = 1.0/N;

	// Allocate everything in advance, in case using OpenMP.
	F.resize( max_k+1 );

	// Simplest OpenMP approach likely okay for our purposes.
	#pragma omp parallel for
	for( int k=0; k<=max_k; k++ )
	{
		Complex acc = 0.0;
		for( size_t i=0; i<N; i++ )
		{
			double arg = tpol * k * x[i];
			acc += Complex( cos(arg), sin(arg) ) * f[i];
		}
		acc *= prefactor;
		F[k] = acc;
	}

	return 1;
}

//
// 1D Fourier synthesis equation:
//
// f(x) = \Sum_k F[k] . exp( -i.arg )
// f(x) = \Sum_k F[k] . ( cos(arg), i.sin(arg) )
//
// Where:
//
// - arg = 2pi.( k.(x[n]/Lx) )
// - f[n], x[n] : N discrete values of function f[] sampled at positions x[], with x on range [0,L)
// - k is an integer wave vector
//
// Notes:
//
// - F[n].exp(-i.arg) rearranged to exp(-i.arg).F[n], for consistency of expression with previous analysis equation
//   implementation. Does not change the result, obviously!
// - Here we determine max_k from the size of the coefficient vector: this vector includes the max_k coefficient from
//   the analysis equation above, so loop bounds of 0 <= k < max_k, NOT 0 <= k <= max_k as in the analysis routine.
//
Complex Synthesis(
	double x,
	double Lx,
	const Coeffs1D& F )
{
	Complex c( 0.0, 0.0 );
	double tpol = 2.0*M_PI/Lx;
	int max_k = (int)F.size(); // F includes max_k, so k < max_k in loop, NOT k <= max_k

	for( int k=0; k<max_k; k++ ) // declare k as int: negation of k as a size_t problematic
	{
		double arg = tpol * (-k) * x;
		c += Complex( cos(arg), sin(arg) ) * F[k];
	}
	return c;
}


//
// 2D Fourier analysis equation:
//
// F(k,l) = 1/N \Sum_{k,l} f[n] . exp( -i.arg )
// F(k,l) = 1/N \Sum_{k,l} f[n] . ( cos(arg), i*sin(arg) )
//
// Where:
//
// - arg = 2pi( k.(x[n]/Lx) + l.(y[n]/Ly) )
// - f[n], x[n], y[n] : N samples of function f at positions x,y, x,y on range [0,Lx|Ly): f[n] = f(x[n],y[n])
// - k, l are integer wave vectors
//
// Notes:
//
// - The prefactor, (1/N) is arbitrary. Fourier descriptions from different sources will have different prefactors.
// - We rearrange f[n].exp(-i.arg) to exp(-i.arg).f[n] to enforce multiplying a complex number by a scalar.
// - k, l coefficients are INCLUSIVE of max_k, max_l
//
int Analysis(
	const std::vector<double>& x,
	const std::vector<double>& y,
	const std::vector<double>& f,
	double Lx, double Ly,
	int max_k, int max_l,
	Coeffs2D& F )
{
	constexpr double twopi = 2.0*M_PI;
	size_t N = x.size();

	if( N<1 || N!=y.size() || N!=f.size() ) return -1;

	double prefactor = 1.0/N;

	// Allocate everything in advance, in case using OpenMP.
	F.resize( max_k+1 );
	for( int k=0; k<=max_k; k++ ) F[k].resize( max_l+1 );

	// Simplest OpenMP approach likely okay for our purposes.
	#pragma omp parallel for
	for( int k=0; k<=max_k; k++ )
	{
		for( int l=0; l<=max_l; l++ )
		{
			Complex acc = 0.0;
			for( size_t i=0; i<N; i++ )
			{
				double arg = twopi*( (x[i]/Lx)*k + (y[i]/Ly)*l );
				acc += Complex( cos(arg), sin(arg) ) * f[i];
			}
			acc *= prefactor;
			F[k][l] = acc;
		}
	}

	return 1;
}


//
// 2D Fourier synthesis equation:
//
// f(x,y) = \Sum_{k,l} F[k,l] . exp( -i.arg )
// f(x,y) = \Sum_{k,l} F[k,l] . ( cos(arg), i*sin(arg) )
//
// Where:
//
// - arg = 2pi( (-k).x[n]/Lx + (-l).y[n]/Ly )
// - f[n], x[n], y[n] : N samples of function f at positions x,y, x,y on range [0,Lx|Ly): f[n] = f(x[n],y[n])
// - k, l are integer wave vectors
//
// Notes:
//
// - F[n].exp(-i.arg) rearranged to exp(-i.arg).F[n], for consistency of expression with previous analysis equation
//   implementation. Does not change the result.
// - We again determine the upper limit of the loops using the size of the coefficient vectors, which are already
//   inclusive of max_k and max_l: therefore use e.g. 0 <= k < max_k and NOT 0 <= k <= max_k for loop bounds.
//
Complex Synthesis(
	double x, double y,
	double Lx, double Ly,
	const Coeffs2D& F )
{
	constexpr double twopi = 2.0*M_PI;
	Complex c( 0.0, 0.0 );

	if( F.size() < 1 || F[0].size() < 1 ) return c;

	//
	// Be careful; the size of these vectors includes the max_k specified in the
	// analysis equations, so use k,l < max, rather than k,l <= max!
	//
	int max_k = F.size();
	int max_l = F[0].size();

	for( int k=0; k<max_k; k++ )
	{
		for( int l=0; l<max_l; l++ )
		{
			double arg = twopi*( (x/Lx)*(-k) + (y/Ly)*(-l) );
			c += Complex( cos(arg), sin(arg) ) * F[k][l];
		}
	}

	return c;
}


}

#endif
