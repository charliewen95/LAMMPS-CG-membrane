#if !defined( STATS_INCLUDED )
#define STATS_INCLUDED

#include <vector>
#include <map>

namespace Stats
{

/*
	Running statistics, based on algorithms of B. P. Welford
	(via Knuth, "The Art of Computer Programming").

	This algorithm provides not only the capability to determine
	variance (and hence, stdev and stderr) from a running input
	with very little storage, it is also robust to catastrophic
	cancellation.
*/

struct Stats
{
	size_t N; // number of samples so far
	double S; // N*sigma^2
	double min, mean, max;
	
	Stats()
	{
		Clear();
	}

	void Clear()
	{
		N = 0;
		S = min = mean = max = 0.0;
	}

	void AddSample( double x )
	{
		N++;
		if( N == 1 )
		{
			min = mean = max = x;
			S = 0.0;
			return;
		}
		
		//
		// Update values (new marked with prime):
		//   mean' = mean + (x-mean)/N
		//   S' = S + (x-mean)*(x-mean')
		//

		double delta = (x-mean);
		mean += delta/N;
		S += delta * (x-mean); // <- note: uses updated value of "mean" as well as old value via "delta".

		if( x < min ) min = x;
		if( x > max ) max = x;
	}

	double Variance() const
	{
		// SAMPLE variance! Note division by N-1 rather than N.
		return (N>1) ? (S/(N-1)) : (0.0);
	}
	double StdDev() const
	{
		// SAMPLE standard deviation, see note in Variance() method.
		return sqrt( Variance() );
	}
	double StdErr() const
	{
		//
		// Estimated standard error of the sample mean:
		//
		// SE = stdev / sqrt(N) : stdev = sample standard deviation
		// SE = sqrt(variance) / sqrt(N) : as stdev = sqrt(variance)
		// SE = sqrt( variance / N ) : as sqrt() is distributive
		//
		// This is a (likely unneccessary) optimization using one sqrt()
		// call rather than two in calculating standard error of the
		// sample mean.
		//
		return (N>1) ? ( sqrt(Variance()/N) ) : (0.0);
	}
	
	//
	// Allows combining separate sets of sample stats using e.g. stats1 += stats2
	// where stats1 and stats2 are both instances of Stats. This obviously assumes
	// both stats1 and stats2 are samples from the same population ...
	//
	Stats & operator += ( const Stats &rhs )
	{
		// Temporary values, in case &rhs == this
		size_t new_N;
		double new_sum, new_mean, new_S;
		
		// Ignore if no data present in rhs
		if( rhs.N < 1 ) return *this;
				
		new_N = N + rhs.N;
		new_sum = (mean*N) + (rhs.mean*rhs.N);
		new_mean = new_sum / new_N; // safe: rhs.N >= 1, so new_N >= 1.
		
		//
		// This is basically the "parallel algorithm" version of the "online" algorithm
		// for calculating variance in one pass when the sample is partitioned into multiple
		// sets. This is attributed to Chan et al, "Updating Formulae and a Pairwise Algorithm for
		// Computing Sample Variances.", Technical Report STAN-CS-79-773, Stanford Comp Sci (1979).
		//
		// http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Parallel_algorithm
		//
		double delta = (mean - rhs.mean);
		new_S = (S + rhs.S) + (delta*delta)*(N*rhs.N)/new_N;
		
		// If no samples in this Stats structure, or if rhs min/max should replace current:
		if( (N<1) || (rhs.min<min) ) min = rhs.min;
		if( (N<1) || (rhs.max>max) ) max = rhs.max;
		
		// Only change N now, as we needed it in the min/max update above.
		N = new_N;
		S = new_S;
		mean = new_mean;
		
		return *this;
	}
};

//
// "Sparse" statistics object for building histograms etc with specified bin width "delta" which acts to
// discretise the sampling domain: this is useful where we don't know the upper and/or lower bounds of the
// sampling domain in advance. This approach also avoids wasted memory where we require a potentially large
// sampling domain but comparatively narrow bins, only some of which may actually be populated.
//
// The "min" value was put in here for the Fourier sampling of bilayer head groups: if minimum wavelength
// is 2pi/L_max, then the first histogram bin can encompass values too small to be physical. That can be
// a problem when plotting etc: the "coordinate" of the data in the bin is at the bin centre, which can
// visually "stretch" the data out in slightly misleading ways, e.g. when plotting on a log scale.
//
class MapStats
{
	protected:

		void save_header( FILE* f )
		{
			if( f == nullptr ) return;
			fprintf( f, "# x, y = coordinate and mean value at coordinate\n" );
			fprintf( f, "#\n" );
			fprintf( f, "# N(y)      = number of samples\n" );
			fprintf( f, "# StdDev(y) = standard deviation: sqrt(variance)\n" );
			fprintf( f, "# StdErr(y) = standard error of the estimated mean: sqrt(variance)/N = StdDev(y)/N\n" );
			fprintf( f, "#\n" );
			fprintf( f, "# Both StdDev(y) and StdErr(y) use the SAMPLE variance\n" );
			fprintf( f, "#\n" );
			fprintf( f, "# %12.12s  %12.12s  %12.12s  %12.12s  %12.12s\n", "x", "y", "N(y)", "StdDev(y)", "StdErr(y)" );
			fprintf( f, "#\n" );
		}

	public:

		std::map<int,Stats> m;
		double delta, min;

		MapStats()                     { Reset( 1.0, 0.0 ); }
		MapStats( double d )           { Reset( d,   0.0 ); }
		MapStats( double d, double m ) { Reset( d,   m   ); }

		void Reset( double delta_, double min_ = 0.0 )
		{
			if( delta_ != 0.0 )  delta = delta_;
			min = min_;
			Clear();
		}
		void Clear() { m.clear(); }

		void AddSampleByBin( int bin, double val ) { m[bin].AddSample( val ); }
		void AddSample( double coord, double val )
		{
			int bin = (int)floor( (coord-min) / delta );
			AddSampleByBin( bin, val );
		}

		int Save( FILE*f, bool print_header = true )
		{
			if( f == nullptr ) return -1;
			if( print_header ) save_header( f );
			for( const auto& it : m )
			{
				int bin_no = it.first;
				const Stats& s = it.second;
				double coord = min + (0.5+bin_no)*delta;
				fprintf( f, "  %e  %e  %12d  %e  %e\n", coord, s.mean, (int)s.N, s.StdDev(), s.StdErr() );
			}
			return 1;
		}
};

//
// Distribution class: designed to build data in a number of passes. For example, if we were interested in an angle distribution
// measured from simulation:
//
// 1. Load trajectory frame.
// 3. Walk angles in trajectory frame, calling AddSample( theta ) for each angle.
// 4. Call Accumulate()
// 5. Repeat from step 1.
//
class Distribution
{
	protected:

		void save_header( FILE* f ) const
		{
			if( f == nullptr ) return;
			fprintf( f, "# x, y = coordinate and mean value at coordinate\n" );
			fprintf( f, "#\n" );
			fprintf( f, "# N(y)      = number of samples\n" );
			fprintf( f, "# StdDev(y) = standard deviation: sqrt(variance)\n" );
			fprintf( f, "# StdErr(y) = standard error of the estimated mean: sqrt(variance)/N = StdDev(y)/N\n" );
			fprintf( f, "# Sum(y) = sum over all y (should be a constant value)\n" );
			fprintf( f, "#\n" );
			fprintf( f, "# Both StdDev(y) and StdErr(y) use the SAMPLE variance\n" );
			fprintf( f, "# P(x) can be plotted as y/Sum(y)\n" );
			fprintf( f, "#\n" );
			fprintf( f, "# %12.12s  %12.12s  %12.12s  %12.12s  %12.12s  %12.12s\n", "x", "y", "N(y)", "StdDev(y)", "StdErr(y)", "Sum(y)" );
			fprintf( f, "#\n" );
		}

	public:

		MapStats pass, total;

		Distribution()                     { Reset( 1.0, 0.0 ); }
		Distribution( double d )           { Reset( d,   0.0 ); }
		Distribution( double d, double m ) { Reset( d,   m   ); }

		void Reset( double delta_, double min_ = 0.0 )
		{
			pass.Reset( delta_, min_ );
			total.Reset( delta_ , min_);
		}

		// NOI CLEAR() METHOD, AS UNCLEAR WHETHER TO CLEAR BOTH PASS AND TOTAL, JUST PASS ETC!

		void AddSampleByBin( int bin ) { pass.m[bin].AddSample( 1.0 ); }
		void AddSample( double coord )
		{
			int bin = (int)floor( (coord-pass.min) / pass.delta );
			AddSampleByBin( bin );
		}
		// Optional prefactor for accumulation.
		void Accumulate( double prefactor = 1.0 )
		{
			// Add total values (i.e. N*mean) from "pass" into "total".
			for( const auto& it : pass.m )
			{
				const auto& key = it.first;
				const auto& s = it.second;
				if( s.N < 1 ) continue; // ignore unsampled data, rather than adding sample of 0 to "total"!
				total.AddSampleByBin( key, prefactor*(s.mean*s.N) );
			}
			pass.Clear();
		}

		int Save( FILE*f, bool print_header = true ) const
		{
			if( f == nullptr ) return -1;
			if( print_header ) save_header( f );

			// Accumulate normalizing value, Sum(y)
			double sum_y = 0.0;
			for( const auto& it : total.m ) { sum_y += it.second.mean; }

			// Write data columns. Note Sum(y) as the final column, so y/Sum(x) = P(x): probability of x.
			for( const auto& it : total.m )
			{
				int bin_no = it.first;
				const Stats& s = it.second;
				double coord = total.min + (0.5+bin_no)*total.delta;
				fprintf( f, "  %e  %e  %12d  %e  %e  %e\n", coord, s.mean, (int)s.N, s.StdDev(), s.StdErr(), sum_y );
			}
			return 1;
		}
};

}

#endif
