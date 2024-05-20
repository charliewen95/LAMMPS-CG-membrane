#include <stdio.h>
#include <stdlib.h>

#include "LAMMPS.h"
#include "Fourier.h"
#include "MembraneGrid.h"
#include "Stats.h"
#include "GridFilter.h"

#include <map>

/*
	Calculate the projected 1-D spectrum from a 2-D discrete Fourier transform of a set of points representing a membrane.

	Building:

		clang++  -Wall -Wextra -pedantic -std=c++11 -O2          Spectrum.cpp -o spectrum
		g++-mp-5 -Wall -Wextra -pedantic -std=c++11 -O2 -fopenmp Spectrum.cpp -o spectrum_openmp

	Use:

		./spectrum  traj=path  max_k=X max_l=X  head_type=X[,X,...] tail_type=X[,X,...]  gx=X gy=X  [delta_q=X] [scale=X] [which=X] [histogram=X] [save_midplane=X] [out_prefix=X] [filter=X] [start=X] [stop=X]

	Where:

		- max_k, max_l : max integer wave numbers on x and y axes respectively.
	 	- head_type, tail_type : LAMMPS atom types for head and terminal tail beads.
	 	- gx, gy  : number of grid cells on x and y axes for midplane calculations.

	 	- delta_q   : OPTIONAL resolution of output spectrum histogram (default: 0.05).
	 	- scale     : OPTIONAL scaling for input->output length units (default: 1.0).
	 	- which     : OPTIONAL switch for which monolayer to use: both, upper, lower, midplane (default: both)
	 	- histogram : OPTIONAL bin width for per-type z density histogram (default: no histograms generated)
	 	- save_midplane : OPTIONAL flag to save the midplane coordinates as xyz (default: no midplane written)
		- out_prefix: OPTIONAL output spectrum file prefix (default: 'spectrum').
		- filter: OPTIONAL grid filter cell size, with contents of isolated cells ignored (default: no filtering).
		- start: OPTIONAL unit-based start frame in trajectory. Negative values ignored (default: -1).
		- stop: OPTIONAL unit-based stop frame in trajectory. Negative values ignored (default: -1).

	 Example:

		 ./spectrum traj=dump.lammpstrj  max_k=10 max_l=10  head_type=1 tail_type=3  gx=20 gy=20  scale=0.1  delta_q=0.05

	 This load dump.lammpstrj, a LAMMPS trajectory of a CG bilayer with lipid head type 1 and tail type 3. Wave vectors of
	 (qx,qy) = 2pi.( k/Lx, l/Ly ) (0 <= k <= 10, 0 <= l <= 10) are used to generate the Fourier coefficients, with the upper and
	 lower monolayers assigned according to head bead z position being above/below the midplanes calculated in a 2-D grid of
	 dimension 20x20 cells on the simulation box.

	 In this case, a scaling factor of 0.1 is applied to all input coordinates and lengths, to e.g. swap from Angstrom (input) to
	 nanometres (output).

	 The histogram bin width is 0.05 x^{-1}, where x is the output length unit (i.e. nm in this case, as the input length scale is
	 Angstrom and we're using a length scale factor of 0.1).

	 The output is two files, "spectrum.standard.txt" and "spectrum.modified.txt", that contain the 1-D spectral information.

	 Analysis, via the structure factor < |u(q)|^2 >:

	 ( e.g. Eqn. 4 and 13 of Brandt et al, Biophys J 100:2104-2111 )

	 	N.< |u(q)|^2 > = kBT / ( APL.kC.q^4 )

	 Where we ignore a surface tension dependent  1/q^2 term, leaving:

		 N = number of lipids in monolayer (or number of samples) fed into discrete Fourier analysis equation
		 u(q) = complex coefficient for q, q = sqrt( qx^2 + qy^2 ), where (qx,qy) = 2pi( k/Lx, l/Ly )
		 APL = area per lipid (or area per sample)
		 kC = bending modulus

	Here, we assume that the coefficients, u(q), are generated as (1\N).[ \Sum_n f(r_n) . exp( -i.2pi.q.r_n) ].
	To avoid the area per lipid value (or area per sample) and simplify the expression, we note that APL = A/N
	and so:

	N.< |u(q)|^2 > = kBT / ( (A/N).(kC.q^4) )
	< |u(q)|^2 > = kBT / ( A.(kC.q^4) )
	A.< |u(q)|^2 > = kBT / ( kC.q^4 )
	< A.|u(q)|^2 > = 1.0 / ( kC.q^4 )

	... with A the instantaneous projected membrane area and kC in units of kBT. This modified structure factor is saved
	in the file "spectrum.modified.txt".

	Notes:

	1. Try not to make the midplane calculation grid too coarse (i.e. small gx,gy), as if the grid squares are large you could
	   have the same monolayer oscillating up and down through the midplane inside a grid square, so pieces of the same monolayer
	   would be assigned to the upper and lower leaflets. This should not affect things where "which=both".

	   Note that the spectrum from the midplane grid ("which=midplane") for extremely small grid squares (i.e. gx,gy so large
	   that we have less than one lipid per grid square on average) should converge onto the spectral results of a single
	   monolayer at small q, as it's basically performing a discrete fourier transform using individual lipid positions, akin
	   to analysis of a single monolayer.

	2. The filtering assigns head group particles to a 3D grid with the specified size. Any cells with no neighbours have
	   their contents ignored.

	3. Start and stop frames are UNIT BASED and INCLUSIVE!
*/

//
// Wrapper for some variable and code to process trajectory data.
//
struct Data
{
	//
	// Input parameters from the command line
	//
	std::string in_path, out_prefix;
	int max_k, max_l;
	int N_gx, N_gy;
	double delta_q, min_q;
	double length_scale;
	std::string which_monolayer;
	double hist_bin_width;
	int save_midplane;
	int start_frame, stop_frame;
	double filter_size;
	int save_raw;

	std::map<int,int> head_types, tail_types;

	//
	// Internal data
	//
	MembraneGrid membrane_grid;
	Stats::MapStats standard_spectrum, modified_spectrum;
	Fourier::Coeffs2D F;

	//
	// In case we want to save the "raw" spectral values. These are
	// regenerated for every trajectory frame.
	//
	std::vector<double> raw_qx, raw_qy, raw_standard, raw_modified;

	bool IsHead( int type ) const
	{
		if( head_types.find(type) != head_types.end() ) return true;
		else return false;
	}
	bool IsTail( int type ) const
	{
		if( tail_types.find(type) != tail_types.end() ) return true;
		else return false;
	}

	Data()
	{
		in_path = "";
		out_prefix = "spectrum";
		max_k = max_l = -1;
		N_gx = N_gy = -1;
		
		delta_q = 0.05;
		min_q = 0.0;

		length_scale = 1.0;

		which_monolayer = "both";

		hist_bin_width = -1.0; // no histograms by default.
		save_midplane = -1; // do not save midplane by default

		filter_size = -1.0; // off by default

		start_frame = stop_frame = -1; // off by default

		save_raw = -1;// off by default
	}

	// Convert string to integer value
	template<typename T> int to_int( const std::string& what, const std::string& str, T& val )
	{
		if( StringUtil::ToInteger(str,val) == -1 )
		{
			printf( "Unable to convert '%s' into integer %s\n", str.c_str(), what.c_str() );
			return -1;
		}
		return 1;
	}

	// Convert string to double precision float value
	int to_dbl( const std::string& what, const std::string& str, double& val )
	{
		if( StringUtil::ToReal(str,val) == -1 )
		{
			printf( "Unable to convert '%s' into real %s\n", str.c_str(), what.c_str() );
			return -1;
		}
		return 1;
	}

	// Extract and sanity check parameters from command line
	int GetParams( int argc, char** argv )
	{
		std::vector<std::string> tokens, subtoks;

		for( int i=1; i<argc; i++ )
		{
			StringUtil::Tokenize( argv[i], tokens, "=" );
			if( tokens.size() < 2 ) continue;

			const std::string& key = tokens[0];
			const std::string& val = tokens[1];

			if( key == "traj" ) in_path = val;
			else if( key == "out_prefix" ) out_prefix = val;
			else if( key == "which" ) which_monolayer = val;
			else if( (key=="max_k") && (to_int(key,val,max_k)==-1) ) return -1;
			else if( (key=="max_l") && (to_int(key,val,max_l)==-1) ) return -1;
			else if( key=="head_type" )
			{
				StringUtil::Tokenize( val.c_str(), subtoks, "," );
				if( subtoks.size() < 1 ) return -1;
				for( const auto& t : subtoks )
				{
					int i;
					if( to_int(key,t,i)==-1 ) return -1;
					head_types[i] = 1;
				}
			}
			else if( key=="tail_type" )
			{
				StringUtil::Tokenize( val.c_str(), subtoks, "," );
				if( subtoks.size() < 1 ) return -1;
				for( const auto& t : subtoks )
				{
					int i;
					if( to_int(key,t,i)==-1 ) return -1;
					tail_types[i] = 1;
				}
			}
			else if( (key=="gx") && (to_int(key,val,N_gx)==-1) ) return -1;
			else if( (key=="gy") && (to_int(key,val,N_gy)==-1) ) return -1;
			else if( (key=="delta_q") && (to_dbl(key,val,delta_q)==-1) ) return -1;
			else if( (key=="min_q") && (to_dbl(key,val,min_q)==-1) ) return -1;
			else if( (key=="scale") && (to_dbl(key,val,length_scale)==-1) ) return -1;
			else if( (key=="histogram") && (to_dbl(key,val,hist_bin_width)==-1) ) return -1;
			else if( (key=="save_midplane") && (to_int(key,val,save_midplane)==-1) ) return -1;
			else if( (key=="filter") && (to_dbl(key,val,filter_size)==-1) ) return -1;
			else if( (key=="start") && (to_int(key,val,start_frame)==-1) ) return -1;
			else if( (key=="stop") && (to_int(key,val,stop_frame)==-1) ) return -1;
			else if( (key=="save_raw") && (to_int(key,val,save_raw)==-1) ) return -1;
		}

		if( in_path == "" ) return -1;
		if( out_prefix == "" ) return -1;

		if( (which_monolayer!="upper") &&
			(which_monolayer!="lower") &&
			(which_monolayer!="both") &&
			(which_monolayer!="midplane") ) return -1;

		if( max_k < 0 || max_l < 0 ) return -1;
		if( head_types.size() < 1 || tail_types.size() < 1 ) return -1;
		if( N_gx < 1 || N_gy < 0 ) return -1;

		if( delta_q <= 0.0 ) return -1;
		if( length_scale <= 0.0 ) return -1;

		if( filter_size == 0.0 ) return -1;

		if( (start_frame>0 && stop_frame>0) && (stop_frame<start_frame) ) return -1;

		standard_spectrum.Reset( delta_q, min_q );
		modified_spectrum.Reset( delta_q, min_q );

		return 1;
	}

	int PrintParams( FILE* f = nullptr )
	{
		if( f == nullptr ) f = stdout;

		fprintf( f, "#\n" );
		fprintf( f, "# Parameters:\n" );
		fprintf( f, "#\n" );
		fprintf( f, "#  Input: %s\n", in_path.c_str() );
		fprintf( f, "#  Output prefix: %s\n", out_prefix.c_str() );
		fprintf( f, "#  max_k, max_l : %d, %d\n", max_k, max_l );

		fprintf( f, "#  head types : " );
		for( const auto&it : head_types ) fprintf( f, "%d ", it.first );
		fprintf( f, "\n" );

		fprintf( f, "#  tail types : " );
		for( const auto&it : tail_types ) fprintf( f, "%d ", it.first );
		fprintf( f, "\n" );

		fprintf( f, "#  gx, gy : %d, %d\n", N_gx, N_gy );
		fprintf( f, "#  delta_q : %g\n", delta_q );
		fprintf( f, "#  scale : %g\n", length_scale );
		fprintf( f, "#  which_monolayer : %s\n", which_monolayer.c_str() );
		fprintf( f, "#  histogram width : %g\n", hist_bin_width );
		fprintf( f, "#  save_midplane   : %d\n", save_midplane );
		fprintf( f, "#  start, stop frames : %d, %d\n", start_frame, stop_frame );
		fprintf( f, "#  save_raw : %d\n", save_raw );
		#if defined( _OPENMP )
			fprintf( f, "#  OpenMP support: omp_get_max_threads() = %d\n", omp_get_max_threads() );
		#else
			fprintf( f, "#  Serial version\n" );
		#endif
		if( filter_size >= 0.0 ) fprintf( f, "#  filter_size : %g\n", filter_size );
		fprintf( f, "#\n" );

		return 1;
	}

	//
	// Determine Fourier coefficients and accumulate spectral data from the
	// specified surface approximation: f[i] = f(x[i],y[i]), with x,y on
	// range [0,Lx) and [0,Ly) and f[i] being the surface "height" (with
	// \Sum_i f[i] = 0 assumed).
	//
	// Note: prefactor should be 1.0 where we're feeding in a single surface/monolayer. If we're
	// using both monolayers of a bilayer, prefactor should be 0.5 to average over the two
	// monolayers in the components of the standard structure factor.
	//
	int Process(
		const std::vector<double>& x,
		const std::vector<double>& y,
		const std::vector<double>& f,
		double Lx, double Ly,
		double prefactor )
	{
		constexpr double twopi = 2.0*M_PI;

		size_t N = f.size();

		if( N < 1 ) return -1;

		//
		// Get Fourier coefficients for the input data => F
		//
		{
			int result = Fourier::Analysis(
				x, y, f,
				Lx, Ly,
				max_k, max_l,
				F );

			if( result == -1 )
			{
				printf( "%s(): Unable to determine Fourier coefficients\n", __func__ );
				return -1;
			}
		}

		//
		// Accumulate 2-D coefficient data into 1-D spectrum.
		// Cutoff = 2pi/L, were L is the largest of Lx and Ly.
		//
		raw_qx.clear();
		raw_qy.clear();
		raw_standard.clear();
		raw_modified.clear();
		{
			double L = (Lx>Ly) ? (Lx) : (Ly);
			double cutoff = twopi/L;

			for( int k=0; k<=max_k; k++ )
			{
				double qx = (twopi/Lx)*k;
				for( int l=0; l<=max_l; l++ )
				{
					double qy = (twopi/Ly)*l;
					double q = sqrt( qx*qx + qy*qy );

					if( q < cutoff ) continue;

					// Target: N < |u(q)|^2 >, where:
					//   - q = sqrt( qx^2 + qy^2 )
					//   - assumes coefficients u(q) were divided by N when generated
					//   - the prefactor should be 1.0 where a single monolayer/surface is
					//     present, else 0.5 to average over the two monolayers/surfaces.
					Fourier::Complex c = std::abs( F[k][l] );
					c *= c;
					standard_spectrum.AddSample( q, c.real()*prefactor*N );

					// Target: < A |u(q)|^2 >, where:
					//   - A = instantaneous projected bilayer area
					//   - q = sqrt( qx^2 + qy^2 )
					//   - assumes coefficients u(q) were divided by N when generated
					//   - we don't need to use the prefactor here, as we treat this the
					//     same regardless of whether there's one or two monolayers (we
					//     have fortuitous cancellation of the N values!).
					modified_spectrum.AddSample( q, (Lx*Ly)*c.real() );

					//
					// Store current raw values, in case we're saving them to file.
					//
					raw_qx.push_back( qx );
					raw_qy.push_back( qy );
					raw_standard.push_back( c.real()*prefactor*N );
					raw_modified.push_back( (Lx*Ly)*c.real() );
				}
			}
		}

		return 1;
	}

};


int main( int argc, char** argv )
{
	FILE *in_traj;

	LAMMPS::TrajectoryFrame frame;
	Data data;

	std::vector<int> heads, tails, upper, lower, combined;
	std::vector<double> x, y, f;

	// Map LAMMPS atom type -> density histogram
	std::map<int,Stats::Distribution*> rho_map;

	std::vector<int> filtered_indices; // only used where filtering enabled.
	GridFilter gf;

	Stats::Stats AreaPerLipid, AreaPerLipidSquared, LipidsPerMonolayer;

	if( data.GetParams(argc,argv) == -1 )
	{
		printf( "\n" );
		printf( "Usage: %s  traj=path  max_k=X max_l=X  head_type=X[,X,...] tail_type=X[,X,...]  gx=X gy=X  [delta_q=X] [scale=X] [which=X] [histogram=X] [out_prefix=X] [filter=X]\n", argv[0] );
		printf( "\n" );
		printf( "Where:\n" );
		printf( "\n" );
		printf( " - max_k, max_l : max integer wave numbers on x and y axes respectively.\n" );
		printf( " - head_type, tail_type : LAMMPS atom types for head and terminal tail beads.\n" );
		printf( " - gx, gy  : grid cell counts on x and y axes for midplane calculations.\n" );
		printf( "\n" );
		printf( " - delta_q   : OPTIONAL resolution of output spectrum histogram (default: 0.05).\n" );
		printf( " - scale     : OPTIONAL scaling for input->output length units (default: 1.0).\n" );
		printf( " - which     : OPTIONAL setting for which monolayer: 'upper', 'lower', 'both', 'midplane' (default: 'both').\n" );
		printf( " - histogram : OPTIONAL per-type histogram bin width in OUTPUT length units (ignored where <= 0.0).\n" );
		printf( " - save_midplane : OPTIONAL flag to save the midplane coordinates as xyz (default: no midplane written).\n" );
		printf( " - out_prefix: OPTIONAL output spectrum file prefix (default: 'spectrum').\n" );
		printf( " - filter: OPTIONAL grid filter cell size, with contents of isolated cells ignored (default: no filtering).\n" );
		printf( " - start: OPTIONAL unit-based start frame in trajectory. Negative values ignored (default: -1).\n" );
		printf( " - stop: OPTIONAL unit-based stop frame in trajectory. Negative values ignored (default: -1).\n" );
		printf( " - save_raw: OPTIONAL flag to save raw (i.e. non-binned) spectral values. Negative values ignored (default: -1).\n" );
		printf( "\n" );
		printf( "Notes:\n" );
		printf( "\n" );
		printf( "Be careful if you use 'which=midplane', as the high frequency components of the spectrum (larger q values)\n" );
		printf( "will be limited by the resolution of the midplane grid (as specified by 'gx' and 'gy' parameters).\n" );
		printf( "\n" );
		printf( "The filtering assigns head group particles to a 3D grid with the specified size. Any cells with no neighbours\n" );
		printf( "have their contents ignored.\n" );
		exit( -1 );
	}
	data.PrintParams();

	if( (in_traj=fopen(data.in_path.c_str(),"r")) == nullptr )
	{
		printf( "Unable to open trajectory '%s'.\n", data.in_path.c_str() );
		exit( -1 );
	}

	//
	// Delete existing midplane trajectory, if present.
	//
	if( data.save_midplane > 0 )
	{
		std::string path = data.out_prefix + ".midplane.xyz";
		FILE* f = fopen( path.c_str(), "w" );
		fclose( f );
	}

	//
	// Delete existing raw values file, if present.
	//
	if( data.save_raw > 0 )
	{
		std::string path = data.out_prefix + ".raw";
		FILE* f = fopen( path.c_str(), "w" );	
		fprintf( f, "# qx, qy, standard, modified. See spectral file headers for more info.\n" );	
		fclose( f );
	}

	int frame_no = 0;
	while( true )
	{
		frame_no++;

		//
		// Try to load a trajectory frame
		//
		{
			int result = frame.Load( in_traj );

			if( result == -1 || frame.x_vec.size() < 3 )
			{
				printf( "Unable to load trajectory frame. Stopping here, frame %d\n", frame_no );
				break;
			}
		}

		if( (data.start_frame>0) && (frame_no<data.start_frame) )
		{
			printf( "Skipping frame %d ...\n", frame_no );
			continue;
		}
		if( (data.stop_frame>0) && (frame_no>data.stop_frame) )
		{
			printf( "Stopping on frame %d.\n", frame_no );
			break;
		}

		//
		// Extract/convert some commonly used info from the frame
		//
		double minx = frame.mins[0] * data.length_scale;
		double miny = frame.mins[1] * data.length_scale;
		double minz = frame.mins[2] * data.length_scale;

		double maxx = frame.maxs[0] * data.length_scale;
		double maxy = frame.maxs[1] * data.length_scale;
		double maxz = frame.maxs[2] * data.length_scale;

		double Lx = maxx - minx;
		double Ly = maxy - miny;
		double Lz = maxz - minz;

		for( size_t i=0, max_i=frame.x_vec.size(); i<max_i; i++ )
		{
			frame.x_vec[i] *= data.length_scale;
			frame.y_vec[i] *= data.length_scale;
			frame.z_vec[i] *= data.length_scale;
		}

		//
		// Determine indices of head & tail particles in the data set, and accumulate histogram information.
		//
		{
			heads.clear();
			tails.clear();

			for( size_t i=0; i<frame.atom_types.size(); i++ )
			{
				int atom_type = frame.atom_types[i];
				double z = frame.z_vec[i];

				// Detect heads/tails
				if( data.IsHead(atom_type) == true ) heads.push_back( i );
				if( data.IsTail(atom_type) == true ) tails.push_back( i );

				// Acumulate histogram data for this frame, if needed.
				if( data.hist_bin_width > 0.0 )
				{
					auto it = rho_map.find( atom_type );
					if( it == rho_map.end() )
					{
						Stats::Distribution* rho = new Stats::Distribution( data.hist_bin_width );
						rho->AddSample( z );
						rho_map[atom_type] = rho;
					}
					else { it->second->AddSample( z ); }
				}
			}

			if( data.hist_bin_width > 0.0 )
			{
				double volume = (Lx*Ly*data.hist_bin_width);
				for( auto& it : rho_map )
				{
					it.second->Accumulate( 1.0/volume ); // convert counts -> densities!
				}
			}
		}

		//
		// Filter data, if needed
		//
		if( data.filter_size > 0.0 )
		{
			int result;

			//
			// Heads ...
			//
			result = gf.Filter(
				frame.x_vec, frame.y_vec, frame.z_vec,
				heads,
				minx,miny,minz, maxx,maxy,maxz,
				data.filter_size,
				frame.PBC[0], frame.PBC[1], frame.PBC[2],
				filtered_indices );
			if( result == -1 )
			{
				printf( "Unable to filter head coordinates!\n" );
				exit( -1 );
			}
			printf( "%d/%d heads passed filtering.\n", result, (int)heads.size() );
			heads = filtered_indices;

			//
			// Tails ...
			//
			result = gf.Filter(
				frame.x_vec, frame.y_vec, frame.z_vec,
				tails,
				minx,miny,minz, maxx,maxy,maxz,
				data.filter_size,
				frame.PBC[0], frame.PBC[1], frame.PBC[2],
				filtered_indices );
			if( result == -1 )
			{
				printf( "Unable to filter tail coordinates!\n" );
				exit( -1 );
			}
			printf( "%d/%d tails passed filtering.\n", result, (int)tails.size() );
			tails = filtered_indices;
		}


		//
		// Assign head particle indices to upper or lower monolayer leaflets.
		// Also generate a combined set of all head particle indices.
		//
		{
			upper.clear();
			lower.clear();
			combined.clear();

			int result = data.membrane_grid.Assign(
				frame.x_vec, frame.y_vec, frame.z_vec,
//				heads, tails,
				heads, heads,
				minx, miny,
				maxx, maxy,
				data.N_gx, data.N_gy,
				upper, lower );

			if( result == -1 )
			{
				printf( "Unable to assign monolayer leaflets. Stopping here, frame %d\n", frame_no );
				break;
			}

			combined = upper;
			for( auto& i : lower ) combined.push_back( i );
		}

		//
		// Transform & accumulate data.
		// We first need to extract the x,y,f data from somewhere ...
		//
		{
			x.clear();
			y.clear();
			f.clear();

			//
			// Get x,y,f data from the appropriate source
			//

			if( data.which_monolayer == "upper" )
			{
				for( auto i : upper )
				{
					double x_ = frame.x_vec[i];
					double y_ = frame.y_vec[i];
					double f_ = frame.z_vec[i];

					x.push_back( x_ );
					y.push_back( y_ );
					f.push_back( f_ );
				}
			}
			else if( data.which_monolayer == "lower" )
			{
				for( auto i : lower )
				{
					double x_ = frame.x_vec[i];
					double y_ = frame.y_vec[i];
					double f_ = frame.z_vec[i];

					x.push_back( x_ );
					y.push_back( y_ );
					f.push_back( f_ );
				}
			}
			else if( data.which_monolayer == "both" )
			{
				for( auto i : lower )
				{
					double x_ = frame.x_vec[i];
					double y_ = frame.y_vec[i];
					double f_ = frame.z_vec[i];

					x.push_back( x_ );
					y.push_back( y_ );
					f.push_back( f_ );
				}

				for( auto i : upper )
				{
					double x_ = frame.x_vec[i];
					double y_ = frame.y_vec[i];
					double f_ = frame.z_vec[i];

					x.push_back( x_ );
					y.push_back( y_ );
					f.push_back( f_ );
				}
			}
			else if( data.which_monolayer == "midplane" )
			{
				int N_gx = data.N_gx;
				int N_gy = data.N_gy;

				for( int xi=0; xi<N_gx; xi++ )
				{
					double x_ = (0.5+xi)*(Lx/N_gx);
					for( int yi=0; yi<N_gy; yi++ )
					{
						int index = (yi*N_gx) + xi;

						double y_ = (0.5+yi)*(Ly/N_gy);
						double f_ = data.membrane_grid.midplane[index]; // midplane generated from scaled particle coords

						if( data.membrane_grid.counts[index] < 1 ) continue;

						x.push_back( x_ );
						y.push_back( y_ );
						f.push_back( f_ );
					}
				}
			}
			else
			{
				printf( "Unknown monolayer input data'%s'!\n", data.which_monolayer.c_str() );
				exit( -1 );
			}

			//
			// Skip empty input data, just in case.
			//
			if( f.size() < 1 ) continue;

			//
			// Sanitise input data:
			// - x,y coords wrapped to range [0,Lx) and [0,Ly).
			// - \Sum_i f[i] = 0.
			//
			{
				double average_f = 0.0;
				for( size_t i=0; i<f.size(); i++ )
				{
					double x_ = x[i];
					double y_ = y[i];
					double f_ = f[i];

					// Ensure x,y ranges [0,Lx), [0,Ly)
					x_ -= minx;
					y_ -= miny;

					// Wrap, if needed.
					if( x_ < 0 ) x_ += Lx;
					else if( x_ >= Lx ) x_ -= Lx;

					if( y_ < 0 ) y_ += Ly;
					else if( y_ >= Ly ) y_ -= Ly;

					x[i] = x_;
					y[i] = y_;
					f[i] = f_;
					average_f += f_;
				}
				average_f /= f.size();
				for( auto& fi : f ) fi -= average_f;
			}

			//
			// Process input data & accumulate derived information
			//
			double prefactor = (data.which_monolayer=="both") ? (0.5) : (1.0);
			int result = data.Process( x, y, f, Lx, Ly, prefactor );
			if( result == -1 )
			{
				printf( "Unable to process frame data. Stopping here, frame %d\n", frame_no );
				break;
			}
		}

		//
		// Let the user know where we're up to
		//
		{
			printf( "Frame %d, dims %g,%g,%g -> %g,%g,%g : %g,%g,%g, %d heads %d tails, %d upper lipids %d lower\n",
				frame_no,
				minx, miny, minz,
				maxx, maxy, maxz,
				Lx, Ly, Lz,
				(int)heads.size(), (int)tails.size(),
				(int)lower.size(), (int)upper.size() );
		}

		//
		// Accumulate some stats
		//
		{
			double Nm = 0.5*(lower.size()+upper.size()); // lipids per monolayer
			double Apl = (Lx*Ly)/Nm; // area per lipid
			LipidsPerMonolayer.AddSample( Nm );
			AreaPerLipid.AddSample( Apl );
			AreaPerLipidSquared.AddSample( Apl*Apl );
		}

		//
		// Save spectral info and histograms:
		//
		// "standard" spectrum = structure factor as per usual.
		// "modified" spectrum = see notes at the start of this file.
		// histograms = one per LAMMPS bead type, if required
		// midplane = if required
		//
		{
			// Standard spectrum
			{
				std::string path = data.out_prefix + ".standard.txt";
				FILE* f = fopen( path.c_str(), "w" );
				if( f == nullptr )
				{
					printf( "Unable to open output file '%s'\n", path.c_str() );
					break;
				}
				data.PrintParams( f );
				fprintf( f, "# Output:\n" );
				fprintf( f, "#\n" );
				fprintf( f, "#  q = projected wavelength: sqrt( qx^2 + qy^2 )\n" );
				fprintf( f, "#  S_q  = structure factor for projected wavelength q: N.< |u(q)|^2 >\n" );
				fprintf( f, "#\n" );
				fprintf( f, "# Fitting a bending modulus, kC:\n" );
				fprintf( f, "#\n" );
				fprintf( f, "# S_q = 1.0/( APL.kC.q^4 ) for small q in the tension-free limit, where\n" );
				fprintf( f, "# APL is the average area per lipid.\n" );
				fprintf( f, "#\n" );
				fprintf( f, "# %12.12s  %12.12s  %12.12s  %12.12s %12.12s\n", "q", "S_q", "N", "StdDev(S_q)", "StdErr(S_q)" );
				fprintf( f, "#\n" );
				data.standard_spectrum.Save( f, false );
				fclose( f );
			}

			// Modified spectrum
			{
				std::string path = data.out_prefix + ".modified.txt";

				FILE* f = fopen( path.c_str(), "w" );
				if( f == nullptr )
				{
					printf( "Unable to open output file '%s'\n", path.c_str() );
					break;
				}
				data.PrintParams( f );
				fprintf( f, "# Output:\n" );
				fprintf( f, "#\n" );
				fprintf( f, "#  q = projected wavelength: sqrt( qx^2 + qy^2 )\n" );
				fprintf( f, "#  S_q  = modified structure factor for projected wavelength q: S_q = < A |u(q)|^2 >\n" );
				fprintf( f, "#\n" );
				fprintf( f, "# Fitting a bending modulus, kC:\n" );
				fprintf( f, "#\n" );
				fprintf( f, "# S_q = 1.0/( kC.q^4 ) for small q in the tension-free limit.\n" );
				fprintf( f, "#\n" );
				fprintf( f, "# %12.12s  %12.12s  %12.12s  %12.12s %12.12s\n", "q", "S_q", "N", "StdDev(S_q)", "StdErr(S_q)" );
				fprintf( f, "#\n" );
				data.modified_spectrum.Save( f, false );
				fclose( f );
			}

			//
			// Histograms
			//
			if( data.hist_bin_width > 0.0 )
			{
				char path[1024];
				for( auto& it : rho_map )
				{
					int atom_type = it.first;
					const Stats::Distribution* dist = it.second;

					sprintf( path, "%s.type_%d.histogram", data.out_prefix.c_str(), atom_type );
					FILE* f = fopen( path, "w" );
					if( f == nullptr )
					{
						printf( "Unable to open output file '%s'\n", path );
						break;
					}
					data.PrintParams( f );
					dist->Save( f );
					fclose( f );
				}
			}

			//
			// Midplane info
			//
			if( data.save_midplane > 0 )
			{
				int ngx = data.N_gx;
				int ngy = data.N_gx;
				double dx = Lx/ngx;
				double dy = Ly/ngy;
				std::string path = data.out_prefix + ".midplane.xyz";

				FILE* f = fopen( path.c_str(), "a" );
				if( f == nullptr )
				{
					printf( "Unable to open output file '%s'\n", path.c_str() );
					break;
				}
				fprintf( f, "%d\n", ngx * ngy );
				fprintf( f, "Midplane coordinates. NaN indicates no samples in grid cell!\n" );
				for( int iy = 0; iy<ngy; iy++ )
				{
					double y = miny + (0.5+iy) * dy;
					for( int ix = 0; ix<ngx; ix++ )
					{
						double x = minx + (0.5+ix) * dx;
						int index = (iy*ngx) + ix;
						if( data.membrane_grid.counts[index] < 1 ) fprintf( f, "%s %g %g NaN\n", "m", x, y );
						else fprintf( f, "%s %g %g %g\n", "m", x/data.length_scale, y/data.length_scale, data.membrane_grid.midplane[index]/data.length_scale ); // original coords
//						else fprintf( f, "%s %g %g %g\n", "m", x, y, data.membrane_grid.midplane[index] ); // scaled coords
					}
				}
				fclose( f );
			}

			//
			// Raw values
			//
			if( data.save_raw > 0 )
			{
				std::string path = data.out_prefix + ".raw";
				FILE* f = fopen( path.c_str(), "a" );
				for( size_t i=0; i<data.raw_qx.size(); i++ )
				{
					fprintf( f, "%e  %e  %e  %e\n", data.raw_qx[i], data.raw_qy[i], data.raw_standard[i], data.raw_modified[i] );
				}
				fclose( f );
			}

			//
			// Filtered data, if needed
			//
			/*
			if( data.filter_size > 0.0 )
			{
				char path[512];
				sprintf( path, "%s.%d.xyz", data.out_prefix.c_str(), frame_no );
				FILE* fptr = fopen( path, "a" );
				if( fptr == nullptr )
				{
					printf( "Unable to open output file '%s'\n", path );
					break;
				}
				fprintf( fptr, "%d\n", (int)x.size() );
				fprintf( fptr, "Filtered coordinates.\n" );
				for( size_t i=0, max_i=x.size(); i<max_i; i++ )
				{
					fprintf( fptr, "%s %g %g %g\n", "f", x[i]/data.length_scale, y[i]/data.length_scale, f[i]/data.length_scale );
				}
				fclose( fptr );
			}
			*/
		}
	}

	fclose( in_traj );

	{
		double Nm = LipidsPerMonolayer.mean;
		double Apl = AreaPerLipid.mean;
		double Apl2 = AreaPerLipidSquared.mean;
		double kA = Apl/( Nm*(Apl2-(Apl*Apl)) );

		printf( "< Nm >    = %g\n", Nm );
		printf( "< Apl >   = %g\n", Apl );
		printf( "< Apl^2 > = %g\n", Apl2 );
		printf( "kA = %g kBT . length^-2 \n", kA );
	}

	return 0;
}

