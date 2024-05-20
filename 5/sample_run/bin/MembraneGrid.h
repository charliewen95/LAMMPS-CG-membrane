#if !defined( MEMBRANEGRID_INCLUDED )

#include <stdio.h>
#include <math.h>

#include <vector>

//
// Generate a Nx*Ny grid on the plane of the membrane, and calculate an average bilayer midplane
// location in each cell from the average z coordinate of tail beads that lie in that cell. We
// may then assign head groups to the upper or lower monolayer by comparing their z coords to the
// midplane of the local x,y cell: above the midplane => upper leflet, below the midplane => lower
// leaflet. THIS IS NOT FOOLPROOF! E.g., too few cells will give a bad estimate of the midplane in
// each cell, but too many can potentially lead to "empty" cells containing a head bead but no tail
// beads (so assignment to upper or lower leaflet fails) etc.
//

class MembraneGrid
{
	private:

		int ToGrid(
			double x, double y,
			double minx, double miny,
			double Lx, double Ly,
			int N_gx, int N_gy )
		{
			double x_hat = (x-minx)/Lx;
			double y_hat = (y-miny)/Ly;

			if( x_hat < 0.0 ) x_hat += 1.0;
			else if( x_hat >= 1.0 ) x_hat -= 1.0;

			if( y_hat < 0.0 ) y_hat += 1.0;
			else if( y_hat >= 1.0 ) y_hat -= 1.0;

			int gx = (int)floor( x_hat * N_gx );
			int gy = (int)floor( y_hat * N_gy );

			if( (gx<0) || (gy<0) || (N_gx<=gx) || (N_gy<=gy) )
			{
				fprintf( stderr, "Bad cell %d,%d: must be 0 <= (x,y) < (%d,%d). Coordinate: %g,%g (cell %g,%g : %g,%g)\n",
					gx,gy, N_gx,N_gy,
					x,y, minx,miny, Lx,Ly );
				return -1;
			}

			return (gy*N_gx) + gx;
		}


	public:

		std::vector<int> counts;      // number of midplane samples in each cell
		std::vector<double> midplane; // midplane coordinate in each cell

		int Assign(
			const std::vector<double>& x_vec,
			const std::vector<double>& y_vec,
			const std::vector<double>& z_vec,
			const std::vector<int>& head_indices,
			const std::vector<int>& tail_indices,
			double minx, double miny,
			double maxx, double maxy,
			int N_gx, int N_gy,
			std::vector<int>& upper,
			std::vector<int>& lower )
		{
			if( (N_gx<1) || (N_gy<1) )
			{
				fprintf( stderr, "Weird membrane grid size: %d x %d\n", N_gx, N_gy );
				return -1;
			}
			if( x_vec.size() < 1 )
			{
				fprintf( stderr, "Weird x data length %d\n", (int)x_vec.size() );
				return -1;
			}
			if( y_vec.size() < 1 || y_vec.size() != x_vec.size() )
			{
				fprintf( stderr, "Weird y data length %d\n", (int)y_vec.size() );
				return -1;
			}
			if( z_vec.size() < 1 || z_vec.size() != x_vec.size() )
			{
				fprintf( stderr, "Weird z data length %d\n", (int)z_vec.size() );
				return -1;
			}

			double Lx = maxx - minx;
			double Ly = maxy - miny;
			int N_atoms = (int)x_vec.size();

			//
			// Clear cell info & monolayer index arrays
			//
			{
				counts.resize( N_gx*N_gy );
				midplane.resize( counts.size() );
				for( size_t i=0; i<counts.size(); i++ )
				{
					counts[i] = 0;
					midplane[i] = 0.0;
				}

				upper.clear();
				lower.clear();
			}

			//
			// Calculate grid cell midplane values from the tail particles.
			//
			for( auto i : tail_indices )
			{
				if( i<0 || i>=N_atoms ) return -1;

				double x = x_vec[i];
				double y = y_vec[i];
				double z = z_vec[i];

				// Get grid index for this coordinate
				int grid_index = ToGrid( x,y, minx,miny, Lx,Ly, N_gx,N_gy );
				if( grid_index == -1 ) return -1;

				// Update midplane information for the current grid cell
				midplane[grid_index] += z;
				counts[grid_index] += 1;
			}

			//
			// Convert midplane acumulation to average, warn re. cells that contain no midplane samples.
			//
			int undefined_midplane = 0;
			for( size_t i=0; i<midplane.size(); i++ )
			{
				if( counts[i] > 0 ) midplane[i] /= counts[i];
				else undefined_midplane++;
			}
			if( undefined_midplane > 0 )
			{
				printf( "WARNING: %d of %d grid squares have undefined midplane coordinate!\n", undefined_midplane, (int)counts.size() );
			}

			//
			// Assign head groups to upper or lower monolayer according to their position relative to midplanes
			//
			for( auto i : head_indices )
			{
				if( i<0 || i>=N_atoms ) return -1;

				double x = x_vec[i];
				double y = y_vec[i];
				double z = z_vec[i];
				
				// Get grid index for this coordinate
				int grid_index = ToGrid( x,y, minx,miny, Lx,Ly, N_gx,N_gy );
				if( grid_index == -1 ) return -1;

				// Check z position vs midplane in this grid cell
				if( z >= midplane[grid_index] ) upper.push_back( i );
				else lower.push_back( i );
			}

			return 1;
		}


};

#define MEMBRANEGRID_INCLUDED
#endif
