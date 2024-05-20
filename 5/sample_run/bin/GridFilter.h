#if !defined( LINKCELLS_CORE_DEFINED )
#define LINKCELLS_CORE_DEFINED

#include <stdio.h>
#include <math.h>

#include <vector>
#include <map>

//
// 3 element vector.
//
template <typename T> struct Vec3
{
	T x, y, z;

	//
	// Define some standard operators, so we can use GridPoint structures as std::map keys etc.
	// All sorting is performed in order: z, y, x (ie z assumed to be "most significant").
	//
	bool operator == (const Vec3<T> &rhs) const
	{
		return ( (x==rhs.x) && (y==rhs.y) && (z==rhs.z) );
	}

	bool operator < (const Vec3<T> &rhs) const
	{
		// sorting rank order: z, y, x.
		if(z != rhs.z) return (z < rhs.z);
		if(y != rhs.y) return (y < rhs.y);
		if(x != rhs.x) return (x < rhs.x);
		return false; // as this == rhs
	}

	bool operator > (const Vec3<T> &rhs) const
	{
		// sorting rank order: z, y, x.
		if(z != rhs.z) return (z > rhs.z);
		if(y != rhs.y) return (y > rhs.y);
		if(x != rhs.x) return (x > rhs.x);
		return false; //  as this == rhs
	}

	bool operator <= (const Vec3<T> &rhs) const
	{
		// sorting rank order: z, y, x.
		if(z != rhs.z) return (z < rhs.z);
		if(y != rhs.y) return (y < rhs.y);
		if(x != rhs.x) return (x < rhs.x);
		return true; // as this == rhs
	}

	bool operator >= (const Vec3<T> &rhs) const
	{
		// sorting rank order: z, y, x.
		if(z != rhs.z) return (z > rhs.z);
		if(y != rhs.y) return (y > rhs.y);
		if(x != rhs.x) return (x > rhs.x);
		return true; // as this == rhs
	}

	struct Vec3<T> operator + (const Vec3<T> &rhs) const
	{
		Vec3<T> gp = { x+rhs.x, y+rhs.y, z+rhs.z };
		return gp;
	}

};

typedef Vec3<int> GridPoint;


class GridFilter
{
	private:
	
		std::vector<GridPoint> neighbour_offsets;
		std::vector<int *> iptrs;

		double normalise( double in_x, double min_x, double Lx, int px ) const
		{
			double x = (in_x-min_x)/Lx;

			if( px == 1 )
			{
				if( x < 0.0 ) x += 1.0;
				else if( x >= 1.0 ) x -= 1.0;
			}

			return x;
		}

	public:

		std::map<GridPoint,int> m;

		GridFilter()
		{
			GridPoint g;
			for( int dz=-1; dz<=1; dz++ )
			{
				for( int dy=-1; dy<=1; dy++ )
				{
					for( int dx=-1; dx<=1; dx++ )
					{
						g.x = dx;
						g.y = dy;
						g.z = dz;
						neighbour_offsets.push_back( g );
					}
				}
			}
		}

		int Filter(
			const std::vector<double>& x_vec,
			const std::vector<double>& y_vec,
			const std::vector<double>& z_vec,
			const std::vector<int>& i_vec,
			double minx, double miny, double minz,
			double maxx, double maxy, double maxz,
			double grid_size,
			int px, int py, int pz,
			std::vector<int>& results
			)
		{
			if( y_vec.size() != x_vec.size() ) return -1;
			if( z_vec.size() != x_vec.size() ) return -1;

			size_t N = i_vec.size();

			double Lx = maxx - minx;
			double Ly = maxy - miny;
			double Lz = maxz - minz;

			int nx = (int) ceil( Lx/grid_size );
			int ny = (int) ceil( Ly/grid_size );
			int nz = (int) ceil( Lz/grid_size );

			m.clear();
			iptrs.resize(N);
			results.clear();

			//
			// Set up map
			//
			for( size_t ii=0; ii<N; ii++ )
			{
				int i = i_vec[ii];

				double x = normalise( x_vec[i], minx, Lx, px );
				double y = normalise( y_vec[i], miny, Ly, py );
				double z = normalise( z_vec[i], minz, Lz, pz );

				int ix = (int) floor( x * nx );
				int iy = (int) floor( y * ny );
				int iz = (int) floor( z * nz );

				GridPoint g;
				g.x = ix;
				g.y = iy;
				g.z = iz;

				if( m.find(g) == m.end() ) m[g] = 0; // assume bad, until we find otherwise!
				auto it = m.find( g );
				iptrs[ii] = &it->second;
			}

			//
			// Check each mapped grid point for neighbours. If none exist, set the value in the map to 0.
			//
			for( auto it = m.begin(); it != m.end(); it++ )
			{
				GridPoint g = it->first;

				//
				// -1 to +1 on each axis.
				//
				bool found = false;
				for( auto& o : neighbour_offsets )
				{
					GridPoint ng = g + o;

					if( ng == g ) continue;

					//
					// Wrap cells, if periodic.
					//
					if( px == 1 )
					{
						if( ng.x < 0 ) ng.x += nx;
						else if( ng.x >= nx ) ng.x -= nx;
					}
					if( py == 1 )
					{
						if( ng.y < 0 ) ng.y += ny;
						else if( ng.y >= ny ) ng.y -= ny;
					}
					if( pz == 1 )
					{
						if( ng.z < 0 ) ng.z += nz;
						else if( ng.z >= nz ) ng.z -= nz;
					}

					if( m.find(ng) != m.end() )
					{
						found = true;
						break;
					}

				}
				// If we found a neighbour, tag contents of this cell as ok!
				if( found == true ) it->second = 1;
			}

			//
			// Per-particle pointers should now be 1 if cell okay, 0 if not ok.
			//
			results.clear();
			for( size_t i=0; i<N; i++ )
			{
				if( *iptrs[i] == 1 ) results.push_back( i_vec[i] );
			}
			return (int)results.size();
		}
};


#endif
