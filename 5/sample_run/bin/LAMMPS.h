#include <stdio.h>

#include <vector>
#include <map>

#include "StringUtil.h"

/*
	Simple handler for LAMMPS trajectory files. The TrajectoryFrame class is in a
	LAMMPS namespace, as this code is actually from another program I wrote.
*/

namespace LAMMPS
{

class TrajectoryFrame
{

	public:
	
		int timestep;
		
		double mins[3], maxs[3];
		int PBC[3];
		
		std::vector<double> x_vec,y_vec,z_vec;
		std::vector<int> atom_ids, atom_types; // per-atom
		
		void Clear()
		{
			timestep = -1;
			for( int i=0; i<3; i++ )
			{
				mins[i] = -0.5;
				maxs[i] = +0.5;
				PBC[i] = 1;
			}

			x_vec.clear();
			y_vec.clear();
			z_vec.clear();
			atom_ids.clear();
			atom_types.clear();
		}
		
		void Wrap( double &x, double &y, double &z ) const
		{
			double Lx = maxs[0]-mins[0];
			double Ly = maxs[1]-mins[1];
			double Lz = maxs[2]-mins[2];
			
			// Wrap deltas
			if( PBC[0] == 1 )
			{
				if(      x <  mins[0] ) x += Lx;
				else if( x >= maxs[0] ) x -= Lx;
			}

			if( PBC[1] == 1 )
			{
				if(      y <  mins[1] ) y += Ly;
				else if( y >= maxs[1] ) y -= Ly;
			}

			if( PBC[2] == 1 )
			{
				if(      z <  mins[2] ) z += Lz;
				else if( z >= maxs[2] ) z -= Lz;
			}
		}
		
		// Basic loading of trajectory frame
		int Load( FILE *f )
		{
			const char *delimiters = " \t\n\r";

			int maxbuf = 1023;
			char buffer[1024];
			
			int ntoks;
			std::vector< std::string > tokens;
			
			fpos_t file_position;

			std::map< std::string, int > token_indices;
			int number_of_atoms = 0; // internal check
			
			token_indices[ "id" ] = -1;
			token_indices[ "type" ] = -1;
			token_indices[ "x" ] = -1;
			token_indices[ "y" ] = -1;
			token_indices[ "z" ] = -1;

			if( f == NULL )
			{
				fprintf( stderr, "LAMMPSTrajectoryFrame::%s(): file pointer is NULL!\n", __func__ );
				return -1;
			}
			
			Clear();

			while( fgets( buffer, maxbuf, f ) != NULL )
			{
				if( (ntoks = StringUtil::Tokenize( buffer, tokens, delimiters )) < 1 || tokens[0] != "ITEM:" ) continue;

				//
				// Trajectory frame timestep number
				//
				if( tokens[1] == "TIMESTEP" )
				{
					if( fgets( buffer, maxbuf, f ) == NULL ) continue;
					if( StringUtil::ToInteger( buffer, timestep ) == -1 )
					{
						fprintf( stderr, "LAMMPSTrajectoryFrame::%s(): unable to convert timestep '%s' into an integer.\n", __func__, buffer );
						return -1;
					}
				}
				//
				// Simulation cell
				//
				else if( tokens[1] == "BOX" )
				{
					// expects ITEM: BOX BOUNDS xx yy zz
					if( tokens.size() != 6 )
					{
						fprintf( stderr, "LAMMPSTrajectoryFrame::%s(): bad box line: '%s'\n", __func__, buffer  );
						return -1;
					}
					
					for( int i=0; i<3; i++ ) PBC[i] = ( tokens[3+i][0] == 'p' ) ? 1 : 0;

					if( fgets( buffer, maxbuf, f ) == NULL ) continue;
					if( (ntoks = StringUtil::Tokenize( buffer, tokens, delimiters)) < 2 ) continue;

					//
					// Get x bounds
					//
					if( StringUtil::ToReal( tokens[0], mins[0] ) == -1 )
					{
						fprintf( stderr, "LAMMPSTrajectoryFrame::%s(): unable to convert min x value '%s' into a number.\n", __func__, tokens[0].c_str() );
						return -1;
					}
					if( StringUtil::ToReal( tokens[1], maxs[0] ) == -1 )
					{
						fprintf( stderr, "LAMMPSTrajectoryFrame::%s(): unable to convert max x value '%s' into a number.\n", __func__, tokens[1].c_str() );
						return -1;
					}

					//
					// Get y bounds
					//
					if( fgets( buffer, maxbuf, f ) == NULL ) continue;
					if( (ntoks = StringUtil::Tokenize( buffer, tokens, delimiters)) < 2 ) continue;

					if( StringUtil::ToReal( tokens[0], mins[1] ) == -1 )
					{
						fprintf( stderr, "LAMMPSTrajectoryFrame::%s(): unable to convert min y value '%s' into a number.\n", __func__, tokens[0].c_str() );
						return -1;
					}
					if( StringUtil::ToReal( tokens[1], maxs[1] ) == -1 )
					{
						fprintf( stderr, "LAMMPSTrajectoryFrame::%s(): unable to convert max y value '%s' into a number.\n", __func__, tokens[1].c_str() );
						return -1;
					}

					//
					// Get z bounds
					//
					if( fgets( buffer, maxbuf, f ) == NULL ) continue;
					if( (ntoks = StringUtil::Tokenize( buffer, tokens, delimiters)) < 2 ) continue;

					if( StringUtil::ToReal( tokens[0], mins[2] ) == -1 )
					{
						fprintf( stderr, "LAMMPSTrajectoryFrame::%s(): unable to convert min z value '%s' into a number.\n", __func__, tokens[0].c_str() );
						return -1;
					}
					if( StringUtil::ToReal( tokens[1], maxs[2] ) == -1 )
					{
						fprintf( stderr, "LAMMPSTrajectoryFrame::%s(): unable to convert max z value '%s' into a number.\n", __func__, tokens[1].c_str() );
						return -1;
					}
				}
				//
				// Number of atoms in frame
				//
				else if( tokens[1] == "NUMBER" && tokens[2] == "OF" && tokens[3] == "ATOMS" )
				{
					if( fgets( buffer, maxbuf, f ) == NULL ) continue;
					if( StringUtil::ToInteger( buffer, number_of_atoms ) == -1 )
					{
						fprintf( stderr, "LAMMPSTrajectoryFrame::%s(): unable to convert number of atoms '%s' into a number.\n", __func__, buffer );
						return -1;
					}
				}
				//
				// Atom coordinates etc
				//
				else if( tokens[1] == "ATOMS" )
				{
					bool xscale = false, yscale = false, zscale = false;

					std::map< std::string, int >::iterator it;
				
					//
					// Clear token indices
					//
					for( it=token_indices.begin(); it!=token_indices.end(); it++ )
					{
						it->second = -1;
					}

					//
					// Determine column indices for the appropriate atomic data
					//
					for( int i=0; i<ntoks; i++ )
					{
						std::string &tok = tokens[i];
						
						if( tok == "xu" ) { tok = "x"; }
						if( tok == "yu" ) { tok = "y"; }
						if( tok == "zu" ) { tok = "z"; }

						if( tok == "xs" ) { tok = "x"; xscale = true; }
						if( tok == "ys" ) { tok = "y"; yscale = true; }
						if( tok == "zs" ) { tok = "z"; zscale = true; }
						
						it = token_indices.find(tok);
						if( it == token_indices.end() ) continue;
						
						token_indices[tok] = i-2; // tokens start with "ITEM:", "ATOMS"
					}

					//
					// Check we've actually defined all the columns for the appropriate atom data
					//
					for( it=token_indices.begin(); it!=token_indices.end(); it++ )
					{
						if( it->second >= 0 ) continue;
						fprintf( stderr, "LAMMPSTrajectoryFrame::%s(): unable to find '%s' info in frame.\n", __func__, it->first.c_str() );
						return -1;
					}
					
					int id_i = token_indices["id"];
					int type_i = token_indices["type"];
					int x_i = token_indices["x"];
					int y_i = token_indices["y"];
					int z_i = token_indices["z"];

					//
					// Store file position, so we can rewind if we hit an ITEM line.
					//
					fgetpos( f, &file_position );

					//
					// Read atom data, using the columns we determined previously
					//
					while( fgets( buffer, maxbuf, f ) != NULL )
					{				
						ntoks = StringUtil::Tokenize( buffer, tokens, delimiters );
						if( ntoks < 1 ) continue;

						//
						// ITEM indicates we've hit the next frame; rewind th start of line and return.
						//
						if( tokens[0] == "ITEM:" )
						{
							fsetpos( f, &file_position );

							if( number_of_atoms != (int)atom_ids.size() )
							{
								fprintf( stderr, "LAMMPSTrajectoryFrame::%s(): WARNING: odd number of atoms (read %d, expected %d)\n", __func__, (int)atom_ids.size(), number_of_atoms );
							}

							return (int)atom_ids.size();
						}
						
						int atom_id, atom_type;
						double x, y, z;

						if( StringUtil::ToInteger( tokens[id_i], atom_id ) == -1 )
						{
							fprintf( stderr, "LAMMPSTrajectoryFrame::%s(): unable to convert id '%s' into a number.\n", __func__, tokens[id_i].c_str() );
							return -1;
						}
						if( StringUtil::ToInteger( tokens[type_i], atom_type ) == -1 )
						{
							fprintf( stderr, "LAMMPSTrajectoryFrame::%s(): unable to convert type '%s' into a number.\n", __func__, tokens[type_i].c_str() );
							return -1;
						}

						if( StringUtil::ToReal( tokens[x_i], x ) == -1 )
						{
							fprintf( stderr, "LAMMPSTrajectoryFrame::%s(): unable to convert x '%s' into a number.\n", __func__, tokens[x_i].c_str() );
							return -1;
						}
						if( StringUtil::ToReal( tokens[y_i], y ) == -1 )
						{
							fprintf( stderr, "LAMMPSTrajectoryFrame::%s(): unable to convert y '%s' into a number.\n", __func__, tokens[y_i].c_str() );
							return -1;
						}
						if( StringUtil::ToReal( tokens[z_i], z ) == -1 )
						{
							fprintf( stderr, "LAMMPSTrajectoryFrame::%s(): unable to convert z '%s' into a number.\n", __func__, tokens[z_i].c_str() );
							return -1;
						}

						atom_ids.push_back( atom_id );
						atom_types.push_back( atom_type );

						if( xscale == true ) x = mins[0] + x*(maxs[0]-mins[0]);
						if( yscale == true ) y = mins[1] + y*(maxs[1]-mins[1]);
						if( zscale == true ) z = mins[2] + z*(maxs[2]-mins[2]);

						x_vec.push_back( x );
						y_vec.push_back( y );
						z_vec.push_back( z );
					}
				}
			}

			return (int)atom_ids.size();
		}		

		//
		// Write xyz file fo the specific particles in the trajectory frame as indicated by "indices".
		//
		void WriteXYZ( FILE *f, const std::vector<int>& indices )
		{
			if( f == nullptr ) return;

			fprintf( f, "%d\n", (int)indices.size() );
			fprintf( f, "\n" );
			for( size_t i=0, max_i=indices.size(); i<max_i; i++ )
			{
				int index = indices[i];

				if( index >= (int)atom_ids.size() ) continue; // should probably be an error!

				double x = x_vec[index];
				double y = y_vec[index];
				double z = z_vec[index];
				fprintf( f, "%s %8.3f %8.3f %8.3f\n", "?", x, y, z );
			}	
		}
};

}