#if !defined( STRINGUTIL_DEFINED )

#include <stdlib.h>
#include <string.h>

#include <string>
#include <vector>

/*
	Various handy string manipulation and conversion routines. These tend to be better than the atoi() / atof()
	etc routines, as it's very simple to know when the conversions have failed!

	Template methods are used wherever possible, so the compiler can deduce what the various data types should be
	on the basis of what parameters you pass in from your own code.
*/

//
// As a class, to avoid unused method warnings that we'd have with static "free" routines in a a namespace.
//
class StringUtil
{
	protected:
	
		//
		// Check if character 'test' is in test_characters (assumes "test_characters" null-terminated!).
		//
		static bool is_in( char test, const char * test_characters )
		{
			if( test_characters == NULL ) return false;

			while( *test_characters != '\0' )
			{
				if( test == *(test_characters++) ) return true;
			}
			return false;
		}

	public:

		//
		// Convert a character sequence into an integer type. Assumes default of base-10.
		// Template method, so should be compatible with all appropriate integer data types.
		//
		template <typename T> static int ToInteger( const char *str, T &result, int base = 10 )
		{
			char *endptr;

			if( str == NULL || base < 2 ) return -1;
	
			result = strtoll( str, &endptr, base );
			if( endptr == str || (*endptr != '\0' && *endptr != '\n' && *endptr != '\r') ) return -1;
			return 1;
		}
		template <typename T> static int ToInteger( const std::string &str, T &result, int base = 10 )
		{
			return ToInteger( str.c_str(), result, base );
		}

		//
		// Extract a set of integers from an appropriately formatted string, using "element_sep"
		// and "range_sep" strings. E.g. with element_sep = "," and range_sep = "-", the following
		// input:
		//
		//  12,13,14,20-22,90
		//
		// should APPEND the list of integers [12,13,14,20,21,22,90] to the "results" vector.
		//
		// NOTES:
		// 1. Ranges are INCLUSIVE of the specified start and end values!
		// 2. The "results" vector IS NOT CLEARED before new values appended (allows progressive generation).
		//
		template <typename T> static int ToIntegers( const char *str, std::vector<T> &results, const char *element_sep, const char *range_sep, int base = 10 )
		{
			std::vector< std::string > tokens, subtokens;
			int start, stop;

			Tokenize( str, tokens, element_sep );
			for( size_t i=0; i<tokens.size(); i++ )
			{
				if( StringUtil::Tokenize( tokens[i].c_str(), subtokens, range_sep ) == 2 )
				{
					for( size_t j=0; j<subtokens.size(); j++ )
					{
						if( StringUtil::ToInteger( subtokens[0].c_str(), start, base ) == -1 ) return -1;
						if( StringUtil::ToInteger( subtokens[1].c_str(), stop, base ) == -1 ) return -1;
						if( stop < start ) return -1;
					}
				}
				else
				{
					if( StringUtil::ToInteger( tokens[i].c_str(), start, base ) == -1 ) return -1;
					stop = start;
				}
				
				for( int value=start; value<=stop; value++ ) results.push_back( value );
			}
			return 1;
		}

		//
		// Convert a character sequence into a floating point type.
		// Template method, so should work fine with both "float" and "double".
		//
		template <typename T> static int ToReal( const char *str, T &result )
		{
			char *endptr;

			if( str == NULL ) return -1;

			result = strtod( str, &endptr );
			if( endptr == str || *endptr != '\0' ) return -1;
			return 1;
		}
		template <typename T> static int ToReal( const std::string &str, T &result )
		{
			return ToReal( str.c_str(), result );
		}

		//
		// This is SLOW, as we repeatedly append to std::string, and then copy it into "results". This approach
		// does avoid using any nasty fixed-size intermediate buffers etc, so it's hopefully safer.
		//
		static int Tokenize( const char * source, std::vector< std::string > &results, const char *delimiters )
		{
			std::string temp;
	
			results.clear(); // BEFORE error return test ...

			if( source == NULL || delimiters == NULL ) return -1;

			size_t src_len = strlen( source );
	
			for( size_t i=0; i<src_len; i++ )
			{
				if( is_in( source[i], delimiters ) == true )
				{
					if( temp.size() > 0 ) results.push_back( temp );
					temp.clear();
				}
				else
				{
					temp += source[i];
				}
			}
			// If there's a trailing token, add to the results vector
			if( temp.size() > 0  ) results.push_back( temp );
	
			return (int)results.size();
		}
		
		//
		// Strip defined whitespace from start and end of string. MODIFIES "str"!
		//
		static int Strip( char * str, const char *whitespace )
		{
			int str_length, start, end;

			if( str == NULL || whitespace == NULL ) return -1;

			str_length = strlen( str );

			// get index of first non-whitespace character in source string
			start = -1;
			for( int i=0; i<str_length && start == -1; i++ )
			{
				if( is_in( str[i], whitespace ) == false ) start = i;
			}

			// get index of last non-whitespace character in source string
			end = -1;
			for( int i=str_length-1; i>=0 && end == -1; i-- )
			{
				if( is_in( str[i], whitespace ) == false ) end = i;
			}

			// note that we could have an empty string, so check for that.
			if( start == -1 || end == -1 || end < start )
			{
				str[0] = '\0';
				return 0;
			}
			else
			{
				int j = 0;
				for( int i=start; i<=end; i++ ) str[j++] = str[i];
				str[j] = '\0';
				return j;
			}

			return 1;
		}

		//
		// Convert arbitrary POD into a character bitstring. Separates groups of "unit_size_bits"
		// with a space for easier reading.
		//
		// WARNING: assumes character buffer is large enough!
		// WARNING: assumes 8-but 'char'! (could change this to CHAR_BIT from <limits.h>)
		//
		template <typename T> static void ToBinary( T data, char *cbuf, int unit_size_bits = 8 )
		{
			T one = 1;
			int nbits = (int)sizeof(T)*8, upto = 0; // e.g. nbits = (int)sizeof(T) * CHAR_BITS
	
			if( cbuf == NULL ) return;

			for( int i=nbits-1; i>=0; i-- )
			{
				cbuf[upto++] = ( data & (one<<i) ) ? '1' : '0';
				if( i > 0 && i%unit_size_bits == 0 ) cbuf[upto++] = ' ';
			}
			cbuf[upto] = '\0';
		}
};


#define STRINGUTIL_DEFINED
#endif
