/**********************************************************************************
                                solve
                                
                            A. Bartholomew 1987, 1992
                            
The programme solve was originally written in ALGOL68 in 1987 as a part of my PhD 
research.  It was re-written in C++ in 1992.

In 2009 it was augmented with the ability to calculate fundamental solutions as well
as extreme fundamental solutions.

In 2010 it was further augmented to support the calculation of a basis for the kernel
to the matching system comprised of tracks, as required by the algorithm in [JSJ].  This
included the introduction of the debug subsystem into solve.

[JSJ] M. J. Dunwoody, Universal JSJ-Decompositions for Finitely Presented Groups,
                    Third Draft, March 2008
                                  
***********************************************************************************/
using namespace std;

#include <fstream>
#include <iostream>
#include <cstdlib>

#include <cstring>
#include <ctype.h>
#include <stdio.h>
#include <iomanip>
#include <ctime>
#include <list>
#include <vector>

extern ofstream debug;
ofstream	output;

#include <util.h>
#include <matrix.h>
//#include <bigint.h>
#include <decomp.h>
#include <n_limit.h>


int decomp_control::DEBUG = decomp_control::OFF;
bool decomp_control::DECOMPOSE_DEBUG = false;
bool decomp_control::SPLIT_DEBUG = false;
bool decomp_control::REDUCE_DEBUG = false;

#include <surd.h>

bool debug_setup(char* argv_ptr);
void check_debug_option_parameters(char* pptr, string option);
void debug_help();

unsigned int surd_control::DEBUG;
unsigned int bigint_control::DEBUG = 0;
unsigned int bigreal_control::DEBUG = 0;
bool rational_control::DEBUG = false;
int bigreal::div_num_places = 9;

matrix<int>* ext_fund_ptr;

bigint w_cuboid_volume_remaining;

bool matrix_control::SINGLE_LINE_OUTPUT = false;
bool		TWISTED_ALLOWED = false; // used by decompose
bool		FLAGS_ON_FIRST_2_TORSION_1_CELL = false; // used by decompose

bool	Greenberg_version = false;
bool 	TRACK_BASIS = true; // the programme default is now to calcualte a track basis for the matching system

/* used for debugging */
bool 		CALCULATE_EXTREME_RAYS = true;
bool		CALCULATE_FUNDAMENTAL_SOLUTIONS = true;
bool 		OPTIMIZED_FUNDAMENTAL_SOLUTION_SEARCH = true;
unsigned long	num_rank_calls;


/********************* Function prototypes ***********************/
void exit_option ();
int non_zero_count(matrix<int>& mat, int row);
int gcdn (const matrix<int>& mat, int row);
int rank (matrix<int>& matrixptr);
bool kbasis(matrix<int>& matrixref, matrix<int>& basis);
int greenberg (matrix<int>& cone_matrix_2, list<ray>& raylist, time_t start_time);
void usage_help();
void initial_help();
//void print(matrix<int> m, ostream& s, int n, string prefix="");
bool satisfies_matching_equations(matrix<int>& matching_system,vector<int> u);
int find_generator(vector<char> generators, char c);
bool split 
(
	vector<int>& 	pattern, 
	K_complex 		K, 
	bool 			short_split, 
	ofstream& 		txt, 
	bool 			return_components = false, 
	matrix<int>**	component_ptr = 0
);

bool check_thresholds(int n, vector<int>& u, n_limit& three_m_limit, list<threshold>::iterator *current_threshold);
//int update_n_limit (n_limit& limit, vector<int> u, n_limit& three_m_limit);
void update_n_limit (n_limit& limit, vector<int> u);
bigint non_excluded_region_size (n_limit& limit, int three_m);
bool reset_current_thresholds(n_limit& limit, n_limit& three_m_limit, vector<int>& u, list<threshold>::iterator *current_threshold);
bigint remaining_volume (n_limit& limit, vector<int> u);
bool next_combination(vector<int>& comb, int n);
void reduce_to_minimal_set (list<vector<int> >& solutions);
unsigned int decompose (vector<int>& track, K_complex& K, bool echo, bool full_decomp, ofstream& txt);
bool compatible(Slice_iter<int> t1,vector<int>& t2, K_complex& K);

/******************* Main Function ************************/
int main (int argc, char* argv[])
{

int			number_of_generators;
int			number_of_relators;

bool		skip_next_char;
bool		equal;
bool		switches = false;
bool    	input_file_provided = false;

int			m, g;
int			number_of_0_cells;
int			occ_count; 		// occurrence count
int			occ_of_cv; //occurrences of centre vertex
int 		number_of_1_cells;
int 		number_of_1_cells_in_rosette;
int			number_of_2_cells;
int			added_0_cell;
int			edge;
int			map_element;
int			being_considered;
int 		disc;
int			number_of_equations;
int			equation_pointer;
int			base_occurrence;
int			working_occurrence;
int			first_positive;
int			second_positive;
int			first_negative;
int			second_negative;
int			repeat_count;
int			column;

string version="23-5-10";
string jobname;
string input_file;

ifstream	input;

/* Determine command line switches and whether input is to 
   come from the screen or a file */

if ( argc > 1 )
{

    if ( *(argv[1]) == '-')
		switches = true;

    if (switches)
    {
		if (strchr(argv[1], 'e'))
			TRACK_BASIS = false; // we'll calculate extreme track instead

		if(strchr(argv[1], 'h') || strchr(argv[1], 'H'))
		{
			if (strchr(argv[1],'#') && strchr(argv[1], 'H') && strchr(argv[1], '!'))
			{
				debug_help();
			}
			else
			{
				usage_help();
				exit(0);
			}
		}

		if (strchr(argv[1], 'g'))
			Greenberg_version = true;

		char* dptr = strchr(argv[1], '#');
		if (dptr != 0)
		{
			debug.open ("solve.dbg"); 

			if (!debug)
			{
				cout << "\nError opening debug file\n";
				exit(0);
			}
			else
				debug << "\ndebug from solve version " << version << endl;

			if (!debug_setup(argv[1]))
			{
				decomp_control::DEBUG = decomp_control::BASIC;
//					decomp_control::DEBUG = true;
				debug << "main: default debug options set" << endl;
			}
			
		}

		if (strchr(argv[1], 'E'))
		{
			CALCULATE_EXTREME_RAYS = false;
			CALCULATE_FUNDAMENTAL_SOLUTIONS = false;
		}

		if (strchr(argv[1], 'F'))
			CALCULATE_FUNDAMENTAL_SOLUTIONS = false;	

		if (strchr(argv[1], 'O'))
			OPTIMIZED_FUNDAMENTAL_SOLUTION_SEARCH = false;

	}

	if ( !switches && argc > 3)
	{
		usage_help();
	}
	if (!switches && argc >= 2)
	{
		input_file = argv[1];
		input_file_provided = true;
	}
	else if (switches && argc >= 3)
	{
		input_file = argv[2];
		input_file_provided = true;
	}
}


if (input_file_provided)
{
	input.open(input_file.c_str());
	if(!input)
   	{
		cout << "\nError opening input file\n";
		exit(0);
   	}	

if (decomp_control::DEBUG >= decomp_control::BASIC)
	debug << "\nmain: input file: " << input_file << "\n";
	
	input >> jobname;	
   	input >> number_of_generators;

if (decomp_control::DEBUG >= decomp_control::BASIC)
	debug << "main: number of generators = " << number_of_generators << endl;

   	input >> number_of_relators;
}
else
{
    cout << "\n                Solve version " << version;
	cout << "\n\nThis program evaluates either: " << endl;
	cout << "\n1. A generating set of tracks in the triangular 2-complex associated with a" << endl;
	cout << "   given finitely presented group, or " << endl;
	cout << "\n2. A basis for the solution space to the matching equations of the triangular" << endl;
	cout << "   2-complex. comprised of a set of tracks." << endl;
	
	cout << "\nPlease enter a jobname for the group you want to analyse, ";
	cout << "or --help for help.\n";
	cout << "\nThe jobname must be a valid file name for your machine";
	cout << "\nJobname ? ";

	cin >> jobname;

	if (jobname.find("--help") != string::npos)
	{
		initial_help();
		exit(0);
	}
	

	cout << "\nPlease enter group presentation:\n";

   	cout << "\nhow many generators? ";
   	cin >> number_of_generators;

if (decomp_control::DEBUG >= decomp_control::BASIC)
	debug << "main: number of generators = " << number_of_generators << endl;

   	cout << "\nhow many relators? ";
   	cin >> number_of_relators;

if (decomp_control::DEBUG >= decomp_control::BASIC)
	debug << "main: number of relators = " << number_of_relators << endl;
}

vector<char>	generator(number_of_generators);
vector<char>	generators_of_order_2(number_of_generators,'#'); 
vector<string> 	relator(number_of_relators);
vector<char> 	generator_corresponding_to_1_cell(2*number_of_generators);
vector<int>		occ_of_order_2_gen(number_of_generators);
vector<int>		length_of_attaching_map(number_of_relators);

if (input_file_provided)
{
	if (decomp_control::DEBUG >= decomp_control::BASIC)
		debug << "main: number of relators = " << number_of_relators << endl;

	for (int i = 0; i < number_of_generators; i++)
    	input >> generator[i];
		
if (decomp_control::DEBUG >= decomp_control::BASIC)
{
	debug << "main: generators = ";
	for (int i = 0; i< number_of_generators; i++)
		debug << generator[i] << ' ';

	debug << endl;
}

	for (int i=0; i < number_of_relators; i++)
    	input >> relator[i];

if (decomp_control::DEBUG >= decomp_control::BASIC)
{
	debug << "main: relators = \n";
	for (int i = 0; i< number_of_relators; i++)
		debug << "main:   " << relator[i] << endl;
}

}
else
{
	cout << "\nEnter the generators to be used, spaces will be ignored.  ";
	cout << "The generators must\nbe distinct single characters, inverses ";
	cout << "not being allowed.";
	cout << "\n\nGenerators? ";

	for (int i = 0; i < number_of_generators; i++)
    	cin >> generator[i];

if (decomp_control::DEBUG >= decomp_control::BASIC)
{
	debug << "main: generators = ";
	for (int i = 0; i< number_of_generators; i++)
		debug << generator[i] << ' ';

	debug << endl;
}


	cout << "\nEnter the relators of the group in the generators '";
	cout << " " << generator[0];

	for (int i=1; i < number_of_generators; i++)
    	cout << "  " << generator[i];

	cout << " '.\n\n";

	cout << "The relators must not contain any spaces but may be separated by ";
	cout << "an arbitrary\nnumber of spaces, or a newline.\n";
	cout << "\nRepresent inverses of generators by preceeding the character ";
	cout <<	 "with '-' eg -a,\nand powers by repetition\n\n";

	cout << "\nWe also require that:\n";

	cout << "1.  relators indicating a generator of finite order should be ";
	cout << "entered in the\n    positive, eg 'aa', not '-a-a';\n";
	cout << "2.  the relators have length > 3 or take one of the forms 'aa'";
	cout << " or 'aaa';\n";
	cout << "3.  the total, absolute, exponent of each generator is at least";
	cout << " two;\n";

	cout << "\nRelators? ";

	for (int i=0; i < number_of_relators; i++)
    	cin >> relator[i];

if (decomp_control::DEBUG >= decomp_control::BASIC)
{
	debug << "main: relators = \n";
	for (int i = 0; i< number_of_relators; i++)
		debug << "main:   " << relator[i] << endl;
}

}

/* Identify and count the generators of order 2 and number of 0-cells,
   rearrange the generators and relators so that any 2-torsion generators 
   and relators appear first
*/
if (decomp_control::DEBUG >= decomp_control::BASIC)
	debug << "main: checking generators and relators for re-ordering" << endl;

vector<char> re_ordered_generators(number_of_generators);
vector<string> re_ordered_relators(number_of_relators);
number_of_0_cells = 1;

for (int i=0; i< number_of_relators; i++)
{
	/* put in the 2-torsion generators and relators first */
	if(relator[i].length() == 2)
	{

if (decomp_control::DEBUG >= decomp_control::BASIC)
	debug << "main:   relator " << i+1 << " " << relator[i] << " has order 2" << endl;

		generators_of_order_2[number_of_0_cells-1] = relator[i][0];
		re_ordered_generators[number_of_0_cells-1] = relator[i][0];
		re_ordered_relators[number_of_0_cells-1] = relator[i];
		number_of_0_cells++;
	}
}

/*  now fill in the rest of the relators */
int place = number_of_0_cells-1;
for (int i=0; i< number_of_relators; i++)
{
	if(relator[i].length() != 2)
	{
		re_ordered_relators[place] = relator[i];
		place++;
	}
}

place = number_of_0_cells-1;
for (int i=0; i< number_of_generators; i++)
{
	/* now the rest of the generators */
	if(find_generator(generators_of_order_2,generator[i]) == NOT_FOUND)
	{
		re_ordered_generators[place] = generator[i];
		place++;
	}
}	

generator = re_ordered_generators;
relator = re_ordered_relators;

if (decomp_control::DEBUG >= decomp_control::BASIC)
{
	debug << "main: re-ordered generators = ";
	for (int i = 0; i< number_of_generators; i++)
		debug << generator[i] << ' ';

	debug << endl;

	debug << "main: re-ordered relators = \n";
	for (int i = 0; i< number_of_relators; i++)
		debug << "main:   " << relator[i] << endl;
}

/* Get time at start of program, after group has been entered. */
time_t start_time = time(0);

/***** Count the number of occurrences in the relators of each order 2 generator **********/

for (int i=0; i < number_of_0_cells-1; i++)
{
    occ_count = 0;
    for (int j = 0; j < number_of_relators; j++ )
    {
		for ( unsigned int k=0; k< relator[j].length(); k++)
		{
		    if ( relator[j][k] == generators_of_order_2[i])
				occ_count += 1;
		}
    }
    occ_of_order_2_gen[i] = occ_count;
}

if (decomp_control::DEBUG >= decomp_control::BASIC)
{
	debug << "main: occurrences of order two generators in relators: ";
	for (int i = 0; i< number_of_0_cells -1; i++)
		debug << occ_of_order_2_gen[i] << ' ';
	debug << endl;
}

/***** Calculate the length of the various attaching maps. **********/

for (int i=0 ; i < number_of_relators ; i++)
    length_of_attaching_map[i] = 0;

for (int i=0 ; i < number_of_relators ; i++)
{
    for (unsigned int j=0; j< relator[i].length(); j++)
	{
		if( relator[i][j] != '-')
		{
		    if (find_generator(generators_of_order_2, relator[i][j]) != NOT_FOUND)
				length_of_attaching_map[i] += 2;
		    else
				length_of_attaching_map[i] += 1;
		}
	}
}

if (decomp_control::DEBUG >= decomp_control::BASIC)
{
	debug << "main: lengths of attaching maps = ";
	for (int i = 0; i< number_of_relators; i++)
		debug << length_of_attaching_map[i] << ' ';	
	debug << endl;
}

/************************************************************************
We now calculate how many occurrences of each 0-cell there are in the
boundary of the relator discs, again, for future use. We already have this
information for the added 0-cells, in 'occ of order 2 gen', so must
calculate the number corresponding to the centre of the rosette.
*************************************************************************/
occ_of_cv = 0;
for ( int i=0; i < number_of_relators ; i++)
    occ_of_cv += length_of_attaching_map[i];

for ( int i=0; i < number_of_0_cells - 1; i++)
    occ_of_cv -= occ_of_order_2_gen[i];

if (decomp_control::DEBUG >= decomp_control::BASIC)
	debug << "main: occurrences of central vertex = " << occ_of_cv << endl;


/******** Calculate the number of 1-cells and 2-cells. ************/

number_of_1_cells = number_of_generators + number_of_0_cells - 1;
number_of_1_cells_in_rosette = number_of_1_cells;
/* '...in rosette' will be used later. */

for ( int i=0 ; i < number_of_relators ; i++)
    number_of_1_cells += length_of_attaching_map[i] - 3;

number_of_2_cells = 0;
for (int i=0; i< number_of_relators; i++)
    number_of_2_cells += length_of_attaching_map[i] - 2;

if (decomp_control::DEBUG >= decomp_control::BASIC)
{
	debug << "main: number of 0-cells = " << number_of_0_cells << endl;
	debug << "main: number of 1-cells = " << number_of_1_cells << endl;
	debug << "main: number of 1-cells in rosette = " << number_of_1_cells_in_rosette << endl;
	debug << "main: number of 2-cells = " << number_of_2_cells << endl;
}

/*************************************************************************
We consider the 0-cells to be labelled 1,2,... according to the rule that
the 0-cell '1' is the 0-cell at the centre of the rosette corresponding to
the group presentation and successively label the added 0-cells 2,... in the
order determined by 'generators_of_order_2', each added 0-cell being
constructed on the corresponding leaf of the rosette.  By moving along the
leaf determined by the first generator we begin to list the 1-cells of the
complex in the order determined by 'generator', moving round each leaf in a
clockwise direction.  Each 1-cell will be recorded as a pair of integers
(m,n) forming a row of a matrix, where m is its initial 0-cell and n is its
terminal 0-cell.  We also record, as a pair of integers, which 1-cells
correspond to each generator (in the case that the generator has order 2 we
store the pair (i,i+1) and otherwise the pair (0,i) for the relevant
value i) and which generator "carries" each 1-cell in the rosette.
*************************************************************************/
matrix<int> one_cell(number_of_1_cells, 2);
matrix<int> one_cells_corresponding_to_generator(number_of_generators, 2);

/*The first added 0-cell will be the 0-cell '2', the first 1-cell ('edge') to
be considered will be 'one_cell[1]'.*/

added_0_cell = 2;
edge = 1;

for( int i=0 ;i< number_of_generators; i++)
{
//    if ( strchr(generators_of_order_2, generator[i] ) )
    if ( find_generator(generators_of_order_2, generator[i]) != NOT_FOUND)
    {
		one_cells_corresponding_to_generator[i][0] = edge;
		one_cells_corresponding_to_generator[i][1] = edge+1;
		
		one_cell[edge-1][0] = 1;
		one_cell[edge-1][1] = added_0_cell;
		
		generator_corresponding_to_1_cell[edge-1] = generator[i];
		edge += 1;

		one_cell[edge-1][0] = added_0_cell;
		one_cell[edge-1][1] = 1;

		generator_corresponding_to_1_cell[edge-1] = generator[i];
		edge += 1;

		added_0_cell += 1;
    }
    else
    {
		one_cells_corresponding_to_generator[i][0] = 0;
		one_cells_corresponding_to_generator[i][1] = edge;

		one_cell[edge-1][0] = 1;
		one_cell[edge-1][1] = 1;

		generator_corresponding_to_1_cell[edge-1] = generator[i];
		edge += 1;
    }
}

if (decomp_control::DEBUG >= decomp_control::BASIC)
{
	debug << "main: one_cells_corresponding_to_generator\n";
	print (one_cells_corresponding_to_generator,debug,0,"main:   ");
}


int interior_1_cell = edge;

/**********************************************************************
We are now in a position to evaluate the attaching maps corresponding to
each relator.  The attaching maps will be recorded as a row of signed
integers corresponding to oriented 1-cells.  The row is obtained by moving
clockwise round the attaching disc determined by each relator (with added
0-cells included) from O(D).  The attaching maps will be stored in a single
matrix, large enough to hold the longest such map.
*********************************************************************/

int max_att_map_len = length_of_attaching_map[0];

for ( int i = 1; i < number_of_relators; i++)
{
    if (length_of_attaching_map[i] > max_att_map_len)
		max_att_map_len = length_of_attaching_map[i];
}

matrix<int> attaching_map (number_of_relators,max_att_map_len);
skip_next_char = false;

for( int i=0; i< number_of_relators; i++)
{
	const char* Relator = relator[i].c_str();
	
    map_element=0;

    for (unsigned int j=0; j < strlen(Relator); j++)
    {
		if (skip_next_char)
		    skip_next_char = false;
		     /*'skip_next_char=true' iff 'Relator[j-1]="-"' in which case
	    	  'Relator[j]' has been dealt with. */
//		else if ((Relator[j] == '-') && strchr(generators_of_order_2, Relator[j+1]))
		else if ((Relator[j] == '-') && find_generator(generators_of_order_2, Relator[j+1])!= NOT_FOUND)
		{
	    	/****** Find Relator[j+1] in generator. ************/
	    	for ( int k = 0; k < number_of_generators; k++)
			if (generator[k] == Relator[j+1])
			{
		    	being_considered = k;
		    	break;
			}

	    	attaching_map[i][map_element] = -1 * one_cells_corresponding_to_generator[being_considered][1];
	    	map_element += 1;
	    	attaching_map[i][map_element] = -1 * one_cells_corresponding_to_generator[being_considered][0];
	    	map_element += 1;
	    	skip_next_char = true;
		}
		else if ( Relator[j] == '-' )
		{
	    	/****** Find relator[i][j+1] in generator. ************/
	    	for ( int k = 0; k < number_of_generators; k++)
			if (generator[k] == Relator[j+1])
			{
		    	being_considered = k;
		    	break;
			}

	    	attaching_map[i][map_element] = -1 * one_cells_corresponding_to_generator[being_considered][1];
	    	map_element += 1;
	    	skip_next_char = true;
    	}
		else if (find_generator(generators_of_order_2, Relator[j])!=NOT_FOUND)
		{
	    	/****** Find relator[i][j+1] in generator. ************/
	    	for ( int k = 0; k < number_of_generators; k++)
			if (generator[k] == Relator[j])
			{
		    	being_considered = k;
		    	break;
			}

	    	attaching_map[i][map_element] = one_cells_corresponding_to_generator[being_considered][0];
	    	map_element += 1;
	    	attaching_map[i][map_element] = one_cells_corresponding_to_generator[being_considered][1];
	    	map_element += 1;
		}
		else
		{
	    	/****** Find relator[i][j+1] in generator. ************/
	    	for ( int k = 0; k < number_of_generators; k++)
			if (generator[k] == Relator[j])
			{
		    	being_considered = k;
		    	break;
			}

	    	attaching_map[i][map_element] =	one_cells_corresponding_to_generator[being_considered][1];
	    	map_element += 1;
		}
    }
}

if (decomp_control::DEBUG >= decomp_control::BASIC)
{
	debug << "main:  attaching_map:" << endl;
	print(attaching_map, debug, 0, "main:   ");
}


/***********************************************************************
We now complete the list of 1-cells, following the canonical order.  For each
relator in turn we move clockwise round the attaching map, beginning at the
first edge determined by 'attaching map[i]'.
***********************************************************************/
for (int i=0; i < number_of_relators; i++)
{
    for (int j=2; j < length_of_attaching_map[i]-1; j++)
    {
		if ( attaching_map[i][j] < 0 )
		{
			 one_cell[edge-1][0] = 1;
			 one_cell[edge-1][1] = one_cell[abs(attaching_map[i][j])-1][1];
			 edge += 1;
		}
		else
		{
			 one_cell[edge-1][0] = 1;
			 one_cell[edge-1][1] = one_cell[attaching_map[i][j]-1][0];
			 edge += 1;
		}
    }
}

if (decomp_control::DEBUG >= decomp_control::BASIC)
{
    debug << "main: one_cell:" << endl;
    print(one_cell,debug,0,"main:   ");
}
    

/*******************************************************************
We now list the 2-cells of the complex in the order determined by
'attaching map[i]' for i=1,...,'number of relators'.  Each 2-cell will be
recorded as a triple of signed integers forming a row of a matrix,
corresponding to the 1-cells forming its boundary (moving clockwise from
O(D) within each 2-cell.  Notice that 'interior 1 cell' is correctly
assigned for this.
************************************************************************/
matrix<int> two_cell(number_of_2_cells,3);
disc = 0;
for ( int i = 0; i < number_of_relators ; i++)
{
    if ( length_of_attaching_map[i] == 3 )
    {
		int j = attaching_map[i][1];
		two_cell[disc][0] = j;
		two_cell[disc][1] = j;
		two_cell[disc][2] = j;
		disc += 1;
    }
    else
    {
		two_cell[disc][0] = attaching_map[i][0];
		two_cell[disc][1] = attaching_map[i][1];
		two_cell[disc][2] = -interior_1_cell;
		disc += 1;
		interior_1_cell += 1;

		for (int j=2; j < length_of_attaching_map[i]-2; j++)
		{
			two_cell[disc][0] = interior_1_cell-1;
			two_cell[disc][1] = attaching_map[i][j];
			two_cell[disc][2] = -interior_1_cell;
			disc += 1;
		    interior_1_cell += 1;
		}

		two_cell[disc][0] = interior_1_cell-1;
		two_cell[disc][1] = attaching_map[i][length_of_attaching_map[i]-2];
		two_cell[disc][2] = attaching_map[i][length_of_attaching_map[i]-1];
		disc += 1;
    }
}

if (decomp_control::DEBUG >= decomp_control::BASIC)
{
    debug << "main: two_cell:" << endl;
    print(two_cell,debug,0,"main:   ");
}    

/***********************************************************************
We now calculate how many equations the matching system for L will contain.
Initially we would have one equation for each triple (p,q,r), where r is a
1-cell of L and p, q are 1-cells in the boundary of the 2-cells P, Q (not
necessarily distinct) such that p and q become identified in L, with p=q=r;
however, this produces redundant information.  By considering 'two cell',
we see that, for each 1-cell i, if i occurs (with a sign) n times in
'two-cell' then we shall need n-1 equations to describe the identifications
at i.  Thus the total number of equations is
3 * number of 2-cells - number of 1-cells.
************************************************************************/
number_of_equations = 3 * number_of_2_cells - number_of_1_cells;

/**********************************************************************
At this stage we calculate the number of occurrences of each 1-cell in
'two cell'.  This information will be used in evaluating decompositions
of G and splitting patterns in the program SEARCH.
**********************************************************************/

int occurrences_of_1_cell[number_of_1_cells];

for( int i = 0; i< number_of_1_cells ; i++)
    occurrences_of_1_cell[i] = 0;

/****************** Original way - 1987 ******************
for (i = 0; i < number_of_2_cells; i++)
    for ( j =0; j < 3; j++)
	*occurrences_of_1_cell.elt(abs(*two_cell.elt(i,j))-1) += 1;
***********************************************************/

/**************** New way  - 1992 *************************/
for (int i=0; i < number_of_relators; i++)
{
    for (int j=0; j < length_of_attaching_map[i]; j++)
		occurrences_of_1_cell[abs(attaching_map[i][j])-1] += 1;
}

for ( int i = number_of_1_cells_in_rosette; i < number_of_1_cells; i++)
    occurrences_of_1_cell[i] = 2;

if (decomp_control::DEBUG >= decomp_control::BASIC)
{
	debug << "main: occurrences of 1 cells: ";
	for (int i=0;i<number_of_1_cells;i++)
    	debug << occurrences_of_1_cell[i] << " ";
	debug << endl;
}

/*************************************************************************
To evaluate the matrix governing the matching system we assign variables to
each 2-cell.  The 2-cells will assigned variables in the order determined by
'two cell', each 2-cell having variables xi,xj,xk (j=i+1,k=j+1) assigned by
the rule that, if the 2-cell has boundary (u,v,w) (as stored in 'two cell'),
xi is assigned to the 'vertex' where u and w meet, xj where u and v meet
and xk where v and w meet.  Thus if we are considering a 1-cell u for which
an occurrence is found as the nth element of two_cell, regarded as a row of
integers (counted from zero), then the corresponding variables are n+1 and n+2, 
if n%3=0 or 1 and n+1 and n-1 otherwise.  This is because the canonical ordering 
of variables is respected by two_cell.  The matrix we construct has the non zero 
entries 1 and -1 only, each row containing nominally four non zero entries, two
positive and two negative.  However, it is possible in the case of a repeated
edge that a positive and negative entry coincide.
************************************************************************/

/* Declare 'matching system'. */
m = number_of_equations;
int n = 3 * number_of_2_cells;
matrix<int> matching_system (m,n);

/* matrices are initialized to zero
for (i = 0 ; i < m; i++)
    for (j =0; j < n; j++)
	*matching_system.elt(i,j) = 0;
*/

equation_pointer = 0;

for ( int i = 0 ; i < number_of_1_cells ; i++ )
{
	bool base_occurrence_found = false;
	
	for (int j=0; j< number_of_2_cells; j++)
	{
		for (int k=0; k< 3; k++)
		{
			if (!base_occurrence_found && abs(two_cell[j][k]) == i+1 )
			{
				base_occurrence = 3*j + k;

if (decomp_control::DEBUG >= decomp_control::BASIC)
{
	debug << "main: base occurrence for one cell " << i+1 << " = " << base_occurrence << endl;
	debug << "main: next occurrences at: ";
}

			    if (base_occurrence % 3 == 0 || base_occurrence % 3 == 1)
			    {
					first_positive = base_occurrence;
					second_positive = base_occurrence + 1;
			    }
			    else
			    {
					first_positive = base_occurrence;
					second_positive = base_occurrence - 2;
			    }
				base_occurrence_found = true;
			}
			else if (base_occurrence_found && abs(two_cell[j][k]) == i+1 )
			{
				working_occurrence = 3*j + k;
if (decomp_control::DEBUG >= decomp_control::BASIC)
	debug << working_occurrence << ' ';

			    if (working_occurrence % 3 == 0 || working_occurrence % 3 == 1)
			    {
					first_negative = working_occurrence;
					second_negative = working_occurrence + 1;
			    }
			    else
			    {
					first_negative = working_occurrence;
					second_negative = working_occurrence - 2;
			    }

			    matching_system[equation_pointer][first_positive] += 1;
			    matching_system[equation_pointer][second_positive] += 1;
			    matching_system[equation_pointer][first_negative] += -1;
			    matching_system[equation_pointer][second_negative] += -1;

			    equation_pointer += 1;

			}
		}
	}
if (decomp_control::DEBUG >= decomp_control::BASIC)
	debug << endl;
}

if (decomp_control::DEBUG >= decomp_control::BASIC)
{
	debug << "main: matching_system: " << m << " x " << n << endl;
	print (matching_system,debug,3, "main: ");
}

/* Declare the cone matrix */
int kdim = n - rank(matching_system);

if (decomp_control::DEBUG >= decomp_control::BASIC)
{
	debug << "main: rank = " << n-kdim << " kernel dimension = " << kdim << endl;
}

matrix<int> cone_matrix(n,kdim);

/* Evaluate basis for kernel of the matching system - ie: cone_matrix */
kbasis(matching_system, cone_matrix);

if (decomp_control::DEBUG >= decomp_control::BASIC)
{
	debug << "main: cone_martrix: " << n << " x " << kdim << endl;
	print (cone_matrix,debug,3, "main: ");
}

/********** Create file "jobname.lis" for transcript of results **********/

//string textfile = jobname + ".out";
string textfile = jobname + "-solve-output.txt";

output.open (textfile.c_str());

if (!output)
{
    cout << "\nError opening " << textfile;
    exit(0);
}

output << "\nGroup presentation:\n\n<";

for (int i=0; i < number_of_generators; i++)
    output << "  " << generator[i];

output << " :";

for (int i=0; i < number_of_relators; i++)
   output << " " << relator[i];

output << " >\n\nJobname: " << jobname << "\n";

output << "\nTriangular 2-complex comprises: " << number_of_0_cells
      << " 0-cells, " << number_of_1_cells << " 1-cells and "
      << number_of_2_cells << " 2-cells."
      << "\n\nMatching system (size " << number_of_equations << " x "
      << 3*number_of_2_cells << ")\n\n";

print(matching_system, output, 2,"");

output << "\n\nCone matrix (size " << n << " x " << kdim << ")\n";

print (cone_matrix, output, 4,"");
output << endl;

/***********************************************************************
Create "jobname.eft" to hold the information needed by 'search'.

NOTE:  the format of this file has changed from the original 1987 version.
       all information will be stored in this file, instead of having
       a .eft file and a .bin file.
	   
	   In the 1987 and 1992 versions this file was in binary format.  
	   Now we have switched to ascii format, since there is no need
	   to have it any other way.

Format of .eft file is

	number_of_generators
	generators
	number_of_relators
	length of relators
	relators
	number_of_0_cells
	number_of_1_cells
	number_of_2_cells
	occ_of_cv
	number_of_1_cells_in_rosette
	length_of_attaching_maps
	attaching_map
	generator_corresponding_to_1_cell
	number of generators of order 2
    [generators_of_order_2] // if any exist
    [occ_of_order_2_gen] // if any exist
	occurrences_of_1_cell
	two_cell
	m // rows in matching_system
	n // columns in matching_system
	matching_system
	s // number of extreme fundamental tracks
	ext_fund_tracks

************************************************************************/

string eftfile = jobname + ".eft";
ofstream eft(eftfile.c_str());

if (!eft)
{
    cout << "\nError opening " << eftfile;
    exit(0);
}

eft << number_of_generators << endl;

for (int i=0; i < number_of_generators; i++)
    eft << generator[i] << ' ';

eft << endl;

eft << number_of_relators << endl;

for (int i=0; i< number_of_relators; i++)
    eft << relator[i] << endl;

eft << number_of_0_cells << endl;
eft << number_of_1_cells << endl; 
eft << number_of_2_cells << endl;
eft << occ_of_cv << endl;
eft << number_of_1_cells_in_rosette << endl; 

for (int i=0; i<number_of_relators; i++)
   eft << length_of_attaching_map[i] << ' ';

eft << attaching_map;

for (int i=0; i<number_of_1_cells_in_rosette; i++)
    eft << generator_corresponding_to_1_cell[i] << ' ';
eft << endl;

int number_of_order_2_generators = number_of_0_cells-1;
eft << number_of_order_2_generators << endl;

if (number_of_order_2_generators > 0)
{
    for (int i=0; i< number_of_order_2_generators; i++)
		eft << generators_of_order_2[i] << ' ';
	eft << endl;

    for (int i=0; i<number_of_order_2_generators; i++)
		eft << occ_of_order_2_gen[i] << ' ';
	eft << endl;
}

for (int i=0; i< number_of_1_cells; i++)
   eft << occurrences_of_1_cell[i] << ' ';

eft << two_cell;

eft << m << ' ' << n;

eft << matching_system;


if (TRACK_BASIS)
{	
	matrix<int> pattern_basis = cone_matrix;
	
	/* exchange the first basis element for u=(1,...,1) */
	for (int i=0; i< n; i++)
		pattern_basis[i][0] = 1;
		
	/* add enough copies of u to every other basis element to produce
	   a basis in the positive quadrant of n-dimensional space
	*/
	for (int i=1;i< kdim; i++)
	{
		int max_negative = 0;
		for (int j=0; j< n; j++)
		{
			if (pattern_basis[j][i] < max_negative)
				max_negative = pattern_basis[j][i];
		}
if (decomp_control::DEBUG >= decomp_control::BASIC)
	debug << "main: basis element " << i << ", max_negative = " << max_negative << endl;
	
	max_negative *= -1;
	
	for (int j=0; j< n; j++)
		pattern_basis[j][i] += max_negative;
	}

	output << "\n\nPattern basis (size " << n << " x " << kdim << ")\n";
	print (pattern_basis, output, 4,"");
	output << endl;

if (decomp_control::DEBUG >= decomp_control::BASIC)
{
	debug << "main: pattern_basis: " << n << " x " << kdim << endl;
	print (pattern_basis,debug,3, "main: ");
}
	
	/* split the patterns in the pattern_basis and select a maximal linearly independent 
	   set of the compeonent tracks.  We only need look for kdim linearly independent tracks
	   since if there were more then the dimension of the kernel would not be kdim.
	*/
	
	
	/***** Prepare a K_complex structure to hold simplicial complex parameters *****/
	K_complex K;
	K.number_of_0_cells = number_of_0_cells;
	K.number_of_1_cells = number_of_1_cells;
	K.number_of_2_cells = number_of_2_cells;
	K.occ_of_cv = occ_of_cv;
	K.generator_corresponding_to_1_cell = generator_corresponding_to_1_cell;
	K.number_of_generators = number_of_generators;
	K.number_of_relators = number_of_relators;
	K.relator = relator;
	K.number_of_1_cells_in_rosette = number_of_1_cells_in_rosette;
	K.generators_of_order_2 = generators_of_order_2;
	K.occ_of_order_2_gen = occ_of_order_2_gen;
	K.occurrences_of_1_cell = occurrences_of_1_cell;
	K.length_of_attaching_map = length_of_attaching_map;
	K.two_cell = &two_cell;
	K.attaching_map = &attaching_map;
	K.matching_system = &matching_system;
	K.s = 0; // There are no extreme fundamental tracks in this case

	//vector<int> w(n);
	
	list<matrix<int> > component_list;
	int num_components = 0;

if (decomp_control::DEBUG >= decomp_control::INTERMEDIATE)
	debug << "main: splitting pattern_basis" << endl;

	for (int i=0; i< kdim; i++)
	{
		matrix<int>* component_ptr;
		vector<int> basis_elt(n);
		
		for (int j=0; j< n; j++)
			basis_elt[j] = pattern_basis[j][i]; // basis elements are in the columns of pattern_basis

if (decomp_control::DEBUG >= decomp_control::INTERMEDIATE)
{
	debug << "main:   basis_elt: ";
	for (int j=0; j<n; j++)
		debug << basis_elt[j] << ' ';
	debug << endl;
}			

		/* split the pattern into components tracks */
		split(basis_elt, K, false,output, true, &component_ptr);
	
if (decomp_control::DEBUG >= decomp_control::INTERMEDIATE)
{
	debug << "main:   components: " << endl;
	print (*component_ptr,debug,3, "main:   ");
}			

		component_list.push_back(*component_ptr);
			
		/* the components are returned as the rows of *component_ptr */
		num_components += component_ptr->numrows();

if (decomp_control::DEBUG >= decomp_control::INTERMEDIATE)
{
	debug << "main:   number of new components = " << component_ptr->numrows() << endl;
	debug << "main:   total number of components = " << num_components << endl;
}			

		/* tidy up memory, we've pushed a copy of component_ptr
		   onto component_list so don't need the copy returned by split */
		delete component_ptr;	
	}
	
	matrix<int> track_basis(kdim,n);
	if (num_components > kdim)
	{
		/* reduce the set of components to a linearly independent set.  We know the
		   pattern_basis elements are linearly independent, so we start with any of 
		   those that are tracks and then add in component tracks of the others as necessary.
		   A set of kdim tracks are linearly independent iff the rank of the matrix formed of
		   that set is kdim.
		*/
		
		list<matrix<int> >::iterator mptr = component_list.begin();
		vector<int> track_flag(kdim);
		int index = 0;	
		
		while (mptr != component_list.end())
		{
			if (mptr->numrows() == 1)
				track_flag[index] = 1;
			index++;
			mptr++;
		}
		
if (decomp_control::DEBUG >= decomp_control::INTERMEDIATE)
{
	debug << "main: track_flags: ";
	for (int i=0; i<kdim; i++)
		debug << track_flag[i] << ' ';
	debug << endl;	
}		
		
		index=0;
		for (int i=0; i< kdim; i++)
		{
			if(track_flag[i] == 1)
			{
				for (int j=0; j< n; j++)
					track_basis[index][j] = pattern_basis[j][i]; // pattern_basis holds elements in columns
				index++;
			}
		}
		
if (decomp_control::DEBUG >= decomp_control::INTERMEDIATE)
{
	debug << "main: initial track_basis: " << endl;
	print (track_basis,debug,3, "main: ");
}			

		/* The rank of track_basis is index, the number of non-zero rows.  We now add rows from 
		   those components in the component_list corresponding to pattern basis elements that are 
		   not tracks.  For each such component we add it to the track basis and check whether the
		   rank increases.  If it does not we select an alternative track component, if it does we 
		   leave the track component as a new track_basis element and look for the next basis element.
		*/
		mptr = component_list.begin();
		
		/* skip over any tracks */
		while (mptr->numrows() == 1)
			mptr++;
			
		unsigned int cpt = 0;
		int old_rank; //only used for debugging
		
		while (index < kdim)
		{
			/* find the next component track to try */
			if (cpt == mptr->numrows())
			{
				mptr++;
				
				/* skip over any tracks */
				while (mptr->numrows() == 1)
					mptr++;
				
				cpt = 0;
			}	
			else
			{
				for (int i=0; i< n; i++)
					track_basis[index][i] = (*mptr)[cpt][i];

if (decomp_control::DEBUG >= decomp_control::INTERMEDIATE)
{
	debug << "main: try adding track: ";
	for (int i=0; i<n; i++)
		debug << track_basis[index][i] << ' ';
	debug << endl;
}			
					
				cpt++;
			}
			
			/* all we need do now is set index to the rank of track_basis.  This will increment index
			   if the new element is independent of the existing rows and will leave index unchanged
			   if it is not.  If index is left unchanged when we look for the next component we shall
			   overwrite the last one we tried, if index is incremented, the last new element is left.
			*/
			old_rank = index;
			index = rank(track_basis); 			

if (decomp_control::DEBUG >= decomp_control::INTERMEDIATE)
{
	debug << "main:   new rank = " << index;
	if (index == old_rank)
		debug << " track is linearly dependent, discard it" << endl;
	else
		debug << " track is linearly independent, keep it" << endl;
	
}			
		}		
	}
	else
	{
		for (int i=0; i< kdim; i++)
		{
			for (int j=0; j< n; j++)
				track_basis[i][j] = pattern_basis[j][i];
		}		
	}

	
	output << "\n\nTrack basis (size " << kdim << " x " << n << ")\n";
	print (track_basis, output, 4, "");
	output << endl;

if (decomp_control::DEBUG >= decomp_control::BASIC)
{
	debug << "main: track_basis: " << endl;
	print (track_basis,debug,3, "main: ");
}			
	
	/* separate the tracks in the track_basis into separating and non-separating tracks */
	list <vector<int> > separating_tracks;

if (decomp_control::DEBUG >= decomp_control::BASIC)
	debug << "main: segregate track_basis into separating and non-separating tracks:" << endl;

	/* For debugging purposes, keep a record of the track basis 
	   element corresponding  to each separating track
	*/
	list<int> track_basis_elt;

	for (int i=0; i< kdim; i++)
	{
		vector<int> track(n);
		for (int j=0; j< n; j++)
			track[j] = track_basis[i][j];		


if (decomp_control::DEBUG >= decomp_control::BASIC)
{
	debug << "main:   track: ";
	for (int j=0; j<n; j++)
		debug << track_basis[i][j] << ' ';
}				
output << "\ncheck track basis element " << i << endl;

		if (decompose(track, K, false, false, output) & DECOMPOSE_NON_SEPARATING)
		{

if (decomp_control::DEBUG >= decomp_control::BASIC)
	debug << "is non-separating" << endl;

		}
		else
		{

if (decomp_control::DEBUG >= decomp_control::BASIC)
	debug << "is separating" << endl;

			separating_tracks.push_back(track);
			track_basis_elt.push_back(i);
		}		
	}
	
if (decomp_control::DEBUG >= decomp_control::BASIC)
	debug << "main: prune list of " << separating_tracks.size() << " separating tracks:" << endl;

output << "\nprune list of " << separating_tracks.size() << " separating tracks:" << endl;

	/* remove any tracks from separating_tracks that give decompositions known to be trivial */
	list<vector<int> >::iterator tptr = separating_tracks.begin();
	list<int>::iterator eptr = track_basis_elt.begin();
	
	while (tptr != separating_tracks.end())
	{
		output << "\n\ntrack basis element " << *eptr << endl;		
		if ((decompose(*tptr, K, false, true, output) & DECOMPOSE_NON_TRIVIAL) == false)
		{

if (decomp_control::DEBUG >= decomp_control::BASIC)
{
	debug << "main:   track basis element " << *eptr << ": ";
	for (int j=0; j<n; j++)
		debug << (*tptr)[j] << ' ';
	debug << "gives a trivial decomposition" << endl;
}				


			tptr = separating_tracks.erase(tptr);
			eptr = track_basis_elt.erase(eptr);
		}
		else
		{
			tptr++;
			eptr++;
		}
	}

if (decomp_control::DEBUG >= decomp_control::BASIC)
	debug << "main: " << separating_tracks.size() << " separating tracks remain after removing those with obvious trivial decompositions" << endl;

output << "\n" << separating_tracks.size() << " separating tracks remain after removing those with obvious trivial decompositions\n" << endl;

	/* remove any tracks from separating_tracks that are either incompatible with 
	   another track in the track_basis 
	
	tptr = separating_tracks.begin();
	eptr = track_basis_elt.begin();
	
	while (tptr != separating_tracks.end())
	{
		bool incompatible = false;
		
		for (int i=0; i< kdim; i++)
		{
			if (!compatible(track_basis[i],*tptr, K))
			{

if (decomp_control::DEBUG >= decomp_control::BASIC)
{
	debug << "main:   track basis element " << *eptr << ": ";
	for (int j=0; j<n; j++)
		debug << (*tptr)[j] << ' ';
	debug << "incompatible with basis element " << i << ": ";
	for (int j=0; j<n; j++)
		debug << track_basis[i][j] << ' ';
	debug << endl;
}				

output << "track basis element " << *eptr << ": ";
for (int j=0; j<n; j++)
	output << (*tptr)[j] << ' ';
output << "incompatible with basis element " << i << ": ";
for (int j=0; j<n; j++)
	output << track_basis[i][j] << ' ';
output << endl;

				incompatible = true;
				break;
			}
		}
		
		if (incompatible)
		{
			tptr = separating_tracks.erase(tptr);
			eptr = track_basis_elt.erase(eptr);
		}
		else
		{
			tptr++;
			eptr++;
		}
	}

if (decomp_control::DEBUG >= decomp_control::BASIC)
{
	debug << "main: there are " << separating_tracks.size() 
	      << " non-trivial separating track_basis elements compatible with all basis elements: " << endl;

	tptr = separating_tracks.begin();
	
	while (tptr != separating_tracks.end())
	{
		debug << "main:   track: ";
		for (int j=0; j<n; j++)
			debug << (*tptr)[j] << ' ';
		debug << endl;
		tptr++;
	}
}				

output << "\n\nthere are " << separating_tracks.size() 
      << " non-trivial separating track_basis elements compatible with all basis elements: " << endl;
*/

	tptr = separating_tracks.begin();
	int track_number = 0;
	while (tptr != separating_tracks.end())
	{
		output << "track " << ++track_number << ": ";
		for (int j=0; j<n; j++)
			output << (*tptr)[j] << ' ';
		output << endl;
		tptr++;
	}

	eft << separating_tracks.size() << endl;

	tptr = separating_tracks.begin();

	while (tptr != separating_tracks.end())
	{
		for (int j=0; j<n; j++)
			eft << (*tptr)[j] << ' ';
		eft << endl;
		tptr++;
	}



	/* Identify the set of 2-cells we are going to fold.  These are the 2-cells of X_1, the 
	   subdivision of X described in [JSJ].  Each 2-cell in two_cell is subdivided into three
	   which we enumerate respecting the canonical order.  Thus, if a 2-cell has vertices 
	   u v w and the boundary in canonical order is uv vw wu then subdividing creates another
	   vertex p and edges pu, pv, pw (oriented from p as indicated by the notation), then the
	   canonical ordering of the 2-cells of the subdivision are uv -pv pu, vw -pw pv, wu -pu pw
	   
	   We shall number the additional vertices in 2-cell order from the current number_of_0_cells
	   and shall number the additional edges from the current number_of_1_cells in 2-cell order
	   and within a 2-cell in the order pu, pv, pw.  Thus the 1-cells in 2-cell i are numbered
	   number_of_1_cells + 3*(i-1) + [1 | 2 | 3]
	
	matrix<int> subdivided_two_cell(3*number_of_2_cells,3);
	
	for (int i=0; i< number_of_2_cells; i++)
	{
		int base_1_cell = number_of_1_cells + 3*i;
		
		subdivided_two_cell[3*i][0] = two_cell[i][0];
		subdivided_two_cell[3*i][1] = -(base_1_cell+2);
		subdivided_two_cell[3*i][2] = (base_1_cell+1);
	
		subdivided_two_cell[3*i+1][0] = two_cell[i][1];
		subdivided_two_cell[3*i+1][1] = -(base_1_cell+3);
		subdivided_two_cell[3*i+1][2] = (base_1_cell+2);
	
		subdivided_two_cell[3*i+2][0] = two_cell[i][2];
		subdivided_two_cell[3*i+2][1] = -(base_1_cell+1);
		subdivided_two_cell[3*i+2][2] = (base_1_cell+3);
	}
	
if (decomp_control::DEBUG >= decomp_control::BASIC)
{
	debug << "main: subdivided_two_cell: " << 3*number_of_2_cells << " x " << 3 << endl;
	print (subdivided_two_cell,debug,4, "main: ");
}
	*/
	
}
else
{
	/**************** Evaluate extreme fundamental tracks *************/
	if ( kdim == 1 )
	{
	    cout << "\n\nThe only fundamental track is (1,...,1).\n";
	
	
	    output << "\n\nThe only fundamental track is (1,...,1)." << endl;
	
	    int s = 1;
	    eft.write( (char*)&s, sizeof(s));
	    for (int i=0; i<n; i++)
			eft.write( (char*)&s, sizeof(s));
	    
		eft.flush();
	
if (decomp_control::DEBUG >= decomp_control::BASIC)
	debug << "\ns = " << s << "\n";
	}
	else
	{
	
		/************************************************
		Check to see if any of the rows of 'cone matrix'
		are duplicated.  If so, note any duplications.
		************************************************/
	
		int repeat[n];
	
		for ( int i = 0; i<n; i++)
	    	repeat[i] = 0;
	
		repeat_count = 0;
	
		/* We shall set 'repeat[i]' to 1 iff row i of 'cone matrix' is a repeat of
		 row j of 'cone matrix'. */
	
		for ( int i = 0; i < n-1 ; i++)
		{
		    if (repeat[i] == 0)
			{
				for (int j = i+1 ; j < n ; j++)
				{
			    	if (repeat[j] == 0)
		    		{
						equal = true;
						int k = 0;
						do
						{
					    	if ( cone_matrix[i][k] != cone_matrix[j][k] )
								equal = false;
						    	k++;
						} while ( equal && k < kdim );
	
						if ( equal )
						{
					    	repeat[j] = 1;
					    	repeat_count++;
						}
		    		}
				}
			}
		}
	
if (decomp_control::DEBUG >= decomp_control::BASIC)
{
	debug << "main: cone matrix row repeat_count = " << repeat_count << endl;
	
	if (repeat_count)
	{
		debug << "main: repeat flags: ";
		for (int i=0; i<n; i++)
			debug << repeat[i] << ' ';			
		debug << endl;
	}
}
	
		/****** Adjust n, the number of rows in the cone_matrix to allow for any repetitions. ********/
		n -= repeat_count;
		
		cout << "\n\nExtreme fundamental tracks determined by a cone in " << kdim
		     << " dimensions,\ndescribed by " << n << " supporting hyperplanes.\n\n";
		
		/************  Declare space for the adjusted 'cone matrix'. See below *******/
		matrix<int> cone_matrix_2(kdim,n);
		
		/***************************************************************
		Let U denote the 'cone matrix'.  We know there is an x0 in real
		'kdim'-space such that Ux0=(1,...,1): suppose r1 and r2 are rows of U such
		that r2=ar1, for some non zero real number a, then r1x0=1 and r2x0=ar1x0=a,
		so a=1.  It follows that there are exactly 'kdim' rows of U containing
		exactly one non zero entry, after repetitions have been removed.  Moreover,
		if r1 contains exactly one non zero entry, it must be positive, r1
		asserting, merely, that xi>=0 (if this non zero entry is in the ith place
		of r1), by the way we evaluated the basis in 'kbasis'.  We work with a
		derivative of the transpose of 'cone matrix', storing this adjusted
		"cone matrix" in 'cone_matrix_2'.  The first 'kdim' columns are assigned to
		be the 'kdim'x'kdim' identity matrix and the remaining n-'kdim' columns are
		assigned to be those distinct rows of 'cone matrix' containing at least two
		non zero entries.
		*******************************************************************/
		for (int i=0; i < kdim; i++)
		{
		    for ( int j = 0 ; j < i-1 ; j++ )
				cone_matrix_2[i][j] = 0;
				
		    cone_matrix_2[i][i] = 1;
		
		    for ( int j = i+1 ; j < kdim ; j++ )
				cone_matrix_2[i][j] = 0;
		}
		
		column = n-1;
		
		for ( unsigned int i = 0; i < cone_matrix.numrows(); i++)
		{
		    if ( repeat[i] == 0	 && non_zero_count(cone_matrix, i) != 1)
		    {
				for (int j = 0; j < kdim; j++)
			    	cone_matrix_2[j][column] = cone_matrix[i][j];
					
				column -= 1;
		    }
		}
		
		output << "\n\nAdjusted cone matrix, size " << kdim << " x " << n << "\n";
		print(cone_matrix_2, output, 4, "");
		output << endl;
		
if (decomp_control::DEBUG >= decomp_control::BASIC)
{
	debug << "main: adjusted cone matrix, size " << kdim << " x " << n << ":\n";
	print(cone_matrix_2, debug, 4, "main: ");
}
		
		
		/* If we're just analysing a cell complex we're done at this point */
		if (CALCULATE_EXTREME_RAYS)
		{
		    /****************** Greenbergs algorithm. ***********************/
		    
		    /*********** Initializations for use in the algorithm. **********/
		    list<ray> raylist;
		    int s = greenberg(cone_matrix_2, raylist, start_time);
		    
		
		    /*************** Evaluate the extreme fundaental tracks *****************/
		    n = cone_matrix.numrows(); // reset n to it's original value.
		    
		    matrix<int> ext_fund_tracks(s,n);
		    
		    list<ray>::iterator rptr = raylist.begin();
		    
		    for (int i=0; i< s; i++)
		    {
		        /* Calculate the contribution to 'ext_fund_tracks' from the raylist entry,
		           which contains an extreme ray of the cone in 'kdim'-space
		           determined by 'cone_matrix' (and 'cone_matrix_2').  Notice it
		           contains some redundant information in the last n-'kdim' places. */
		        for (int j=0; j<n; j++)
		    	{
		    		for (int k=0; k<kdim; k++)
		    		    ext_fund_tracks[i][j] += cone_matrix[j][k] * rptr->r[k];
		    	}
		    	rptr++;
		    }
		    
		    /* Divide down by gcds */
		    for (int i=0; i<s; i++)
		    {
		        g = gcdn(ext_fund_tracks, i);
		        if ( g !=0 )
		    	for (int j=0; j<n; j++)
		    	    ext_fund_tracks[i][j] /= g;
		    }
		    
			output << "\nVertex solutions (extreme fundamental tracks), n=" << n << " s=" << s << endl;
			int count=1;
			for (int i=0;i< s;i++)
			{
				output << setw(3) << count++ << ". ";
				for (int j=0;j<n;j++) 
					output << setw(3) << ext_fund_tracks[i][j];		
				output << endl;
			}
		
if (decomp_control::DEBUG >= decomp_control::BASIC)
{
	debug << "main: found " << s << " vertex solutions (extreme fundamental tracks):" << endl;
//    debug << "\ns = " << s << endl;
	print(ext_fund_tracks, debug, 3, "main: ");
}
		    
		
		    /* If we're just calculating vertex solutions we can skip the fundamental solution calculation */
		    if (CALCULATE_FUNDAMENTAL_SOLUTIONS)
		    {
		    	
		        /*                               Added September 2009 
		           Work out the set of minimal fundamental solutions to the matching equations.
		           First, we determine the point z in n-space given by z[i] = max(t[i]) for all 
		           extreme fundamental tracks t
		        */
		
				/***** Prepare a K_complex structure to hold simplicial complex parameters *****/
		        K_complex K;
		        K.number_of_0_cells = number_of_0_cells;
		        K.number_of_1_cells = number_of_1_cells;
		        K.number_of_2_cells = number_of_2_cells;
		        K.occ_of_cv = occ_of_cv;
		        K.generator_corresponding_to_1_cell = generator_corresponding_to_1_cell;
		        K.number_of_generators = number_of_generators;
		        K.number_of_relators = number_of_relators;
		        K.relator = relator;
		        K.number_of_1_cells_in_rosette = number_of_1_cells_in_rosette;
		        K.generators_of_order_2 = generators_of_order_2;
		        K.occ_of_order_2_gen = occ_of_order_2_gen;
		        K.occurrences_of_1_cell = occurrences_of_1_cell;
		        K.length_of_attaching_map = length_of_attaching_map;
		        K.two_cell = &two_cell;
		        K.attaching_map = &attaching_map;
		        K.matching_system = &matching_system;
		        K.s = s;
		  
		        vector<int> w(n);
		
if (decomp_control::DEBUG >= decomp_control::BASIC)        
	debug << "main: calculate upper bound w for fundamental solutions" << endl;
			
		        for (int c=0; c< n; c++)
		        for (int r=0; r< s; r++)
		        {
					w[c] += ext_fund_tracks[r][c];
				}
		
				if (OPTIMIZED_FUNDAMENTAL_SOLUTION_SEARCH)
				{		/* all fundamental solutions are <= w but we want to make an n-limit
						   out of w, and n-limit inequalities are strict, so we need to increment
						   each component of w by one.
						*/
						for (int c=0; c< n; c++)
							w[c] += 1;
				}
		
		
			       
if (decomp_control::DEBUG >= decomp_control::BASIC)
{
	debug << "main: w = (" << w[0];
	for (int i=1; i< n; i++)
		debug << ',' << w[i];
	debug << ')' << endl;
}
		        
		        cout << "main: w = (" << w[0];
		      	for (int i=1; i< n; i++)
		      		cout << ',' << w[i];
		       	cout << ')' << endl;
		        
		        /* check all integral points in the n-cuboid determined by the origin and w
		           to see if they are solutions to the matching equations.  These will be the
		           set of fundamental solutions described in the 2009 AB/MJD paper.
		        
		           By considering the start of the first relator disc for a relation of
		           length > 3 it can be seen that one of the matching equations is 
		           
		                  x1+x3 = x4+x                               (1)
		                           
		           There is also a relation x1+x2 = xi+xj but i and j are determined by the 
		           relation and therefore vary.  Now at the end of the last relator disc we have 
		           a corresponding situation.  Suppose for example there is a single relator disc
		           of length 5, so there are three triangular 2-cells and variables x1,...,x9.  
		           The last internal 1-cell gives the relation 
		           
		                 x4+x6 = x7+x8                               (2)
		                 
		           which does not involve x9 whereas our earlier relation does involve x1.  Thus if 
		           we enumerate the integral points of the n-cuboid in the natural way (from the right) 
		           as if we were increasing an n-digit number, we will have situations where the relation 
		           corresponding to (2) is not satisfied but we repeatedly increment the last digit to its 
		           maximum value and potentially have to check all the matching equations until we discover
		           that (2) fails.  We shall therefore enumerate the n-cuboid from the left.      
		        
		           To check the matching system efficiently, we establish an sx4 matrix whose first two 
		           columns contian the two positive indices from a row of the matching system and the second 
		           two columns contian the two negative indices.  Since the matching system sometimes contains
		           a row with only one positive and one negative entry (when a relator starts aa... for example)
		           we need a method of detecting this and so store indices from 1, allowing 0 to indicate an 
		           unused entry.
		        */
		        
		        matrix<int> matching_indices(m,4);
		        for (int i=0; i< m; i++) // m = matching_system.numrows()
		        {
		        	for (int j=0; j< n; j++)
		        	{
		        		if (matching_system[i][j] > 0)
		        		{
		        			if (matching_indices[i][0] == 0)
		        				matching_indices[i][0] = j+1;
		        			else
		        				matching_indices[i][1] = j+1;
		        		}
		        		else if (matching_system[i][j] < 0)
		        		{
		        			if (matching_indices[i][2] == 0)
		        				matching_indices[i][2] = j+1;
		        			else
		        				matching_indices[i][3] = j+1;
		        		}
		        	}
		        }
		        
if (decomp_control::DEBUG >= decomp_control::BASIC)
{
	debug << "main: matching_indices = " << endl;
	print(matching_indices,debug,0,"main:  ");
}
		        
		        vector<int> u(n);
		        list<vector<int> > fundamental_solutions;
		
		       
				if (OPTIMIZED_FUNDAMENTAL_SOLUTION_SEARCH)
				{
						/* The search algorithm used to determine the minimal fundamental solutions is described in the
						   notes "fundamental-solutions" that should accompany the software distribution, they are also 
						   available at www.layer8.co.uk/maths/tracks
						*/
						
				  		/* create an n_limit initialized to the w-cuboid */
				  		n_limit three_m_limit(w);
				  		
if (decomp_control::DEBUG >= decomp_control::BASIC)
{
	debug << "main: starting enumeration of w-cuboid" << endl;
//	three_m_limit.print(debug,0,"main#: ");
}	
				
				        bigint w_cuboid_volume=1;
				        for (int i=0; i< n; i++)
				        	w_cuboid_volume *= w[i]+1;
				        	
						output << "\nvolume of w-cuboid = " << w_cuboid_volume << endl;
						cout << "volume of w-cuboid = " << w_cuboid_volume << endl;
				
if (decomp_control::DEBUG >= decomp_control::BASIC)
	debug << "main: volume of w-cuboid = " << w_cuboid_volume << endl;
					
						/* we know the extreme fundamental tracks are minimal fundamental solutions, so add these 
						   to the list and update three_m_limit with each one.
						*/
					
						for (int i=0; i< s; i++)
						{
					
							for (int j=0; j< n; j++)
								u[j] = ext_fund_tracks[i][j];       	
							
				       		fundamental_solutions.push_back(u);		
				
if (decomp_control::DEBUG >= decomp_control::BASIC)
{
	debug << "main: updating three_m_limit with extreme fundamental track " << i+1 << ": ";
	for (int j=0; j< n; j++)
		debug << u[j] << " ";
	debug << endl;	
}    			
				      			
				       		// update the n-limit with u and adjust the current thresholds
				       		
				       		update_n_limit (three_m_limit, u);
				
				       		bigint volume_remaining = non_excluded_region_size (three_m_limit, n);
				       		cout << "\n i = " << i << " volume remaining = " << volume_remaining << endl; //*bigint(100)/w_cuboid_volume << "%" << endl;
				       		output << "\n i = " << i << " volume remaining = " << volume_remaining << endl; //*bigint(100)/w_cuboid_volume << "%" << endl;
						}
				
						w_cuboid_volume_remaining = non_excluded_region_size (three_m_limit, n);
				       	cout << "\nvolume remaining after adding extreme fundamental solutions = " << w_cuboid_volume_remaining << endl;
				       	output << "\nvolume remaining after adding extreme fundamental solutions = " << w_cuboid_volume_remaining << endl;
				
			if (decomp_control::DEBUG >= decomp_control::BASIC)
				debug << "main: remaining volume after adding extreme fundamental solutions = " << w_cuboid_volume_remaining << endl;
				
						/* reset u to the zero vector */
						for (int i=0; i< n; i++)
							u[i] = 0;
							
				        bool cuboid_enumeration_complete = false;
				
if (decomp_control::DEBUG >= decomp_control::BASIC)
{
	debug << "main: three_m_limit after adding extreme fundamental tracks:" << endl;
	three_m_limit.print(debug,0,"main: ");
	debug << endl;
}
					
						/* create an array of threshold iterators to hold pointers to the current thresholds
						   and initialize the array to point to the first threshold for each variable.
						*/	
						list<threshold>::iterator current_threshold[n];
						current_threshold[n-1] = three_m_limit.thresholds.begin();
						for (int i=n-2; i >= 0; i--) // note n = 3m > 3
							current_threshold[i] = current_threshold[i+1]->second->thresholds.begin();
				
if (decomp_control::DEBUG >= decomp_control::BASIC)
{
	debug << "main: current thresholds: ";
	for (int i=0; i< n; i++)
		debug << *current_threshold[i] << ' ';
	debug << endl;
}	
				
					
						//cout << "main: three_m_limit after adding extreme fundamental tracks:" << endl;
						//three_m_limit.print(cout,0,"");
						
				        do
				        {
							if (++u[0] >= current_threshold[0]->first)
							{
/*if (decomp_control::DEBUG >= decomp_control::BASIC)
{
	debug << "main: u[0] = " << u[0] << " has reached the threshold value " 
	      << current_threshold[0]->first << " checking thresholds" << endl;
}
*/
								cuboid_enumeration_complete = check_thresholds(1,u,three_m_limit,current_threshold);
							}
				
				        	if (!cuboid_enumeration_complete)
				        	{       
				        		/* check if u satisfies the matching equations and if so add 
				        		   it to the list of fundamental solutions
				        		*/
/*if (decomp_control::DEBUG >= decomp_control::BASIC)
{
	debug << "main: checking u = ";
	for (int i=0; i< n; i++)
		debug << u[i] << ' ';
	debug << endl;
}*/
				
				        		if (satisfies_matching_equations(matching_indices,u))
				        		{
if (decomp_control::DEBUG >= decomp_control::BASIC)
{
	debug << "main: found fundamental solution ";
	for (int i=0; i< n; i++)
		debug << u[i] << ' ';
	debug << endl;
	debug << "main: updating 3m-limit" << endl;
}
									vector<int> saved_u(u);
									
				        			fundamental_solutions.push_back(u);		
				        			
				        			/* update the n-limit with u and adjust the current thresholds */
				        			update_n_limit (three_m_limit, u);
				        			cuboid_enumeration_complete = reset_current_thresholds(three_m_limit,three_m_limit,u,current_threshold);
				        			
				        			for (int i=0; i<n; i++)
				        			{
										if (u[i] != saved_u[i])
										{
if (decomp_control::DEBUG >= decomp_control::BASIC)
	debug << "main: reset_current_thresholds has adjusted u, decrementing u[0]" << endl;
											u[0]--;
											break;
										}
									}
				        			
				        			/* calculate how much more of the w-cuboid we have left to check
				        			bigint volume_remaining = non_excluded_region_size (three_m_limit, n);        			
				        			cout << "\nvolume remaining after updating with u = " << volume_remaining << endl;
									*/
								}
				        	}	
				        } while (!cuboid_enumeration_complete);
				
if (decomp_control::DEBUG >= decomp_control::BASIC)
	debug << "main: optimized cuboid enumeration complete" << endl;	
				}
				else
				{        
				        bool cuboid_enumeration_complete = false;
				        bigint total_num_points=1;
				        for (int i=0; i< n; i++)
				        	total_num_points *= w[i]+1;
				        	
if (decomp_control::DEBUG >= decomp_control::BASIC)
	debug << "main: number of integral points in w n-cuboid = " << total_num_points << endl;
					
				        bigint percentile = total_num_points/bigint(100);
				        bigint bigcount = 1;
				        bigint threshold = percentile;
				        
				        cout << "\nEvaluating fundamental solutions" << endl;
				        cout << "Total number of points to consider = " << total_num_points << endl;
				        do
				        {
				        	/* display progress information */
				        	if (++bigcount == threshold)
				        	{
				        		cout << "\r" << threshold/percentile << "%" << flush;
				        		threshold += percentile;
				        	}
				        	
				        	/* increment u from the left as an n-digit number */
				        	int place=0;
				        	bool done=false;
				        	do
				        	{
				        		if (++u[place] > w[place]) 
				        		{
				        			u[place]=0;
				        			place++;
				        		}
				        		else
				        			done=true;			
				        	} while (!done && place < n);
				        	
				        	if (place == n)
				        		cuboid_enumeration_complete = true;
				        		
				        	if (!cuboid_enumeration_complete)
				        	{
				        		/* check if u satisfies the matching equations and if so add 
				        		   it to the list of fundamental solutions
				        		*/
				        		if (satisfies_matching_equations(matching_indices,u))
				        			fundamental_solutions.push_back(u);		
				        	}	
				        } while (!cuboid_enumeration_complete);
				        
if (decomp_control::DEBUG >= decomp_control::BASIC)
	debug << "main: original cuboid enumeration complete" << endl;	
				
				} // end of simple search for fundamental solutions
		
				output << "\n\nCuboid search list of fundamental solutions, n=" << n << " s=" << fundamental_solutions.size() << endl;
				int count = 1;
				list<vector<int> >::iterator fptr = fundamental_solutions.begin();	
				while (fptr != fundamental_solutions.end())
				{
					output << setw(3) << count++ << ". ";
					for (int j=0; j< n; j++)
						output << setw(3) << (*fptr)[j];
					output << endl;			
					fptr++;
				}
		        
if (decomp_control::DEBUG >= decomp_control::BASIC)
{
	debug << "main: list of fundamental solutions found, n=" << n << " s=" << fundamental_solutions.size() << endl;
    list<vector<int> >::iterator fptr = fundamental_solutions.begin();
        	
    while (fptr != fundamental_solutions.end())
    {
		debug << "main:   ";
    	for (int j=0; j< n; j++)
    		debug << setw(3) << (*fptr)[j];
    	debug << endl;
        			
    	fptr++;
    }
}
		        
				reduce_to_minimal_set(fundamental_solutions);
		    
if (decomp_control::DEBUG >= decomp_control::BASIC)
{
   	debug << "main: found " << fundamental_solutions.size() << " minimal solutions:" << endl;
   	list<vector<int> >::iterator fptr = fundamental_solutions.begin();
        	
  	while (fptr != fundamental_solutions.end())
   	{
		debug << "main:   ";
        for (int j=0; j< n; j++)
        	debug << setw(3) << (*fptr)[j];
        debug << endl;
        			
        fptr++;
  	}
}
		        
				cout << "\n\nMinimal solutions, n=" << n << " s=" << fundamental_solutions.size() << endl;
				count = 1;
				fptr = fundamental_solutions.begin();	
				while (fptr != fundamental_solutions.end())
				{
					cout << setw(3) << count++ << ". ";
					for (int j=0; j< n; j++)
						cout << setw(3) << (*fptr)[j];
					cout << endl;
					fptr++;
				}     
		
		//    	eft << fundamental_solutions.size() << endl;
		    	   	
				output << "\n\nMinimal solutions, n=" << n << " s=" << fundamental_solutions.size() << endl;
				count = 1;
				fptr = fundamental_solutions.begin();	
				while (fptr != fundamental_solutions.end())
				{
					output << setw(3) << count++ << ". ";
					for (int j=0; j< n; j++)
						output << setw(3) << (*fptr)[j];
					output << endl;
		//			for (int j=0; j< n; j++)
		//				eft << setw(3) << (*fptr)[j];
		//			eft << endl;
					
					fptr++;
				}
				
				/* Calculate the minimal untwisted solutions from the minimal solutions.  We construct a list
				   of the following:
				   
				   1. separating minimal solutions
				   2. untwisted non-separating mimimal solutions
				   3. doubled twisted minimal solutions
				   4. all combinations of distinct non-separating minimal solutions that sum to give tracks
				   
				   we then reduce this list to a minimal set.
				*/
		        list<vector<int> > minimal_untwisted_solutions;
		        
if (decomp_control::DEBUG >= decomp_control::BASIC)
	debug << "main: checking for minimal twisted fundamental solutions:" << endl;
		        
		        fptr = fundamental_solutions.begin();
		        
		        /* Add the first three types of minimal untwisted solutions to the list */
		        while (fptr != fundamental_solutions.end())
		        {
		        
if (decomp_control::DEBUG >= decomp_control::BASIC)
{
  	debug << "main:   fptr ";
	for (int i=0; i< n; i++)
		debug << setw(3) << (*fptr)[i];
}
		        
		        	vector<int> doubled_track(n);
		        	for (int i=0; i< n; i++)
		        		doubled_track[i] = 2* (*fptr)[i];
		        
		        	bool twisted_track = split(doubled_track,K,true, output); // nothing will be written to output because short_split is true
		        	
		        	if (twisted_track)
		        	{
if (decomp_control::DEBUG >= decomp_control::BASIC)
   	debug << " twisted, adding doubled track" << endl;
		        	
		        		minimal_untwisted_solutions.push_back(doubled_track);
		        	}
		        	else
		        	{
if (decomp_control::DEBUG >= decomp_control::BASIC)
  	debug << " untwisted, adding track" << endl;
		        
		        		minimal_untwisted_solutions.push_back(*fptr);
		        	}
		        	
		      		fptr++;
		        }
		        
		        /* Evaluate all combinations of distinct non-separating minimal solutions that sum to give tracks.
		           First we isolate in a matrix those tracks that are non-separating minimal solutions.
		        */
		        matrix<int>nsm_solutions(fundamental_solutions.size(),n); //non-separating minimal solutions
		        int nsm_count=0;
		        
				fptr = fundamental_solutions.begin();	
				while (fptr != fundamental_solutions.end())
				{
					if (decompose(*fptr, K, false, false, output) & DECOMPOSE_NON_SEPARATING) // no echo no full_decomp
					{
						for (int i=0; i< n; i++)
							nsm_solutions[nsm_count][i] = (*fptr)[i];
						nsm_count++;
					}
					
					fptr++;
				}        
		
if (decomp_control::DEBUG >= decomp_control::BASIC)
{
   	debug << "main: " << nsm_count << " non-separating minimal solutions:" << endl;
   	for (int i=0; i< nsm_count; i++)
   	{
		debug << "main:   " << i+1 << ". ";
		for (int j=0; j< n; j++)
			debug << nsm_solutions[i][j] << ' ';
		debug << endl;
	}
}        
		
				/* consider all combinations of distinct non-separating minimal solutions to see if they give tracks */
				for (int k = 2; k <= nsm_count; k++)
				{
					vector<int> comb(k);
					for (int i=0; i<k; i++)
						comb[i]=i;
						
					do
					{
						/* test the combination of non-separating minimal tracks determined by comb */
if (decomp_control::DEBUG >= decomp_control::BASIC)
{
	debug << "main: testing combination ";
	for (int i=0; i< k; i++)
		debug << comb[i] << ' ';
	debug << endl;
}
						vector<int> pattern(n);
						
						for (int i=0; i<k; i++)
						{
							for (int j=0; j< n; j++)
								pattern[j] += nsm_solutions[comb[i]][j];
						}
						
if (decomp_control::DEBUG >= decomp_control::BASIC)
{
	debug << "main:   resulting pattern ";
	for (int i=0; i< n; i++)
		debug << pattern[i] << ' ';
	debug << endl;
}
						bool track = split(pattern,K,true, output);
						
						if (track)
						{
if (decomp_control::DEBUG >= decomp_control::BASIC)
	debug << "main:   is a track, adding to list of minimal untwisted solutions" << endl;
			
							minimal_untwisted_solutions.push_back(pattern);
						}
						else
						{
if (decomp_control::DEBUG >= decomp_control::BASIC)
	debug << "main:   is a not track" << endl;
						}
						
					} while (next_combination(comb, nsm_count));
				}
		
if (decomp_control::DEBUG >= decomp_control::BASIC)
{
   	debug << "main: " << minimal_untwisted_solutions.size() << " minimal untwisted solutions before minimizing:" << endl;
   	list<vector<int> >::iterator fptr = minimal_untwisted_solutions.begin();
    int count=0;
   	while (fptr != minimal_untwisted_solutions.end())
   	{
   		debug << "main: " << ++count << ". ";
   		for (int j=0; j< n; j++)
   			debug << setw(3) << (*fptr)[j];
   		debug << endl;
        			
   		fptr++;
   	}
}		
				/* now reduce the list of minimal untwisted solutions to a minimal set */
				reduce_to_minimal_set(minimal_untwisted_solutions);
		        
if (decomp_control::DEBUG >= decomp_control::BASIC)
{
   	debug << "main: found " << minimal_untwisted_solutions.size() << " minimal untwisted solutions:" << endl;
   	list<vector<int> >::iterator fptr = minimal_untwisted_solutions.begin();
        	
    int count=0;
   	while (fptr != minimal_untwisted_solutions.end())
   	{
       	debug << "main: " << ++count << ". ";
   		for (int j=0; j< n; j++)
   			debug << setw(3) << (*fptr)[j];
   		debug << endl;
        			
   		fptr++;
   	}
}
		        	
				cout << "\n\nMinimal untwisted solutions, n=" << n << " s=" << minimal_untwisted_solutions.size() << endl;
				count = 1;
				fptr = minimal_untwisted_solutions.begin();	
				while (fptr != minimal_untwisted_solutions.end())
				{
					cout << setw(3) << count++ << ". ";
					for (int j=0; j< n; j++)
						cout << setw(3) << (*fptr)[j];
					cout << endl;
					fptr++;
				}     
		
		    	eft << minimal_untwisted_solutions.size() << endl;
		    	   	
				output << "\n\nMinimal untwisted solutions, n=" << n << " s=" << minimal_untwisted_solutions.size() << endl;
				count = 1;
				fptr = minimal_untwisted_solutions.begin();	
				while (fptr != minimal_untwisted_solutions.end())
				{
					output << setw(3) << count++ << ". ";
					for (int j=0; j< n; j++)
						output << setw(3) << (*fptr)[j];
					output << endl;
					for (int j=0; j< n; j++)
						eft << setw(3) << (*fptr)[j];
					eft << endl;
					
					fptr++;
				}
		    }
		    else
		    {
				/* not CALCULATE_FUNDAMENTAL_SOLUTIONS */
				eft << s;
				eft << ext_fund_tracks << endl;
				cout << "\n\n" << s << " vertex solutions (extreme fundamental tracks)" << endl;
				output << "\n\nVertex solutions (extreme fundamental tracks), n=" << n << " s=" << s << endl;
				int count=1;
				for (int i=0;i< s;i++)
				{
					output << setw(3) << count++ << ". ";
					for (int j=0;j<n;j++) 
						output << setw(3) << ext_fund_tracks[i][j];		
					output << endl;
				}
			}
		} // end of if (CALCULATE_EXTREME_RAYS)
	
	}  // end of else clause following if ( kdim == 1 )
}

cout << "\nTranscript written to " << textfile << ", search data to " << eftfile << endl;


output.close();
eft.close();

if (decomp_control::DEBUG >= decomp_control::BASIC)
    debug.close();

return 0;

} // End of function main


/***************  Functions required for setting up solve specific debug **************/

void set_main_debug_option_parameter(char* pptr, string option);
void set_main_debug_option(char* start, char* end)
{
	char  loc_buf[end-start+2];
	char* c1 = start;
	char* c2 = loc_buf;

	/* if both start and end are zero, display debug help information.
	   this has been included here in this manner so that each time a debug option
	   is added, the help will be updated (hopefully!)
	*/
	if (start == 0 && end == 0)
	{
		set_main_debug_option_parameter(0,"decomp");
		cout << "\t\tdecompose, boolean" << endl;
		cout << "\t\tsplit, boolean" << endl;
		cout << "\t\treduce, boolean" << endl;
		return;
	}

	do
	{
		*c2++ = *c1++;
	} while (c1 <= end);
	
	*c2 = '\0';

	char* pptr = strchr(loc_buf,'{');
	if (pptr)
	{
		*pptr++ = '\0';
	}

	/* now, even if there are parameters, the first part of loc_buf 
	   is a C-string that identifies the option */
	
	if (!strcmp(loc_buf,"decomp"))
	{
		decomp_control::DEBUG = decomp_control::SUMMARY;
		debug << "main::set_main_debug_option: setting debug option decomp_control::DEBUG = decomp_control::SUMMARY\n";		
		
		if (pptr)
		{
			/* check for any parameters */
			check_debug_option_parameters(pptr, "decomp");
		}
	}
	else if (!strcmp(loc_buf,"decompose"))
	{
		decomp_control::DECOMPOSE_DEBUG = true;
		debug << "main::set_main_debug_option: setting debug option decomp_control::DECOMPOSE_DEBUG\n";		
	}
	else if (!strcmp(loc_buf,"split"))
	{
		decomp_control::SPLIT_DEBUG = true;
		debug << "main::set_main_debug_option: setting debug option decomp_control::SPLIT_DEBUG\n";		
	}
	else if (!strcmp(loc_buf,"reduce"))
	{
		decomp_control::REDUCE_DEBUG = true;
		debug << "main::set_main_debug_option: setting debug option decomp_control::REDUCE_DEBUG\n";		
	}

	debug.flush();

}

void set_main_debug_option_parameter(char* pptr, string option)
{
	if (option == "decomp")
	{
		if (!pptr)
		{
			cout << "\t\tdecomp{summary|1:basic|2:intermediate|3:detail|4:exhaustive|5}, integer: default 0=off, no parameters sets basic" << endl;
		}
		else
		{
			if (!strcmp(pptr,"summary") || !strcmp(pptr,"1") )
			{
				decomp_control::DEBUG = decomp_control::SUMMARY;
				debug << "main::set_main_debug_option_parameter: setting debug option decomp_control::DEBUG = decomp_control::SUMMARY\n";		
			}
			if (!strcmp(pptr,"basic") || !strcmp(pptr,"2") )
			{
				decomp_control::DEBUG = decomp_control::BASIC;
				debug << "main::set_main_debug_option_parameter: setting debug option decomp_control::DEBUG = decomp_control::BASIC\n";		
			}
			else if (!strcmp(pptr,"intermediate") || !strcmp(pptr,"3"))
			{
				decomp_control::DEBUG = decomp_control::INTERMEDIATE;
				debug << "main::set_main_debug_option_parameter: setting debug option decomp_control::DEBUG = decomp_control::INTERMEDIATE\n";		
			}
			else if (!strcmp(pptr,"detail") || !strcmp(pptr,"4"))
			{
				decomp_control::DEBUG = decomp_control::DETAIL;
				debug << "main::set_main_debug_option_parameter: setting debug option decomp_control::DEBUG = decomp_control::DETAIL\n";		
			}
			else if (!strcmp(pptr,"exhaustive") || !strcmp(pptr,"5"))
			{
				decomp_control::DEBUG = decomp_control::DETAIL;
				debug << "main::set_main_debug_option_parameter: setting debug option decomp_control::DEBUG = decomp_control::EXHAUSTIVE\n";		
			}
		}
	}
}

void main_display_default_options() 
{
	cout << "\t\t[decomp{basic}]" << endl;
}


