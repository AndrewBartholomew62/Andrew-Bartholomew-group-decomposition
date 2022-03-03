using namespace std;

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <list>
#include <cstring>
#include <ctype.h>
#include <stdio.h>
#include <iomanip>
#include <vector>

extern ofstream	debug;

#include <util.h>
#include <matrix.h>
#include <decomp.h>


bool		TWISTED_ALLOWED = false;
bool		FLAGS_ON_FIRST_2_TORSION_1_CELL = false;

int decomp_control::DEBUG = decomp_control::OFF;
bool decomp_control::DECOMPOSE_DEBUG = false;
bool decomp_control::SPLIT_DEBUG = false;
bool decomp_control::REDUCE_DEBUG = false;

bool matrix_control::SINGLE_LINE_OUTPUT = false;


/********************* Function prototypes ***********************/
void decompose_all(string eftfile, ofstream& txt, vector<ray>& raylist, K_complex K);
void single(string eftfile, ofstream& txt, vector<ray>& raylist, K_complex K);
void range(string eftfile, ofstream& txt, vector<ray>& raylist, K_complex K);
void add(string eftfile, ofstream& txt, vector<ray>& raylist,K_complex K);
void enquire(ofstream& txt,K_complex K);
void decompose_given(ofstream& txt, K_complex& K);
bool match_ok (matrix<int> matching_system, vector<int>& p);
unsigned int decompose (vector<int>& track, K_complex& K, bool echo, bool full_decomp, ofstream& txt);
void compatibility_check(string eftfile, ofstream& txt, vector<ray>& raylist, K_complex K);
bool split 
(
	vector<int>& 	pattern, 
	K_complex 		K, 
	bool 			short_split, 
	ofstream& 		txt, 
	bool 			return_components = false, 
	matrix<int>**	component_ptr = 0
);
void initial_help();
//void print(matrix<int> m, ostream& s, int n, string prefix="");
bool compatible (ray& t1, ray& t2, K_complex K);
list<ray>::iterator find_cpt(list<ray>& cpt_list, ray cpt);
void display_characteristics(string eftfile, ofstream& txt, vector<ray>& raylist, K_complex K, bool test_all_tracks);
void maximal_sageev_n_cubes(string eftfile, ofstream& txt, vector<ray>& raylist, K_complex K);
void incompatible_clases(string eftfile, ofstream& txt, vector<ray>& raylist, K_complex K);

template <class T, class St> istream& operator >> (istream& s, matrix<T,St>& m)
{
	for (unsigned int i=0; i< m.numrows(); i++)
	for (unsigned int j=0; j< m.numcols(); j++)
		s >> m[i][j];
	return s;
}

bool debug_setup(char* argv_ptr);
void check_debug_option_parameters(char* pptr, string option);
void debug_help();

/******************* Main Function ************************/
int main (int argc, char* argv[])
{
   
bool		switches = false;
string version="23-5-10";


/* Determine command line switches, check command 
   line for a valid number of arguments */
if ( argc > 1 )
{
    if ( *(argv[1]) == '-')
	switches = true;

    if (switches)
    {
		if (strchr(argv[1], 'h') || strchr(argv[1], 'H'))
		{
			if (strchr(argv[1],'#') && strchr(argv[1], 'H') && strchr(argv[1], '!'))
			{
				debug_help();
			}
			else
			{
				initial_help();
				exit(0);
			}
		}

		if (strchr(argv[1], 'f'))
			FLAGS_ON_FIRST_2_TORSION_1_CELL = true;

		if (strchr(argv[1], 't'))
			TWISTED_ALLOWED = true;

		char* dptr = strchr(argv[1], '#');
		if (dptr != 0)
		{
			debug.open ("search.dbg"); 

			if (!debug)
			{
				cout << "\nError opening debug file\n";
				exit(0);
			}
			else
				debug << "\ndebug from search version " << version << endl;

			if (!debug_setup(argv[1]))
			{
				decomp_control::DEBUG = decomp_control::BASIC;
//					decomp_control::DEBUG = true;
				debug << "main: default debug options set" << endl;
			}
			
		}
	}

	if ( (!switches && argc > 2) || (switches && argc > 3) )
    {
		cout << "\n\nUsage: search [-#[1|2|3]th] [<input_file>]" << endl;
		cout << "      #[1|2|3]: generate debug output" << endl;
		cout << "      t: decompose twisted tracks t, not the untwisted 2t " << endl;
		cout << "      h: show help information" << endl;
		exit(0);
    }
}

cout << "\n               Search version " << version;
cout << "\n\nThis program provides a number of tools  that may be \n"
     << "applied to a set of fundamental tracks corresponding to a\n"
     << "finite group presentation.  The program 'solve' should have\n"
     << "been run to produce an .eft file containing the tracks before\n"
     << "this program is used.  The tracks are considered to be numbered\n"
     << "as in the corresponding .out file.\n\n";

string jobname;

if (switches && argc == 3)
    jobname = argv[2];
else
if (!switches && argc == 2)
    jobname = argv[1];
else
{
    cout << "Please enter the jobname of the group you want to analyse.\n";
    cin >> jobname;
}
   
   
/***** Open corresponding .eft file *******/
string eftfile = jobname+".eft";

ifstream eft(eftfile.c_str());
if (eft)
{

if (decomp_control::DEBUG >= decomp_control::BASIC)
	debug << "\nmain: input file: " << eftfile << "\n";

}
else
{
    cout << "\nError opening " << eftfile;
    exit(0);
}

/***********************************************************************
 "jobname.eft" holds the information created by 'solve'.

NOTE:  the format of this file has changed from the original version.
       all information has been stored in this file, instead of having
       a .eft file and a .bin file.

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
int number_of_generators;
eft >> number_of_generators;

char generator[number_of_generators];
for (int i=0; i < number_of_generators; i++)
    eft >> generator[i];

int number_of_relators;
eft >> number_of_relators;


vector<string> relator(number_of_relators);
for (int i=0; i< number_of_relators; i++)
    eft >> relator[i];
	
	
int number_of_0_cells;
eft >> number_of_0_cells;

int number_of_1_cells;
eft >> number_of_1_cells; 

int number_of_2_cells;
eft >> number_of_2_cells;

int occ_of_cv;
eft >> occ_of_cv;

int number_of_1_cells_in_rosette;
eft >> number_of_1_cells_in_rosette; 

vector<int> length_of_attaching_map(number_of_relators);

for (int i=0; i<number_of_relators; i++)
   eft >> length_of_attaching_map[i];

int max_att_map_len = length_of_attaching_map[0];
for ( int i = 1; i < number_of_relators; i++)
{
    if (length_of_attaching_map[i] > max_att_map_len)
		max_att_map_len = length_of_attaching_map[i];
}

matrix<int> attaching_map(number_of_relators, max_att_map_len);
eft >> attaching_map;

vector<char> generator_corresponding_to_1_cell(number_of_1_cells_in_rosette);
for (int i=0; i<number_of_1_cells_in_rosette; i++)
    eft >> generator_corresponding_to_1_cell[i];

int number_of_order_2_generators;
eft >> number_of_order_2_generators;

vector<char> generators_of_order_2(number_of_order_2_generators);
vector<int> occ_of_order_2_gen(number_of_order_2_generators);

if (number_of_order_2_generators > 0)
{
    for (int i=0; i< number_of_order_2_generators; i++)
		eft >> generators_of_order_2[i];

    for (int i=0; i<number_of_order_2_generators; i++)
		eft >> occ_of_order_2_gen[i];
}

int occurrences_of_1_cell[number_of_1_cells];
for (int i=0; i< number_of_1_cells; i++)
   eft >> occurrences_of_1_cell[i];

matrix<int> two_cell(number_of_2_cells,3);
eft >> two_cell;

int m;
eft >> m;

int n;
eft >> n;

matrix<int> matching_system(m,n);
eft >> matching_system;

int s;
eft >>s;

vector<ray> raylist;

/* Read the ray list into raylist */
for ( int i=0; i< s; i++)
{
	ray next_ray(n);
	for (int j=0; j<n; j++)
		eft >> next_ray.r[j];
	next_ray.index = i+1;
	raylist.push_back(next_ray);
}

if(decomp_control::DEBUG >= decomp_control::BASIC)
{
    debug << "main: read from " << eftfile << endl;
    debug << "main:   number_of_generators = " << number_of_generators << endl;
    debug << "main:   generators: ";

    for (int i=0; i<number_of_generators; i++)
		debug << generator[i] << ' ';
	debug << endl;

    debug << "main:   number_of_relators = " << number_of_relators << endl;

    debug << "main:   relator lengths: ";
    for (int i=0; i< number_of_relators; i++)
		debug << relator[i].length() << ' ';
	debug << endl;

    debug << "main:   relators:\n";
    for (int i=0; i< number_of_relators; i++)
		debug << "main:     " << relator[i] << endl;

    debug << "main:   number of 0-cells = " << number_of_0_cells << endl;
    debug << "main:   number of 1-cells = " << number_of_1_cells << endl;
    debug << "main:   number of 2-cells = " << number_of_2_cells << endl;
    debug << "main:   occ of cv = " << occ_of_cv << endl;
    debug << "main:   number of 1-cells in rosette = "  << number_of_1_cells_in_rosette << endl;

    debug << "main:   length of attaching maps: ";
    for (int i=0; i<number_of_relators; i++)
		debug << length_of_attaching_map[i] << ' ';
	debug << endl;

    debug << "main:   attaching maps:\n";
	print(attaching_map, debug, 3, "main:   ");

    debug << "main:   generator_corresponding_to_1_cell: ";
    for (int i=0; i<number_of_1_cells_in_rosette; i++)
		debug << generator_corresponding_to_1_cell[i] << ' ';
	debug << endl;

    debug << "main:   length of generators_of_order_2 = " << number_of_order_2_generators << endl;

    if (number_of_order_2_generators > 0)
    {
		debug << "main:   generators_of_order_2: ";
		for (int i=0; i < number_of_order_2_generators; i++)
			debug  << generators_of_order_2[i] << ' ';
		debug << endl;
		
		debug << "main:   occ_of_order_2_gen: ";
		for (int i=0; i<number_of_order_2_generators; i++)
	    	debug << occ_of_order_2_gen[i] <<  ' ';
	    debug << endl;
    }

    debug << "main:   occurrences_of_1_cell: ";
	for (int i=0; i< number_of_1_cells; i++)
		debug << occurrences_of_1_cell[i] << ' ';
    debug << endl;
    
	debug << "main:   two_cell\n";
	print(two_cell, debug, 3, "main:   ");

    debug << "main:   m = " << m << " n = " << n << endl;

    debug << "main:   matching_system\n";
	print(matching_system, debug, 3, "main:   ");

    debug << "main:   s = " << s << endl; 
    debug << "main: end of " << eftfile << " data" << endl;

} 

/***** Prepare the K_complex structure to hold simplicial complex parameters *****/
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

/******* Establish a file to log the output of this program ******/
string txtfile = jobname + "-search-output.txt";
ofstream txt(txtfile.c_str());

if (!txt)
{
    cout << "\nError opening " << txtfile;
    exit(0);
}

txt << "\nInformation read from " << eftfile 
    << ".\n\nGroup presentation is\n\n<";

for (int i=0; i < number_of_generators; i++)
    txt << "  " << generator[i];

txt << " :";

for (int i=0; i < number_of_relators; i++)
   txt << " " << relator[i];

txt << " >\n" << endl;


cout << "\n" << s << " fundamental tracks read from " << eftfile << endl;

char selection;
//cin >> selection;
//if (selection == 'y' || selection == 'Y')
{
	bool repeat;
    do
    {
		repeat = true;
		cout <<"\n"; // \n\t\tOPTIONS AVAILABLE\n" << endl;
		cout << "\n\tA. Evaluate the decomposition given by every track in " << eftfile << endl;
		cout << "\n\tB. Evaluate decomposition given by a single track in " << eftfile << endl;
		cout << "\n\tC. Evaluate decomposition given by a collection of tracks in " << eftfile << endl;
		cout << "\n\tD. Add tracks together" << endl;
		cout << "\n\tE. Evaluate decomposition from a given track"  << endl;
		cout << "\n\tF. Decide whether a given N-tuple determines a track" << endl;
		cout << "\n\tG. Split a given pattern into component tracks" << endl;
		cout << "\n\tH. Evaluate the equivalence classes of incompatible tracks in " << eftfile << endl;	
		cout << "\n\tI. Determine whether a pair of patterns are compatible" << endl;	
		cout << "\n\tJ. Determine the characteristics of a given track" << endl;	
		cout << "\n\tK. Determine the characteristics of every track in " << eftfile << endl;	
		cout << "\n\tL. Determine the maximal Sageev n-cubes corresponding to the tracks in " << eftfile << endl;	
		cout << "\n\tM. Display help information" << endl;	
		cout << "\n\tQ. Quit menu" << endl;
		cout << "\n\nEnter selection letter: ";			       
		cin >> selection;


		switch (selection)
		{
	    	case 'a':
	    	case 'A':    
						decompose_all(eftfile,txt,raylist,K); repeat = false; break;
	    	case 'b':
	    	case 'B':   
						single(eftfile,txt,raylist,K); break;
	    	case 'c':
	    	case 'C':   
						range(eftfile,txt,raylist,K); break;
	    	case 'd':                      
	    	case 'D':   
						add(eftfile,txt,raylist,K); break;
	    	case 'e':
	    	case 'E':   
						decompose_given(txt, K); break;
	    	case 'f':
	    	case 'F':   
						enquire(txt, K); break;
			case 'g':                      		 
	    	case 'G':
						enquire(txt, K); break; /* This is essentially the same as the previous option. */
	    	case 'h':							 
	    	case 'H':   
						incompatible_clases(eftfile,txt,raylist,K); break;
	    	case 'i':							 
	    	case 'I':   
						compatibility_check(eftfile,txt,raylist,K); break;
	    	case 'j':							 
	    	case 'J':   
						display_characteristics(eftfile,txt,raylist,K,false); break;
	    	case 'k':							 
	    	case 'K':   
						display_characteristics(eftfile,txt,raylist,K,true); break;
	    	case 'l':							 
	    	case 'L':   
						maximal_sageev_n_cubes(eftfile,txt,raylist,K); break;
	    	case 'm':							 
	    	case 'M':   
						initial_help(); break;


	    	case 'q':							 
	    	case 'Q':   
						repeat=false; break;

    	    default:    cout << "\nSelection MUST be one of the index letters of the menu\n";;
		}

		if (repeat)
		{
			cout << "\nPress m to return to the menu, q to quit\n> ";
			cin >> selection;
			if (selection == 'q' || selection =='Q')
		    	repeat = false;
		}

    } while (repeat);

	cout << "\n\nTranscript written on " << txtfile << endl;
}

txt << endl;
txt.close();
eft.close();
if (decomp_control::DEBUG >= decomp_control::BASIC)
	debug.close();
	
	return 0;
}


/* decompose_all produces decompositions for all the ext fund tracks */
void decompose_all(string eftfile, ofstream& txt, vector<ray>& raylist, K_complex K)
{
    cout << "\nProcessing " << K.s << " tracks\n";
    txt << "\nDecomposing all tracks from " << eftfile << ".\n" << endl;

	bool non_trivial_decomp_found = false;
	
    for (int i=0; i<K.s; i++)
    {
		cout << "\r" << i+1 << flush;
		if (decompose(raylist[i].r, K, false, true, txt) & DECOMPOSE_NON_TRIVIAL)
			non_trivial_decomp_found = true;
    }

    if ( !non_trivial_decomp_found )
    {
		txt << "\n\nNo non-trivial decomposition found.";
		cout << "\n\nNo non-trivial decomposition found.";
    }
}

/*********************************************************************
single evaluates the decomposition of a single track from the eft file
*********************************************************************/
void single(string eftfile, ofstream& txt, vector<ray>& raylist, K_complex K)
{
    int number;
	int s = K.s;
    char c;

    cout << "\nFor which track do you want the decomposition?\nEnter number: ";
    cin >> number; /* Tracks numbered in .eft file from 1 */

    if (number > s || number < 1 )
	{
		cout << "\nSorry, the tracks in " << eftfile << " are numbered 1 to " << s << ".";
	}
    else
    {
	
		txt << "\nTrack number " << number << " selected from " << eftfile << ".\n\n";
		decompose(raylist[number-1].r, K, true, true, txt);
    }

    cout << "\nDo you want this option again [y/n]? ";
    cin >> c;
	
    if ( c == 'y' || c == 'Y' )
	single(eftfile,txt,raylist,K);
}

/******************************************************************
 range function produces the decompositions from a collection of
 tracks, either contiguously or a general collection.
 ******************************************************************/
void range(string eftfile, ofstream& txt, vector<ray>& raylist, K_complex K)
{
    cout << "\nDo you want a contiguous range of the tracks in the .eft file or a discrete collection?";
	cout << "\nEnter c for a contiguous range or d for a discrete collection [c/d]: ";

	int s = raylist.size();
	bool contiguous;
	bool repeat;
    do
    {
		repeat = false;
		char c;
		cin >> c;
		if ( c == 'c' || c == 'C' )
			contiguous = true;
		else if ( c == 'd' || c == 'D' )
	    	contiguous = false;
		else
		{
	    	cout << "\nPlease enter c or d: ";
	    	repeat = true;
		}
    } while (repeat);

    if (contiguous)
    {
		cout << "\nTo decompose tracks row1 to row2 inclusive, enter row1 and row2 with row1 less than row2.";
		int row1;
		int row2;
		do
		{
	    	repeat = false;
	    	cout << "\nrow1: ";
	    	cin >> row1;
	    	cout << "\nrow2: ";
	    	cin >> row2;

	    	if ( row2 < row1 || row1 < 1 || row2 > s )
	    	{
				cout << "\nrow1 MUST be less than row2 and both in the range 1 to " << s;
				repeat = true;
	    	}
		} while (repeat);

		txt << "\nTracks " << row1 << " to " << row2 << "selected.";

		for ( int i = row1-1; i< row2; i++)
			decompose(raylist[i].r, K, true, true, txt);

    }
    else
    {
		cout << "\nEnter the number of tracks in your collection: ";
		int num;
		cin >> num;

		int input_col[num];
		do
		{
	    	repeat = false;
	    	cout << "\nEnter collection: ";
			for (int i=0; i< num; i++)
			{
				cin >> input_col[i];
				if (input_col[i] < 1 || input_col[i] > s)
				{
					cout <<"\nvalue entered is outside allowable range of 1 to " << s << endl;
					repeat = true;
					break;
				}
			}
			
			if (repeat)
				continue;
			else
			{
				for ( int i = 0; i < num; i++)
					decompose(raylist[input_col[i]-1].r, K,true, true, txt);
			}	
		} while (repeat);
    }
    cout << "\nDo you want this option again [y/n]? ";
	char c;
    cin >> c;
    if ( c == 'y' || c == 'Y' )
	range(eftfile,txt,raylist,K);

}

/******************************************************
 add will simply add a collection of patterns together.
 *****************************************************/
void add(string eftfile, ofstream& txt, vector<ray>& raylist, K_complex K)
{
	matrix<int>& matching_system = *K.matching_system;
	int n = matching_system.numcols();
    vector<int> sum(n);

	int s = raylist.size();
    bool include_eft = true;
	bool entering_patterns = true;	

    txt << "\n\nAdding patterns\n\n";
    cout <<"\nDo you want to include any of the tracks from " << eftfile << "? [y/n] ";
	char c;	
	cin >> c;
    if (c == 'n' || c == 'N')
		include_eft = false;

    if (include_eft)
    {
		cout << "\nEnter the number of tracks you want to include ";
		int eftnum;
		cin >> eftnum;
		int input_col[eftnum];
		bool repeat;
		
    	cout << "\nEnter the number(s) of the tracks\n";
		for (int i=0; i< eftnum; i++)
		{
			do
			{
	    		repeat = false;    			
				cin >> input_col[i];
				if (input_col[i] < 1 || input_col[i] > s)
				{
					cout << "\ntrack numbers MUST be in the range 1 to " << s;
					repeat = true;
				}
			} while(repeat);
		}
			
		txt << "\nTracks from .eft file:\n";
		for (int i = 0; i< eftnum; i++)
		{
	    	ray track = raylist[input_col[i]-1];
	    	txt << "\nTrack " << input_col[i] <<". ";
			for (int j=0; j< n; j++)
			{
	    		txt << setw(4) << track.r[j];
				sum[j] += track.r[j];
			}
				
	    	txt << endl;
		}

		cout << "\nDo you want to add enter any more patterns? [y/n] ";
		cin >> c;
		if ( c == 'n' || c == 'N' )
	    	entering_patterns = false;
    }

    if (entering_patterns)
    {
	    vector<int> pbuf(n);
		txt << "\n\nPatterns entered from keyboard:\n";
		bool repeat;
		
		do
		{
	    	repeat = true;
	    	cout << "\nEnter a " << n <<"-tuple corresponding to a pattern.\n";
			for (int i=0; i< n; i++)
		    	cin >> pbuf[i];;

	    	if (match_ok(matching_system, pbuf))
	    	{
				txt << "\n";

				for (int j=0; j< n; j++)
				{
		    		txt << setw(4) << pbuf[j];
					sum[j] += pbuf[j];
				}			

				cout << "\nAny more? [y/n] ";
				cin >> c;
				if (c == 'n' || c== 'N')
		    		repeat = false;
	    	}
	    	else
	    	{
				cout << "\nThe " << n <<"-tuple entered does not satisfy the "
		    		 << "matching equations;\nit is therefore not a "
			    	 << "pattern.\n";
	    	}
		} while(repeat);
    }

    cout << "\nSum of patterns is the pattern:";
    txt <<  "\nSum of patterns is the pattern:";

    split(sum,K,false,txt);

    cout << "\nDo you want this option again [y/n]? ";
    cin >> c;
    if ( c == 'y' || c == 'Y' )
		add(eftfile, txt, raylist, K);
}

/**************************************************
enquire reports whether a given pattern is a track; 
if not it delivers the component tracks.
***************************************************/
void enquire(ofstream& txt,K_complex K)
{
	matrix<int>& matching_system = *K.matching_system;
 	int n = matching_system.numcols();
    vector<int> pbuf(n);

    txt << "\nEnquiry regarding N-tuple.\n";

    bool repeat;
    do
    {
		repeat = false;
		cout << "\nEnter a " << n <<"-tuple you want to check.\n";
		for (int i=0; i< n; i++)
	    	cin >> pbuf[i];;

		if (match_ok(matching_system, pbuf))
		{
	    	txt << "\n" << n <<"-tuple entered is:";
	    	cout << "\n" << n <<"-tuple entered is:";
	    	split (pbuf,K,false,txt);
		}
		else
		{
	    	cout << "\nThe " << n <<"-tuple entered does not satisfy the "
			 << "matching equations;\nit is therefore not a "
			 << "pattern.\n";
	    	repeat = true;
		}
    } while(repeat);

    char c;
    cout << "\nDo you want this option again [y/n]? ";
    cin >> c;
    if ( c == 'y' || c == 'Y' )
		enquire(txt,K);
}

void decompose_given(ofstream& txt, K_complex& K)
{
	matrix<int>& matching_system = *K.matching_system;
	int n = matching_system.numcols();
    vector<int> pbuf(n);

    txt << "\nDecomposing a given track.\n";

    bool repeat;
    do
    {
		repeat = false;
		cout << "\nEnter the " << n <<"-tuple corresponding to the track "
	    	 << "you want to decompose.\n";
		for (int i=0; i< n; i++)
	    	cin >> pbuf[i];;

		if (match_ok(matching_system, pbuf))
		{
	    	if (split (pbuf,K,true,txt))
			{
				ray track(pbuf);
				decompose (track.r, K, true, true, txt);
			}
	    	else
	    	{
				cout << "\nError!  You have entered the pattern:";
				txt << "\nError!  You have entered the pattern:";
				split(pbuf,K,false,txt);
	    	}
		}
		else
		{
	    	cout << "\nThe N-tuple entered does not satisfy the "
			 << "matching equations;\nit is therefore not a "
			 << "pattern.\n";
	    	repeat = true;
		}
    } while(repeat);

    cout << "\nDo you want this option again [y/n]? ";
    char c;
    cin >> c;
    if ( c == 'y' || c == 'Y' )
		decompose_given(txt, K);
}

void compatibility_check(string eftfile, ofstream& txt, vector<ray>& raylist, K_complex K)
{

	matrix<int>& matching_system = *K.matching_system;
	int n = matching_system.numcols();
	ray raypair[2] = {raylist[0],raylist[0]}; //initialized since there is no default constructor

	int s = raylist.size();
    bool include_eft = true;
	bool entering_patterns = true;	

    txt << "\n\nChecking track compatibility\n\n";
    cout <<"\nDo you want to use a track from " << eftfile << "? [y/n] ";
	char c;	
	cin >> c;
    if (c == 'n' || c == 'N')
		include_eft = false;

	int eftnum=0;

    if (include_eft)
    {

		do 
		{
			cout << "\nEnter the number of the tracks you want to include (1 or 2): ";
			cin >> eftnum;
		} while (eftnum != 1 && eftnum != 2);
		
		int input_col[eftnum];
		bool repeat;
		
    	cout << "\nEnter the number(s) of the tracks\n";
		for (int i=0; i< eftnum; i++)
		{
			do
			{
	    		repeat = false;    			
				cin >> input_col[i];
				if (input_col[i] < 1 || input_col[i] > s)
				{
					cout << "\ntrack numbers MUST be in the range 1 to " << s;
					repeat = true;
				}
			} while(repeat);
		}
			
		txt << "\n\ntracks from eftfile:\n";

		for (int i = 0; i< eftnum; i++)
		{
	    	raypair[i] = raylist[input_col[i]-1];

	    	txt << "\nTrack " << input_col[i] <<". ";
			for (int j=0; j< n; j++)
			{
	    		txt << setw(4) << raypair[i].r[j];
			}
		}

    	if (eftnum == 2)
			entering_patterns = false;
		else
			cout << "\nEnter the " << n <<"-tuple corresponding to the second pattern.\n";		
    }
	else
	{
		cout << "\nEnter the " << n <<"-tuples corresponding to a pair of patterns.\n";
	}
	

    if (entering_patterns)
    {
	    vector<int> pbuf(n);
		txt << "\n\npatterns entered from keyboard:\n";
		
		for (int i= eftnum; i<2; i++)
		{
	    	cout << "\nPattern " << n <<"-tuple:\n";
			for (int j=0; j< n; j++)
		    	cin >> pbuf[j];

	    	if (match_ok(matching_system, pbuf))
	    	{
				txt << "\n";

				for (int j=0; j< n; j++)
		    		txt << setw(4) << pbuf[j];

				raypair[i] = ray(pbuf);
	    	}
	    	else
	    	{
				cout << "\nThe " << n <<"-tuple entered does not satisfy the "
		    		 << "matching equations;\nit is therefore not a "
			    	 << "pattern.\n";
				i--;
	    	}
		}
    }

	if (compatible(raypair[0],raypair[1],K))
	{
	    cout << "\nPatterns entered are compatible" << endl;
	    txt << "\nPatterns entered are compatible" << endl;
	}
	else
	{
	    cout << "\nPatterns entered are incompatible" << endl;
	    txt << "\nPatterns entered are incompatible" << endl;
	}


    cout << "\nDo you want this option again [y/n]? ";
    cin >> c;
    if ( c == 'y' || c == 'Y' )
		compatibility_check(eftfile,txt,raylist,K);
}

/* display_characteristics currently just identifies whether a given pattern is twisted
   and if untwisted whether it is separating or not.  To do this it calls decompose
*/
void display_characteristics(string eftfile, ofstream& txt, vector<ray>& raylist, K_complex K, bool test_all_tracks)
{

	matrix<int>& matching_system = *K.matching_system;
	int n = matching_system.numcols();
	ray track(n);

	int s = raylist.size();
    bool include_eft = true;
	bool entering_patterns = true;	

    txt << "\n\nDetermining track characteristics\n\n";

	if (test_all_tracks)
	{
   		for (int i = 0; i< s; i++)
   		{
			if (split(raylist[i].r,K,true,txt)) //it's a track
				decompose(raylist[i].r, K, false, false, txt); //no echo no full_decomp
			else
			{
				cout << "The entry (";
				for(int j=0; j< n-1; j++)
					cout << raylist[i].r[j] << ", ";       			
				cout << raylist[i].r[n-1] << ")  is not a track" << endl;

				txt << "The entry (";
				for(int j=0; j< n-1; j++)
					txt << raylist[i].r[j] << ", ";       			
				txt << raylist[i].r[n-1] << ")  is not a track" << endl;
			}
		}
 	}
	else
	{
		cout <<"\nDo you want to test a track from " << eftfile << "? [y/n] ";
		char c;	
		cin >> c;
		if (c == 'n' || c == 'N')
			include_eft = false;

        if (include_eft)
        {
        	cout << "\nEnter the number of the track: ";
			int track_num;
			bool repeat;
    		do
    		{
        		repeat = false;    			
    			cin >> track_num;
    			if (track_num < 1 || track_num > s)
    			{
    				cout << "\ntrack numbers MUST be in the range 1 to " << s;
    				repeat = true;
    			}
    		} while(repeat);
    			
    		txt << "\n\ntrack from eftfile:\n";
    
   	    	track = raylist[track_num-1];
    
   	    	txt << "\nTrack " << track_num <<". ";
   			for (int j=0; j< n; j++)
   			{
   	    		txt << setw(4) << track.r[j];
   			}
    
    		entering_patterns = false;
        }
    	else
    	{
    		cout << "\nEnter the " << n <<"-tuple corresponding to the track.\n";
    	}
    	
    
        if (entering_patterns)
        {
    	    vector<int> tbuf(n);
    		txt << "\n\ntrack entered from keyboard:\n";
    		
    		bool satisfies_matching_equations;
    		do
    		{
				satisfies_matching_equations = true;
    	    	cout << "\nTrack " << n <<"-tuple:\n";
    			for (int j=0; j< n; j++)
    		    	cin >> tbuf[j];
    
    	    	if (match_ok(matching_system, tbuf))
    	    	{
    				txt << "\n";
    
    				for (int j=0; j< n; j++)
    		    		txt << setw(4) << tbuf[j];
    
    				track = ray(tbuf);
    	    	}
    	    	else
    	    	{
    				cout << "\nThe " << n <<"-tuple entered does not satisfy the "
    		    		 << "matching equations;\nit is therefore not a "
    			    	 << "pattern.\n";
    				satisfies_matching_equations = false;
    	    	}
    		} while (!satisfies_matching_equations);
        }
		
		if (split(track.r,K,true,txt)) //it's a track
			decompose(track.r, K, false, false, txt); // no echo no full_decomp
		else
		{
			cout << "\nThe " << n <<"-tuple entered is not a track" << endl;
			txt << "\nThe " << n <<"-tuple entered is not a track" << endl;
		}

		cout << "\nDo you want this option again [y/n]? ";
		cin >> c;
		if ( c == 'y' || c == 'Y' )
			display_characteristics(eftfile,txt,raylist,K,false);
	}
}


/************************************************************
 match_ok checks whether the pattern supplied 
 satisfies the matching equations and returns true if it does
 ************************************************************/
bool match_ok(matrix<int> matching_system, vector<int>& pattern)
{
    bool matches = true;

    for (unsigned int i = 0; i< matching_system.numrows(); i++)
    {
		int test = 0;
		for (unsigned int j = 0; j<matching_system.numcols(); j++)
	    	test += matching_system[i][j] * pattern[j];

		if (test != 0)
		{
	    	matches = false;
	    	break;
		}
    }

    return (matches);
}


/***************  Functions required for setting up search specific debug **************/

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
			cout << "\t\tsearch{summary|1:basic|2:intermediate|3:detail|4:exhaustive|5}, integer: default 0=off, no parameters sets basic" << endl;
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


