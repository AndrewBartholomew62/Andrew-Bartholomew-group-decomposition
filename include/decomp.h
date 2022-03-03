
/* the following structure is used to pass a group of parameters relating to the 
   simplicial complex K to various functions without having huge parameter lists
*/
struct K_complex
{
	int number_of_0_cells;
	int number_of_1_cells;
	int number_of_2_cells;
	int occ_of_cv;
	vector<char> generator_corresponding_to_1_cell;
	int number_of_generators;
	int number_of_relators;
	vector<string> relator;
	int number_of_1_cells_in_rosette;
	vector<char> generators_of_order_2;
	vector<int> occ_of_order_2_gen;
	int* occurrences_of_1_cell;
	vector<int> length_of_attaching_map;
	matrix<int>* two_cell;
	matrix<int>* attaching_map;
	matrix<int>* matching_system;
	int s;
};

struct ray
{
	vector<int> r;
	int index;  // used to hold indices into a list of rays for output information

	ray(int n):r(n){}
	ray(vector<int> v):r(v){}
	bool operator == (const ray&);
};

struct decomp_control 
{
	enum level {OFF, SUMMARY, BASIC, INTERMEDIATE, DETAIL, EXHAUSTIVE};
	static int DEBUG;
	static bool DECOMPOSE_DEBUG;
	static bool SPLIT_DEBUG;
	static bool REDUCE_DEBUG;
};


#define NOT_FOUND -1

#define DECOMPOSE_NON_TRIVIAL 		1
#define DECOMPOSE_NON_SEPARATING	2
#define DECOMPOSE_TWISTED 			4
