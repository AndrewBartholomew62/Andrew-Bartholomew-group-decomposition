using namespace std;

#include <fstream>
#include <iostream>
#include <vector>
#include <cstring>
#include <iomanip>
#include <list>
#include <vector>

extern unsigned int	DEBUG;
extern bool		TWISTED_ALLOWED;
extern ofstream	debug;

#include <util.h>
#include <matrix.h>
#include <decomp.h>


bool compatible (ray& t1, ray& t2, K_complex K);
list<ray>::iterator find_cpt(list<ray>& cpt_list, ray cpt);
void find(int label,int number_of_1_cells_in_rosette, matrix<int>** labels_for_1_cell, vector<int>& position);
int bvertex(int label, int number_of_0_cells, int numbcols, vector<int>** blabels_for_0_cell);


bool ray::operator == (const ray& R)
{
	bool result = true;
	
	if (r.size() != R.r.size())
		result = false;
	else
	{
		for (unsigned int i; i < r.size(); i++)
		{
			if (r[i] != R.r[i])
			{
				result = false;
				break;
			}
		}
	}
	
	return result;
}

void initial_help()
{
	cout << "\n\nUsage: search [-t][-h] [<infile>]\n\n";
	cout << "This program may be run be typing 'search' at the command line\n"
	     << "and following the instructions given, or the stub of an .eft file may\n"
	     << "be entered, as in 'search fred', where fred.eft must exist in the\n"
	     << "current directory (fred.eft created by running 'solve').\n"
	     << "Optionally, the switch '-t' may be used to tell the program not to\n"
	     << "double up twisted tracks and decompose the double track automatically.\n"
	     << "The -h option simply repeats this information.\n" << endl;
}


/* incompatible_classes determines the equivalence classes of the relation determined by incompatibility.
   Write s~~t if tracks s and t are incompatible, then define an equivalence relation ~ by 
   u~v if there is a sequence of tracks such that u~~t_1~~ ... ~~t_n~~v.
   
   To evaluate these equivalence clases we use a tree based approach, constructing a tree for
   each class.  For each set of tree vertices k edges from the root we check whether there are any
   tracks remaining in the list that are incompatible with a tree vertex, and if so add that new track 
   to the tree joining it by an edge to an incompatible track.  Thus at each stage the tracks added to 
   the tree are k+1 edges from the root.
   
   Since we do not know how many classes there are, and we need to know which extreme fundamental tracks
   have been assigned to a class, we set up an array of s (number of extreme fundamental tracks) integers, 
   initialized to 0 and set to i if extreme fundamental track i is in the ith class we find.
   
   The tree we create for each class is a tree of track indices in the raylist, though we shall refer to a 
   vertex of the tree as being a track for the purposes of description.  In fact, we do not need to record
   the full tree structure, only the level of each vertex to aid our breadth-first search.  We therefore, 
   store the tree as an array of s integers, initialized to -1, recording the level of each track.  The root's 
   level is set to zero and we may therefore use non-negative levels to determine tree membership.
*/
void incompatible_clases(string eftfile, ofstream& txt, vector<ray>& raylist, K_complex K)
{

int s = raylist.size();
int equivalence_class[s];

for (int i=0; i< s; i++)
	equivalence_class[i] = 0;

int eq_class_start = 0; // we shall use this variable to indicate the first track not yet assigned to a class
int eq_class = 0; // used to count classes and flag tracks in equivalence_class


cout << "\nProcessing " << s << " extreme fundamental tracks from " << eftfile << endl;

if (decomp_control::DEBUG >= decomp_control::INTERMEDIATE)
	debug << "main: evaluating incompatibility classes:" << endl;


int tracknum = 0;

bool equivalence_classes_incomplete;

do
{
	equivalence_classes_incomplete = false;
	
	/* look for another equivalence class */
	for (int i = eq_class_start; i< s; i++)
	{
		if (equivalence_class[i] != 0)
		{
			eq_class_start++; // don't need to check this track again
		}
		else		
		{
			equivalence_classes_incomplete = true;
			break;
		}		
	}
	
	if (equivalence_classes_incomplete)
	{

if (decomp_control::DEBUG >= decomp_control::INTERMEDIATE)
	debug << "main: found equivalence class with eq_class_start = " << eq_class_start << endl;
	
		eq_class++; // there's another class
		
		/* build the next equivalence class */
		int tree[s];
		for (int i=0; i< s; i++)
			tree[i] = -1;

		/* the root of this equivalence class tree will be the track
		   referenced by eq_class_start */
		tree[eq_class_start] = 0;  // note it's parent is left as -1
		 
		/* this track is in the new equivalence class */
		equivalence_class[eq_class_start] = eq_class;		

		if (s > 100)
			cout << "\r" << ++tracknum << flush;

		bool incompatible_tracks_found;
		int level = -1;
		do
		{
			incompatible_tracks_found = false;
			level++;
					
			/* look for all tree vertices at the current level, and for each one
			   check whether there are any incompatible tracks not already assigned
			   to an equivalence class */
			for (int i=0; i < s; i++)
			{
				if (tree[i] == level)
				{
					/* this track is in the tree and at the right level, check 
					   the remaining tracks indicated by equivalence_class.
					*/
					for (int j= eq_class_start+1; j < s; j++)
					{						
						if (equivalence_class[j] == 0 && !compatible(raylist[i], raylist[j], K))
						{					
							if (s > 100)
								cout << "\r" << ++tracknum << flush;
							
							/* add track j to this class */
							equivalence_class[j] = eq_class;
							
							/* record the level of track j in the tree */
							tree[j] = level+1;
if (decomp_control::DEBUG >= decomp_control::INTERMEDIATE)
{
	debug << "main:   adding track " << j << " to class at level " << level+1 << " via track " << i << "\ntree:\n\t";
	for (int i=0; i< s; i++)
		debug << setw(3) << tree[i];
	debug << "main: equivalence_class: ";
	for (int i=0; i< s; i++)
		debug << setw(3) << equivalence_class[i];
	debug << endl;
}
							incompatible_tracks_found = true;
						}
					}
				}
			}
		} while (incompatible_tracks_found);

if (decomp_control::DEBUG >= decomp_control::BASIC)
{
	debug << "main: equivalence classes of incompatible tracks: " << endl;
	for (int i=0; i< s; i++)
		debug << equivalence_class[i] << ' ';
	debug << endl;
}

	}
} while(equivalence_classes_incomplete);

cout << "\nThere " << (eq_class == 1? "is ": "are ") << eq_class << " equivalence " << (eq_class == 1? "class": "classes")
	 << " of incompatible tracks" << endl;
txt << "\nThere " << (eq_class == 1? "is ": "are ") << eq_class << " equivalence " << (eq_class == 1? "class": "classes")
	 << " of incompatible tracks" << endl;
}

/*
A Sageev n-cube corresponds to a set of vertex solutions (extreme fundamental tracks) 
any pair of which are incompatible.  To find the maximal n-cubes we take the first track
in raylist then find the first subsequent track with which it is incompatible. Then look 
for one that is incompatible with those two and then carry on until one gets to the end of the list.

This produces a maximal sized cube containing the first track as a vertex.

We now repeat this process but sequentially look for the second, then third track incompatible with 
the first, and so on.  Each time having found the kth incompatible track we proceed to build up a list of
tracks that are pairwise incompatible.

Each such list must be checked to see if it is a sub-cube of the first one, if not it corresponds to a 
different maximal n-cube.
*/
void maximal_sageev_n_cubes(string eftfile, ofstream& txt, vector<ray>& raylist, K_complex K)
{
	/* each maximal cube will be a list of vertex solutions but we don't know how many 
	   there will be, so we use a list of lists.
	*/
	list<list<ray> > maximal_cubes;
	list<ray> candidate; // a candidate maximal cube
	
	bool incompatible_track_found;
	
	/* we'll execute the following look for each maximal n-cube but we don't want to check 
	   tracks against the first track repeatedly, so we maintain a starting offset for 
	   checking raylist, this will automatically ensure we're skipping the correct number 
	   of tracks.
	 */
	unsigned int raylist_starting_offset = 1;
	do
	{
if (decomp_control::DEBUG >= decomp_control::BASIC)
	debug << "maximal_sageev_n_cubes: maximal_cubes.size()=" << maximal_cubes.size() << endl;

if (decomp_control::DEBUG >= decomp_control::BASIC)
	debug << "maximal_sageev_n_cubes: starting to look for maximal cube from raylist offset " << raylist_starting_offset << endl;
	
		incompatible_track_found = false;
		candidate.clear();
		candidate.push_back(raylist[0]);
			
		for (unsigned int i=raylist_starting_offset; i< raylist.size(); i++)
		{
			if (!compatible(raylist[i], raylist[0], K))
			{
				incompatible_track_found = true;
				candidate.push_back(raylist[i]);
				raylist_starting_offset = i;
if (decomp_control::DEBUG >= decomp_control::BASIC)
	debug << "maximal_sageev_n_cubes: incompatible track found at raylist offset " << i << endl;

				break;
			}
		}
		
		if (incompatible_track_found)
		{
			/* start looking for pairwise incompatible tracks from raylist offset 1,
			   skipping over the first incompatible track we've found (identifiable from 
			   raylist_starting_offset).
			*/
			for (unsigned int i = 1; i< raylist.size(); i++)
			{
				if (i != raylist_starting_offset)
				{
					/* if raylist[i] is pairwise compatible with the candidate list, 
					   add it to that list
					*/
					list<ray>::iterator rptr = candidate.begin();
					bool pairwise_incompatible = true;
					while (rptr != candidate.end())
					{
						if (compatible(raylist[i],*rptr, K))
						{
							pairwise_incompatible = false;
							break;
						}
						rptr++;
					}
				
					if (pairwise_incompatible)
					{
						candidate.push_back(raylist[i]);
if (decomp_control::DEBUG >= decomp_control::BASIC)
	debug << "maximal_sageev_n_cubes: added raylist offset " << i << " to candidate n-cube" << endl;
					}
				}
			}
			
			/* candidate now holds a candidate maximal cube, if this is not the first
			   maximal cube, check it's not a sub-cube of the first one.  The boolean
			   sub_cube is initialized false and set true if the maximal_cubes list is non-empty
			   so that sub-cube is false for the first maximal n-cube we find.
			*/
			bool sub_cube = false;
			
			if (maximal_cubes.size() !=0)
			{
				sub_cube = true;
				list<ray>::iterator rptr = candidate.begin();
				while (rptr != candidate.end())
				{
					/* check to see if *rptr is in the first maximal n-cube 
					   we can use find_cpt for this purpose
					*/
					list<ray>& first_maximal_cube = *maximal_cubes.begin();
					if (find_cpt(first_maximal_cube,*rptr) == first_maximal_cube.end())
					{
						/* *rptr is not a vertex of the first maximal subcube therefore
						   candidate is not a sub-cube of the first maximal subcube */
						sub_cube = false;
						break;
					}
					rptr++;
				}
			}
			
			if (!sub_cube)
			{
if (decomp_control::DEBUG >= decomp_control::BASIC)
{
	if (maximal_cubes.size() != 0)
		debug << "maximal_sageev_n_cubes: candidate is not a sub-cube of the first maximal n-cube, adding to maximal_cubes list" << endl;
	else
		debug << "maximal_sageev_n_cubes: adding candidate to maximal_cubes list as the first maximal n-cube" << endl;
}
				maximal_cubes.push_back(candidate);
				cout << "\rFound " << maximal_cubes.size() << " maximal Sageev n-cubes" << flush;
			}
			else
			{
if (decomp_control::DEBUG >= decomp_control::BASIC)
	debug << "maximal_sageev_n_cubes: candidate is a sub-cube of the first maximal n-cube, discarding candidate" << endl;
			}
			
			raylist_starting_offset++;  //ready for the next candidate
		}
		
	} while (incompatible_track_found);
	
	txt << "\nFound " << maximal_cubes.size() << " maximal Sageev n-cubes\n" << endl;
	list<list<ray> >::iterator cptr = maximal_cubes.begin();
	int cube_count = 0;
	
	while (cptr != maximal_cubes.end())
	{
		++cube_count;
//		cout << "\ncube " << cube_count << " dimension " << (*cptr).size() << endl;
		txt << "\ncube " << cube_count << " dimension " << (*cptr).size() << endl;
		list<ray>::iterator rptr = cptr->begin();
		while (rptr != cptr->end())
		{
//				cout << setw(4) << rptr->index << ". ";
				txt << setw(4) << rptr->index << ". ";
				for (unsigned int i=0; i< rptr->r.size(); i++)
				{
//					cout << rptr->r[i] << ' ';
					txt << rptr->r[i] << ' ';
				}
//				cout << endl;
				txt << endl;
				rptr++;
		}
		cptr++;
	}
}

