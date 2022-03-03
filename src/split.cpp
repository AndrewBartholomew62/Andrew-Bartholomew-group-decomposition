using namespace std;

#include <fstream>
#include <iostream>
//#include <cstdlib>
//#include <list>
//#include <cstring>
//#include <ctype.h>
//#include <stdio.h>
#include <iomanip>
#include <vector>

//#define		MAX_NUM_GENERATORS	10
//#define		MAX_NUM_RELATORS    10


extern unsigned int DEBUG;
//bool		TWISTED_ALLOWED = false;

extern ofstream	debug;

//#include <util.h>
#include <debug.h>
#include <matrix.h>
//#include <ray_struct.h>
#include <K_complex.h>


/*************************************************************************
This procedure decomposes the pattern P (:=P('pattern')) into its connected
components and returns the N-tuples corresponding to each component.  The
strategy for the decomposition is as follows.  For each relator disc we
label the points of intersection of P with the 1-skeleton of the disc with
consecutive integers (without repetitions) starting at O and labelling
clockwise around the boundary, then for each interior 1-cell ("clockwise")
starting at O.  We then use these labels to build up equivalence classes of
labels on each disc for the relator lbl1~lbl2 iff lbl1 and lbl2 are joined
by an edge of P.  This enables us to discover how the pattern crosses each
disc. We then build classes of "end" labels of the arcs across each disc by
considering how the discs are identified under the attaching maps (we shall
have one class for each component of P), finally returning to the crossing
arcs to determine the tuples corresponding to each component.
************************************************************************/
bool split (valarray<int>& pattern, K_complex K, bool short_split, ofstream& txt, bool return_components, matrix<int>** component_ptr)     
{

if (DEBUG >= INTERMEDIATE)
{
	debug << "split: pattern provided: ";
	for (unsigned int i=0; i<pattern.size(); i++)
		debug << pattern[i] << ' ';
	debug << endl;
}

    /**********  read variables required by split from the K_complex **********/
	int number_of_1_cells = K.number_of_1_cells;
	int number_of_2_cells = K.number_of_2_cells;
	int number_of_relators = K.number_of_relators;
	int number_of_1_cells_in_rosette = K.number_of_1_cells_in_rosette;
	int* occurrences_of_1_cell = K.occurrences_of_1_cell;
	int* length_of_attaching_map = K.length_of_attaching_map;
	matrix<int>& two_cell = *K.two_cell;
	matrix<int>& attaching_map = *K.attaching_map;
	


    /***********************************************************************
    For the 1-cells of the rosette, the total number of occurrences of 1-cell
    i in 'two cell' is the number of occurrences of i in the boundary of some
    relator disc.  We therefore set out, for each such 1-cell, a matrix with
    r rows and c columns, where r is the number of occurrences of the 1-cell
    in the boundary of the relator discs and c is the number of points of P
    intersecting that 1-cell.  For the interior 1-cells, which occur exactly
    once in a unique relator disc, we set out a 1xc matrix, where c is
    as above.  We store in these matrices the labels attached to the 1-cells
    corresponding to their rows (in accordance with the canonical order
    through the relator discs, ignoring identifications) in such a way that
    the columns of each matrix corresponding to a 1-cell in the rosette
    determines how those labels become identified in L.
    ***********************************************************************/

    matrix<int>* labels_for_1_cell[number_of_1_cells];
    /* labels_for_1_cell[i] will reference the "label matrix" 
       for the ith 1-cell.
    */

    int n_pts_on_1_cell[number_of_1_cells];
    /* "number of points on 1-cell" will be used later. */

    /*  First deal with the 1-cells in the rosette. */
    for (int i = 0; i< number_of_1_cells_in_rosette; i++)
    {
		bool one_cell_found = false;

		for (int j =0; j< number_of_2_cells; j++)
		{       
	    	for ( int k = 0; k<3; k++)
	    	{
				if ( abs(two_cell[j][k]) == i+1)
				{
		    		if (k == 2)	    
						n_pts_on_1_cell[i] = pattern[3*j+k] + pattern[3*j];
		    		else 
						n_pts_on_1_cell[i] = pattern[3*j+k] + pattern[3*j+k+1];

		    		one_cell_found = true;
		    		break;
				}  	     
	    	} 
			
	    	if (one_cell_found)
			break;
		}

		labels_for_1_cell[i] = new matrix<int> (occurrences_of_1_cell[i],n_pts_on_1_cell[i]);	      
    }

    /* Now deal with the interior 1-cells, for which purpose we shall need a
    different marker for 'labels for 1 cell'.  We need to know the first
    variable on each relator disc, so evaluate these first. */

    int start[number_of_relators];
    start[0] = 0;   

    for (int i = 1; i< number_of_relators; i++)
		start[i] = start[i-1] + 3 * (length_of_attaching_map[i-1] - 2);  

    int interior_1_cell = number_of_1_cells_in_rosette;
    for (int i =0; i< number_of_relators; i++)
    {
		int s = start[i];
		for (int j = s + 3 ; j < s + 3 * (length_of_attaching_map[i] - 2) - 2 ; j += 3 )
		{
	    	n_pts_on_1_cell[interior_1_cell] = pattern[j] + pattern[j+1];

		   labels_for_1_cell[interior_1_cell] = new matrix<int> (1,n_pts_on_1_cell[interior_1_cell]);   
		   interior_1_cell ++;
		}
    }

if (DEBUG >= INTERMEDIATE)
{
	debug << "split: number of points on 1-cells: ";
	for (int i=0; i<number_of_1_cells; i++)
		debug << n_pts_on_1_cell[i] << ' ';
	debug << endl;
}

    /********************************************************************
    Assign the labels to the space just allocated.  We take the convention
    that, for the 1-cells in the rosette, the orientation of the labels in
    the label matrix is left to right if the 1-cell is "positive in the
    attaching map" and right to left otherwise.  This ensures the
    identification matching mentioned above.  We also record maximum and
    minimum label on each disc and the maximum label on the boundary of
    each disc, all for future use.
    ********************************************************************/
    int marker[number_of_1_cells_in_rosette];
    /* These markers will mark the rows in the label matrices. */

    for (int i=0; i< number_of_1_cells_in_rosette; i++)
		marker[i] = 0;

    int label = 1;
    /* 'one cell' will denote a 1-cell of the rosette; label will 
      be used to allocate labels. */

    interior_1_cell = number_of_1_cells_in_rosette;
	
    int minlb[number_of_relators];
    int maxlb[number_of_relators];   
    int mbl[number_of_relators];
    /*"minimum label","maximum label","maximum boundary label"*/

    for ( int i = 0; i < number_of_relators ; i++)
    {

		/* First assign to 'minlb[i]'. */
		minlb[i] = label;
		for (int j = 0; j < length_of_attaching_map[i]; j++)
		{
	    	int one_cell = attaching_map[i][j];
	    	if ( one_cell > 0 )
	    	{
				Slice_iter<int> eltptr = labels_for_1_cell[one_cell-1]->row(marker[one_cell-1]);
								 
				for (int k=0; k< n_pts_on_1_cell[one_cell-1]; k++)
		    		*eltptr++ = label++;

				marker[one_cell-1]++;
	    	}
	    	else
	    	{
				one_cell *= -1;
				
				Slice_iter<int> eltptr = labels_for_1_cell[one_cell-1]->row(marker[one_cell-1]);

				/* new code to allocate labels in reverse order 18/2/06 */
				label += n_pts_on_1_cell[one_cell-1]-1;
				
				for (int k=0; k< n_pts_on_1_cell[one_cell-1]; k++)
	    			*eltptr++ = label--;
				
				label += n_pts_on_1_cell[one_cell-1]+1;
				marker[one_cell-1]++;
	    	}
		}

		/* Now assign the maximum label on the boundary of
		 this disc to 'mbl[i]'.*/
		mbl[i] = label-1;

		for (int j=0; j< length_of_attaching_map[i]-3; j++)
		{
			Slice_iter<int> eltptr = labels_for_1_cell[interior_1_cell]->row(0);

	    	for( int k =0; k< n_pts_on_1_cell[interior_1_cell]; k++)
			*eltptr++ = label++;
			
	    	interior_1_cell++;
		}

		/* Now assign to 'maxlb[i]'. */
		maxlb[i] = label-1;

    }

    if (DEBUG >= INTERMEDIATE)
    {
		debug << "split: minlb = ";
		for (int i=0; i< number_of_relators; i++)
			debug << setw(4) <<minlb[i] << ' ';
		debug << endl;
			
		debug << "split: maxlb = ";
		for (int i=0; i< number_of_relators; i++)
			debug << setw(4) <<maxlb[i] << ' ';
		debug << endl;

		debug << "split: mbl = ";
		for (int i=0; i< number_of_relators; i++)
			debug << setw(4) <<mbl[i] << ' ';
		debug << endl;

		debug << "split: label matracies\n";
		for (int i = 0; i< number_of_1_cells; i++)
		{
	    	debug << "split: one cell " << i+1 << endl;

	    	debug << *labels_for_1_cell[i];
	    	debug << endl;
		}

    }

    /*****************************************************************
    We next determine which labels are joined by an edge of P on each
    relator disc; we store the information as pairs, forming the rows of a
    matrix.  We have one pair for each edge of P.
                                                    
         Given a variable and the two sets of labels the edges of P
    determined by that variable "join", we make the pairing by moving "into
    the vertex", ie if 'pattern[variable]'=x we must have x labels in each
    list paired by the corresponding edges of t, thus the pairings must be of
    the x labels closest to that vertex on the relator disc.

         Each variable lies in a two cell of one of three types, A, B, C.
    Type A 2-cells are those that appear first in a relator disc (with the
    canonical ordering), type B 2-cells are the remaining 2-cells of that
    disc, except the last, and type C 2-cells are those that appear last in a
    relator disc.  Within any 2-cell a variable appears in place 1, 2, or 3,
    determined, again, by the canonical order.
   
    We store the pairs corresponding to each relator disc separately, and
    with each pair shall store the corresponding variable,
    to assist later calculations.
    ********************************************************************/

    matrix<int>* pair[number_of_relators];
    int partial_sum[number_of_relators];

    /* 'partial sum' will indicate how many 
       pairs arise from each relator disc. */

    for (int i = 0 ; i < number_of_relators ; i++)
		partial_sum[i] = 0;
   
    /* Assign partial sums. */
    int variable = 0;

    for ( int i = 0; i < number_of_relators ; i++)
	for (int j = 0; j < 3*(length_of_attaching_map[i]-2) ; j++)
	    partial_sum[i] += pattern[variable++];
   
    /* Declare space for 'pair'. */
    for (int i=0; i<number_of_relators; i++)
		pair[i] = new matrix<int>(partial_sum[i], 3);   

    /****************************************************************
    We traverse the relator discs in the canonical order, thereby moving
    through 'pattern'.  We make use of 'marker' to determine which set of
    labels we should consider for the boundary 1-cells and 'interior 1 cell'
    to take care of the rest.
    ***************************************************************/

    for (int i=0; i< number_of_1_cells_in_rosette; i++)
		marker[i] = 0;

    interior_1_cell = number_of_1_cells_in_rosette;
    variable = 0; /* Reset 'variable' for assignation to 'pair'.*/
       
    for ( int i = 0; i< number_of_relators; i++)
    {
		int pair_marker = 0; /* reset pair_marker for next relator disk */

		if (length_of_attaching_map[i] == 3)
		{
	    	/* The relator must be of the form 'aaa' for some generator. */
	    	int one_cell = attaching_map[i][0]; 

	    	int mark = marker[one_cell-1];
	    	int pv = pattern[variable];      

	    	Slice_iter<int> r1 = labels_for_1_cell[one_cell-1]->row(mark);
	    	Slice_iter<int> r2 = labels_for_1_cell[one_cell-1]->row(mark+1);
	    	Slice_iter<int> r3 = labels_for_1_cell[one_cell-1]->row(mark+2);

	    	matrix<int>& matrixptr = *pair[i]; /* matrixptr now points at the target
						 of the asignation */
	    	/* Place 1. */
	    	r1 += pv-1;
	    	r3 += pv;

	    	for (int j=0; j< pv; j++)
			{
				matrixptr[pair_marker][0] = *r1--;
				matrixptr[pair_marker][1] = *r3++;
				matrixptr[pair_marker][2] = variable;
				pair_marker++;
			}

            variable++;

	    	/* Place 2. */
	    	r1 += pv+1;
	    	r2 += pv-1;

	    	for (int j=0; j< pv; j++)
			{
				matrixptr[pair_marker][0] = *r1++;
				matrixptr[pair_marker][1] = *r2--;
				matrixptr[pair_marker][2] = variable;
				pair_marker++;
			}

	    	variable++;

	    	/* Place 3. */
	    	r2 += pv+1;
	    	r3 -= pv+1; /* i.e. the pv th label in the row */

	    	for (int j=0; j< pv; j++)
			{
				matrixptr[pair_marker][0] = *r2++;
				matrixptr[pair_marker][1] = *r3--;
				matrixptr[pair_marker][2] = variable;
				pair_marker++;
			}

	    	variable++;

	    	marker[one_cell-1] += 3;
		}
		else
		{
	    	matrix<int>& matrixptr = *pair[i]; /* matrixptr now points at the target
						 of the asignation */
						 
	    	/* Type A 2-cell. */
	    	int one_cell1 = attaching_map[i][0];
	    	int one_cell2 = attaching_map[i][1];
			bool s1;
			bool s2;

	    	if (one_cell1 < 0 )
	    	{
				s1 = false;
				one_cell1 *= -1;
	    	}
	    	else
				s1 = true;

	    	if (one_cell2 < 0 )
	    	{
				s2 = false;
				one_cell2 *= -1;
	    	}
	    	else
				s2 = true;

	    	Slice_iter<int> r1 = labels_for_1_cell[one_cell1-1]->row(marker[one_cell1-1]);

			Slice_iter<int> r2 = r1; // assignation just for the declaration, since there's no default constructor 
			
	    	if (one_cell1 == one_cell2)
				r2 = labels_for_1_cell[one_cell1-1]->row(marker[one_cell1-1]+1); 						
	    	else
				r2 = labels_for_1_cell[one_cell2-1]->row(marker[one_cell2-1]); 

	    	Slice_iter<int> r3 = labels_for_1_cell[interior_1_cell]->row(0); 

	    	int n1 = labels_for_1_cell[one_cell1-1]->numcols();
	    	int n2 = labels_for_1_cell[one_cell2-1]->numcols();
	    	int n3 = labels_for_1_cell[interior_1_cell]->numcols();    

	    	/* Place 1. */
	    	int pv = pattern[variable];      

	    	if ( s1 ) r1 += pv-1; else r1 += n1-pv;
	    	r3 += pv-1;

	    	for (int j=0; j< pv; j++)
	    	{

				matrixptr[pair_marker][0] = *r1;
				matrixptr[pair_marker][1] = *r3--;
				matrixptr[pair_marker][2] = variable;
				pair_marker++;
				
				if ( s1 ) r1--; else r1++;
	    	}
	    	variable++;

	    	/* Place 2. */
	    	pv = pattern[variable];      

	    	if ( s1 ) r1 += n1-pv+1; else r1 -= n1-pv+1;
	    	if ( s2 ) r2 += pv-1; else r2 += n2-pv;

	    	for (int j=0; j< pv; j++)
	    	{
				matrixptr[pair_marker][0] = *r1;
				matrixptr[pair_marker][1] = *r2;
				matrixptr[pair_marker][2] = variable;
				pair_marker++;

				if ( s1 ) r1++; else r1--;
				if ( s2 ) r2--; else r2++;
	    	}
	    	variable++;

	    	/* Place 3. */
	    	pv = pattern[variable];      

	    	if ( s2 ) r2 += n2-pv+1; else r2 -= n2-pv+1;
	    	r3 += n3-pv+1; 

	    	for (int j=0; j< pv; j++)
	    	{
				matrixptr[pair_marker][0] = *r2;
				matrixptr[pair_marker][1] = *r3++;
				matrixptr[pair_marker][2] = variable;
				pair_marker++;

				if ( s2 ) r2++; else r2--;
	    	}
	    	variable++;

	    	marker[one_cell1-1] += 1;
	    	marker[one_cell2-1] += 1;
	    	/* Notice that if one_cell1=one_cell2 then this code is
	    	 still correct. */

	    	/* Type B 2-cells. */
	    	int a_map_place = 2;
	    	for ( int j =0; j< length_of_attaching_map[i]-4; j++)
	    	{
				int in_rosette = attaching_map[i][a_map_place];

				if (in_rosette < 0 )
				{
		    		s2 = false;
		    		in_rosette *= -1;
				}
				else
		    		s2 = true;

				Slice_iter<int> r1 = labels_for_1_cell[interior_1_cell]->row(0); 
				Slice_iter<int> r2 = labels_for_1_cell[in_rosette-1]->row(marker[in_rosette-1]);
				Slice_iter<int> r3 = labels_for_1_cell[interior_1_cell+1]->row(0); 
				
				int n1 = labels_for_1_cell[interior_1_cell]->numcols();
				int n2 = labels_for_1_cell[in_rosette-1]->numcols();
				int n3 = labels_for_1_cell[interior_1_cell+1]->numcols();

				/* Place 1. */
				pv = pattern[variable];      
				r1 += pv-1;
				r3 += pv-1;
				
				for (int k=0; k< pv; k++)
				{
					matrixptr[pair_marker][0] = *r1--;
					matrixptr[pair_marker][1] = *r3--;
					matrixptr[pair_marker][2] = variable;
					pair_marker++;
				}
					
				variable++;

				/* Place 2. */
				pv = pattern[variable];      
				r1 += n1-pv+1;
				if ( s2 ) r2 += pv-1; else r2 += n2-pv;

				for (int k=0; k< pv; k++)		
				{
					matrixptr[pair_marker][0] = *r1++;
					matrixptr[pair_marker][1] = *r2;
					matrixptr[pair_marker][2] = variable;
					pair_marker++;

		    		if (s2) r2--; else r2++;
				}
				variable++;

				/* Place 3. */
				pv = pattern[variable];      
				if (s2) r2 += n2-pv+1; else r2 -= n2-pv+1;
				r3 += n3-pv+1;

				for (int k=0; k< pv; k++)		
				{
					matrixptr[pair_marker][0] = *r2;
					matrixptr[pair_marker][1] = *r3++;
					matrixptr[pair_marker][2] = variable;
					pair_marker++;
		    		if (s2) r2++; else r2--;
				}
				variable++;
				interior_1_cell++;
				marker[in_rosette-1] += 1;
				a_map_place ++;
	    	}

	    	/* Type C 2-cell. */
	    	one_cell1 = attaching_map[i][a_map_place];
	    	one_cell2 = attaching_map[i][a_map_place+1];

			bool s3;
	    	if (one_cell1 < 0 )
	    	{
				s2 = false;
				one_cell1 *= -1;
	    	}
	    	else
				s2 = true;

	    	if (one_cell2 < 0 )
	    	{
				s3 = false;
				one_cell2 *= -1;
	    	}
	    	else
				s3 = true;

	    	r1 = labels_for_1_cell[interior_1_cell]->row(0); 
	    	r2 = labels_for_1_cell[one_cell1-1]->row(marker[one_cell1-1]); 
		
	    	if (one_cell1 == one_cell2)
				r3 = labels_for_1_cell[one_cell1-1]->row(marker[one_cell1-1]+1); 
	    	else
				r3 = labels_for_1_cell[one_cell2-1]->row(marker[one_cell2-1]); 

           	n1 = labels_for_1_cell[interior_1_cell]->numcols();
	    	n2 = labels_for_1_cell[one_cell1-1]->numcols();
	    	n3 = labels_for_1_cell[one_cell2-1]->numcols();

	    	/* Place 1. */
	    	pv = pattern[variable];          
	    	r1 += pv-1;
	    	if ( s3 ) r3 += n3-pv; else r3 += pv-1;

	    	for (int j=0; j< pv; j++)
	    	{
				matrixptr[pair_marker][0] = *r1--;
				matrixptr[pair_marker][1] = *r3;
				matrixptr[pair_marker][2] = variable;
				pair_marker++;
				if ( s3 ) r3++; else r3--;
	    	}
	    	variable++;

	    	/* Place 2. */
	    	pv = pattern[variable];             
	    	r1 += n1-pv+1;
	    	if (s2) r2 += pv-1; else r2 += n2-pv;

	    	for (int j=0; j< pv; j++)
	    	{
				matrixptr[pair_marker][0] = *r1++;
				matrixptr[pair_marker][1] = *r2;
				matrixptr[pair_marker][2] = variable;
				pair_marker++;
				if ( s2 ) r2--; else r2++;
	    	}
	    	variable++;

	    	/* Place 3. */
	    	pv = pattern[variable];             
	    	if (s2) r2 += n2-pv+1; else r2 -= n2-pv+1;
	    	if (s3) r3 -= n3-pv+1; else r3 += n3-pv+1;

	    	for (int j=0; j< pv; j++)
	    	{
				matrixptr[pair_marker][0] = *r2;
				matrixptr[pair_marker][1] = *r3;
				matrixptr[pair_marker][2] = variable;
				pair_marker++;
				if (s2) r2++; else r2--;
				if (s3) r3--; else r3++;
	    	}
	    	variable++;
	    	interior_1_cell += 1;
	    	marker[one_cell1-1] += 1;
	    	marker[one_cell2-1] += 1;
	    	/* Again, this code is correct for all cases. */
		}
    }


if (DEBUG >= INTERMEDIATE)
{
	debug << "split: pair matracies\n";
	for (int i=0; i< number_of_relators; i++)
	{
	    debug << "split: relator " << i+1;
		debug << *pair[i];
	    debug << endl;
	}
}

    /*******************************************************************
    We now determine the component arcs on each relator disc.
   
    Declare space for the disc permutations.  For each disc, i, we set up
    a row with maxlb[i]-minlb[i] elements.  Note that each label occurs
    in "pair[i]" exactly once iff the label is on the boundary of the disc 
    and exactly twice iff it is on an interior 1-cell.
    *********************************************************************/

    valarray<int>* dperm[number_of_relators]; /*"disc permutations" */
    for( int i =0; i< number_of_relators ; i++)
		dperm[i] = new valarray<int> (maxlb[i] - minlb[i] + 1);
	
    /*********************************************************************
    Assign disc permutations.  By moving successively through the rows of 
    "pair[i]" we deal with the "variable corners" in the
    canonical order and, by nature of our construction of "pair[i]", have the 
    following properties satisfied by each row:
   
    1.  If a label in the first column is <='mbl[i]', it is the initial
	label of a component arc on the disc.  If it is
	>'mbl[i]', the pair in that row determines an edge of P whose
	initial point lies on an interior 1-cell.  In this case the label
	at this interior point will have occured exactly once, in the
	second column of a row of "pair[i]", BEFORE the current pair.
   
    2.  If a label in the second column is <='mbl[i]', it is the terminal 
	point of a component arc on the disc.  If it is >'mbl[i]', the 
	pair in that row determines an edge of P whose terminal point
	lies on an interior 1-cell.  In this case the label
	at this terminal point will occur exactly once, in the first
	column of a row of "pair[i]", AFTER the current pair.
   
        The disc permutation will allow us to find the component arcs.  A
    label in the first column of "pair[i]" determines the initial point of
    some edge of P and therefore shall be made to index the label at the
    terminal point of that edge in the permutation.  If a label in the second
    column is <='mbl[i]' it terminates a component arc and, so, shall be made
    to index zero, to signify this.  If the label in the second column is
    >'mbl[i]' it will occur, also, in the first column in a subsequent row and
    so no assignation will be made at this first encounter.  By properties 1
    and 2 above, it can be seen that this procedure completely determines the 
    permutation.
    *********************************************************************/

    for (int i=0;i<number_of_relators; i++)
    {
		matrix<int>& pairs = *pair[i];
		valarray<int>& perm = *dperm[i];
		int min = minlb[i];
		int bdry_max = mbl[i];

		for ( int j = 0; j< partial_sum[i]; j++)
		{
		    perm[pairs[j][0] - min] = pairs[j][1] - min;

		    if ( pairs[j][1] <= bdry_max )
				perm[pairs[j][1] - min] = 0;
		}
    }

if (DEBUG >= INTERMEDIATE)
{
	debug << "split: disc permutations\n";
	for (int i=0; i< number_of_relators; i++)
	{
	    debug << "split: relator " << i+1 << ": ";
		valarray<int>& perm = *dperm[i];

		for (unsigned int j=0; j < perm.size(); j++)
			debug << perm[j] << ' ';
	    debug << endl;
	}
}
    /********************************************************************
    Evaluate the number of components on each disc.  This is simply a matter 
    of counting the terminating points of arcs on the disc, which are
    identified by a unique zero amongst the first 'mbl[i]' elements of the 
    corresponding disc permutation.
    *************************************************************/
    int ndc[number_of_relators]; /* "number of disc components" */
    for ( int i=0; i< number_of_relators ;i++)
		ndc[i] = 0;

    for ( int i =0; i< number_of_relators; i++)
    {
		valarray<int>& perm = *dperm[i];
		for ( int j = 0; j < mbl[i] - minlb[i] + 1 ; j++)
		{
		    if ( perm[j] == 0 )
				ndc[i] += 1;
		}
    }
   
    /***************************************************************
    Declare space for disc arcs.  The maximum number of edges an arc may 
    contain is the number, N, of 2-cells in the disc under consideration. 
    The number of labels on an arc is thus <=N+1 and, by storing the 
    terminating zero, we must record a maximum of N+2 integers for each arc.  These
    integers will be held as the rows of a matrix.
    *******************************************************************/
    matrix<int>* darcs[number_of_relators]; /*"disc arcs"*/
    for (int i=0; i < number_of_relators ; i++)
		darcs[i] = new matrix<int>(ndc[i], length_of_attaching_map[i]);
   
    /********************************************************************
    Assign disc arcs.  The initial and terminal labels of each arc may be
    recovered from the indexes of the first 'mbl[i]' components of the 
    corresponding permutation (this number must therefore be even).  
    The terminal labels are indicated by zero and so the labels on each arc may be obtained by 
    "tracing through" the permutation until a zero is encountered.  
    This zero is stored to indicate the termination of that arc.
    *******************************************************************/
    
    for (int i=0; i< number_of_relators; i++)
    {
		valarray<int>& perm = *dperm[i];
		matrix<int>& arcs = *darcs[i];
		int min = minlb[i];
		int arcn = 0; /* arc number */

		for ( int j  = 0 ; j < mbl[i] - minlb[i] + 1; j++)
		{
		    if ( perm[j] != 0)
		    {
				arcs[arcn][0] = j + min;
				int arc_cpt = 1; /* "arc component"*/
				bool arc_not_complete = true;
				do
				{
				    arcs[arcn][arc_cpt] = perm[arcs[arcn][arc_cpt-1]-min] + min;
				    arc_cpt++;
	
				    if (arcs[arcn][arc_cpt-1] == min)
				    {
						arcs[arcn][arc_cpt-1] = 0;
						arc_not_complete = false;
				    }
	
				} while (arc_not_complete);
				arcn++;
		    }
		}
	}

if (DEBUG >= INTERMEDIATE)
{
	debug << "split: number of disc components\n";
	for (int i=0; i< number_of_relators; i++)
		debug << ndc[i] << ' ';
	debug << endl;

	debug << "split: disc arcs\n";
	for (int i=0; i< number_of_relators; i++)
	{
	    debug << "split: relator " << i+1 << ":" << endl;
	    	debug << *darcs[i];
	    debug << endl;
	}
}

       
    /********************************************************************
    It remains to determine how the disc arc components are identified when 
    the relator discs are attached to the rosette.  To accomplish this we 
    utilize a matrix whose rows are made up of the end labels of each disc 
    component.
    ********************************************************************/
    int tndc = 0; //"total number of disc components" 
    for ( int i =0; i<number_of_relators ; i++)
		tndc += ndc[i];

    matrix<int> ael(tndc,2); /* "arc end labels" */
    int elm = 0; /* "end label marker" */

    for ( int i =0; i< number_of_relators ; i++)
    {
		matrix<int>& arcs = *darcs[i];
		for (int j=0; j< ndc[i]; j++)
		{
		    ael[elm][0] = arcs[j][0];
		    int k = 1;
		    do  
				k++; 
		    while (arcs[j][k] != 0);
			
		    /* We know there are at least two labels on each arc.  When 
		    this loop terminates, k indicates the terminating zero of the
		    current  arc. */

		   ael[elm][1] = arcs[j][k-1];
		   elm++;
		}
    }
   

    /********************************************************************
    We build up components of P as classes of boundary labels, by copying 
    required columns from the label matrices into a row of a component
    matrix.  The  labels in these columns are identified by the attaching
    maps, as discussed above. Given one column, we may find the other 
    columns that contain labels, identified in L, belonging to the same 
    component of P by moving  along an arc in some relator disc from each 
    label in turn, ie examining the corresponding place in 'ael'.  We shall 
    set elements of 'ael' to be negative to prevent us obtaining the same 
    column twice thus enabling a more efficient method of constructing the 
    classes.
   
	 We have 'tndc' as an upper bound for the number of components and 
    the number of boundary labels as an upper bound for the size of a 
    component.  Note the labels on interior 1-cells are irrelevant to this 
    procedure, since no identification of  discs occurs along interior 
    1-cells.
    *********************************************************************/
    /* calculate the number of boundary labels */
    label = 0;
    for (int i=0; i< number_of_relators; i++)
    	label += mbl[i] - minlb[i] + 1;

    matrix<int> cptm(tndc,label+1); /* "component matrix" We need to add 1
				    to this total to ensure we have a zero
				    at the end of each row of cptm, for 
				    later use */
 
    /*Build up component classes.*/
    int cpt = 0;
  	bool all_components_found;
	
    do /* while !all_components_found */
    {
		/* Find a boundary label not already recorded. */
		all_components_found = true;
		int new_label;
		for (int i = 0; i< tndc; i++)
		{
			for (int j = 0; j < 2; j++)
			{
		    	if (ael[i][j] > 0)
		    	{
					all_components_found = false;
					new_label = ael[i][j];
					goto new_label_found;
	    		}
			}
		}
		new_label_found:
		
		if (!all_components_found)
		{
	    	int celt = 0; // component element
	    	Slice_iter<int> pointer = cptm[cpt]; /* pointer will reference a label in the current comonent */
	    	bool component_complete = false;
	    	do /* While !component_complete */
	    	{
				int one_cell;
				int column;
				
				/* Look for this new label in the label matrices. */
				for (int i=0; i< number_of_1_cells_in_rosette; i++)
				{
		    		matrix<int>& matrixptr = *labels_for_1_cell[i];
					
		    		for (unsigned int j=0; j < matrixptr.numrows(); j++)
		    		{
						for (unsigned int k=0; k< matrixptr.numcols(); k++)
						{
				    		if ( matrixptr[j][k] == new_label)
				    		{
								one_cell = i;
								column = k;
								goto found_in_label_matrices;
				    		}
						}
		    		}
				}
				found_in_label_matrices:

				 /************************************************************************				   
				 Set those integers in 'ael' that appear in the column just found to be 
				 negative, so we do not look for this column again.  These integers will 
				 currently be positive, since '*new_label' was positive and each label 
				 occurs in a unique place amongst the label matrices.  The integers may be 
				 found by the symmetry of 'jba'; at the same time we copy this column 
				 into the current component.
				 ************************************************************************/

	    		matrix<int>& matrixptr = *labels_for_1_cell[one_cell];
				for (unsigned int i=0; i< matrixptr.numrows(); i++)
				{
		    		int label = matrixptr[i][column];
					for (int i = 0; i< tndc; i++)
					{
						for (int j = 0; j < 2; j++)
						{
					    	if (ael[i][j] == label)
					    	{
								ael[i][j] *= -1;
								goto found_in_arc_end_labels;
				    		}
						}
					}
					found_in_arc_end_labels:
		    		
					cptm[cpt][celt++] = label;
				}

				/* Now find the next label to look for. */
				bool label_found = false;
				do
				{
		    		if ( *pointer == 0 )
					component_complete = true;
		    		else
		    		{
						/* Look for *pointer in ael.  When we find it,
						   we check to see if we already have the labels
						   at the other end of this arc.  If we have them,
						   increment pointer and look again */
						   
						for (int i=0; i< tndc; i++)
						{
							for (int j=0; j< 2; j++)				
							{
				    			if ( abs(ael[i][j]) == *pointer )
					    		{
									if ( j == 0 && ael[i][1] > 0)
									{
							    		new_label = ael[i][1];
							    		label_found = true;
									}
									else
									if ( j == 1 && ael[i][0] > 0) 
									{
							    		new_label = ael[i][0];
							    		label_found = true;
									}
									break;
					    		}
							}
						}

						if (!label_found)
			    			pointer++;
		    		}
				} while(!component_complete && !label_found);
	    	} while (!component_complete);
	    	cpt++;
		}
    } while (!all_components_found);

if (DEBUG >= INTERMEDIATE)
{
	debug << "split: total number of disc arcs = " << tndc;
	debug << "split: arc end labels:\n";
	debug << ael << endl;
	debug << "split: component matrix:\n";
	debug << cptm << endl;
}
	
	/* we can now do some clearing up of memory*/
	for (int i=0; i< number_of_1_cells; i++)
		delete labels_for_1_cell[i];

    /*********************************************************************
    We may now use the rows of 'cptm' to determine the tuples corresponding to 
    each component of P.  Each label of a disc arc component occurs exactly 
    once in 'cptm', clearly with both end labels of a given arc in the same 
    component.  We therefore determine which component each "initial" label 
    lies in and then use that disc arc to build up the tuple for the   
    corresponding component.  Notice 'cpt' now indicates how many components 
    P has.
    **********************************************************************/

    /* Check whether the components are required by call to function. */
    if ( short_split )
		return (cpt == 1? true: false);
		 
	if (!return_components)
	{
	    cout << "\n";
		for (unsigned int i=0; i<pattern.size(); i++)
			cout << pattern[i] << ' ';
		cout << endl;

    	txt << "\n";
		for (unsigned int i=0; i<pattern.size(); i++)
			txt << pattern[i] << ' ';
		txt << endl;
	}

	matrix<int> ctracks(cpt, pattern.size()); /* "component tracks" */

	int b = cptm.numcols();
	for (int i=0; i< number_of_relators; i++)
	{
		matrix<int>& arcs = *darcs[i];
		matrix<int>& pairs = *pair[i];

		for (int j=0;  j < ndc[i]; j++)
		{
			/* Determine in which component the initial label of 'arcs[j,]' lies. */
			bool label_found = false;
			int mark = arcs[j][0];
			int disc_cpt; // "disc_component"
			int pattern_cpt; // "pattern component"

			for (int k=0; k < cpt; k++)
			{
				for (int l=0; l < b; l++)
		    	{
					label = cptm[k][l];
					if ( label == mark )
					{
						disc_cpt = j;
						pattern_cpt = k;
						label_found = true;
						break;
					}
					else
					if ( label == 0 )
						break;
				} 

				if (label_found)
					break;
			}

			 /**********************************************************************
			 We now traverse the arc 'arcs[disc_cpt,]', at each stage determining with which 
			 variable the corresponding edge is associated. In the loop below, 
			 'arcs[disc_cpt,k]' is the initial label of some edge on disc i and so it appears 
			 exactly once in the first column of 'pairs'.  We now determine in which 
			 row, then, the third element of that row determines the variable 
			 associated with the pair, so we add 1 to the corresponding variable in 
			 'ctracks[in cpt,]'.
			 ***********************************************************************/

			for (unsigned int k = 0; k < arcs.numcols(); k++)
			{
				label = arcs[disc_cpt][k];
				if (label == 0)
					break;

		    	label_found = false;
				for (int l = 0; l< partial_sum[i]; l++)
				{
					if ( pairs[l][0] == label )
					{
						ctracks[pattern_cpt][pairs[l][2]] += 1;
						label_found = true;
						break;
					}
				}
			}
	    }
	}
		
if (DEBUG >= INTERMEDIATE)
{
	debug << "split: component tracks\n";           
	debug << ctracks << endl;
}

	if (return_components)
		*component_ptr = new matrix<int>(ctracks);
	else
	{
	    if ( cpt == 1 )
		{
			cout << "which is a track\n";
			txt << "which is a track\n";
		}
		else
		{
			cout <<	"which has components\n";           
			txt <<	"which has components\n";           
			cout << ctracks << endl;
			txt << ctracks << endl;
		}
	}

    /* This case is 'short_split = false', we just need to return 
       a bool which will be discarded at the call. */
	
	for (int i=0; i< number_of_relators; i++)
	{
		delete pair[i];
		delete dperm[i];
		delete darcs[i];
	}

    return true;
}

