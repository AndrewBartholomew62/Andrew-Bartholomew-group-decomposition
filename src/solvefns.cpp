/* The two functions
	
   bool check_thresholds(int n, vector<int>& u, n_limit& three_m_limit, list<threshold>::iterator *current_threshold);
   int update_n_limit (n_limit& limit, vector<int> u);
   
   defined herein form part of the search algorithm used to determine minimal fundamental solutions to the matching
   equations.  The search algorithm is described in the notes "fundamental-solutions" that should accompany the 
   software distribution, they are also available at www.layer8.co.uk/maths/tracks
*/

using namespace std;

#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <ctype.h>
#include <stdio.h>
#include <iomanip>
#include <ctime>
#include <list>
#include <vector>

/******************* External variables ************************/
extern unsigned long	num_rank_calls;
extern bool		Greenberg_version;
//extern unsigned int DEBUG;
extern ofstream		output;
extern ofstream		debug;

#include <util.h>
#include <matrix.h>
#include <bigint.h>
#include <n_limit.h>
#include <decomp.h>

extern matrix<int>* ext_fund_ptr;

extern bigint w_cuboid_volume_remaining;

//void print(matrix<int> m, ostream& s, int n, string prefix="");
bigint remaining_volume (n_limit& limit, vector<int> u);
bigint non_excluded_region_size (n_limit& limit, int three_m);
bool compatible (ray& P1, ray& P2, K_complex K);

void exit_option ()
{

    char	ch;

	cout << "\nContinue [y/n]";
	cin >> ch;

	if (ch == 'y')
		return;
	else
	{
	    cout << "Compile time parameters need to be changed before ";
	    cout << "this group may be analyzed.";
	    exit (0);
	}
}

/*******************************************************************
'non_zero count' counts the number of non zero components of the 
row in mat
*******************************************************************/
int non_zero_count(matrix<int>& mat, int row)
{
   int count = 0;
   for ( unsigned int i = 0; i < mat.numcols(); i++)
       if ( mat[row][i] != 0)
	   count += 1;
   return (count);
}

int gcdn (const matrix<int>& mat, int row)
{
    int g = abs(mat[row][0]);

    for (unsigned int i = 1; i<mat.numcols(); i++)
       g = gcd(g, abs(mat[row][i]));

    return (g);
}

int gcdn (int* iptr, int n)
{
    int g = *iptr++;

    for (int i = 1; i<n; i++)
       g = gcd(g, abs(*iptr++));

    return (g);
}

int gcdn (vector<int> v, int n)
{
    int g = v[0];

    for (int i = 1; i<n; i++)
       g = gcd(g, abs(v[i]));

    return (g);
}

/***********************************************************************
Old version:

If 'full echelon' is true at the call, 'echelon' reduces 'matrix' to
echelon form and assigns this echelon form to 'matrix'. If 'full echelon'
is false at the call then the function sends a row permutation to
'reference to indep rows' the first 'RANK matrix' of which determine a
maximal set of indep' rows of 'matrix'.  In this case 'matrix' will not
have been assigned to from 'echelon'.  'reference to indep rows' must
therefore be a reference to a [UPB matrix]INT multiple.

New version: 

echelon reduces mat to echelon form
***********************************************************************/
//void echelon (matrix<int>& matrixref, bool full_echelon)
void echelon (matrix<int>& matrixref)
{
	matrix<int> mat(matrixref);
    int m = mat.numrows();
	int n = mat.numcols();

    int perm [m]; //Creates row permutation vector.

    for (int i = 0 ; i< m; i++)
		perm[i] = i;

    int r = 0;
    int c = 0;
    bool complete = false;

    if (non_zero_count(mat,0) > 1 )
    {
		int g = gcdn (mat,0);
		for (int i = 0; i < n ; i++)
		   mat[0][i] /= g;

if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << "echelon: row 0 scaled by " << g << endl;
    }


    do
    {
		int a = r;
		int p = r;
		bool found = false;

		do
		{
	    	do
	    	{
				if (mat[perm[p]][c] != 0)
				{
		    		int t = perm[a];
		    		perm[a] = perm[p];
		    		perm[p] = t;
		    		found = true;
if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << "echelon: found non-zero entry at " << p << "," << c << endl;
				}
				else
		    		p += 1;
	    	} while ( !found && p < m);

	    	if (!found)
	    	{
				p = a;
				c += 1;
	    	}
		} while (!found && c < n );

		if (!found)
	    	complete = true;
		else
		{
	    	if ( non_zero_count(mat,perm[r]) > 1)
	    	{
				int g = gcdn (mat,perm[r]);
				for (int i= c; i< n ; i++)
			    	mat[perm[r]][i] /= g;

if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << "echelon: row " << r << " scaled by " << g << endl;
	    	}
			
	    	for (int i = r+1 ; i < m; i++)
			{
				if (mat[perm[i]][c] != 0)
				{
		    		if(non_zero_count(mat,perm[i]) > 1)
		    		{
						int g = gcdn(mat,perm[i]);
						for (int j = c ;j < n ; j++)
				    		mat[perm[i]][j] /= g;
if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << "echelon: row " << i << " scaled by " << g << endl;
		    		}

		    		a = mat[perm[r]][c];
		    		p = mat[perm[i]][c];
		    		for (int j = c ; j< n; j++)
		    		   mat[perm[i]][j] = a * mat[perm[i]][j] - p * mat[perm[r]][j];
				}
			}
if (decomp_control::DEBUG >= decomp_control::DETAIL)
{
	debug << "echelon: matrix row reduced to" << endl;
	print (mat,debug,3,"echelon:   ");
}
			
		}

		if (c == n-1)
	    	complete = true;
		else
		{
	    	c += 1;
	    	r += 1;
		}
    } while (!complete && r < m);


	for (int i=0 ; i< m ; i++)
	{
	    for (int j = 0; j<n; j++)
			matrixref[i][j] = mat[perm[i]][j];
	}

}

int rank (matrix<int>& matrixref)
{
    matrix<int> mat(matrixref);
    echelon(mat);

if (decomp_control::DEBUG >= decomp_control::DETAIL)
{
	debug << "rank: echelon form:" << endl;
	print (mat, debug, 3, "rank:   ");
}
	
    int m = mat.numrows();
	int n = mat.numcols();
	
    int rank = n;

    bool free_found = false;

    int r = 0;
	int c = 0;

    do
    {
		if(mat[r][c] != 0)
		{
	    	r += 1;
	    	c += 1;
	    	if (c == n) 
				free_found = true;
		}
		else
		{
if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << "rank: free variable identified at " << r << " " << c << endl;
	    	rank -= 1;
	    	c += 1;
	    	if (c == n) 
				free_found = true;

	    	if (!free_found)
	    	{
				bool non_zero = false;
				do
				{
		    		if ( mat[r][c] == 0)
		    		{
if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << "rank: free variable identified at " << r << " " << c << endl;
						rank -= 1;
						c += 1;
						if ( c == n ) 
							free_found = true;
		    		}
		    		else
		    		{
						non_zero = true;
						r += 1;
						c += 1;
						if ( c == n ) 
							free_found = true;
		    		}
				} while (!non_zero && !free_found);
	    	}
		}
    } while (!free_found && r < m );

   	if (!free_found)
   	{
   		for (int i = c ; i < n ; i++ ) 
		{
if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << "rank: free variable identified at " << r << " " << c << endl;
			rank -= 1;
		}
	}
	
	return rank;
}

/*************************************************************************
The matrix basis must be ready to accept a basis for the kernel of matrixref
in its COLUMNS (should the basis be non trivial).  If  the kernel is trivial, 
the function returns the bool value false.  If the kernel is non trivial, 
a basis is stored in the matrix basis and the bool value true is returned.
The basis is calculated in the standard way, isolating "free" variables and 
assigning 1 to each in turn, the others remaining zero, then back substituting.  
The only modifications are that the arithmetic is entirely integral and, should
the ith component of a basis element, e say, be negative after the back
substitution (where the ith variable is free and initially assigned 1 to
calculate e), we multiply e through by -1.  This simplifies matters later.
************************************************************************/
bool kbasis(matrix<int>& matrixref, matrix<int>& basis)
{
    matrix<int> mat(matrixref);
    int m = mat.numrows();
	int n = mat.numcols();
	
    echelon(mat);

    int r = 0;
	int c = 0;

    /************ Isolate free variables. *********************/
    int free[n];
    for (int i=0; i < n ; i++ )
		free[i] = 0;

    int kdim = 0;
    bool free_found = false;
    do
    {
		if( mat[r][c] != 0)
		{
	    	r += 1;
	    	c += 1;
	    	if (c == n) 
				free_found = true;
		}
		else
		{
	    	kdim += 1;
	    	free[c] = 1;
	    	c += 1;
	    	if (c == n) 
				free_found = true;

	    	if (!free_found)
	    	{
				bool non_zero = false;
				do
				{
		    		if ( mat[r][c] == 0)
		    		{
						kdim += 1;
						free[c] = 1;
						c += 1;
						if ( c == n )
							free_found = true;
		    		}
		    		else
		    		{
						non_zero = true;
						r += 1;
						c += 1;
						if ( c == n ) 
							free_found = true;
		    		}
				} while (!non_zero && !free_found);
	    	}
		}
    } while (!free_found && r < m );

    if (!free_found)
	{
		for (int i = c ; i < n ; i++ )
		{
		    kdim += 1;
		    free[i] = 1;
		}
	}

    if ( kdim == 0 )
		return false;
    else // Calculate basis elements.
    {
		/*********************************************************
		First determine which variable will be determined (after
		back substituting) by each of the first n-kdim (ie the
		non zero) rows of 'mat'.
		************************************************************/
		int variable[n-kdim];
		
		int variable_number = 0;
		for ( int i = 0; i < n ; i++)
		{
	    	if (free[i] == 0)
	    	{
				variable[variable_number] = i;
				variable_number += 1;
	    	}
		}
		
		/*******************************************
		Now we know row i of 'mat' (i=1,...,n-kdim)
		controls 'variable[i]'.
		*******************************************/

		int eltnumber = 0;
	    int element[n];
		for (int i = 0; i < n; i++)
		{
	    	if (free[i] == 1)
	    	{
				/************************************************
				If the ith variable is free we shall get a basis
				element for the kernel.
				*************************************************/

				/********** Reset element. ****************/
				for (int j=0 ; j < n ; j++)
		    		element[j] = 0;

				element[i] = 1;

				/********* Back substitute. **********/
				for (int j = n-kdim-1 ; j>=0 ; j-- )
				{
		    		int velt = variable[j];
		    		for (int k = variable[j]+1 ; k < n ; k++)
						element[velt] -= mat[j][k] * element[k];

		    		int melt = mat[j][velt];
		    		for (int k = variable[j]+1 ; k< n; k++)
						element[k] *= melt;

				}

				int g = gcdn(element, n);
				
				if (element[i] < 0)
		    		g *= -1; /*This has the effect of multiplying through by -1.*/

				for (int j = 0; j < n; j ++)
		    		element[j] /= g;

				for (int j = 0; j < n; j++)
		    		basis[j][eltnumber] = element[j];
					
				eltnumber += 1;
	    	}
		}

		return true;
    }
}

/*******************************************************************
bool extreme_ray(row* zrecord, int zcount)
This procedure considers those rows of 'cone matrix 2' indicated by
'zrecord' and returns the bool true if the collection contains a
linearly independent subcollection of kdim-1 elements, false otherwise.
The integer 'zcount' should indicate how many rows 'zrecord' determines,
we assume that this is >=kdim-1.
*******************************************************************/
bool extreme_ray(matrix<int>& cone_matrix_2, int* zrecord, int zcount)
{
    int zind[zcount];
	
    int n = cone_matrix_2.numcols();
    int kdim = cone_matrix_2.numrows();
	

	bool not_extreme = true; 
	bool not_finished = true;


    matrix<int> test_matrix(kdim-1,kdim);
 
    /* isolate the indices of the indicated planes. */
    int zc = 0;

    for (int i=0; i < n; i++)
    {
		if ( zrecord[i] == 1 )
		    zind[zc++] = i;
    }

    if (decomp_control::DEBUG >= decomp_control::DETAIL)
    {
		debug << "\nextreme_ray: zind: ";
		for (int i=0; i<zcount; i++)
			debug << zind[i] << ' ';
		debug << endl;
		debug << "extreme_ray: zcount = " << zcount << endl;
    }

    int block[kdim-1];
	int max[kdim-1]; //[kdim-1]INT block,max;

    /* Initialize */
    for ( int i = 0 ; i< kdim-1 ; i++)
    {
		block[i] = i;
		max[i] = zcount-kdim+1+i;
    }

    if (decomp_control::DEBUG >= decomp_control::DETAIL)
    {
		debug << "extreme_ray: max: " << endl;
		debug << "extreme_ray:   ";
		for (int i=0; i<kdim-1; i++)
			debug << max[i] << ' ';
		debug << endl;
    }

    do
    {
		do
		{
	    	/* Test the planes determined by 'block'. */

if (decomp_control::DEBUG >= decomp_control::DETAIL)
{
	debug << "extreme_ray: test block: " << endl;
	debug << "extreme_ray:   ";
	for (int i=0; i<kdim-1; i++)
		debug << block[i] << ' ';
	debug << endl;
}

	    	cout << "\r" << ++num_rank_calls;

	    	for (int i = 0 ; i < kdim-1 ; i++)
			for (int j = 0 ; j < kdim ; j++)
		    	test_matrix[i][j] = cone_matrix_2[j][zind[block[i]]];

	    	if ( rank(test_matrix) == kdim-1 )
				not_finished = not_extreme = false;

	    	block[kdim-2] += 1;

		} while ( not_finished && block[kdim-2] < zcount );

		if (not_finished)
		{
if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << "extreme_ray: reset test" << endl;
	
	    	int test = kdim-3;
if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << "\ntest = " << test << endl;
	
	    	bool test_not_complete = true;
	    	do
	    	{
if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << "extreme_ray: block[" << test << "] = " << block[test] << endl;
	
				if ( block[test] == max[test] )
				{
		    		if ( --test == -1 )
		    		{
						not_finished = test_not_complete = false;
if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << "extreme_ray: test complete" << endl;
		    		}
if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << "extreme_ray: test decremented to " << test << endl;
				}
				else
				{
		    		block[test] += 1;
					
		    		for ( int i = 1; i < kdim-1-test; i++)
						block[test+i] = block[test]+i;

		    		test_not_complete = false;
				}

	    	} while (test_not_complete);
		}
    } while (not_finished);

    return !not_extreme;
}


/* the greenberg algorithm function returns the number of rows in emat when the algorithm completes */
int greenberg (matrix<int>& cone_matrix_2, list<ray>& raylist, time_t start_time)
{

int n = cone_matrix_2.numcols();
int kdim = cone_matrix_2.numrows();

/* 'nr ineq[i]' will be set to zero iff the inequality determined by
the ith column of 'cone matrix 2' becomes redundant. */
int nr_ineq[n];
for (int i=0; i<n ; i++)
    nr_ineq[i] = 1;

int new_ineq = kdim+1;
int s = kdim;

/* Initialize the ray list to the rows of cone_matrix_2 */
for ( int i=0; i< s; i++)
{
	vector<int> ray(n);
	for (int j=0; j<n; j++)
		ray[j] = cone_matrix_2[i][j];
	raylist.push_back(ray);
}


bool initial_step = true;

/******************** main loop clause ***********************/
while (new_ineq <= n) 
{
	/*Declare space for, and assign, "Ev".  Recall this, in fact, is
	  the transpose of "Ev", as defined by Greenberg and we also carry
	  all the rows atEv (sic, see Greenberg's algorithm) forming Fv as 
	  the first part.  This makes alpha easy to calculate.*/

	matrix<int> ematrix(s,n);


	if (initial_step) 	       /*We simply require 'cone matrix 2'.*/
	{
		ematrix = cone_matrix_2;
	}
	else /* We must read ematrix from raylist */
	{

		list<ray>::iterator rptr = raylist.begin();
		int i = 0;
		do
		{
			for (int j=0; j< n; j++)
				ematrix[i][j] = rptr->r[j];
			i++;
			rptr++;
		} while (rptr != raylist.end());

	}

if (decomp_control::DEBUG >= decomp_control::DETAIL)
{
    debug << "greenberg: ematrix" << endl;
   	print(ematrix, debug, 3, "greenberg:   ");
}

	/*Isolate 'alpha' as the 'new_ineq'th column of ematrix.*/
	/* Used to be REF[]INT alpha=ematrix[,new ineq]; */

	int alpha = new_ineq-1;

if (decomp_control::DEBUG >= decomp_control::DETAIL)
{
    debug <<"greenberg: new_ineq = " << new_ineq << "; alpha =" << endl;
	debug << "greenberg:   ";
    for (int i=0; i<s;i++)
		debug << ematrix[i][alpha] << " ";
	debug << endl;
}

	/*Determine the size of N+ and N-.*/
	int pcount = 0;
	int ncount = 0;
	for (int i =0; i< s; i++)
	{
	    if (ematrix[i][alpha] > 0 )
			pcount += 1;
	    else
		if (ematrix[i][alpha] < 0 )
			ncount += 1;
	}

if (decomp_control::DEBUG >= decomp_control::DETAIL)
    debug << "greenberg: pcount = " << pcount << " ncount = " << ncount << endl;

	int new_rays = 0;
	/*'new rays' will count how many rays we add to "Ev+1".*/

	if ( ncount != 0)
	{
	    /*Determine N+ and N-.*/
	    int pind[pcount];
	    int nind[ncount];
	    pcount = ncount = 0;

	    for (int i =0 ; i< s; i++)
		if ( ematrix[i][alpha] > 0 )
		    pind[pcount++] = i;
		else
		if ( ematrix[i][alpha] < 0 )
		    nind[ncount++] = i;

if (decomp_control::DEBUG >= decomp_control::DETAIL)
{
	debug << "greenberg: pind: ";
	for (int i=0; i< pcount; i++)
		debug << pind[i] << ' ';
	debug << endl;

	debug << "greenberg: nind: ";
	for (int i=0; i< ncount; i++)
		debug << nind[i] << ' ';
	debug << endl;
}

	    /*Update 'nr_ineq'; that is, determine Jv+1.*/
	    for (int i = 0 ; i < new_ineq-1 ; i++)
		{
			if ( nr_ineq[i] == 1 )
			{
			    bool dependent = true;
			    for (int j = 0; j < s; j++)
				{
					if ( ematrix[j][i] == 0 // beta[j]
							    && ematrix[j][alpha] > 0 )
					{
					    dependent = false;
					    break;
					}
				}
				
			    if ( dependent )
					nr_ineq[i] = 0;
			}
		}

	    if (decomp_control::DEBUG >= decomp_control::DETAIL)
	    {
			debug << "greenberg: nr_ineq" << endl;
			
			debug << "greenberg:   ";
			for (int i=0; i< n; i++)
				debug << nr_ineq[i] << ' ';
			debug << endl;
	    }

	    /* remove from raylist the rows of 'ematrix' we should loose,
	      ie those indexed by 'nind'; that is, columns of Ev+1
	      determined by N0 and N+.*/

		list<ray>::iterator rptr = raylist.begin();
	    int row_counter = 0;
	    for (int i=0; i < ncount ; i++)
		{
			/* advance rptr til it points at the nind[i]-th element of the list */
			while (row_counter <  nind[i] )
	    	{
				rptr++;
				row_counter++;

	    	}

			/* now remove that element of the list */
if (decomp_control::DEBUG >= decomp_control::DETAIL)
debug << "greenberg: removing list element " << row_counter << endl;

			rptr = raylist.erase(rptr);
			row_counter++;
		}


		/*We now determine the required e(j,i) for j in N+ and i in N-, 'nr ineq'
		 indicates the non redundant set of inequalities to be considered for this
		 purpose.  We compute any such rays in 'ray buffer'.*/

	    if ( pcount != 0 )
	    {
			int common_zero = 0;
			for (int i = 0; i< ncount; i++)
			{
		    	int ni = nind[i];
if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << "greenberg: starting nind " << ni << endl;

		    	for (int j = 0; j<pcount; j++)
		    	{
					int zrecord[n];
					int pj = pind[j];
if(decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << "greenberg: starting pind " << pj << endl;

					if (!initial_step)
					{
			    		common_zero = 1;

			    		for (int k=0; k < n ; k++)
							zrecord[k] = 0;

			    		zrecord[new_ineq-1] = 1;

			    		for (int k=0; k< new_ineq-1; k++)
						{
							if (nr_ineq[k] == 1 && ematrix[pj][k] == 0 && ematrix[ni][k] == 0 )
							{
					    		common_zero++;
					    		zrecord[k] = 1;
							}
						}
					}

if (decomp_control::DEBUG >= decomp_control::DETAIL)
{
	if (initial_step)
	{
		debug << "greenberg: initial step" << endl;
	}
	else
	{
		debug << "greenberg: checking zrecord ";
		for (int i=0; i<n; i++)
			debug << zrecord[i] << ' ';
		debug << endl;

		debug << "greenberg:   common_zero = " << common_zero;
		if ( common_zero < kdim-1 )
			debug << "greenberg:   not a candidate" << endl;
	}
}

					if ( common_zero >= kdim-1 || initial_step )
					{
			    		if ( initial_step || Greenberg_version || extreme_ray(cone_matrix_2,zrecord,common_zero) )
			    		{
if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << "greenberg:   extreme" << endl;
							new_rays += 1;
							int a1 = ematrix[pj][alpha];
							int a2= -ematrix[ni][alpha];

							vector<int> ray_buffer(n);

							for (int k=0; k < n ; k++)
				    			ray_buffer[k] = a2*ematrix[pj][k]+ a1*ematrix[ni][k];
							int g = gcdn(ray_buffer, n);

							for (int k =0; k <n; k++)
				    			 ray_buffer[k] /= g;

							ray vector(ray_buffer);
							raylist.push_back(vector);

if (decomp_control::DEBUG >= decomp_control::DETAIL)
{
	debug << "greenberg:   new ray ";
	for (int i=0; i< n; i++)
		debug << ray_buffer[i] << ' ';
	debug << endl;	
}

			    		}
			    		else
						{
if (decomp_control::DEBUG >= decomp_control::DETAIL)
debug << "greenberg:   not extreme" << endl;
						}
					}
		    	}
			}
	    }

	    /*Adjust s to denote the number of rays determined by "Ev+1".*/
	    s += (new_rays - ncount);
	    initial_step=false;

	}
	else /* We must adjust nr_ineq to accommodate this dependance */
	{
	    nr_ineq[alpha] = 0;
	}

if (false)
{
	output << "\nProcessed hyperplane " << new_ineq-kdim << " of "
		<< n-kdim << ", "<< ncount << " ray(s) removed " << new_rays
		<< " added." << "  s = " << s;
}
	cout << "\nProcessed hyperplane " << new_ineq-kdim << " of "
		<< n-kdim << ", "<< ncount << " ray(s) removed " << new_rays
		<< " added." << "  s = " << s;

if (decomp_control::DEBUG >= decomp_control::DETAIL)
{
	debug << "greenberg: processed hyperplane " << new_ineq-kdim << " of " << n-kdim << ", "<< ncount 
  		  << " ray(s) removed " << new_rays << " added." << "  s = " << s << endl;
}



	/* Evaluate time taken so far */
	time_t current_time = time(0);
	int time_exp = static_cast<int>(difftime(current_time, start_time));
	int hh = time_exp/3600;
	int mm = (time_exp - hh * 3600)/60;
	int ss = time_exp - hh * 3600 - mm * 60;
if (false)
{
	output << "\nTotal time expired " << hh << " hours " << mm
	      << " minutes " << ss << " seconds\n";
}
	cout << "\nTotal time expired " << hh << " hours " << mm
		<< " minutes " << ss << " seconds\n";

	if ( new_ineq < n )
	    cout << "\nRay test\n" << 0;

	new_ineq +=1;

	/* Reset num_rank_calls */
	num_rank_calls = 0;

} /************ end of main clause ****************/


return s;

}

void usage_help()
{
	cout << "\n\nUsage: solve [-#[1|2|3]bEFgh] [<input_file>] [<outfile>]" << endl;
	cout << "      #[1|2|3]: generate debug output" << endl;
	cout << "      b: only calculate a track basis for the matching system kernel." << endl;
	cout << "      E: do not calculate extreme fundamental solutions." << endl;
	cout << "      F: do not calculate fundamental solutions" << endl;
	cout << "      g: use original Greenberg algorithm" << endl;
	cout << "      h: show help information" << endl;
	cout << "\n\nThe E, F, and g options are only appicable if the b option is not provided." << endl;
	cout << "In this case the programme will, by default, evaluate a set of fundamental" << endl;
	cout << "solutions to the matching system.  The F option then allows the calculation" << endl;
	cout << "to stop after the extreme fundamental solutions have been identified and the" << endl;
	cout << "E option permits the analysis of a groups cell-complex via the debug capabilitiy" << endl;
	cout << "without any calculation of extreme fundamental solutions." << endl;
	exit(0);
}

void initial_help()
{
	cout << "\n\nUsage: solve [-#[1|2|3]bEFGgh] [<input_file>] [<outfile>]\n\n";
	cout << "This programme may be run by typing 'solve' at the command line or by\n"
	     << "specifying a parameter file <infile> by typing 'solve <infile>'.\n"
	     << "A parameter file must contain the following information, on separate lines\n"
		 << "as shown.\n";
	cout << "\n<jobname>";
	cout << "\n<number of generators>";
	cout << "\n<number of relators>";
	cout << "\n<generator 1> <generator 2> ...";
	cout << "\n<relator 1> <relator 2> ...\n\n";

	cout << "The relators must not contain any spaces but may be separated by an arbitrary\n";
	cout << "number of spaces, or a newline.\n";
	cout << "\nRepresent inverses of generators by preceeding the character with '-' eg -a,\n";
	cout <<	 "and powers by repetition\n";

	cout << "\nWe also require that:\n";

	cout << "1.  relators indicating a generator of finite order should be ";
	cout << "entered in the\n    positive, eg 'aa', not '-a-a';\n";
	cout << "2.  the relators have length > 3 or take one of the forms 'aa'";
	cout << " or 'aaa';\n";
	cout << "3.  the total, absolute, exponent of each generator is at least";
	cout << " two;\n";

	cout << "\nRequirements 1 to 3 are ESSENTIAL for the running of the program.\n\n",

	cout << "The b option is used to calculate a track basis for the kernel of the \n"
		 << "matching system.\n\n";

	cout << "The #, E and F options are intended for debugging, which writes information \n"
		 << "in a file called solve.dbg.  These options can usually be ignored." << endl;

	cout << "\nThe E, F, and g options are only appicable if the b option is not provided." << endl;
	cout << "In this case the programme will, by default, evaluate a set of fundamental" << endl;
	cout << "solutions to the matching system.  The F option then allows the calculation" << endl;
	cout << "to stop after the extreme fundamental solutions have been identified and the" << endl;
	cout << "E option permits the analysis of a groups cell-complex via the debug capabilitiy" << endl;
	cout << "without any calculation of extreme fundamental solutions." << endl;
		 
	cout << "The '-g' switch causes the original Greenberg algorithm to be used, this \n"
		 << "creates LOTS more tracks but could be quicker in some instances - no promises\n"
		 << "though.  The -h switch displays basic usage information\n" << endl;
}


/* satisfies_matching_system requires that the matching_indices matrix has been
   set up so that it is an rx4 matrix whose first two columns contian the two positive 
   idices from a row of the matching system and the second two columns contian the two 
   negative indices.  Since the matching system sometimes contains
   a row with only one positive and one negative entry (when a relator starts aa... for example)
   we need a method of detecting this and so store indices from 1, allowing 0 to indicate an 
   unused entry.
*/
bool satisfies_matching_equations(matrix<int>& matching_indices,vector<int> u)
{
	for (unsigned int i=0; i< matching_indices.numrows(); i++)
	{
		int sum = u[matching_indices[i][0]-1];
		if (matching_indices[i][1] != 0)
			sum += u[matching_indices[i][1]-1];
		sum -= u[matching_indices[i][2]-1];
		if (matching_indices[i][3] != 0)
			sum -= u[matching_indices[i][3]-1];
					
		if (sum != 0)
			return false;
	}
	return true;
}


/* in the function check_thresholds, the integer n indicates the variable n not an offset into u
   the function returns true if the enumaration is complete, and false otherwise.
*/
bool check_thresholds(int n, vector<int>& u, n_limit& three_m_limit, list<threshold>::iterator *current_threshold)
{	    
	int three_m = u.size();
	string prefix = "check_thresholds: ";
	for (int i=0; i< n-1; i++)
		prefix += "    ";

 	if (n == 8*three_m/10)
	{
		bigint volume_remaining = remaining_volume (three_m_limit, u);
		cout << "\r"//volume remaining = " << volume_remaining << endl;
		     << bigint(100) - volume_remaining*bigint(100)/w_cuboid_volume_remaining << "%" << flush;
	}

	/* check whether there is another threshold for this variable in the current n_limit 
	   the current threshold is indicated by the appropriate current_threshold element
	*/
	list<threshold>::iterator& tptr = current_threshold[n-1];

if (decomp_control::DEBUG >= decomp_control::BASIC)
{
	if (decomp_control::DEBUG >= decomp_control::EXHAUSTIVE && n==1)
	{
		debug << prefix << "current three_m_limit:" << endl;
		three_m_limit.print(debug,0,"check_thresholds: ");
		debug << prefix << "current_thresholds: ";
		for (int i=0; i < three_m; i++)
			debug << *current_threshold[i] << ' ';
		debug << endl;
	}
	
	debug << prefix << "n = " << n << endl;
	debug << prefix << "current threshold pair for variable " << n << ": " << *tptr << endl;
	debug << prefix << "current u = ";
	for (int i=0; i< three_m; i++)
		debug << u[i] << ' ';
	debug << endl;
}
	
	if (n == three_m)
	{
		/* we check if there is another threshold in the 3m_limit */
		if ( ++tptr == three_m_limit.thresholds.end())
		{

if (decomp_control::DEBUG >= decomp_control::BASIC)
	debug << prefix << "no further thresholds for x_3m, enumeration complete" << endl;
	
			/* we've reached the end of the enumeration */
			return true;
		}
		else
		{

if (decomp_control::DEBUG >= decomp_control::BASIC)
	debug << prefix << "advanced current threshold pair for 3m-variable " << n+1 << "to " << *tptr << endl;

			/* re-set the current_threshold pointers to reflect the new (3m-1)-limit indicated by tptr 
			   we start at the beginning of each n-limit for n=3m-1,...,1
			*/
			for (int i=three_m-2; i >= 0; i--) // note n = 3m > 3
				current_threshold[i] = current_threshold[i+1]->second->thresholds.begin();

			/* set the variables u_1,...,u_{3m-1} to zero */
			for (int i=0; i< three_m-1; i++)
				u[i] = 0;			

if (decomp_control::DEBUG >= decomp_control::BASIC)
{	debug << prefix << "current_thresholds adjusted to: " << endl;
	for (int i=0; i < three_m; i++)
		debug << *current_threshold[i] << ' ';
	debug << endl;
}
if (decomp_control::DEBUG >= decomp_control::BASIC)
{
	debug << prefix << "updated u = ";
	for (int i=0; i< three_m; i++)
		debug << u[i] << ' ';
	debug << endl;
}
		}
	}
	else
	{
		/* we have to determine the contextural n-limit for tptr, which is the n-limit indicated by 
		   the current_threshold element for the next variable.
		*/
		list<threshold>::iterator& context_ptr = current_threshold[n];

if (decomp_control::DEBUG >= decomp_control::BASIC)
	debug << prefix << "contextural " << n << "-limit for variable " << n << " is indicated by threshold " << *context_ptr << endl;
	
		list<threshold>& context = context_ptr->second->thresholds;
		
		if ( ++tptr == context.end())
		{

if (decomp_control::DEBUG >= decomp_control::BASIC)
	debug << prefix << "no further thresholds for variable " << n << ", checking variable " << n+1 << ", currently = " << u[n] << endl;

			/* we've come to the end of the thesholds for the current variable (x_{n+1} given that 
			   here n counts variables from zero).  Advance x_{n+1} and see whether the current
			   threshold for x_{n+1} has been breached.
			*/
			if (++u[n] >= current_threshold[n]->first)
			{			  
				/* recurse and look to see if there's another threshold for the next variable 
				   (x_{n+1}) by calling check_thresholds with n+1
				*/

if (decomp_control::DEBUG >= decomp_control::BASIC)
{
	debug << prefix << "incrementing variable " << n+1 << " breaches current threshold " << *current_threshold[n] << endl;
	debug << prefix << "recursing with n = " << n+1 << endl;
}

				return check_thresholds(n+1, u,three_m_limit,current_threshold);
			}
			else
			{
				/* re-set the current_threshold pointers to reflect the current n-limit 
				   indicated by current_threshold[n], since the same n-limit is applied
				   for each value of u[n] (variable n+1) until the threshold indicated
				   by current_threshold[n] is reached. Thus we have to start at the 
				   beginning of each k-limit for k=n,...,1
				*/
				for (int i=n-1; i >= 0; i--) 
					current_threshold[i] = current_threshold[i+1]->second->thresholds.begin();
				
				/* set the variables u_1,...,u_n to zero */
				for (int i=0; i<n; i++)
					u[i]=0;
if (decomp_control::DEBUG >= decomp_control::BASIC)
{
	debug << prefix << "incrementing variable " << n+1 << " stays within current threshold " << *current_threshold[n] << endl;
	debug << prefix << "reset u to ";
	for (int i=0; i< three_m; i++)
		debug << u[i] << ' ';
	debug << endl;
}
			}
		}
		else
		{

if (decomp_control::DEBUG >= decomp_control::BASIC)
	debug << prefix << "advanced current threshold pair for variable " << n << " to " << *tptr << endl;

			/* re-set the current_threshold pointers to reflect the new (n-1)-limit indicated by tptr 
			   we start at the beginning of each k-limit for k=n-1,...,1
			*/
			for (int i=n-2; i >= 0; i--) // note n = 3m > 3
				current_threshold[i] = current_threshold[i+1]->second->thresholds.begin();
				
			/* set the variables u_1,...,u_{n-1} to zero */
			for (int i=0; i< n-1; i++)
				u[i] = 0;						

if (decomp_control::DEBUG >= decomp_control::BASIC)
{
	debug << prefix << "current_thresholds adjusted to: ";
	for (int i=0; i < three_m; i++)
		debug << *current_threshold[i] << ' ';
	debug << endl;
}

if (decomp_control::DEBUG >= decomp_control::BASIC)
{
	debug << prefix << "updated u = ";
	for (int i=0; i< three_m; i++)
		debug << u[i] << ' ';
	debug << endl;
}
		}
	}
	
	return false;
}

void update_n_limit (n_limit& limit, vector<int> u)
{
	int n = limit.n;
	int three_m = u.size();
	ostringstream oss;
	oss << "update_n_limit(" << n << "): ";
	string prefix = oss.str();
	for (int i=n; i< three_m; i++)
		prefix += "  ";

if (decomp_control::DEBUG >= decomp_control::DETAIL)
{
	debug << prefix << "initial limit thresholds: " << limit.thresholds << endl;
//	debug << prefix << "n-limit:" << endl;
//	limit.print(debug,three_m-n,"update_n_limit#: ");
	
	debug << prefix << "updating with u = ";
	for (int i=0; i< three_m; i++)
		debug << u[i] << ' ';
	debug << endl;
}

if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << prefix << "limit thresholds: " << limit.thresholds << endl;


	/* find the first threshold for the nth variable that is >= u[n-1] (n counts variables from 1) */
	list<threshold>::iterator tptr = limit.thresholds.begin();
	while (tptr != limit.thresholds.end())
	{
		if (tptr->first >= u[n-1]) // n counts variables from 1
			break;
		
		tptr++;
	}
	

if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << prefix << "first threshold >= u[" << n-1 << "] (x_" << n << ") = " << u[n-1] << " is " << *tptr << endl;
	
		
	/* If (u_0,...,u_{n-2}) = (0,...0) we can remove all the thresholds after u_{n-1}, since there
	   are no points to check in the (n-1)-limit, which would be used to update any other thresholds
	   in the contextural list, since L_k > (0,...,0) for any limit L_k (see algorithm notes).
	   
	   If we do need to update the remainder of the contextural list of thresholds we remove duplicates
	   from the list, so establish a compare iterator to point at threshold we may end up removing. 	   
	   	   
	*/
	list<threshold>::iterator cmp_ptr;
	
	bool zero_n_minus_one_limit = true;
	for (int i=0; i< n-1; i++)
	{
		if (u[i] !=0)
		{
			zero_n_minus_one_limit = false;
			break;
		}
	}

	if (u[n-1] < tptr->first)
	{
		/* if u_n > 0 (i.e u[n-1]) create a new threshold (u_n, L_{n-1}^i) 
		   in front of tptr (see algorithm notes) 
		*/	
		if (u[n-1] > 0)	
		{	
			threshold new_threshold;
			new_threshold.first = u[n-1];	
			new_threshold.second = new n_limit(*tptr->second);	
	
if (decomp_control::DEBUG >= decomp_control::DETAIL)
{
	debug << prefix << "new_threshold.threshold = " << new_threshold.first << endl;
	debug << prefix << "new_threshold.n_limit at " << new_threshold.second << ":" << endl;
//	new_threshold.second->print(debug,0,prefix);
	debug << prefix << "new_threshold.nlimit thresholds: " << new_threshold.second->thresholds << endl;
}

			limit.thresholds.insert(tptr,1,new_threshold); //insert 1 copy in front of tptr

if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << prefix << "revised threshold list: " << limit.thresholds << endl;

		}				
		else
		{
if (decomp_control::DEBUG >= decomp_control::DETAIL)
		debug << prefix << "no new_threshold due to zero value for u[" << n-1 << "]" << endl;
		}
	
		if (zero_n_minus_one_limit)
		{
if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << prefix << "zero (n-1)-limit, removing thresholds from " << *tptr << " onwards"<< endl;
	
			limit.thresholds.erase(tptr,limit.thresholds.end());
			tptr = limit.thresholds.end();

if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << prefix << "threshold list reduced to: " << limit.thresholds << endl;
	
		}
		else
		{
			/* replace the (n-1)-limit indicated by tptr->second with the (n-1)-limit updated with u. */
			if (n == 2)
			{		
				list<threshold>& one_limit_thresholds = tptr->second->thresholds;
				/* one_limit_thresholds contains just a single threshold, since it is
		           a 1-limit, it will be of the form (t_1,0), so we only need update
			       the first component with u[0]
				*/
				one_limit_thresholds.begin()->first = u[0];
		
if (decomp_control::DEBUG >= decomp_control::DETAIL)
{
//	debug << prefix << "updated 1-limit at " << tptr->first << " to become:" << endl;
//	tptr->second->print(debug,0,prefix);
	debug << prefix << "updated 1-limit at " << *tptr << " to become:" << tptr->second->thresholds << endl;
}

			}
			else
			{
				/* update the (n-1)-limit indicated by tptr->second and record the maximum
			       threshold to which it was updated as the copare_threshold.
				*/

if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << prefix << "recursing to update " << n-1 << "-limit at " << tptr->second << endl;

				update_n_limit(*tptr->second,u);

if (decomp_control::DEBUG >= decomp_control::DETAIL)
{
//	debug << prefix << "updated " << n-1 << "-limit at " << tptr->first << " to become:" << endl;
//	tptr->second->print(debug,0,prefix);
	debug << prefix << "updated " << n-1 << "-limit at " << *tptr << " to become:" << tptr->second->thresholds << endl;
}
			}

			/* tptr is the first threshold against which we compare any subsequent thresholds for duplicates */
			cmp_ptr = tptr;

if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << prefix << " compare_threshold_ptr indicates " << *cmp_ptr << endl;	
		}
	}
	else
	{
		/* u[n-1] = tptr->first.  This case will only occur during the initialization of the 3m-limit 
		   since otherwise (during the enumeration) we will be updating the current thresholds when a 
		   solution is found.  This means that the current threshold for x_n will always be greater 
		   than u[n-1].  For comments on this clause's code, see above,
		   since it is the same after tptr has been advanced.
		*/
if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << prefix << "initialization case for u[" << n-1 << "] = " << tptr->first << endl;

		tptr++;
		
		if (tptr != limit.thresholds.end())
		{
			if (zero_n_minus_one_limit)
			{
if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << prefix << "zero (n-1)-limit, removing thresholds from " << *tptr << " onwards"<< endl;
	
				limit.thresholds.erase(tptr,limit.thresholds.end());
				tptr = limit.thresholds.end();

if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << prefix << "threshold list reduced to: " << limit.thresholds << endl;

			}
			else
			{
				if (n == 2)
				{		
					list<threshold>& one_limit_thresholds = tptr->second->thresholds;
					one_limit_thresholds.begin()->first = u[0];
		
if (decomp_control::DEBUG >= decomp_control::DETAIL)
{
//	debug << prefix << "updated 1-limit at " << tptr->first << " to become:" << endl;
//	tptr->second->print(debug,0,prefix);
	debug << prefix << "updated 1-limit at " << *tptr << " to become:" << tptr->second->thresholds << endl;
}

				}
				else
				{

if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << prefix << "recursing to update " << n-1 << "-limit at " << tptr->second << endl;

					update_n_limit(*tptr->second,u);
					cmp_ptr = tptr;

if (decomp_control::DEBUG >= decomp_control::DETAIL)
{
//	debug << prefix << "updated " << n-1 << "-limit at " << tptr->first << " to become:" << endl;
//	tptr->second->print(debug,0,prefix);
	debug << prefix << "updated " << n-1 << "-limit at " << *tptr << " to become:" << tptr->second->thresholds << endl;
}
				}
			
				cmp_ptr = tptr;

if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << prefix << " compare_threshold_ptr indicates " << *cmp_ptr << endl;	
			}
		}
		else
		{
if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << prefix << "no additional thresholds to consider" << endl;
		}
	}
	
	/* now advance along the contextural list after tptr updating the (n-1)-limits with u and removing
	   any duplicate thresholds from the list.  This is step 5 of the process in the algorithm notes.
	*/
	if (tptr != limit.thresholds.end())
	{
		tptr++;

		if (tptr == limit.thresholds.end())
		{
if (decomp_control::DEBUG >= decomp_control::DETAIL)		
	debug << prefix << "no further thresholds to update" << endl;
		}
		
		while (tptr != limit.thresholds.end())
		{		
			if (n==2)
			{
				/* check whether the 1-limit indicated by tptr->second is bigger than u[0], 
				   if it is we update it and can remove the preceeding threshold from the list,
				   as indicated by cmp_ptr
				*/
				list<threshold> one_limit_thresholds = tptr->second->thresholds;
				if (one_limit_thresholds.begin()->first > u[0])
				{
if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << prefix << "1-limit at threshold " << tptr->first << " greater than " << u[0] << ", removing previous threshold ";
					one_limit_thresholds.begin()->first = u[0];

if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << *cmp_ptr << " to give context" << endl;
	
	
					tptr = limit.thresholds.erase(cmp_ptr); // erase previous threshold
					cmp_ptr = tptr;
					
if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << prefix << limit.thresholds << endl;

				}
				else
				{
if (decomp_control::DEBUG >= decomp_control::DETAIL)
{
	debug << prefix << "1-limit at threshold " << *tptr << " is " << one_limit_thresholds.begin()->first 
	      << ", which is not greater than " << u[0] << endl;
}	      
					break;
				}				
			}
			else
			{
				/* check whether the (n-1)-limit indicated by tptr is greater than (u_0,...u_{n-1}),  
				   if it is we update it with u and check whether can remove the preceeding threshold
	    		   from the list.
			    */		
			
				if (*tptr->second > u)
				{

if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << prefix << n-1 << "-limit at threshold " << *tptr << " greater than (u[0],...,u[n-1])" << endl;

					update_n_limit(*tptr->second, u);

					if (*tptr->second == *cmp_ptr->second)
					{
						
if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << prefix << "  removing previous threshold ";

if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << *cmp_ptr << " to give context ";

						tptr = limit.thresholds.erase(cmp_ptr); // erase previous threshold
						cmp_ptr = tptr;

if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << limit.thresholds << endl;
					
					}
					else
					{
						break;
					}					
				}
				else
				{
					break;
				}		
			}		
		
			tptr++;
		}		
	}

if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << prefix << "final limit thresholds: " << limit.thresholds << endl;
}

bool reset_current_thresholds(n_limit& limit, n_limit& three_m_limit, vector<int>& u, list<threshold>::iterator *current_threshold)
{
	int n = limit.n;
	int three_m = u.size();
	string prefix = "reset_current_thresholds: ";
	for (int i=n; i< three_m; i++)
		prefix += "    ";

	bool cuboid_enumeration_complete = false;
	
if (decomp_control::DEBUG >= decomp_control::DETAIL)
{
	debug << prefix << "resetting current thresholds for u = ";
	for (int i=0; i< n; i++)
		debug << u[i] << ' ';
	debug << endl;
	debug << prefix << "limit thresholds: " << limit.thresholds << endl;
}

	/* find the first threshold for the nth variable that is >= u[n-1] (n counts variables from 1) */
	list<threshold>::iterator tptr = limit.thresholds.begin();
	while (tptr != limit.thresholds.end())
	{
		if (tptr->first > u[n-1]) // n counts variables from 1
			break;
		
		tptr++;
	}	

	if (tptr != limit.thresholds.end())
	{
if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << prefix << "first threshold > u[" << n-1 << "] (x_" << n << ") = " << u[n-1] << " is " << *tptr << endl;
	
		/* set current threshold for x_n to be tptr */
		current_threshold[n-1] = tptr;

if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << prefix << "reset current_threshold[" << n-1 << "] to" << *tptr << endl;
	
		/* tptr->second is the contextural limit for x_{n-1} so recurse to set the current threshold for x_{n-1} */
		cuboid_enumeration_complete = reset_current_thresholds(*tptr->second,three_m_limit,u,current_threshold);
	}
	else if (n<three_m)
	{
if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << prefix << "no " << n << " threshold greater than u[" << n-1 << "] call check_thresholds for variable x_" << n << endl;
		
		/* set current_threshold[n-1] to the last threshold in the list for 
		   variable n, ready for the call to check_thresholds
		*/
		tptr--;
		current_threshold[n-1] = tptr;
		cuboid_enumeration_complete = check_thresholds(n, u, three_m_limit, current_threshold);
	}
	else
	{
		cuboid_enumeration_complete = true;
	}

if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << prefix << "returning cuboid_enumeration_complete = " << (cuboid_enumeration_complete? "true": "false") << endl;
	
	return cuboid_enumeration_complete;
}

/* excluded_region_size calculates the non-excluded volume of the limit */
bigint non_excluded_region_size (n_limit& limit, int three_m)
{
	int n = limit.n;
	string prefix = "non_excluded_region_size: ";
	for (int i=n; i< three_m; i++)
		prefix += "    ";

if (decomp_control::DEBUG >= decomp_control::DETAIL)
{
//	debug << prefix << "n-limit:" << endl;
	debug << prefix << "limit thresholds: " << limit.thresholds << endl;
//	limit.print(debug,three_m-n,"non_excluded_region_size: ");
}
   
	list<threshold>::iterator tptr = limit.thresholds.begin();

	bigint lower_bound;
	bigint upper_bound;
	bigint volume;
	
	while (tptr != limit.thresholds.end())
	{
		bigint upper_bound(tptr->first);
		bigint limiting_value;

if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << prefix << "tptr = " << *tptr << endl;

		if (n==2)
		{
			/* tptr->second is a 1-limit pointer to a list of one threshold of the form (t,0) */
			limiting_value = bigint(tptr->second->thresholds.begin()->first);

if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << prefix << "n = 2, limiting_value = " << limiting_value << endl;
		}
		else
		{

if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << prefix << "n = " << n << "> 2, recursing to calculate limiting_value" << endl;

			limiting_value = non_excluded_region_size(*tptr->second,three_m);

if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << prefix << "recursion returned limiting_value = " << limiting_value << endl;

		}
		
if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << prefix << "upper_bound = " << upper_bound << " lower_bound = " << lower_bound << endl;
	
		volume += (upper_bound-lower_bound)*limiting_value;

if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << prefix << "added " << (upper_bound-lower_bound)*limiting_value << " to volume" << endl;

		lower_bound = upper_bound;
		tptr++;
	}

if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << prefix << "returning volume = " << volume << endl;

	return volume;
}

/* remaining_volume calculates the volume of the non-excluded region of limit from u 
   during the course of the enumeration.
*/
bigint remaining_volume (n_limit& limit, vector<int> u)
{
	int n = limit.n;
	int three_m = u.size();
	string prefix = "remaining_volume: ";
	for (int i=n; i< three_m; i++)
		prefix += "    ";

if (decomp_control::DEBUG >= decomp_control::DETAIL)
{
	debug << prefix << "u = ";
	for (int i=0; i< three_m; i++)
		debug << u[i] << ' ';
	debug << endl;
	debug << prefix << "n-limit:" << endl;
//	debug << prefix << "limit thresholds: " << limit.thresholds << endl;
	limit.print(debug,three_m-n,"remaining_volume: ");
}

	/* find first threshold >= u[n-1] */
	list<threshold>::iterator tptr = limit.thresholds.begin();
	while (tptr != limit.thresholds.end())
	{
		if (tptr->first >= u[n-1]) // n counts variables from 1
			break;
		
		tptr++;
	}

if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << prefix << "first threshold >= u[" << n-1 << "] (x_" << n << ") = " << u[n-1] << " is " << *tptr << endl;
	
	bigint volume=0;
	
	/* first calculate the remainder of the current (n-1)-limit */
if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << prefix << "calculating volume of the remainder of the current " << n-1 << "-limit at " << tptr->second << endl;

	if (tptr->first > u[n-1])
	{
		if (n==2) 
			volume =  bigint(tptr->second->thresholds.begin()->first- u[0]);
		else
			volume = remaining_volume (*tptr->second,u);
	}
	else
	{
		tptr++;

		if (tptr != limit.thresholds.end())
		{
			if (n==2) 
				volume =  bigint(tptr->second->thresholds.begin()->first- u[0]);
			else
				volume =  remaining_volume (*tptr->second,u);
		}
		else
			volume = 0; //u[n-1] is outside the 2-limit in this case

		tptr--; // have to step back for the lower bound of the first slice below.
	}
	
if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << prefix << "volume of the remainder of the current " << n-1 <<"-limit = " << volume << endl;
    
    
	bigint upper_bound;
	bigint lower_bound;
	bigint limiting_value;
	
	/* The location of the first "slice" of volume depends on whether tptr->first > u[n-1] or not 
	   Note also that we have already calculated the remaining volume of the curent (n-1)-limit
	   (corresponding to u[n-1]) so in the first slice the lower bound is incremented by one
	*/
	if (tptr->first > u[n-1])
	{
		upper_bound=tptr->first;
		lower_bound=u[n-1]+1;
	}
	else
	{
		lower_bound = tptr->first + 1;
		tptr++;
		if (tptr != limit.thresholds.end())
			upper_bound = tptr->first;
	}

	if (tptr != limit.thresholds.end())
	{
		if (n==2)
		{
			/* tptr->second is a 1-limit pointer to a list of one threshold of the form (t,0) */
			limiting_value = bigint(tptr->second->thresholds.begin()->first);

if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << prefix << "n = 2, limiting_value = " << limiting_value << endl;
		}
		else
		{
			limiting_value = non_excluded_region_size(*tptr->second,three_m);

if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << prefix << "n > 2, limiting_value = " << limiting_value << endl;
		}

if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << prefix << "upper_bound = " << upper_bound << " lower_bound = " << lower_bound << endl;

		volume += (upper_bound-lower_bound)*limiting_value;

if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << prefix << "added " << (upper_bound-lower_bound)*limiting_value << " to volume" << endl;
	
		lower_bound = upper_bound;
		tptr++;
	}
	
	while (tptr != limit.thresholds.end())
	{
		upper_bound = tptr->first;

if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << prefix << "tptr = " << *tptr << endl;

		if (n==2)
		{
			/* tptr->second is a 1-limit pointer to a list of one threshold of the form (t,0) */
			limiting_value = bigint(tptr->second->thresholds.begin()->first);

if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << prefix << "n = 2, limiting_value = " << limiting_value << endl;
		}
		else
		{
			limiting_value = non_excluded_region_size(*tptr->second,three_m);

if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << prefix << "n > 2, limiting_value = " << limiting_value << endl;
		}

if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << prefix << "upper_bound = " << upper_bound << " lower_bound = " << lower_bound << endl;

		volume += (upper_bound-lower_bound)*limiting_value;

if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << prefix << "added " << (upper_bound-lower_bound)*limiting_value << " to volume" << endl;
	
		lower_bound = upper_bound;
		tptr++;
	}

if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << prefix << "returning volume = " << volume << endl;
	
	return volume;
}

/* next_combination sets comb to the next combination of k integers in the range
   0,...,n-1 where k is the size comb.  The next combination is determined by the 
   natural enumeration of the binnomial coefficient (n,k) number of combinations
   of k objects from a collection of n.  This enumeration starts with the integers
   0,...,k-1 and ends with n-k,...,n-1.  At each stage of the enumeration the upper
   bound for the ith element of the combination (i=0,...k-1) is n-k+i. We check from 
   the (k-1)-th place to see if the upper bound has been reached.  If the upper bound
   has not been reached for the ith place, we increment the ith place by one and set
   all places > i one greater than the place to their left.
      
   The function returns true if another combination exists, in which case comb is 
   set to the next combination and false otherwise.  If no next combination exists
   comb is unaltered by the function.
*/
bool next_combination(vector<int>& comb, int n)
{
	int k = comb.size();
	bool found = false;
	
	for (int i=k-1; i >= 0; i--)
	{
		if (comb[i] < n-k+i)
		{
			comb[i]++;
			for (int j=i+1;j<k;j++)
				comb[j] = comb[j-1]+1;
				
			found=true;
			break;
		}
	}
	
	return found;
}


/* reduce_to_minimal_set reduces the given list of solutions to a minimal set, 
   it removes from the list any u for which there is a w in the set for which 
   w_i <= u_i for all i.
*/
void reduce_to_minimal_set (list<vector<int> >& solutions)
{
    list<vector<int> >::iterator sptr = solutions.begin();
    int n = sptr->size();
    
    while (sptr != solutions.end())
    {	
        	
if (decomp_control::DEBUG >= decomp_control::INTERMEDIATE)
{
	debug << "reduce_to_minimal_set:   sptr ";
    for (int i=0; i< n; i++)
		debug << setw(3) << (*sptr)[i];
    debug << endl;
}
        
		/* we assume sptr is minimal and that we can't remove it.  If it is not
           minimal there is a w with w < sptr, that is w_i < sptr_i for each i.
        */
        bool minimal = true;
        list<vector<int> >::iterator chk_ptr = solutions.begin();
        while (chk_ptr != solutions.end())
        {
			
if (decomp_control::DEBUG >= decomp_control::INTERMEDIATE)
{
	debug << "reduce_to_minimal_set:        ";
	for (int i=0; i< n; i++)
		debug << setw(3) << (*chk_ptr)[i];
}
			if (chk_ptr != sptr)
			{
				int i;
				/* assume we've found a w with w < sptr */
				bool w_found = true;
				for (i=0; i< n; i++)
				{
					if ((*chk_ptr)[i] > (*sptr)[i])
					{
//if (decomp_control::DEBUG >= decomp_control::INTERMEDIATE)
//	debug << " (*chk_ptr)[" << i << "] = " << (*chk_ptr)[i] << " (*sptr)[" << i << "] = " << (*sptr)[i];
						w_found = false;				
						break;
					}
				}			
        		
				if (w_found)
				{
if (decomp_control::DEBUG >= decomp_control::INTERMEDIATE)
  	debug << " w found" << endl;
					minimal = false;
        			break;
        		}
        		else
        		{
if (decomp_control::DEBUG >= decomp_control::INTERMEDIATE)
   	debug << " element " << i << " > sptr" << endl;
        		}
        	}
        	else
        	{
if (decomp_control::DEBUG >= decomp_control::INTERMEDIATE)
  	debug << " is sptr" << endl;
        	}
        	chk_ptr++;
        }
        	
        if (minimal)
        {
if (decomp_control::DEBUG >= decomp_control::INTERMEDIATE)
  	debug << "reduce_to_minimal_set:   minimal, advancing" << endl;
        		sptr++;
        }
        else
        {
if (decomp_control::DEBUG >= decomp_control::INTERMEDIATE)
  	debug << "reduce_to_minimal_set:   not minimal, erasing" << endl;
        		sptr = solutions.erase(sptr);
       	}
    }
}

bool compatible(Slice_iter<int> s,vector<int>& t, K_complex& K)
{
	ray r1(t.size());
	
	for (unsigned int i=0; i< t.size(); i++)
		r1.r[i] = s[i];

	ray r2(t);

if (decomp_control::DEBUG >= decomp_control::DETAIL)
{
	debug << "compatible(Slice_iter<int> s,vector<int>& t, K_complex& K):" << endl;
	debug << "compatible:   r1 = ";
	for (unsigned int i=0; i< t.size(); i++)
		debug << r1.r[i] << ' ';
	debug << endl;
	debug << "compatible:   r2 = ";
	for (unsigned int i=0; i< t.size(); i++)
		debug << r2.r[i] << ' ';
	debug << endl;
}

	return (compatible(r1,r2,K));
	
}
