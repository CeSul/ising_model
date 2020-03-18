#include"hdf5.h"
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<omp.h>


void initialState(int N, int ** lattice_ptr);

void mcmove(int N, int ** lattice, double temp, double external_field);

void printLattice(int N, int ** lattice_ptr, double * temp_ptr,int time_step, hid_t * file_id_ptr);

int rand_range(int min_n, int max_n);

int * site_neighbors (int row, int column, int lattice_size, int** lattice_ptr);
