/*----------------------------------------------------
Functions used in the main code
/----------------------------------------------------*/

#include"ising_functions.h"

int rand_range(int min_n, int max_n){
    // https://stackoverflow.com/a/2999087
    // How to generate random number within range
    return rand() % (max_n - min_n +1) + min_n;
}


void initialState(int N, int ** lattice_ptr){
    *lattice_ptr  = (int*)malloc(N*N*sizeof(int));
    int * lattice = *lattice_ptr;

    int i=0;
    int j=0;
    int rand;
    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
            lattice[i*N + j] = i*N+j;
            /*
            if( i>(j+4)){
                lattice[i*N + j] = +1;
            }else{
                lattice[i*N + j] = -1;
            */
            rand = rand_range(1,100);
            if( rand >= 70){
                lattice[i*N+j] = -1;
            }else{
                lattice[i*N+j] = +1;
            }
        }
    }
    return;
}
int * site_neighbors (int row, int column, int lattice_size, int** lattice_ptr){
        int * lattice = *lattice_ptr;

        int* neighbors = (int*)malloc(4*sizeof(int));

        int pbc_top_offset = 0;
        int pbc_bottom_offset = 0;
        int pbc_left_offset = 0;
        int pbc_right_offset = 0;

        int top_neighbor;
        int bottom_neighbor;
        int left_neighbor;
        int right_neighbor;
        int site_spin;
        // If we pick a site on the top or bottom,
        // enforce periodic boundary conditions
        if( row==0 ){
            // If top, move to down N rows
            pbc_top_offset = lattice_size;
        }

        if( row == lattice_size-1 ){
            // If bottom, move back up 3N rows
            pbc_bottom_offset = -lattice_size;
        }

        if( column == 0 ){
            // If left edge,
            pbc_left_offset = lattice_size;
        }

        if( column == lattice_size - 1 ){
            // If right edge,
            pbc_right_offset = -lattice_size;
        }

        top_neighbor    = lattice[(row - 1 + pbc_top_offset)*lattice_size
                        + (column + 0)];
        bottom_neighbor = lattice[(row + 1 + pbc_bottom_offset)*lattice_size
                        + (column + 0)];
        left_neighbor   = lattice[(row + 0 )*lattice_size
                        + (column - 1 + pbc_left_offset  )];
        right_neighbor  = lattice[(row + 0)*lattice_size
                        + (column + 1 + pbc_right_offset )];
        site_spin       = lattice[row*lattice_size+column];

        /*
        printf("Site and neighbor map for column=%d,row=%d\n",column,row);
        printf("xxx %03d xxx\n",top_neighbor );
        printf("%03d %03d %03d\n",left_neighbor, site_spin,right_neighbor );
        printf("xxx %03d xxx\n",bottom_neighbor );
        */

        neighbors[0] = top_neighbor;
        neighbors[1] = bottom_neighbor;
        neighbors[2] = left_neighbor;
        neighbors[3] = right_neighbor;



    return neighbors;

}
int total_spin(int lattice_size, int** lattice_ptr){
    int * lattice= *lattice_ptr;
    int i=0;
    int j=0;
    int spin_sum=0;
    for(i=0;i<lattice_size;i++){
        for(j=0;j<lattice_size;j++){
            spin_sum+=lattice[i*lattice_size+j];
            }
        }
    return spin_sum;
}

void flip_sites(int** lattice_ptr, int ** toFlip_ptr, int lattice_size){
    int* lattice = *lattice_ptr;
    int* toFlip = *toFlip_ptr;
    int index,column,row;
    int N_sqr=lattice_size*lattice_size;

    // Parallelizeable
    #pragma omp parallel
    {
        #pragma omp for
            for(index=0; index< N_sqr; index++){
                if(toFlip[index] != 0){
                    lattice[index] *= -1;
                }
            }
    }
}


void mcmove(int lattice_size, int** lattice_ptr,double temp, double external_field){
    int* lattice = *lattice_ptr;
    int i;
    // Spin of neighbor
    int * neighbors;
    int site_spin;
    int column, row;
    double energy;
    double energy_if_flip;
    double dE;

    double beta=1/(temp);
    double spin_interaction=1.0;
    double flip_prob;

    int id,nthreads;

    int* toFlip = (int*)calloc(lattice_size*lattice_size, sizeof(int));

    double randCompare;
    int shouldFlip;

    //printf("\n\nStarting mcmove for loop...\n\n");
    #pragma omp parallel
    {
    #pragma omp for collapse (2) private(column,row,neighbors,site_spin,energy,energy_if_flip,dE,flip_prob,shouldFlip,randCompare)
    for(column=0; column < lattice_size; column++ )
        for(row=0; row < lattice_size; row++){

        if (rand_range(1,100) < 80)
            continue;

        // Get neighbors of site
        neighbors=site_neighbors(column,row,lattice_size,lattice_ptr);


        site_spin = lattice[column*lattice_size + row];

        energy = (-1*spin_interaction*site_spin*(neighbors[0]+ neighbors[1] +
               + neighbors[2] + neighbors[3])+ site_spin*external_field);

        energy_if_flip = (-1*spin_interaction*-1*site_spin*(neighbors[0]+ neighbors[1] +
               + neighbors[2] + neighbors[3])+ -1*site_spin*external_field);
        dE = energy_if_flip-energy;
        // If flipping the spin is a lower energy state, change site_spin
        // OR
        // If flipping the spin increases the energy, flip randomly
        // Randomness based on exp(-cost*beta), note that beta is 1/T, where
        // T is the actual Temp, not the most intuitive
        randCompare = (double) rand()/(double)RAND_MAX;
        flip_prob   = exp(-dE * beta);
        shouldFlip  = randCompare < flip_prob;
/*
        if(dE >0){

        printf("dE=%f, randCompare=%f flip_prob=%f shouldFlip=%d\n",dE,randCompare,flip_prob,shouldFlip );
        }
*/
        if(dE < 0 || shouldFlip ){

            // Track flip status
            toFlip[column*lattice_size + row]=1;

            //printf("Logging site [%03d,%03d]\n", column,row);
        }

        free(neighbors);

        }
    }

    flip_sites(lattice_ptr, &toFlip, lattice_size);
    free(toFlip);
    return;
}

void printLattice(int lattice_size, int ** lattice_ptr, double * temp_ptr,
                  int time_step, hid_t * file_id_ptr){
    int * lattice=*lattice_ptr;
    int i,j,site_spin;

    hid_t file_id = *file_id_ptr;
    // Dataset properties
    hid_t dataset_id, dataspace_id;
    hsize_t dset_dims[2];

    // Dataset attribute properties
    hid_t attribute_id, a_dataspace_id;
    hsize_t a_dims[1];

    herr_t status;

    char dataset_name[20];

    sprintf(dataset_name,"ising_step_%04i",time_step);
    dset_dims[0]=lattice_size;
    dset_dims[1]=lattice_size;
    a_dims[0] = 1;



    // Setup dataspace and datasets
    //printf("Setting up dataspce and datasets\n");
    dataspace_id = H5Screate_simple(2,dset_dims,NULL);

    dataset_id = H5Dcreate(file_id, dataset_name, H5T_NATIVE_INT, dataspace_id,
                           H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);


    a_dataspace_id = H5Screate_simple(1,a_dims,NULL);

    attribute_id = H5Acreate2(dataset_id,"Temperature", H5T_NATIVE_DOUBLE,
                              a_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);

    // Write dataset and attribute
    //printf("Writing dataset and attribute\n");
    status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
                      H5P_DEFAULT,lattice );

    //printf("Status of H5Dwrite is %d\n",status );

    status = H5Awrite(attribute_id, H5T_NATIVE_DOUBLE,temp_ptr);

    //printf("Status of H5Awrite is %d\n",status );
    //

    //printf("\nClosing datasets, dataspaces, ......\n");

    status = H5Dclose(dataset_id);
    //printf("Status of H5Dclose is %d\n",status );

    status = H5Sclose(dataspace_id);
    //printf("Status of H5Sclose is %d\n",status );

    status = H5Sclose(a_dataspace_id);
    //printf("Status of H5Sclose is %d\n",status );

    status = H5Aclose(attribute_id);
    //printf("Status of H5Aclose is %d\n",status );

    //printf("Leaving printLattice()\n");
    return;
}
