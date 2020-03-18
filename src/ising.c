/*************************************************
*************************************************/

#include "ising_functions.h"

void usage(char* program_name){
   printf("Usage is %s [options]\n", program_name);
   printf("Options:\n");
   printf("-N  <lattice_size> Sets N for NxN ising lattice\n");
   printf("-t  <time_steps>   Sets number of time steps\n");
   printf("-i <temp_initial>  Sets initial temp of lattice\n");
   printf("-f <temp_final>    Sets final temp of lattice\n");
   printf("-B  <field_mag>    Sets external field strength\n");
   printf("-h                 Displays this message       \n\n");
}

void set_args(int argc, char* argv[], int* nt, unsigned int* N,
              double* external_field, double* initial_temp, double* final_temp){
   int ntIsSet=0;
   int NIsSet=0;
   int external_fieldIsSet=0;
   int initial_tempIsSet=0;
   int final_tempIsSet=0;
   char* name=argv[0];
   // Parse arguments
   while( (argc>1) && (argv[1][0] =='-')){
   //argv[1][1] is the option character
      switch(argv[1][1]){
         case 'N':
            // Set lattice size
            *N=atoi(argv[2]);
             NIsSet=1;
            break;
         case 't':
            // Set num time steps
            *nt=atoi(argv[2]);
             ntIsSet=1;
            break;
         case 'i':
            // Set initial lattice temp
            *initial_temp=atoi(argv[2]);
             initial_tempIsSet=1;
            break;
         case 'f':
            // Set final lattice temp
            *final_temp=atoi(argv[2]);
             final_tempIsSet=1;
            break;
         case 'B':
            // Set N
            *external_field=atoi(argv[2]);
             external_fieldIsSet=1;
            break;
         case 'h':
            // Display help message
            usage(name);
            exit(0);
         default:
            printf("Bad option %s\n", argv[1]);
            usage(name);
            exit(0);
            break;
      }
      argv+=2;
      argc-=2;

   }


 // Check if defaults were overridden, if not use default values
   if(!ntIsSet){
      *nt=1000;
   }
   if(!NIsSet){
      *N=128;
   }
   if(!external_fieldIsSet){
      *external_field=0.02;
   }
   if(!initial_tempIsSet){
      *initial_temp=0.5;
   }
   if(!final_tempIsSet){
      *final_temp=0.0;
   }
}
// Set problem size
int main(int argc, char* argv[]){
    int nt;       // Number of temp steps
    unsigned int N ;       //Size of lattice, N x N
    int id, nthreads;


    double  initial_temp;
    double  final_temp;
    double  external_field;

   set_args(argc, argv, &nt, &N,
              &external_field, &initial_temp, &final_temp);
 // Diagnostic info to check that args are set
   printf("N is %d\n", N);
   printf("nt is %d\n", nt);
   printf("initial_temp is %f\n", initial_temp);
   printf("final_temp is %f\n", final_temp);
   printf("external_field is %f\n", external_field);


    double  current_temp = initial_temp;
    int seed = time(NULL);
    int *lattice = NULL;
    double  dt;

    // HDF5 objects
    hid_t file_id;
    herr_t status;

    file_id = H5Fcreate("data.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    dt = (final_temp - initial_temp)/nt;
    srand(seed);
    initialState(N,&lattice);
    printLattice(N,&lattice,&current_temp,0,&file_id);
    int i;
    #pragma omp parallel private(id)
    {
    id = omp_get_thread_num();
    printf("Thread %d starting...\n", id);

    if(id==0) {
            nthreads = omp_get_num_threads();
            printf("Number of threads = %d \n", nthreads);
        }
    }

    for(i=1; i<nt; i++){
        mcmove(N,&lattice, current_temp,external_field);
        current_temp += dt;

        if( i%(2) || i==(nt-1)){
            //printf("Writing timestep %i\n", i);
            printLattice(N,&lattice,&current_temp,i,&file_id);
        }
    }
    free(lattice);
    status = H5Fclose(file_id);
    return 0;
}
