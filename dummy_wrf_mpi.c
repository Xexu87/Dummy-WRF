#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <mpi.h>

#define MASTER 0
#define EXTERN -1

#define UP 0
#define DOWN 1
#define LEFT 2
#define RIGHT 3

#define VERBOSE 1

int my_rank, nprocs;

typedef struct point {
    float x;
    float y;
    float z;
    float t;
} point;

/*
    This function decomposes sets the process grid given the number of processes.
    @input param num: Number of processes.
    @output para y: The number of Rows. Cols can be computed outside with num/ROWS.
 */
int partition(int num){
    int a,b;
    int x = 1;
    int y = num;
    int i;

    for (i = 1; i < num; i++){
        if (num % i == 0){
            a = i;
            b = num / i;
            if(  (b-a)*(b-a) < (x-y)*(x-y) ){
                x = a;
                y = b;
            }
        }
    }
    return y;
}

/*
    This function prints the matrix size.
 */
void print_size(unsigned long long int s){
    int i = 0;
    unsigned long long int size = s;

    printf("\nSize = %d Bytes.", size);

    while ( size >= 1024){
        size = size / 1024;
        i = i + 1;
    }
    switch(i){
        case 0:
            printf("\nSize = %d Bytes.", size); break;
        case 1:
            printf("\nSize = %d KB.", size);break;
        case 2:
            printf("\nSize = %d MB.", size);break;
        case 3:
            printf("\nSize = %d GB.", size);break;
        case 4:
            printf("\nSize = %d TB.", size);break;
        default:
            printf("\nSomething went wrong. Over TB ?");
    }
    printf("\n\n");
}

/*
    This function sets the neighbours of a given process.
    V [  UP,  DOWN,  LEFT,  RIGHT  ]
    -1    : means NO neigh --extern border-- no computable
    value : means process neigh
 */
void set_neighs(int* V, int pid, int cols, int rows){
    if(pid < cols)            {V[UP]  = EXTERN;} else {V[UP]    = pid - cols;}
    if(pid >= (cols)*(rows-1)){V[DOWN]= EXTERN;} else {V[DOWN]  = pid + cols;}
    if(pid % cols == 0)       {V[LEFT]= EXTERN;} else {V[LEFT]  = pid - 1;}
    if((pid+1) % cols == 0)   {V[RIGHT]=EXTERN;} else {V[RIGHT] = pid + 1;}
}

// DEGUB ONLY
void print_neighs(int pid, int* V){
    printf("\nProcess %d, UP = %d, DOWN = %d, LEFT = %d, RIGHT = %d\n", pid, V[UP], V[DOWN], V[LEFT], V[RIGHT]);
}

point* generate_matrix(int W, int H){
    point* aux;
    aux = (point*)malloc(sizeof(point) * W * H);
    return aux;
}

void init_matrix(point* matrix, int W, int H) {
    int L = 10;
    for (int i = 0; i < H; i++){
        for (int j = 0; j < W; j++){
            matrix[(i*W)+j].x = lrand48() % L;
            matrix[(i*W)+j].y = lrand48() % L;
            matrix[(i*W)+j].z = lrand48() % L;
            matrix[(i*W)+j].t = lrand48() % L;
        }
    }
}

// DEGUB ONLY
void print_matrix(point* m, int W, int H){
    for(int i = 0; i < H; i++){
        for(int j = 0; j < W; j++){
            printf("\033[0;31m %.2f", i,j,m[i*W + j].x);
        }
        printf("\n");
    }
}

// DEGUB ONLY
void print_work_rank(point* m, int rank, int W, int H){
    if (rank == 0){printf("\033[0;32m\nMASTER  \n");}
    if (rank == 1){printf("\033[0;37m\nWORKER_1\n");}
    if (rank == 2){printf("\033[0;35m\nWORKER_2\n");}
    if (rank == 3){printf("\033[0;34m\nWORKER_3\n");}

    for(int i = 0; i < H ; i++){
        for (int j = 0; j < W ; j++){
            if (rank == MASTER){
                printf("\033[0;32mMA[%d][%d]%.2f  - ", i,j,m[i*W + j].x);
            }
            if (rank == 1){
                printf("\033[0;37mW1[%d][%d]%.2f  - ", i,j,m[i*W + j].x);
            }
            if (rank == 2){
                printf("\033[0;35mW2[%d][%d]%.2f  - ", i,j,m[i*W + j].x);
            }
            if (rank == 3){
                printf("\033[0;34mW3[%d][%d]%.2f  - ", i,j,m[i*W + j].x);
            }
        }
        printf("\n");
    }
}

void init_work(point* matrix, int W, int H) {
    int i,j;
    for (i = 0; i < H; i++){
        for (j = 0; j < W; j++){
            matrix[(i*W)+j].x = 0;
            matrix[(i*W)+j].y = 0;
            matrix[(i*W)+j].z = 0;
            matrix[(i*W)+j].t = 0;
        }
    }
}


/*
    Stencil of average of 4-neigh
    -  X  -
    X  -  X
    -  X  -
 */
float stencil_x(point* u, point* d, point* l, point* r){
    float tmp = ((u->x + d->x) + (l->x + r->x)) / 4.0;
    return tmp;
}

/*
    Stencil of average of 8-neigh and cell
    X  X  X
    X  X  X
    X  X  X
 */
float stencil_y(point* ul, point* uc, point* ur, 
                point* cl, point*  c, point* cr, 
                point* dl, point* dc, point* dr){
    float tmp = ( (ul->y + uc->y + ur->y)+
                  (cl->y +  c->y + cr->y)+
                  (dl->y + dc->y + dr->y));
    return tmp/9.0;
}

/*
    Stencil of average of 4-neigh
    -  X  -
    X  X  X
    -  X  -
 */
float stencil_z(point* c, point* u, point* d, point* l, point* r){
    float tmp = ( c->z + u->z + d->z + l->z + r->z) / 5.0;
    return tmp;
}

/*
    Average of the 3 performed stencils in x,y,z.
 */
float stencil_t(point* c){
    float average = (c->x + c->y + c->z) / 3;
    return average;
}


/*
    Each step, {x,y,z,t} from each cell are being updated/computed.
    x -> 4-stencil
    y -> 9-stencil
    z -> 5-stencil
    t -> average of previous 3.
 */
point* step(point* m, point* aux, int R, int W, int H){
    for (int i = 1; i < (H-1); i++){
        for (int j = 1; j < (W-1); j++){
            int up    = ((i-1)*W) + j;
            int down  = ((i+1)*W) + j;
            int left  =  (i * W) + j - 1;
            int right =  (i * W) + j + 1;
            aux[(i*W)+j].x = stencil_x(&m[up],  &m[down],  &m[left],  &m[right]);
        }
    }
}

/*
    Calcula el promig dels elements de la diagonal X
 */
float checksum(point* m, int W, int H){
    int lower = W;
    if (H < W){
        lower = H;
    }

    float sum = 0.0;

    for (int i = 0; i < lower; i++){
        sum = sum + (m[((i)*W)+i].x/lower);
    }
    return sum;
}

void print_welcome(){
    printf("\033[0;34m");
    printf("===============================================================\n");
    printf("\033[0;32m");
    printf("      _                                                     __ \n");
    printf("   __| |_   _ _ __ ___  _ __ ___  _   _     __      ___ __ / _|\n");
    printf("  / _` | | | | '_ ` _ \\| '_ ` _ \\| | | |____\\ \\ /\\ / / '__| |_ \n");
    printf(" | (_| | |_| | | | | | | | | | | | |_| |_____\\ V  V /| |  |  _|\n");
    printf("  \\__,_|\\__,_|_| |_| |_|_| |_| |_|\\__, |      \\_/\\_/ |_|  |_|  \n");
    printf("                                  |___/\n\n");
    printf("\033[0;34m");
    printf("===============================================================\n");
    printf("\033[0m");
}

void main(int argc, char **argv){
    
    // INIT MPI
    // ----------------------------------------------------------------------
    MPI_Init (&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    MPI_Request request;

    if(nprocs < 2){
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if (my_rank == MASTER){
        print_welcome();
        printf("\nUSAGE: wrf.exe [W] [H] [iters] [seed]\n");
    }

    // INIT DATA
    // ----------------------------------------------------------------------
    int W=2048, H=1024, Z=50, S=1;
    int verbose = VERBOSE;
    unsigned long long int size;

    if (argc>1) {  W = atoi(argv[1]); }
    if (argc>2) {  H = atoi(argv[2]); }
    if (argc>3) {  Z = atoi(argv[3]); }
    if (argc>4) {  S = atoi(argv[4]); }
    srand48((long int) S);

    // Create time variables
    struct timeval start, end;

    // Compute how many processes in X,Y direcctions. GRID PROCESSES.
    int n_rows = partition(nprocs);
    int n_cols = nprocs / n_rows;

    // Compute length of work in X,Y for each process [will be Work matrix size]
    int w_size = W  / n_cols;
    int h_size = H  / n_rows;

    // Vector to store processes indexes (MASTER)
    int I[nprocs];

    // Compute offsets for each process 
    int start_w, start_h, end_w, end_h;
    int this_process = 0;

    for (int i = 0; i < n_rows; i++){
        for (int j = 0; j < n_cols; j++){
            if (this_process == my_rank){
                start_w = j * w_size;
                start_h = i * h_size;
                end_w   = (j + 1) * w_size;
                end_h   = (i + 1) * h_size;
            }
            if(my_rank == MASTER){
                start_w = j * w_size;
                start_h = i * h_size;
                end_w   = (j + 1) * w_size;
                end_h   = (i + 1) * h_size;
                I[this_process] = start_h * W + start_w; // STORE index offset of this process
            }
            this_process += 1;
        }        
    }
    if(my_rank == MASTER){
        start_w = 0;
        start_h = 0;
        end_w   = w_size;
        end_h   = h_size;
    }

    // VERBOSE INFO
    // ----------------------------------------------------------------------
    if (verbose && my_rank == MASTER){
        printf("\nMatrix size = %d x %d.", W, H);
        size = W * H * sizeof(point);
        print_size(size);
        printf("\nComputing %d iterations", Z);fflush(stdin);
        printf("\nPartition MPI is: %d rows with %d columns per row.\n", n_rows, n_cols);
    }

    // Allocate memory for matrix and work. 
    /*
        Improvement --> Only master needs to have matrix allocated.
        Slaves only require work.
        Since matrix size is not big enough, it won't cause any performance lack.
     */
    point* matrix;
    point* work;
    point* aux;
    point* tmp;
   
    matrix = generate_matrix(W,H);
    work   = generate_matrix(2 + w_size, 2 + h_size);
    aux    = generate_matrix(2 + w_size, 2 + h_size);

    init_matrix(matrix, W, H);
    init_work(work, w_size+2, h_size+2);

    int* neighs = (int*)malloc(4 * sizeof(int));
    set_neighs(neighs, my_rank, n_cols, n_rows);

    if (1 && my_rank == MASTER){
        //print_neighs(my_rank, neighs);
    }

    if (my_rank == MASTER){    
        gettimeofday(&start, NULL);
    }

    // Create the point struct class MPI
    // https://stackoverflow.com/questions/9864510/struct-serialization-in-c-and-transfer-over-mpi
    const int nitems = 4;
    int blocklenghts[4] = {1,1,1,1};
    MPI_Datatype types[4] = {MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT};
    MPI_Datatype mpi_point_type;
    MPI_Aint offsets[4];

    offsets[0] = offsetof(point, x);
    offsets[1] = offsetof(point, y);
    offsets[2] = offsetof(point, z);
    offsets[3] = offsetof(point, t);

    MPI_Type_create_struct(nitems, blocklenghts, offsets, types, &mpi_point_type);
    MPI_Type_commit(&mpi_point_type);

    // Create the VECTOR COLUMN datatype to send HALO columns.
    MPI_Datatype column_type;
    int row = w_size + 2;
    MPI_Type_vector(h_size, 1, row, mpi_point_type, &column_type);
    MPI_Type_commit(&column_type);

    // Set initial conditions
    int offset = start_h * W + start_w;

    for (int i = 0; i < h_size; i++){
        for (int j = 0; j < w_size; j++){
            work[ ((1 + i)*(w_size + 2))+ 1 + j ] = matrix[offset + j];
        }
        offset = offset + W;
    }

    MPI_Barrier(MPI_COMM_WORLD);  

    // MAIN WORK
    // =========================================================
    for (int t = 0; t < Z; t++){

        MPI_Barrier(MPI_COMM_WORLD);

        // UPDATE HALO
        if ( neighs[UP] != EXTERN){
            offset = w_size + 2;
            MPI_Isend(&work[offset], row, mpi_point_type, neighs[UP],  0, MPI_COMM_WORLD, &request);

            offset = 0;
            MPI_Recv(&work[offset], row, mpi_point_type, neighs[UP], 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        if ( neighs[DOWN] != EXTERN){
            offset = (w_size + 2) * (h_size );
            MPI_Isend(&work[offset], row, mpi_point_type, neighs[DOWN], 1, MPI_COMM_WORLD, &request);

            offset = (w_size + 2) * (h_size + 1);
            MPI_Recv(&work[offset], row, mpi_point_type, neighs[DOWN],  0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
        }

        if ( neighs[LEFT] != EXTERN){
            offset = 1 + row;
            MPI_Isend(&work[offset], 1, column_type, neighs[LEFT],  2, MPI_COMM_WORLD, &request );

            offset = row;
            MPI_Recv(&work[offset], 1, column_type, neighs[LEFT], 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
        }

        if ( neighs[RIGHT] != EXTERN){
            offset = (w_size + 2) + w_size;
            MPI_Isend(&work[offset], 1, column_type, neighs[RIGHT], 3, MPI_COMM_WORLD, &request );

            offset = (w_size + 2) + w_size + 1;
            MPI_Recv(&work[offset], 1, column_type, neighs[RIGHT],  2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        MPI_Barrier(MPI_COMM_WORLD);

        step(work, aux, my_rank, 2 + w_size, 2 + h_size);

        // Aux matrix becomes work matrix. Avoid Copy.
        tmp  = work;
        work = aux;
        aux  = tmp;
        
        // SYNC
        MPI_Barrier(MPI_COMM_WORLD);
    }

    //  ========================================================= END ====
    
    // RETRIEVE FINAL MATRIX
    int block = (2 + w_size)*(2 + h_size);
    
    // ALL PROCESSES SEND THEIR WORK AND INDEX TO MASTER
    // Index is required to set the location of the work inside the matrix.
    if(my_rank != MASTER){       
        MPI_Send(&work[0], block, mpi_point_type, MASTER, my_rank, MPI_COMM_WORLD);
    }

    // MASTER RECEIVE WORK AND FULLFIL FINAL MATRIX 
    else{        
        // MASTER'S CHUNK  
        offset = 0;
        for (int i = 0; i < h_size; i++){
            for (int j = 0; j < w_size; j++){
                matrix[offset + j] = work[ ((1 + i)*(w_size + 2))+ 1 + j ];
            }
            offset = offset + W;
        }

        // SLAVES' CHUNK          
        for(int p = 1; p < nprocs; p++){
            offset = I[p]; // GET INDEX OF THIS PROCESS (offset of where his work start relative to matrix)
            MPI_Recv(&work[0],  block, mpi_point_type, p, p, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int i = 0; i < h_size; i++){
                for (int j = 0; j < w_size; j++){
                    matrix[offset + j] = work[ ((1 + i)*(w_size + 2))+ 1 + j ];
                }
                offset = offset + W;
            }
        }            
    }
    // ========================================================= 

    MPI_Barrier(MPI_COMM_WORLD);

    if (my_rank == MASTER){
        gettimeofday(&end, NULL);
        end.tv_sec = end.tv_sec - start.tv_sec;
        end.tv_usec = end.tv_usec - start.tv_usec;
        if (end.tv_usec < 0) {  end.tv_usec += 1000000; end.tv_sec--; }
        printf("\nTotal Time is %ld seconds and %ld microseconds\n", end.tv_sec, end.tv_usec);    
        printf("\n\nCheksum is %f\n", checksum(matrix, W, H));
    }

    if (my_rank == MASTER){
        printf("\nFreeing Memory.");
        printf("\nDone.\n\n");
    }

    free(matrix);
    free(work);
    free(aux);
    MPI_Type_free(&mpi_point_type);
    MPI_Finalize();
}
