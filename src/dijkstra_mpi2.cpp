#include <iostream> 
#include <cassert>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstdarg>
#include <cstring>
#include <iomanip>
#include <mpi.h>

using namespace std;

const int INF = 0xffff; 
int my_rank, processes;
MPI_Comm comm;

void log(bool every_process, const char *fmt, ...);
void panic(const string &text);
MPI_Datatype make_blk_col_type(int n, int my_n);
void Dijkstra(double *my_mat, double *my_dist_to, int *my_point_to, int my_n, int g_n, MPI_Comm comm);

int main(int argc, char **argv) {
    double *my_mat, *my_dist_to, *g_dist_to;
    int *my_point_to, *g_point_to;
    int my_n, g_n;
    ifstream fis;
    // MPI_Datatype blk_col_mpi_t;

    // 
    // Init the MPI system
    // 
    MPI_Init(nullptr, nullptr);

    // 
    // the comm for for all the processes, 
    // so it is MPI_COMM_WORLD
    // 
    comm = MPI_COMM_WORLD;

    // 
    // Get how many processes we start
    // 
    MPI_Comm_size(comm, &processes);
    // 
    // Get my rank (the id from 0 ~ processes - 1)
    // 
    MPI_Comm_rank(comm, &my_rank);

    //
    // Input style
    //
    if (argc != 3) {
        panic("Please use ./dijkstra_mpi [resources name] [out_file_name]");
    }

    std::string out_file_name(argv[2]);
    // 
    // Read n:
    // Process 0 should read the n
    // and boardcast to every other process. 
    // 
    if (my_rank == 0) {
        fis.open(argv[1]);
        fis >> g_n;
    }
    MPI_Bcast(&g_n, 1, MPI_INT, 0, comm);
    log(true, "-- My rank is %d, the processes is %d\b", my_rank, processes);

    // 
    // my_n is the number of data this process should deal with
    //
    if (g_n < processes) {
        panic("processes number is too large");
    } else if (g_n % processes != 0) {
        panic("g_n % process is not 0");
    }
    my_n = g_n / processes;
    log(true, "my_n is %d, g_n is %d", my_n, g_n);

    // 
    // Alloc the space for local data 
    // 
    my_mat = (double *) calloc(g_n * my_n, sizeof(double));
    my_dist_to = (double *) calloc(my_n, sizeof(double));
    my_point_to = (int *) calloc(my_n, sizeof(int));

    //
    // make the mpi block col type
    //
    // blk_col_mpi_t = make_blk_col_type(g_n, my_n);
    
    // 
    // Allocate the global array
    // 
    if (my_rank == 0) {
        g_dist_to = (double *) calloc(g_n, sizeof(double));
        g_point_to = (int *) calloc(g_n, sizeof(int));
    }

    //
    // Read the metrix and send to every processes..
    //
    double *mat = NULL;
    if (my_rank == 0) {
        mat = (double *) calloc(g_n * g_n, sizeof(double));
        for (int i = 0; i < g_n * g_n; ++i) {
            mat[i] = INF;
        }
        for (int i = 0; i < g_n; ++i) {
            mat[i * g_n + i] = 0;
        }

        int edges;
        fis >> edges;
        // 
        // read the edges 
        // 
        for (int i = 0; i < edges; ++i) {
            int from, to;
            double distance;
            fis >> from >> to >> distance;  

            mat[from * g_n + to] = distance;
        }
        {
            ostringstream oss;
            int i, j;
            oss << "\n";
            oss << "  Distance matrix:\n";
            oss << "\n";
            oss << "\t";
            for (i = 0; i < g_n; i++)
                oss << "\t[" << i << "]";
            oss << endl;
            for (i = 0; i < g_n; i++) {
                oss << "\t[" << i << "]";
                for (j = 0; j < g_n; j++) {
                    if (mat[(i) * g_n + j] == INF) {
                        oss << "\tInf";
                    } else {
                        oss << "\t" << mat[(i) * g_n + j];
                    }
                }
                oss << "\n";
            }
            // std::cout << oss.str() << std::endl;
        }
    }

    log(false, "successfully reading data");


    //
    // Sends data from one process to all other processes in a communicator
    //
    // int MPI_Scatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
    //          void *recvbuf, int recvcount, MPI_Datatype recvtype, int root,
    //          MPI_Comm comm)
    //
    MPI_Scatter(mat, my_n * g_n, MPI_DOUBLE, my_mat, g_n * my_n, MPI_DOUBLE, 0, comm);
    log(false, "successfully scatter data");
    //
    // Log the mat recv
    //
    // {
    //     ostringstream oss;
    //     int i, j;
    //     oss << "\n";
    //     oss << "  Distance matrix:\n";
    //     oss << "\n";
    //     oss << "\t ";
    //     for (i = 0; i < g_n; i++)
    //         oss << "\t[" << i << "]";
    //     oss << endl;
    //     for (i = 0; i < my_n; i++) {
    //         oss << "\t[" << my_rank + i << "]";
    //         for (j = 0; j < g_n; j++) {
    //             if (my_mat[i * g_n + j] == -1) {
    //                 oss << "\tInf";
    //             } else {
    //                 oss << "\t" << my_mat[i * g_n + j];
    //             }
    //         }
    //         oss << "\n";
    //     }

    //     log(true, "%s", oss.str().c_str());
    // }
    if (my_rank == 0) {
        free(mat);
    }

    log(false, "successfully print data!");

    // Start dijkstra
    Dijkstra(my_mat, my_dist_to, my_point_to, my_n, g_n, comm);

    // Getter the data together...
    MPI_Gather(my_dist_to, my_n, MPI_DOUBLE, g_dist_to, my_n, MPI_DOUBLE, 0, comm);

    // Show the result
    if (my_rank == 0) {
        cout << "\n";
        cout << "  Minimum distances from node 0:\n";
        cout << "\n";
        for (int i = 0; i < g_n; i++) {
          cout << "  " << setw(2) << i << "  " << setw(2) << g_dist_to[i] << "\n";
        }
    }

    //
    // write to file... 
    //
    if (my_rank == 0) {
        std::ofstream fos(out_file_name);
        for (int i = 0; i < g_n; ++i) {
            fos << i << ": " << g_dist_to[i] << endl;
        }
    }

    
    //
    // Free the type we make 
    //
    // MPI_Type_free(&blk_col_mpi_t);
    //
    // Quit the mpi system
    //
    MPI_Finalize();
    return 0;
}

//
// This funcion should called before 
// MPI_init
//
void panic(const string &text) {
    log(false, "=== Error: %s", text.c_str());
    MPI_Finalize();
    exit(1);
}

void log(bool every_process, const char *fmt, ...) {
    char buf[4096] = { 0 };
    va_list args;
    va_start(args, fmt);
    vsprintf(buf, fmt, args);
    
    if (every_process) {
        if (my_rank == 0) {
            cout << "[0] " << buf << endl;
            for (int i = 1; i < processes; ++i) {
                char buf_from_other_process[4096] = { 0 };
                MPI_Recv(buf_from_other_process, 4096, MPI_CHAR, i, 0, comm, MPI_STATUS_IGNORE);
                printf("[%d] %s\n", i, buf_from_other_process);
            }
        } else {
            MPI_Send(buf, 4096, MPI_CHAR, 0, 0, comm);
        }
    } else {
        if (my_rank == 0) {
            cout << "[singal] " << buf << endl;
        }
    }
    va_end(args);
}

MPI_Datatype make_blk_col_type(int n, int my_n) {
    MPI_Aint lb, extent;
    MPI_Datatype block_mpi_t;
    MPI_Datatype first_bc_mpi_t;
    MPI_Datatype blk_col_mpi_t;

    // 
    // Creates a contiguous datatype
    // Here I think it is a int array, named block_mpi_t. 1
    //
    MPI_Type_contiguous(my_n, MPI_INT, &block_mpi_t);

    //
    // Get the lower bound and extent for a Datatype
    // Example:
    //  if my_n is 100, and the bytes of int is 4,
    //  so extent is 400.
    //
    MPI_Type_get_extent(block_mpi_t, &lb, &extent);
    log(true, "the lower bound is %ld, the extent is %ld", lb, extent);

    // 
    // Create a vector of block_mpi_t.
    // 
    /* MPI_Type_vector(numblocks, elts_per_block, stride, oldtype, *newtype) */
    MPI_Type_vector(n, my_n, n, MPI_INT, &first_bc_mpi_t);

    /* This call is needed to get the right extent of the new datatype */
    MPI_Type_create_resized(first_bc_mpi_t, lb, extent, &blk_col_mpi_t);

    // 
    // commit the type we make 
    //
    MPI_Type_commit(&blk_col_mpi_t);

    MPI_Type_free(&block_mpi_t);
    MPI_Type_free(&first_bc_mpi_t);

    return blk_col_mpi_t;
}

// void Dijkstra_Init(double loc_mat[], int loc_pred[], double loc_dist[], int loc_known[],
//                    int my_rank, int loc_n) {
//     int loc_v;
// 
//     if (my_rank == 0)
//         loc_known[0] = 1;
//     else
//         loc_known[0] = 0;
// 
//     for (loc_v = 1; loc_v < loc_n; loc_v++)
//         loc_known[loc_v] = 0;
// 
//     for (loc_v = 0; loc_v < loc_n; loc_v++) {
//         loc_dist[loc_v] = loc_mat[0 * loc_n + loc_v];
//     }
// }
// 
// int Find_min_dist(double loc_dist[], int loc_known[], int loc_n) {
//     int loc_u = -1, loc_v;
//     double shortest_dist = INF;
// 
//     for (loc_v = 0; loc_v < loc_n; loc_v++) {
//         if (!loc_known[loc_v]) {
//             if (loc_dist[loc_v] < shortest_dist) {
//                 shortest_dist = loc_dist[loc_v];
//                 loc_u = loc_v;
//             }
//         }
//     }
//     return loc_u;
// }
// 
// 
// 
// 
// void Dijkstra(double loc_mat[], double loc_dist[], int loc_pred[], int loc_n, int n,
//               MPI_Comm comm) {
// 
//     int i, loc_v, loc_u, glbl_u, new_dist, my_rank;
//     double dist_glbl_u;
//     int *loc_known;
//     char my_min[12];
//     char glbl_min[12];
//     double my_min0, glbl_min0;
//     int my_min1, glbl_min1;
//     // int my_min[2];
//     // int glbl_min[2];
// 
//     loc_known = (int *) malloc(loc_n * sizeof(int));
// 
//     Dijkstra_Init(loc_mat, loc_pred, loc_dist, loc_known, my_rank, loc_n);
// 
//     /* Run loop n - 1 times since we already know the shortest path to global
//        vertex 0 from global vertex 0 */
//     for (i = 0; i < n - 1; i++) {
//         loc_u = Find_min_dist(loc_dist, loc_known, loc_n);
// 
//         if (loc_u != -1) {
//             // my_min[0] = loc_dist[loc_u];
//             my_min0 = loc_dist[loc_u];
//             // my_min[1] = loc_u + my_rank * loc_n;
//             my_min1 = loc_u + my_rank * loc_n;
//         }
//         else {
//             my_min0 = INF;
//             my_min1 = -1;
//         }
// 
//         /* Get the minimum distance found by the processes and store that
//            distance and the global vertex in glbl_min
//         */
//         memcpy(my_min, &my_min0, sizeof(double));
//         memcpy(my_min + sizeof(double), &my_min1, sizeof(int));
//         // int MPI_Allreduce(const void *sendbuf, void *recvbuf, int count,
//         //          MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
//         MPI_Allreduce(my_min, glbl_min, 1, MPI_DOUBLE_INT, MPI_MINLOC, comm);
//         memcpy(&glbl_min0, glbl_min, sizeof(double));
//         memcpy(&glbl_min1, glbl_min + sizeof(double), sizeof(int));
//         
//         // MPI_Allreduce(my_min, glbl_min, 1, MPI_2INT, MPI_MINLOC, comm);        
//         log(true, "my_min0 is %lf and my_min1 is %d", my_min0, my_min1);
//         log(false, "glbl_min0 is %lf and glbl_min1 is %d", glbl_min0, glbl_min1);
// 
//         dist_glbl_u = glbl_min[0];
//         glbl_u = glbl_min[1];
// 
//         /* This test is to assure that loc_known is not accessed with -1 */
//         if (glbl_u == -1)
//             break;
// 
//         /* Check if global u belongs to process, and if so update loc_known */
//         if ((glbl_u / loc_n) == my_rank) {
//             loc_u = glbl_u % loc_n;
//             loc_known[loc_u] = 1;
//         }
// 
//         // log for the myself connect node
//         {
//             std::ostringstream oss;
//             for (int i = 0; i < loc_n; ++i) {
//                 oss << loc_known[i] << " ";
//             }
//             log(true, "connect: %s", oss.str().c_str());
//         }
// 
//         /* For each local vertex (global vertex = loc_v + my_rank * loc_n)
//            Update the distances from source vertex (0) to loc_v. If vertex
//            is unmarked check if the distance from source to the global u + the
//            distance from global u to local v is smaller than the distance
//            from the source to local v
//          */
//         for (loc_v = 0; loc_v < loc_n; loc_v++) {
//             if (!loc_known[loc_v]) {
//                 new_dist = dist_glbl_u + loc_mat[glbl_u * loc_n + loc_v];
//                 if (new_dist < loc_dist[loc_v]) {
//                     loc_dist[loc_v] = new_dist;
//                 }
//             }
//         }
//     }
//     free(loc_known);
// }



void Dijkstra(double *my_mat, double *my_dist_to, int *my_point_to, int my_n, int g_n, MPI_Comm comm) {
    double my_md, g_md;
    int my_mv, g_mv;
    int *my_connected; 

    // init my connect
    my_connected = (int *) calloc(my_n, sizeof(int));
    if (my_rank == 0) {
        my_connected[0] = 1;
    } 
    
    //
    // Init my_dist_to to one step distance
    //
    for (int i = 0; i < my_n; ++i) {
        my_dist_to[i] = my_mat[0 * my_n + i];
    }

    //
    // dijkstra
    //
    for (int i = 0; i < g_n - 1; ++i) {
        //     
        // Find min dist
        //     
        my_md = INF; // shortest_dist
        my_mv = -1; // loc_u
        for (int i = 0; i < my_n; ++i) {
            if (!my_connected[i] && my_dist_to[i] < my_md) {
                my_md = my_dist_to[i];
                my_mv = i + my_rank * my_n;
            }
        }

        //
        // make an Allreduce operation 
        //
        char my_data[12], g_data[12];
        memcpy(my_data, &my_md, sizeof(double));
        memcpy(my_data + sizeof(double), &my_mv, sizeof(int));
        // int MPI_Allreduce(const void *sendbuf, void *recvbuf, int count,
        //          MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
        MPI_Allreduce(my_data, g_data, 1, MPI_DOUBLE_INT, MPI_MINLOC, comm);
        memcpy(&g_md, g_data, sizeof(double));
        memcpy(&g_mv, g_data + sizeof(double), sizeof(int));

        // log(true, "my_md is %lf and my_mv is %d", my_md, my_mv);
        // log(false, "g_md is %lf and g_mv is %d", g_md, g_mv);

        if (g_mv == -1) {
            panic("Exit because g_mv == 0");
        }

        bool update = (g_mv / my_n == my_rank);
        if (update) {
            int idx = g_mv % my_n;
            my_connected[idx] = 1;
        }
        // log(true, "I will %s my_connected", update ? "update" : "not update");

        // log for the myself connect node
        // {
        //     std::ostringstream oss;
        //     for (int i = 0; i < my_n; ++i) {
        //         oss << my_connected[i] << " ";
        //     }
        //     log(true, "connect: %s", oss.str().c_str());
        // }
        // log for the myself dist node
        // {
        //     std::ostringstream oss;
        //     for (int i = 0; i < my_n; ++i) {
        //         oss << my_dist_to[i] << " ";
        //     }
        //     log(true, "my_dist_to: %s", oss.str().c_str());
        // }

        for (int i = 0; i < my_n; ++i) {
            if (!my_connected[i]) {
                double new_dist = g_md + my_mat[g_mv * my_n + i];
                if (new_dist < my_dist_to[i]) {
                    my_dist_to[i] = new_dist;
                }
            }
        }

    }
    free(my_connected);
}

