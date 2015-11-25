#include <stdlib.h>
#include <mpi.h>
#include <utility>
#include <math.h>
#include <iostream>



int main(int argc, char[] argv[]) {
    long long int total;
    int commSz;
    int myRank;
    long long int numIterations;
    long long int numHits = 0;
    long long int totalHits = 0;
    MPI_Init(NULL, NULL);

    MPI_Comm_size(MPI_COMM_WORLD, &commSz); 
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);


    // First tell each worker how many iterations they need to do.
    if (myRank == 0) {
        std::cout << "Please enter the number of iterations you want to run per process" << std::endl;
        std::cin >> numIterations;
        for (int i = 1; i < commSz; i++) 
            MPI_Send(&numIterations, 1, MPI_LONG_LONG_INT, i, 0, MPI_COMM_WORLD);

    } 
    else {
        MPI_Recv(&numIterations, 1, MPI_LONG_LONG_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (int i = 0; i < numIterations; i++) {
            if (distSquared(throwDart()) <= 1.0) 
                numHits += 1;
        }
    }
    MPI_Reduce(&numHits, &totalHits, 1, MPI_LONG_LONG_INT, MPI_Sum, 0, MPI_COMM_WORLD);
    MPI_Finalize();
    std::cout << "Our predicted value for PI is the following: " << 4.0 * ((double) totalHits / \
                                                                           (double) (numIterations * (commSz - 1))) << std::endl;
}


std::pair<double, double> throwDart() {
    double x = (((double) rand() / RAND_MAX) * 2) - 1;
    double y = (((double) rand() / RAND_MAX) * 2) - 1;
    return std::make_pair(x,y);
}
double distSquared(std::pair<double, double> coords) {
    return coords.first * coords.first + coords.second * coords.second;
}
