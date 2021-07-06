#include "front.h"

#include <iostream>
#include <deal.II/base/mpi.h>


constexpr int p_dim = 2;

int main(int argc, char** argv)
{
    std::string paramsPath;
    if (argc < 2) {
        std::cout << "ERROR: no parameters file has been specified" << std::endl;
        return -1;
    } else {
        paramsPath = argv[1];
        std::cout << "NOTICE: parameters file provided: " << paramsPath << std::endl;
    }
    
    dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);
    std::cout << "    Number of cores       = " << dealii::MultithreadInfo::n_cores() << std::endl;
    std::cout << "    Number of threads     = " << dealii::MultithreadInfo::n_threads() << std::endl;
    std::cout << "    Is single thread mode = " << (dealii::MultithreadInfo::is_running_single_threaded() ? "yes" : "no") << std::endl;

    Front<p_dim> problem(paramsPath);

    // First call
    // returns a vector<double>
    problem.run();

    // Second call just for demonstration can be here
}
