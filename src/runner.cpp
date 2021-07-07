#include "front.h"

#include <iostream>
#include <deal.II/base/mpi.h>

int main(int argc, char** argv)
{
    if (argc < 3) {
        std::cout << "ERROR: please specify required arguments [params_path] [geometry]" << std::endl;
        std::cout << "    [params_path] - parameters path" << std::endl;
        std::cout << "    [geometry] - which geometry to analyze:" << std::endl;
        std::cout << "        2d - 2D case" << std::endl;
        std::cout << "        3d - 3D case" << std::endl;
        return -1;
    }

    // Params path
    std::string paramsPath(argv[1]);

    // Which case to analyze: simple or not
    std::string geometry(argv[2]);
    bool is3D = geometry == "3d";
    
    dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);
    std::cout << "    Number of cores       = " << dealii::MultithreadInfo::n_cores() << std::endl;
    std::cout << "    Number of threads     = " << dealii::MultithreadInfo::n_threads() << std::endl;
    std::cout << "    Is single thread mode = " << (dealii::MultithreadInfo::is_running_single_threaded() ? "yes" : "no") << std::endl;

    if (is3D) {
        Front<3> problem(paramsPath);
        problem.run();
    } else {
        Front<2> problem(paramsPath);
        problem.run();
    }
}
