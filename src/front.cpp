#include "spheres.h"
#include "front.h"
#include "parameter_reader.h"
#include "tools/global.h"

#include <deal.II/base/parameter_handler.h>

// Main model itself
template <int dim>
Front<dim>::Front(const std::string& paramsPath)
    : _spheres(nullptr)
{
    std::cout << "Building Spheres frontend ..." << std::endl;
    TOOLS::initGlobals();

    dealii::ParameterHandler prm;
    ParameterReader reader(prm);
    reader.read_parameters(paramsPath);

    _spheres = std::make_unique<Spheres<dim>>(prm);
}

// Main model itself
template <int dim>
Front<dim>::~Front()
{
    TOOLS::timer::instance().print_summary();
}

template <int dim>
void Front<dim>::run()
{
    _spheres->run();
}

template class Front<2>;
template class Front<3>;
