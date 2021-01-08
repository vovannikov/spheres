#include "bone.h"
#include "bone_front.h"
#include "parameter_reader.h"
#include "tools/global.h"

#include <deal.II/base/parameter_handler.h>

// Main model itself
template <int dim>
BoneFront<dim>::BoneFront(const std::string& paramsPath)
    : _bone(nullptr)
{
    std::cout << "Building Bone frontend ..." << std::endl;
    TOOLS::initGlobals();

    dealii::ParameterHandler prm;
    ParameterReader reader(prm);
    reader.read_parameters(paramsPath);

    _bone = std::make_unique<Bone<dim>>(prm);
}

// Main model itself
template <int dim>
BoneFront<dim>::~BoneFront()
{
    TOOLS::timer::instance().print_summary();
}

template <int dim>
void BoneFront<dim>::run()
{
    _bone->run();
}

template class BoneFront<2>;
//template class BoneFront<3>;
