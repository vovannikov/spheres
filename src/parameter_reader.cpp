
#include "parameter_reader.h"

ParameterReader::ParameterReader(dealii::ParameterHandler &paramhandler)
    : prm(paramhandler)
{}

void ParameterReader::declare_parameters()
{
    prm.enter_subsection("Geometry");
    {
        prm.declare_entry("mesh", "", dealii::Patterns::Anything(), "Mesh file");
        prm.declare_entry("x1", "0.0", dealii::Patterns::Double(0), "Sphere 1");
        prm.declare_entry("x2", "0.0", dealii::Patterns::Double(0), "Sphere 2");
    }
    prm.leave_subsection();
    
    prm.enter_subsection("Physical constants");
    {
        prm.declare_entry("E", "0", dealii::Patterns::Double(0), "Young's modulus");
        prm.declare_entry("nu", "0", dealii::Patterns::Double(0), "Poisson ratio");
        prm.declare_entry("tensor", "true", dealii::Patterns::Bool(), "Use tensor material");
    }
    prm.leave_subsection();

    prm.enter_subsection("Output");
    {
        prm.declare_entry("path_vtk", ".", dealii::Patterns::Anything(), "Save path for vtk");
        prm.declare_entry("path_log", ".", dealii::Patterns::Anything(), "Save path for log");
        prm.declare_entry("save_vtk", "true", dealii::Patterns::Bool(), "Save vtk output");
        prm.declare_entry("save_log", "true", dealii::Patterns::Bool(), "Save log output");
    }
    prm.leave_subsection();
}

void ParameterReader::read_parameters(const std::string &parameter_file)
{
    declare_parameters();
    prm.parse_input(parameter_file);
}