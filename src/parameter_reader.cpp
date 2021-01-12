
#include "parameter_reader.h"

ParameterReader::ParameterReader(dealii::ParameterHandler &paramhandler)
    : prm(paramhandler)
{}

void ParameterReader::declare_parameters()
{
    prm.enter_subsection("Geometry");
    {
        prm.declare_entry("mesh", "", dealii::Patterns::Anything(), "Mesh file");
    }
    prm.leave_subsection();
    
    prm.enter_subsection("Physical constants");
    {
        prm.declare_entry("D_mM1", "0", dealii::Patterns::Double(0), "DM1::D_mM1");
        prm.declare_entry("K_mM1_D", "0", dealii::Patterns::Double(0), "DM1::K_mM1_D");

        prm.declare_entry("D_mM2", "0", dealii::Patterns::Double(0), "DM2::D_mM2");
        prm.declare_entry("K_mM2_D", "0", dealii::Patterns::Double(0), "DM2::K_mM2_D");

        prm.declare_entry("C_TNFM1", "0", dealii::Patterns::Double(0), "CM1::C_TNFM1");
        prm.declare_entry("K_TNFM1_C", "0", dealii::Patterns::Double(0), "CM1::K_TNFM1_C");
        prm.declare_entry("C_mM1", "0", dealii::Patterns::Double(0), "CM1::C_mM1");
        prm.declare_entry("K_mM1_C", "0", dealii::Patterns::Double(0), "CM1::K_mM1_C");

        prm.declare_entry("A_TNFM1", "0", dealii::Patterns::Double(0), "AM1::A_TNFM1");
        prm.declare_entry("A_sM1", "0", dealii::Patterns::Double(0), "AM1::A_sM1");
        prm.declare_entry("K_TNFM1_A", "0", dealii::Patterns::Double(0), "AM1::K_TNFM1_A");
        prm.declare_entry("K_sM1_A", "0", dealii::Patterns::Double(0), "AM1::K_sM1_A");

        prm.declare_entry("A_IL4M2", "0", dealii::Patterns::Double(0), "AM2::A_IL4M2");
        prm.declare_entry("A_sM2", "0", dealii::Patterns::Double(0), "AM2::A_sM2");
        prm.declare_entry("K_IL4M2_A", "0", dealii::Patterns::Double(0), "AM2::K_IL4M2_A");
        prm.declare_entry("K_sM2_A", "0", dealii::Patterns::Double(0), "AM2::K_sM2_A");

        prm.declare_entry("I_IL4M1", "0", dealii::Patterns::Double(0), "F5::I_IL4M1");
        prm.declare_entry("I_TGFM1", "0", dealii::Patterns::Double(0), "F5::I_TGFM1");
        prm.declare_entry("K_IL4M1_I", "0", dealii::Patterns::Double(0), "F5::K_IL4M1_I");
        prm.declare_entry("K_TGFM1_I", "0", dealii::Patterns::Double(0), "F5::K_TGFM1_I");

        prm.declare_entry("Y_M1_h", "0", dealii::Patterns::Double(0), "YM1::Y_M1_h");
        prm.declare_entry("Y_M1_L", "0", dealii::Patterns::Double(0), "YM1::Y_M1_L");
        prm.declare_entry("K_TNFM1_Y", "0", dealii::Patterns::Double(0), "YM1::K_TNFM1_Y");

        prm.declare_entry("Y_M2_h", "0", dealii::Patterns::Double(0), "YM2::Y_M2_h");

        prm.declare_entry("E_TNF_h", "0", dealii::Patterns::Double(0), "ETNF::E_TNF_h");
        prm.declare_entry("K_TNF", "0", dealii::Patterns::Double(0), "ETNF::K_TNF");

        prm.declare_entry("E_IL4_h", "0", dealii::Patterns::Double(0), "EIL4::E_IL4_h");
        prm.declare_entry("K_IL4", "0", dealii::Patterns::Double(0), "EIL4::K_IL4");

        prm.declare_entry("D_TNF_h", "0", dealii::Patterns::Double(0), "DTNF::D_TNF_h");

        prm.declare_entry("D_IL4_h", "0", dealii::Patterns::Double(0), "DIL4::D_IL4_h");

        prm.declare_entry("d_TNF", "0", dealii::Patterns::Double(0), "HTNF::d_TNF");

        prm.declare_entry("d_IL4", "0", dealii::Patterns::Double(0), "HIL4::d_IL4");

        prm.declare_entry("D_d", "0", dealii::Patterns::Double(0), "Dd::D_d");
        prm.declare_entry("d_d", "0", dealii::Patterns::Double(0), "Dd::d_d");
    }
    prm.leave_subsection();

    prm.enter_subsection("Initial conditions");
    {
        prm.declare_entry("cm1", "0.1", dealii::Patterns::Double(0), "cm1");
        prm.declare_entry("cm2", "0.0", dealii::Patterns::Double(0), "cm2");
        prm.declare_entry("gTNF", "0.1", dealii::Patterns::Double(0), "gTNF");
        prm.declare_entry("gIL4", "0.0", dealii::Patterns::Double(0), "gIL4");
        prm.declare_entry("md", "1.0", dealii::Patterns::Double(0), "md");
        prm.declare_entry("mf", "0.1", dealii::Patterns::Double(0), "mf");
    }
    prm.leave_subsection();

    prm.enter_subsection("Time integration");
    {
        prm.declare_entry("time_start", "0", dealii::Patterns::Double(0), "Start time, h");
        prm.declare_entry("time_end", "1", dealii::Patterns::Double(0), "End time, h");
        prm.declare_entry("time_step_init", "1e-2", dealii::Patterns::Double(0), "Time step initial, h");
        prm.declare_entry("time_step_min", "1e-3", dealii::Patterns::Double(0), "Time step min, h");
        prm.declare_entry("time_step_max", "1e-1", dealii::Patterns::Double(0), "Time step max, h");
    }
    prm.leave_subsection();

    prm.enter_subsection("Output");
    {
        prm.declare_entry("save_path", ".", dealii::Patterns::Anything(), "Save path");
    }
    prm.leave_subsection();
}

void ParameterReader::read_parameters(const std::string &parameter_file)
{
    declare_parameters();
    prm.parse_input(parameter_file);
}