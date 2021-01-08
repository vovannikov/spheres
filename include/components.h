#include <deal.II/base/point.h>
#include <deal.II/base/function.h>

#include "pde/material_object.h"

// Forward declarations
namespace PDE {
    template <int dim>
    class VariableValues;

    template <int dim>
    class Variable;
}

// Component 1
template <int dim>
class Di : public PDE::MaterialObject<double, dim>
{
private:
    std::shared_ptr<const PDE::Variable<dim>> _mVar;

    double _D_mi;
    double _K_mi_D;

    double _K_mi_D2;

public:
    Di(std::shared_ptr<const PDE::Variable<dim>> mVar, double D_mi, double K_mi_D);

    virtual double
    f(const PDE::VariableValues<dim>& variablesValues, unsigned int qp) const override;

    virtual double
    dfVal(const PDE::VariableValues<dim>& variablesValues, std::shared_ptr<const PDE::Variable<dim>> ivar, unsigned int qp) const override;

    virtual double
    dfGrad(const PDE::VariableValues<dim>& variablesValues, std::shared_ptr<const PDE::Variable<dim>> ivar, unsigned int qp) const override;
};

// Component 2
template <int dim>
class CM1 : public PDE::MaterialObject<double, dim>
{
private:
    std::shared_ptr<const PDE::Variable<dim>> _gTNFVar;
    std::shared_ptr<const PDE::Variable<dim>> _mVar;
    std::shared_ptr<const PDE::Variable<dim>> _cm1Var;

    const double _C_TNFM1;
    const double _K_TNFM1_C;
    const double _C_mM1;
    const double _K_mM1_C;

    const double _K_mM1_C2;
    const double _K_TNFM1_C2;

public:
    CM1(std::shared_ptr<const PDE::Variable<dim>> gTNFVar,
        std::shared_ptr<const PDE::Variable<dim>> mVar, 
        std::shared_ptr<const PDE::Variable<dim>> cm1Var,
        double C_TNFM1, double K_TNFM1_C, double C_mM1, double K_mM1_C);

    virtual double
    f(const PDE::VariableValues<dim>& variablesValues, unsigned int qp) const override;

    virtual double
    dfVal(const PDE::VariableValues<dim>& variablesValues, std::shared_ptr<const PDE::Variable<dim>> ivar, unsigned int qp) const override;

    virtual double
    dfGrad(const PDE::VariableValues<dim>& variablesValues, std::shared_ptr<const PDE::Variable<dim>> ivar, unsigned int qp) const override;
};

// Component 3
template <int dim>
class Ai : public PDE::MaterialObject<double, dim>
{
private:
    std::shared_ptr<const PDE::Variable<dim>> _mVar;
    std::shared_ptr<const PDE::Variable<dim>> _cm1Var;
    std::shared_ptr<const PDE::Variable<dim>> _cm2Var;
    std::shared_ptr<const PDE::Variable<dim>> _gjVar;

    const double _A_ji;
    const double _A_si;
    const double _K_ji_A;
    const double _K_si_A;

    const double _K_ji_A2;
    const double _K_si_A2;

public:
    Ai(std::shared_ptr<const PDE::Variable<dim>> mVar,
       std::shared_ptr<const PDE::Variable<dim>> cm1Var,
       std::shared_ptr<const PDE::Variable<dim>> cm2Var,
       std::shared_ptr<const PDE::Variable<dim>> gjVar,
       double A_ji, double A_si, double K_ji_A, double K_si_A);

    virtual double
    f(const PDE::VariableValues<dim>& variablesValues, unsigned int qp) const override;

    virtual double
    dfVal(const PDE::VariableValues<dim>& variablesValues, std::shared_ptr<const PDE::Variable<dim>> ivar, unsigned int qp) const override;

    virtual double
    dfGrad(const PDE::VariableValues<dim>& variablesValues, std::shared_ptr<const PDE::Variable<dim>> ivar, unsigned int qp) const override;
};

// Component 4
template <int dim>
class F5 : public PDE::MaterialObject<double, dim>
{
private:
    std::shared_ptr<const PDE::Variable<dim>> _gIL4Var;
    std::shared_ptr<const PDE::Variable<dim>> _gTGFVar;

    const double _I_IL4M1;
    const double _I_TGFM1;
    const double _K_IL4M1_I;
    const double _K_TGFM1_I;

    const double _K_IL4M1_I2;
    const double _K_TGFM1_I2;

public:
    F5(std::shared_ptr<const PDE::Variable<dim>> gIL4Var,
       std::shared_ptr<const PDE::Variable<dim>> gTGFVar,
       double I_IL4M1, double I_TGFM1, double K_IL4M1_I, double K_TGFM1_I);

    virtual double
    f(const PDE::VariableValues<dim>& variablesValues, unsigned int qp) const override;

    virtual double
    dfVal(const PDE::VariableValues<dim>& variablesValues, std::shared_ptr<const PDE::Variable<dim>> ivar, unsigned int qp) const override;

    virtual double
    dfGrad(const PDE::VariableValues<dim>& variablesValues, std::shared_ptr<const PDE::Variable<dim>> ivar, unsigned int qp) const override;
};

// Component 5
template <int dim>
class YM1 : public PDE::MaterialObject<double, dim>
{
private:
    std::shared_ptr<const PDE::Variable<dim>> _gTNFVar;

    const double _Y_M1_h;
    const double _Y_M1_L;

    const double _K_TNFM1_Y;
    const double _K_TNFM1_Y2;

public:
    YM1(std::shared_ptr<const PDE::Variable<dim>> gTNFVar,
       double Y_M1_h, double Y_M1_L, double K_TNFM1_Y);

    virtual double
    f(const PDE::VariableValues<dim>& variablesValues, unsigned int qp) const override;

    virtual double
    dfVal(const PDE::VariableValues<dim>& variablesValues, std::shared_ptr<const PDE::Variable<dim>> ivar, unsigned int qp) const override;

    virtual double
    dfGrad(const PDE::VariableValues<dim>& variablesValues, std::shared_ptr<const PDE::Variable<dim>> ivar, unsigned int qp) const override;
};

// Component 6
template <int dim>
class Dconst : public PDE::MaterialObject<double, dim>
{
private:
    const double _value;

public:
    Dconst(double value);

    virtual double
    f(const PDE::VariableValues<dim>& variablesValues, unsigned int qp) const override;

    virtual double
    dfVal(const PDE::VariableValues<dim>& variablesValues, std::shared_ptr<const PDE::Variable<dim>> ivar, unsigned int qp) const override;

    virtual double
    dfGrad(const PDE::VariableValues<dim>& variablesValues, std::shared_ptr<const PDE::Variable<dim>> ivar, unsigned int qp) const override;
};

// Component 7
template <int dim>
class ETNF : public PDE::MaterialObject<double, dim>
{
private:
    std::shared_ptr<const PDE::Variable<dim>> _mdVar;
    std::shared_ptr<const PDE::Variable<dim>> _gTNFVar;

    const double _E_TNF_h;
    const double _K_TNF;
    const double _K_TNF2;

public:
    ETNF(std::shared_ptr<const PDE::Variable<dim>> mdVar,
        std::shared_ptr<const PDE::Variable<dim>> gTNFVar,
        double E_TNF_h, double K_TNF);

    virtual double
    f(const PDE::VariableValues<dim>& variablesValues, unsigned int qp) const override;

    virtual double
    dfVal(const PDE::VariableValues<dim>& variablesValues, std::shared_ptr<const PDE::Variable<dim>> ivar, unsigned int qp) const override;

    virtual double
    dfGrad(const PDE::VariableValues<dim>& variablesValues, std::shared_ptr<const PDE::Variable<dim>> ivar, unsigned int qp) const override;
};

// Component 8
template <int dim>
class EIL4 : public PDE::MaterialObject<double, dim>
{
private:
    std::shared_ptr<const PDE::Variable<dim>> _gIL4Var;

    const double _E_IL4_h;
    const double _K_IL4;
    const double _K_IL42;

public:
    EIL4(std::shared_ptr<const PDE::Variable<dim>> gIL4Var,
        double E_IL4_h, double K_IL4);

    virtual double
    f(const PDE::VariableValues<dim>& variablesValues, unsigned int qp) const override;

    virtual double
    dfVal(const PDE::VariableValues<dim>& variablesValues, std::shared_ptr<const PDE::Variable<dim>> ivar, unsigned int qp) const override;

    virtual double
    dfGrad(const PDE::VariableValues<dim>& variablesValues, std::shared_ptr<const PDE::Variable<dim>> ivar, unsigned int qp) const override;
};

// Component 9
template <int dim>
class H : public PDE::MaterialObject<double, dim>
{
private:
    std::shared_ptr<const PDE::Variable<dim>> _cm1Var;
    std::shared_ptr<const PDE::Variable<dim>> _cm2Var;

    const double _dj;

public:
    H(std::shared_ptr<const PDE::Variable<dim>> cm1Var, std::shared_ptr<const PDE::Variable<dim>> cm2Var, double dj);

    virtual double
    f(const PDE::VariableValues<dim>& variablesValues, unsigned int qp) const override;

    virtual double
    dfVal(const PDE::VariableValues<dim>& variablesValues, std::shared_ptr<const PDE::Variable<dim>> ivar, unsigned int qp) const override;

    virtual double
    dfGrad(const PDE::VariableValues<dim>& variablesValues, std::shared_ptr<const PDE::Variable<dim>> ivar, unsigned int qp) const override;
};