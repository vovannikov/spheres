#include "components.h"

#include "pde/common.h"
#include "pde/variable_values.h"

// Component 1
template <int dim>
Di<dim>::Di(std::shared_ptr<const PDE::Variable<dim>> mVar, double D_mi, double K_mi_D)
    : _mVar(mVar), _D_mi(D_mi), _K_mi_D(K_mi_D)
    , _K_mi_D2(_K_mi_D*_K_mi_D)
{}

template <int dim>
double Di<dim>::f(const PDE::VariableValues<dim>& variablesValues, unsigned int qp) const
{
    const double m = PDE::getVarValue(variablesValues, _mVar, qp);

    return _D_mi * m / (_K_mi_D2 + m*m);
}

template <int dim>
double Di<dim>::dfVal(const PDE::VariableValues<dim>& variablesValues, std::shared_ptr<const PDE::Variable<dim>> ivar,
    unsigned int qp) const
{
    if (ivar == _mVar) {
        const double m = PDE::getVarValue(variablesValues, _mVar, qp);

        return -2*_D_mi*m*m/std::pow(_K_mi_D2 + m*m, 2) + _D_mi/(_K_mi_D2 + m*m);
    } else {
        return 0;
    }
}

template <int dim>
double Di<dim>::dfGrad(const PDE::VariableValues<dim>& /*variablesValues*/, std::shared_ptr<const PDE::Variable<dim>> /*ivar*/,
    unsigned int /*qp*/) const
{
    return 0;
}


// Component 2 - changed sing after moving term to the left
template <int dim>
CM1<dim>::CM1(std::shared_ptr<const PDE::Variable<dim>> gTNFVar, 
         std::shared_ptr<const PDE::Variable<dim>> mVar, 
         std::shared_ptr<const PDE::Variable<dim>> cm1Var,
        double C_TNFM1, double K_TNFM1_C, double C_mM1, double K_mM1_C)
    : _gTNFVar(gTNFVar), _mVar(mVar), _cm1Var(cm1Var)
    , _C_TNFM1(C_TNFM1), _K_TNFM1_C(K_TNFM1_C), _C_mM1(C_mM1), _K_mM1_C(K_mM1_C)
    , _K_mM1_C2(_K_mM1_C*_K_mM1_C)
    , _K_TNFM1_C2(_K_TNFM1_C*_K_TNFM1_C)
{}

template <int dim>
double CM1<dim>::f(const PDE::VariableValues<dim>& variablesValues, unsigned int qp) const
{
    const double gTNF = PDE::getVarValue(variablesValues, _gTNFVar, qp);
    const double m = PDE::getVarValue(variablesValues, _mVar, qp);
    const double cm1 = PDE::getVarValue(variablesValues, _cm1Var, qp);

    const double gTNF2 = gTNF*gTNF;
    const double m2 = m*m;

    return - _C_TNFM1 * gTNF / (_K_TNFM1_C2 + gTNF2) * _C_mM1 * m / (_K_mM1_C2 + m2) * cm1;
}

template <int dim>
double CM1<dim>::dfVal(const PDE::VariableValues<dim>& variablesValues, std::shared_ptr<const PDE::Variable<dim>> ivar,
    unsigned int qp) const
{
    const double gTNF = PDE::getVarValue(variablesValues, _gTNFVar, qp);
    const double m = PDE::getVarValue(variablesValues, _mVar, qp);
    const double cm1 = PDE::getVarValue(variablesValues, _cm1Var, qp);

    const double gTNF2 = gTNF*gTNF;
    const double m2 = m*m;

    if (ivar == _gTNFVar) {
        return 2*_C_TNFM1*_C_mM1*cm1*gTNF2*m/(std::pow(_K_TNFM1_C2 + gTNF2, 2)*(_K_mM1_C2 + m2)) + _C_TNFM1*_C_mM1*cm1*m/((_K_TNFM1_C2 + gTNF2)*(_K_mM1_C2 + m2));
    } else if (ivar == _mVar) {
        return 2*_C_TNFM1*_C_mM1*cm1*gTNF*m2/((_K_TNFM1_C2 + gTNF2)*std::pow(_K_mM1_C2 + m2, 2)) + _C_TNFM1*_C_mM1*cm1*gTNF/((_K_TNFM1_C2 + gTNF2)*(_K_mM1_C2 + m2));
    } else if (ivar == _cm1Var) {
        return - _C_TNFM1*_C_mM1*gTNF*m/((_K_TNFM1_C2 + gTNF2)*(_K_mM1_C2 + m2));
    } else {
        return 0;
    }
}

template <int dim>
double CM1<dim>::dfGrad(const PDE::VariableValues<dim>& /*variablesValues*/, std::shared_ptr<const PDE::Variable<dim>> /*ivar*/,
    unsigned int /*qp*/) const
{
    return 0;
}



// Component 3 - changed sing after moving term to the left
template <int dim>
Ai<dim>::Ai(std::shared_ptr<const PDE::Variable<dim>> mVar,
       std::shared_ptr<const PDE::Variable<dim>> cm1Var,
       std::shared_ptr<const PDE::Variable<dim>> cm2Var,
       std::shared_ptr<const PDE::Variable<dim>> gjVar,
       double A_ji, double A_si, double K_ji_A, double K_si_A)
    : _mVar(mVar), _cm1Var(cm1Var), _cm2Var(cm2Var), _gjVar(gjVar)
    , _A_ji(A_ji), _A_si(A_si), _K_ji_A(K_ji_A), _K_si_A(K_si_A)
    , _K_ji_A2(_K_ji_A*_K_ji_A)
    , _K_si_A2(_K_si_A*_K_si_A)
{}

template <int dim>
double Ai<dim>::f(const PDE::VariableValues<dim>& variablesValues, unsigned int qp) const
{
    const double m = PDE::getVarValue(variablesValues, _mVar, qp);
    const double cm1 = PDE::getVarValue(variablesValues, _cm1Var, qp);
    const double cm2 = PDE::getVarValue(variablesValues, _cm2Var, qp);
    const double gj = PDE::getVarValue(variablesValues, _gjVar, qp);

    const double s = m + cm1 + cm2;

    const double gj2 = gj*gj;
    const double s2 = s*s;

    return - _A_ji * gj / (_K_ji_A2 + gj2) * _A_si * s / (_K_si_A2 + s2);
}

template <int dim>
double Ai<dim>::dfVal(const PDE::VariableValues<dim>& variablesValues, std::shared_ptr<const PDE::Variable<dim>> ivar,
    unsigned int qp) const
{
    const double m = PDE::getVarValue(variablesValues, _mVar, qp);
    const double cm1 = PDE::getVarValue(variablesValues, _cm1Var, qp);
    const double cm2 = PDE::getVarValue(variablesValues, _cm2Var, qp);
    const double gj = PDE::getVarValue(variablesValues, _gjVar, qp);

    const double s = m + cm1 + cm2;

    const double gj2 = gj*gj;
    const double s2 = s*s;

    if (ivar == _mVar || ivar == _cm1Var || ivar == _cm2Var) {
        return 2*_A_ji*_A_si*gj*s2/((_K_ji_A2 + gj2)*std::pow(_K_si_A2 + s2, 2)) + _A_ji*_A_si*gj/((_K_ji_A2 + gj2)*(_K_si_A2 + s2));
    } else if (ivar == _gjVar) {
        return 2*_A_ji*_A_si*gj2*s/(std::pow(_K_ji_A2 + gj2, 2)*(_K_si_A2 + s2)) + _A_ji*_A_si*s/((_K_ji_A2 + gj2)*(_K_si_A2 + s2));
    } else {
        return 0;
    }
}

template <int dim>
double Ai<dim>::dfGrad(const PDE::VariableValues<dim>& /*variablesValues*/, std::shared_ptr<const PDE::Variable<dim>> /*ivar*/,
    unsigned int /*qp*/) const
{
    return 0;
}




// Component 4
template <int dim>
F5<dim>::F5(std::shared_ptr<const PDE::Variable<dim>> gIL4Var,
       std::shared_ptr<const PDE::Variable<dim>> gTGFVar,
       double I_IL4M1, double I_TGFM1, double K_IL4M1_I, double K_TGFM1_I)
    : _gIL4Var(gIL4Var), _gTGFVar(gTGFVar)
    , _I_IL4M1(I_IL4M1), _I_TGFM1(I_TGFM1)
    , _K_IL4M1_I(K_IL4M1_I), _K_TGFM1_I(K_TGFM1_I)
    , _K_IL4M1_I2(_K_IL4M1_I*_K_IL4M1_I)
    , _K_TGFM1_I2(_K_TGFM1_I*_K_TGFM1_I)
{}

template <int dim>
double F5<dim>::f(const PDE::VariableValues<dim>& variablesValues, unsigned int qp) const
{
    const double gIL4 = PDE::getVarValue(variablesValues, _gIL4Var, qp);
    const double gTGF = PDE::getVarValue(variablesValues, _gTGFVar, qp);

    const double gIL42 = gIL4*gIL4;
    const double gTGF2 = gTGF*gTGF;

    return _I_IL4M1 * gIL42 / (_K_IL4M1_I2 + gIL42) * _I_TGFM1 * gTGF2 / (_K_TGFM1_I2 + gTGF2);
}

template <int dim>
double F5<dim>::dfVal(const PDE::VariableValues<dim>& variablesValues, std::shared_ptr<const PDE::Variable<dim>> ivar,
    unsigned int qp) const
{
    const double gIL4 = PDE::getVarValue(variablesValues, _gIL4Var, qp);
    const double gTGF = PDE::getVarValue(variablesValues, _gTGFVar, qp);

    const double gIL42 = gIL4*gIL4;
    const double gIL43 = gIL4*gIL4*gIL4;
    const double gTGF2 = gTGF*gTGF;
    const double gTGF3 = gTGF*gTGF*gTGF;

    if (ivar == _gIL4Var) {
        return -2*_I_IL4M1*_I_TGFM1*gIL43*gTGF2/(std::pow(_K_IL4M1_I2 + gIL42,2)*(_K_TGFM1_I2 + gTGF2)) + 2*_I_IL4M1*_I_TGFM1*gIL4*gTGF2/((_K_IL4M1_I2 + gIL42)*(_K_TGFM1_I2 + gTGF2));
    } else if (ivar == _gTGFVar) {
        return -2*_I_IL4M1*_I_TGFM1*gIL42*gTGF3/((_K_IL4M1_I2 + gIL42)*std::pow(_K_TGFM1_I2 + gTGF2,2)) + 2*_I_IL4M1*_I_TGFM1*gIL42*gTGF/((_K_IL4M1_I2 + gIL42)*(_K_TGFM1_I2 + gTGF2));
    } else {
        return 0;
    }
}

template <int dim>
double F5<dim>::dfGrad(const PDE::VariableValues<dim>& /*variablesValues*/, std::shared_ptr<const PDE::Variable<dim>> /*ivar*/,
    unsigned int /*qp*/) const
{
    return 0;
}


// Component 5
template <int dim>
YM1<dim>::YM1(std::shared_ptr<const PDE::Variable<dim>> gTNFVar,
       double Y_M1_h, double Y_M1_L, double K_TNFM1_Y)
    : _gTNFVar(gTNFVar)
    , _Y_M1_h(Y_M1_h), _Y_M1_L(Y_M1_L)
    , _K_TNFM1_Y(K_TNFM1_Y)
    , _K_TNFM1_Y2(_K_TNFM1_Y*_K_TNFM1_Y)
{}

template <int dim>
double YM1<dim>::f(const PDE::VariableValues<dim>& variablesValues, unsigned int qp) const
{
    const double gTNF = PDE::getVarValue(variablesValues, _gTNFVar, qp);

    const double gTNF2 = gTNF*gTNF;

    return _Y_M1_h - _Y_M1_L * gTNF / (_K_TNFM1_Y2 + gTNF2);
}

template <int dim>
double YM1<dim>::dfVal(const PDE::VariableValues<dim>& variablesValues, std::shared_ptr<const PDE::Variable<dim>> ivar,
    unsigned int qp) const
{
    const double gTNF = PDE::getVarValue(variablesValues, _gTNFVar, qp);

    const double gTNF2 = gTNF*gTNF;

    if (ivar == _gTNFVar) {
        return 2*_Y_M1_L*gTNF2/std::pow(_K_TNFM1_Y2 + gTNF2, 2) - _Y_M1_L/(_K_TNFM1_Y2 + gTNF2);
    } else {
        return 0;
    }
}

template <int dim>
double YM1<dim>::dfGrad(const PDE::VariableValues<dim>& /*variablesValues*/, std::shared_ptr<const PDE::Variable<dim>> /*ivar*/,
    unsigned int /*qp*/) const
{
    return 0;
}



// Component 6
template <int dim>
Dconst<dim>::Dconst(double value)
    : _value(value)
{}

template <int dim>
double Dconst<dim>::f(const PDE::VariableValues<dim>& /*variablesValues*/, unsigned int /*qp*/) const
{
    return _value;
}

template <int dim>
double Dconst<dim>::dfVal(const PDE::VariableValues<dim>& /*variablesValues*/, std::shared_ptr<const PDE::Variable<dim>> /*ivar*/,
    unsigned int /*qp*/) const
{
    return 0;
}

template <int dim>
double Dconst<dim>::dfGrad(const PDE::VariableValues<dim>& /*variablesValues*/, std::shared_ptr<const PDE::Variable<dim>> /*ivar*/,
    unsigned int /*qp*/) const
{
    return 0;
}



// Component 7 - changed sing after moving term to the left
template <int dim>
ETNF<dim>::ETNF(std::shared_ptr<const PDE::Variable<dim>> mdVar,
    std::shared_ptr<const PDE::Variable<dim>> gTNFVar,
    double E_TNF_h, double K_TNF)
    : _mdVar(mdVar), _gTNFVar(gTNFVar)
    , _E_TNF_h(E_TNF_h), _K_TNF(K_TNF)
    , _K_TNF2(_K_TNF*_K_TNF)
{}

template <int dim>
double ETNF<dim>::f(const PDE::VariableValues<dim>& variablesValues, unsigned int qp) const
{
    const double md = PDE::getVarValue(variablesValues, _mdVar, qp);
    const double gTNF = PDE::getVarValue(variablesValues, _gTNFVar, qp);

    const double gTNF2 = gTNF*gTNF;

    return - _E_TNF_h * md / (_K_TNF2 + gTNF2);
}

template <int dim>
double ETNF<dim>::dfVal(const PDE::VariableValues<dim>& variablesValues, std::shared_ptr<const PDE::Variable<dim>> ivar,
    unsigned int qp) const
{
    const double md = PDE::getVarValue(variablesValues, _mdVar, qp);
    const double gTNF = PDE::getVarValue(variablesValues, _gTNFVar, qp);

    const double gTNF2 = gTNF*gTNF;

    if (ivar == _mdVar) {
        return - _E_TNF_h/(_K_TNF2 + gTNF2);
    } else if (ivar == _gTNFVar) {
        return 2*_E_TNF_h*gTNF*md / std::pow(_K_TNF2 + gTNF2, 2);
    } else {
        return 0;
    }
}

template <int dim>
double ETNF<dim>::dfGrad(const PDE::VariableValues<dim>& /*variablesValues*/, std::shared_ptr<const PDE::Variable<dim>> /*ivar*/,
    unsigned int /*qp*/) const
{
    return 0;
}



// Component 8 - changed sing after moving term to the left
template <int dim>
EIL4<dim>::EIL4(std::shared_ptr<const PDE::Variable<dim>> gIL4Var,
        double E_IL4_h, double K_IL4)
    : _gIL4Var(gIL4Var)
    , _E_IL4_h(E_IL4_h), _K_IL4(K_IL4)
    , _K_IL42(_K_IL4*_K_IL4)
{}

template <int dim>
double EIL4<dim>::f(const PDE::VariableValues<dim>& variablesValues, unsigned int qp) const
{
    const double gIL4 = PDE::getVarValue(variablesValues, _gIL4Var, qp);

    const double gIL42 = gIL4*gIL4;

    return - _E_IL4_h / (_K_IL42 + gIL42);
}

template <int dim>
double EIL4<dim>::dfVal(const PDE::VariableValues<dim>& variablesValues, std::shared_ptr<const PDE::Variable<dim>> ivar,
    unsigned int qp) const
{
    const double gIL4 = PDE::getVarValue(variablesValues, _gIL4Var, qp);

    const double gIL42 = gIL4*gIL4;

    if (ivar == _gIL4Var) {
        return 2*_E_IL4_h*gIL4 / std::pow(_K_IL42 + gIL42, 2);
    } else {
        return 0;
    }
}

template <int dim>
double EIL4<dim>::dfGrad(const PDE::VariableValues<dim>& /*variablesValues*/, std::shared_ptr<const PDE::Variable<dim>> /*ivar*/,
    unsigned int /*qp*/) const
{
    return 0;
}


// Component 9
template <int dim>
H<dim>::H(std::shared_ptr<const PDE::Variable<dim>> cm1Var, std::shared_ptr<const PDE::Variable<dim>> cm2Var, double dj)
    : _cm1Var(cm1Var)
    , _cm2Var(cm2Var)
    , _dj(dj)
{}

template <int dim>
double H<dim>::f(const PDE::VariableValues<dim>& variablesValues, unsigned int qp) const
{
    const double cm1 = PDE::getVarValue(variablesValues, _cm1Var, qp);
    const double cm2 = PDE::getVarValue(variablesValues, _cm2Var, qp);

    return _dj * (1 + cm1 + cm2);
}

template <int dim>
double H<dim>::dfVal(const PDE::VariableValues<dim>& /*variablesValues*/, std::shared_ptr<const PDE::Variable<dim>> ivar,
    unsigned int /*qp*/) const
{

    if (ivar == _cm1Var || ivar == _cm2Var) {
        return _dj;
    } else {
        return 0;
    }
}

template <int dim>
double H<dim>::dfGrad(const PDE::VariableValues<dim>& /*variablesValues*/, std::shared_ptr<const PDE::Variable<dim>> /*ivar*/,
    unsigned int /*qp*/) const
{
    return 0;
}


template class Di<2>;
template class CM1<2>;
template class Ai<2>;
template class F5<2>;
template class YM1<2>;
template class Dconst<2>;
template class ETNF<2>;
template class EIL4<2>;
template class H<2>;