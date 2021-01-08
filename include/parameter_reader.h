#include <deal.II/base/subscriptor.h>
#include <deal.II/base/parameter_handler.h>

class ParameterReader : public dealii::Subscriptor
{
public:
    ParameterReader(dealii::ParameterHandler &);
    void read_parameters(const std::string &);

private:
    void declare_parameters();
    dealii::ParameterHandler &prm;
};