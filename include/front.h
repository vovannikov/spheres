#include <vector>
#include <string>
#include <memory>

// This is a pimpl wrapper

// Forward declaration of the Spheres model
template <int dim>
class Spheres;

template <int dim>
class Front
{
private:
    std::unique_ptr<Spheres<dim>> _spheres;

public:
    Front(const std::string& paramsPath);

    ~Front();

    void run();
};
