#include <vector>
#include <string>
#include <memory>

// This is a pimpl wrapper

// Forward declaration of the Bone model
template <int dim>
class Bone;

template <int dim>
class BoneFront
{
private:
    std::unique_ptr<Bone<dim>> _bone;

public:
    BoneFront(const std::string& paramsPath);

    ~BoneFront();

    void run();
};
