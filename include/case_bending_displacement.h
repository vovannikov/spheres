#ifndef CASE_BENDING_DISPLACEMENT_H
#define CASE_BENDING_DISPLACEMENT_H

#include "case_two_sections.h"

template <int dim, typename VectorType, typename MatrixType, typename ModelType>
class CaseBendingDisplacement : public CaseTwoSections<dim, VectorType, MatrixType, ModelType>
{
private:
    double _magnitude;

    unsigned int _activeId;
    unsigned int _fixedId;

public:
    CaseBendingDisplacement(std::shared_ptr<CORE::ProblemBase<dim, VectorType>> problem, 
        std::shared_ptr<const ModelType> model,
        std::shared_ptr<CORE::LinearSOE<VectorType, MatrixType>> linearSOE,
        const dealii::Point<dim>& O1, const dealii::Point<dim>& O2,
        bool isRightActive)
        : CaseTwoSections<dim, VectorType, MatrixType, ModelType>(problem, model, linearSOE, O1, O2)
        //, _magnitude(isRightActive ? -1.0 : 1.0)
        , _magnitude(-1.0)
        , _activeId(isRightActive ? 1 : 0)
        , _fixedId(isRightActive ? 0 : 1)
    {
    }

    virtual void imposeBoundaryConditions(std::shared_ptr<CORE::Model<dim, VectorType>> model) override
    {
        // Fix left section
        for (const auto& vertex : this->getSections()[_fixedId].vertices) {
            const auto& point = vertex.second;

            for(unsigned int id = 0; id < dim; id++) {
                model->addVertexConstraint(id, point, 0);
            }
        }

        // Impose restraints of the right section
        /*
        const auto upperVertex = std::max_element(
            this->getSections()[1].vertices.begin(), 
            this->getSections()[1].vertices.end(), 
            [](const auto& a, const auto& b) {
                return a.second(1) < b.second(1);
            }
        );
        const auto& upperPoint = upperVertex->second;
        */

        for (const auto& vertex : this->getSections()[_activeId].vertices) {
            const auto& point = vertex.second;

            // ux = 0
            model->addVertexConstraint(0, point, 0);

            // uy = 1
            model->addVertexConstraint(1, point, _magnitude);
/*
            double distTol = 1e-6;
            if (point.distance(upperPoint) > distTol) {
                model->addGeneralizedConstraint([upperPoint, point] (dealii::DoFHandler<dim>& dofHandler, dealii::AffineConstraints<double>& constraints) {

                    double yMax = upperPoint(1);
                    double y = point(1);

                    // ux proportional to radius
                    double ratio = y / yMax;

                    const unsigned int componentNum = 0;

                    const int globalDofIndexUpper = CORE::findGlobalDofIndex(dofHandler, upperPoint, componentNum);
                    const int globalDofIndex = CORE::findGlobalDofIndex(dofHandler, point, componentNum);

                    if (globalDofIndex >= 0) {
                        constraints.add_line(globalDofIndex);
                        constraints.add_entry(globalDofIndex, globalDofIndexUpper, ratio);
                    } else {
                        TOOLS::cout::instance() << "WARNING: Vertex dependency constraint for point " << point 
                            << " was not imposed since the point does not seem to be a vertex. " << std::endl;
                    }
                });
            }
*/
        }
    }

    virtual void addReactions(std::shared_ptr<CORE::TableOutput> table) override
    {
        unsigned int componentQy = 1; // shear force Qy
        unsigned int componentMz = 5; // bending moment Mz

        for (const auto& c : this->getSections()) {
            table->addEntry(this->createTableEntry(c.center, componentQy));
            table->addEntry(this->createTableEntry(c.center, componentMz));
        }
    }

    virtual double loadMagnitude() const override
    {
        return _magnitude;
    }
};

#endif // CASE_BENDING_DISPLACEMENT_H