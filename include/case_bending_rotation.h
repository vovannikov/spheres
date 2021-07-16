#ifndef CASE_BENDING_ROTATION_H
#define CASE_BENDING_ROTATION_H

#include "case_two_sections.h"

template <int dim, typename VectorType, typename MatrixType, typename ModelType>
class CaseBendingRotation : public CaseTwoSections<dim, VectorType, MatrixType, ModelType>
{
private:
    double _magnitude = 0.01;

    unsigned int _activeId;
    unsigned int _fixedId;

public:
    CaseBendingRotation(std::shared_ptr<CORE::ProblemBase<dim, VectorType>> problem, 
        std::shared_ptr<const ModelType> model,
        std::shared_ptr<CORE::LinearSOE<VectorType, MatrixType>> linearSOE,
        const dealii::Point<dim>& O1, const dealii::Point<dim>& O2,
        bool isRightActive)
        : CaseTwoSections<dim, VectorType, MatrixType, ModelType>(problem, model, linearSOE, O1, O2)
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
        for (const auto& vertex : this->getSections()[_activeId].vertices) {
            const auto& point = vertex.second;

            double y = point(1);

            // ux = theta*y
            model->addVertexConstraint(0, point, _magnitude * y);

            // uy = 0
            model->addVertexConstraint(1, point, 0);
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

#endif // CASE_BENDING_ROTATION_H