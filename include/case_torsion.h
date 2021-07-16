#ifndef CASE_TORSION_H
#define CASE_TORSION_H

#include "case_two_sections.h"

template <int dim, typename VectorType, typename MatrixType, typename ModelType>
class CaseTorsion : public CaseTwoSections<dim, VectorType, MatrixType, ModelType>
{
private:
    std::map<const Section<dim>*, dealii::Point<dim>> _values;

public:
    CaseTorsion(std::shared_ptr<CORE::ProblemBase<dim, VectorType>> problem, 
        std::shared_ptr<const ModelType> model,
        std::shared_ptr<CORE::LinearSOE<VectorType, MatrixType>> linearSOE,
        const dealii::Point<dim>& O1, const dealii::Point<dim>& O2,
        bool isRightActive)
        : CaseTwoSections<dim, VectorType, MatrixType, ModelType>(problem, model, linearSOE, O1, O2)
    {
        unsigned int activeId = isRightActive ? 1 : 0;
        unsigned int fixedId = isRightActive ? 0 : 1;

        dealii::Point<dim> values;

        _values[&this->getSections()[fixedId]] = values;

        values(0) = 0.01;
        _values[&this->getSections()[activeId]] = values;
    }

    virtual void imposeBoundaryConditions(std::shared_ptr<CORE::Model<dim, VectorType>> model) override
    {
        for (const auto& c : this->getSections()) {

            auto valuesAtSection = _values[&c];

            for (const auto& vertex : c.vertices) {

                const auto& point = vertex.second;

                auto radius = point - c.center;
                auto displ = dealii::cross_product_3d(valuesAtSection, radius);

                for(unsigned int id = 0; id < dim; id++) {
                    model->addVertexConstraint(id, point, displ[id]);
                }
            }
        }
    }

    virtual void addReactions(std::shared_ptr<CORE::TableOutput> table) override
    {
        unsigned int componentMx = 3; // torsion Mx

        for (const auto& c : this->getSections()) {
            table->addEntry(this->createTableEntry(c.center, componentMx));
        }
    }
};

#endif // CASE_TORSION_H