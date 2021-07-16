#ifndef CASE_TENSION_H
#define CASE_TENSION_H

#include "case_two_sections.h"

template <int dim, typename VectorType, typename MatrixType, typename ModelType>
class CaseTension : public CaseTwoSections<dim, VectorType, MatrixType, ModelType>
{
private:
    std::map<const Section<dim>*, std::vector<double>> _values;

public:
    CaseTension(std::shared_ptr<CORE::ProblemBase<dim, VectorType>> problem, 
        std::shared_ptr<const ModelType> model,
        std::shared_ptr<CORE::LinearSOE<VectorType, MatrixType>> linearSOE,
        const dealii::Point<dim>& O1, const dealii::Point<dim>& O2)
        : CaseTwoSections<dim, VectorType, MatrixType, ModelType>(problem, model, linearSOE, O1, O2)
    {
        _values[&this->getSections()[0]] = std::vector<double>({0, 0, 0});
        _values[&this->getSections()[1]] = std::vector<double>({1, 0, 0});
    }

    virtual void imposeBoundaryConditions(std::shared_ptr<CORE::Model<dim, VectorType>> model) override
    {
        for (const auto& c : this->getSections()) {

            auto valuesAtSection = _values[&c];

            for (const auto& vertex : c.vertices) {
                const auto& point = vertex.second;

                for(unsigned int id = 0; id < dim; id++) {
                    model->addVertexConstraint(id, point, valuesAtSection[id]);
                }
            }
        }
    }

    virtual void addReactions(std::shared_ptr<CORE::TableOutput> table) override
    {
        unsigned int componentFx = 0; // tension Fx

        for (const auto& c : this->getSections()) {
            table->addEntry(this->createTableEntry(c.center, componentFx));
        }
    }
};

#endif // CASE_TENSION_H