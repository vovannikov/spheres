#ifndef CASE_TENSION_POINTWISE_H
#define CASE_TENSION_POINTWISE_H

#include "study_case.h"

#include "core/linear_soe.h"
#include "core/problem_base.h"
#include "core/scalar_quantity_reaction.h"
#include "core/table_entry_scalar.h"
#include "core/table_output.h"

template <int dim, typename VectorType, typename MatrixType, typename ModelType>
class CaseTensionPointwise : public StudyCase<dim, VectorType>
{
private:
    double _magnitude = 1.0;

    std::shared_ptr<CORE::ProblemBase<dim, VectorType>> _problem;
    std::shared_ptr<const ModelType> _model;
    std::shared_ptr<CORE::LinearSOE<VectorType, MatrixType>> _linearSOE;

    struct Constraint {
        dealii::Point<dim> point;
        std::vector<double> values;
    };

    std::vector<Constraint> _constraints;

public:
    CaseTensionPointwise(std::shared_ptr<CORE::ProblemBase<dim, VectorType>> problem, 
    std::shared_ptr<const ModelType> model,
        std::shared_ptr<CORE::LinearSOE<VectorType, MatrixType>> linearSOE,
        const dealii::Point<dim>& O1, const dealii::Point<dim>& O2)
        : _problem(problem)
        , _model(model)
        , _linearSOE(linearSOE)
    {
        std::vector<double> values(dim);

        // Fix left sphere
        _constraints.push_back(Constraint{O1, values});

        // Define unit tension at the right one
        values[0] = _magnitude;
        _constraints.push_back(Constraint{O2, values});
    }

    virtual void imposeBoundaryConditions(std::shared_ptr<CORE::Model<dim, VectorType>> model) override
    {
        for (const auto& c : _constraints) {
            for(unsigned int id = 0; id < c.values.size(); id++) {
                model->addVertexConstraint(id, c.point, c.values[id]);
            }
        }
    }

    virtual void addReactions(std::shared_ptr<CORE::TableOutput> table) override
    {
        dealii::types::global_dof_index component_num = 0;

        unsigned int precision = 10;

        for (const auto& c : _constraints) {
            auto scalarReaction = std::make_shared<CORE::ScalarQuantityReaction<dim, VectorType, MatrixType, ModelType>>(
                _problem, _model, _linearSOE, component_num, c.point);
            table->addEntry(std::make_shared<CORE::TableEntryScalar<double>>(scalarReaction, precision));
        }
    }

    virtual double loadMagnitude() const override
    {
        return _magnitude;
    }
};

#endif // CASE_TENSION_POINTWISE_H