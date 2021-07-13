#ifndef CASE_TORSION_H
#define CASE_TORSION_H

#include "study_case.h"
#include "scalar_quantity_section_resultant.h"
#include "scalar_quantity_torsion.h"

#include "core/linear_soe.h"
#include "core/problem_base.h"
#include "core/table_entry_scalar.h"
#include "core/table_output.h"

template <int dim, typename VectorType, typename MatrixType, typename ModelType>
class CaseTorsion : public StudyCase<dim, VectorType>
{
private:
    double _magnitude = 0.01;

    std::shared_ptr<CORE::ProblemBase<dim, VectorType>> _problem;
    std::shared_ptr<const ModelType> _model;
    std::shared_ptr<CORE::LinearSOE<VectorType, MatrixType>> _linearSOE;

    struct Constraint {
        std::map<unsigned int, dealii::Point<dim>> vertices;
        dealii::Point<dim> values;
        dealii::Point<dim> center;
    };

    std::vector<Constraint> _constraints;

public:
    CaseTorsion(std::shared_ptr<CORE::ProblemBase<dim, VectorType>> problem, 
    std::shared_ptr<const ModelType> model,
        std::shared_ptr<CORE::LinearSOE<VectorType, MatrixType>> linearSOE,
        const dealii::Point<dim>& O1, const dealii::Point<dim>& O2)
        : _problem(problem)
        , _model(model)
        , _linearSOE(linearSOE)
    {
        dealii::Point<dim> values;

        _constraints.push_back(Constraint{{}, values, O1});

        values(0) = _magnitude;
        _constraints.push_back(Constraint{{}, values, O2});

        double posTol = 1e-6;

        for (const auto &cell : _problem->getApproximation().getDoFHandler().active_cell_iterators())
        //for (const auto &tria_cell : triangulation->active_cell_iterators())
        {
            if (cell->is_locally_owned())
            {
                /*
                auto cell = typename dealii::DoFHandler<dim>::cell_iterator(*tria_cell, 
                    &_problem->getApproximation().getDoFHandler());
                */

                for (unsigned int v = 0; v < cell->n_vertices(); ++v)
                {
                    auto& point = cell->vertex(v);

                    if (std::abs(point(0) - O1(0)) < posTol) {
                        _constraints[0].vertices[cell->vertex_index(v)] = point;
                    } else if (std::abs(point(0) - O2(0)) < posTol) {
                        _constraints[1].vertices[cell->vertex_index(v)] = point;
                    }
                }
            }
        }

/*
        for (auto ip : cLeft) {
            std::cout << ip->index() << std::endl;
        }
*/
        std::cout << "# of constraint sets: " << _constraints.size() << std::endl;
        std::cout << "# of left nodes:  " << _constraints[0].vertices.size() << std::endl;
        std::cout << "# of right nodes: " << _constraints[1].vertices.size() << std::endl;
    }

    virtual void imposeBoundaryConditions(std::shared_ptr<CORE::Model<dim, VectorType>> model) override
    {
        for (const auto& c : _constraints) {
            for (const auto& vertex : c.vertices) {
                const auto& point = vertex.second;

                auto radius = point - c.center;
                auto displ = dealii::cross_product_3d(c.values, radius);

                for(unsigned int id = 0; id < dim; id++) {
                    model->addVertexConstraint(id, point, displ[id]);
                }
            }
        }
    }

    virtual void addReactions(std::shared_ptr<CORE::TableOutput> table) override
    {
        unsigned int precision = 10;
        unsigned int componentNum = 3; // torsion Mx

        for (const auto& c : _constraints) {
            
            auto scalarTorsionMoment = std::make_shared<ScalarQuantitySectionResultant<dim, VectorType, MatrixType, ModelType>>(
                _problem, _model, _linearSOE, c.center, componentNum);
                
            table->addEntry(std::make_shared<CORE::TableEntryScalar<double>>(scalarTorsionMoment, precision));
        }
    }
};

#endif // CASE_TORSION_H