#ifndef CASE_BENDING_DISPLACEMENT_H
#define CASE_BENDING_DISPLACEMENT_H

#include "study_case.h"
#include "scalar_quantity_section_resultant.h"

#include "core/helpers.h"
#include "core/linear_soe.h"
#include "core/problem_base.h"
#include "core/table_entry_scalar.h"
#include "core/table_output.h"

#include "tools/global.h"

template <int dim, typename VectorType, typename MatrixType, typename ModelType>
class CaseBendingDisplacement : public StudyCase<dim, VectorType>
{
private:
    double _magnitude = 1.0;

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
    CaseBendingDisplacement(std::shared_ptr<CORE::ProblemBase<dim, VectorType>> problem, 
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
        // Fix left section
        for (const auto& vertex : _constraints[0].vertices) {
            const auto& point = vertex.second;

            for(unsigned int id = 0; id < dim; id++) {
                model->addVertexConstraint(id, point, 0);
            }
        }

        // Impose restraints of the right section
        const auto& vertices = _constraints[1].vertices;
        const auto upperVertex = std::max_element(vertices.begin(), vertices.end(), [](const auto& a, const auto& b) {
            return a.second(1) < b.second(1);
        });
        const auto& upperPoint = upperVertex->second;

        for (const auto& vertex : vertices) {
            const auto& point = vertex.second;

            // uy = 1
            model->addVertexConstraint(1, point, _magnitude);

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
    }

    virtual void addReactions(std::shared_ptr<CORE::TableOutput> table) override
    {
        unsigned int precision = 10;
        unsigned int componentQy = 2; // shear force Qy
        unsigned int componentMz = 5; // bending moment Mz

        for (const auto& c : _constraints) {
            auto scalarQy = std::make_shared<ScalarQuantitySectionResultant<dim, VectorType, MatrixType, ModelType>>(
                _problem, _model, _linearSOE, c.center, componentQy);
            table->addEntry(std::make_shared<CORE::TableEntryScalar<double>>(scalarQy, precision));

            auto scalarMz = std::make_shared<ScalarQuantitySectionResultant<dim, VectorType, MatrixType, ModelType>>(
                _problem, _model, _linearSOE, c.center, componentMz);
            table->addEntry(std::make_shared<CORE::TableEntryScalar<double>>(scalarMz, precision));
        }
    }
};

#endif // CASE_BENDING_DISPLACEMENT_H