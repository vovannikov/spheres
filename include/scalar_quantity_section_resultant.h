#ifndef SCALAR_QUANTITY_SECTION_RESULTANT_H
#define SCALAR_QUANTITY_SECTION_RESULTANT_H

#include "core/linear_soe.h"
#include "core/problem_base.h"
#include "core/scalar_quantity.h"

#include <deal.II/base/point.h>

template <int dim, typename VectorType, typename MatrixType, typename ModelType>
class ScalarQuantitySectionResultant : public CORE::ScalarQuantity<double>
{
private:
    std::shared_ptr<CORE::ProblemBase<dim, VectorType>> _problem;
    std::shared_ptr<const ModelType> _model;
    mutable std::shared_ptr<CORE::LinearSOE<VectorType, MatrixType>> _linearSOE;
    dealii::Point<dim> _center;
    unsigned int _component_num;

    struct VertexData {
        dealii::Point<dim> point;
        std::vector<unsigned int> dofs;
    };

    mutable std::map<unsigned int, VertexData> _vertices;

    std::string _name;

public:
    ScalarQuantitySectionResultant(
        std::shared_ptr<CORE::ProblemBase<dim, VectorType>> problem,
        std::shared_ptr<const ModelType> model,
        std::shared_ptr<CORE::LinearSOE<VectorType, MatrixType>> linearSOE,
        dealii::Point<dim> center,
        unsigned int component_num)
        : _problem(problem)
        , _model(model)
        , _linearSOE(linearSOE)
        , _center(center)
        , _component_num(component_num)
    {
        std::stringstream sstr;
        sstr << "sec_force_" <<  _component_num << "_pt_" << _center;
        _name = sstr.str();

        std::transform(_name.begin(), _name.end(), _name.begin(), [](char ch) {
            return ch == ' ' ? '_' : ch;
        });
    }

    virtual double evaluate() const
    {
        buildVertices();

        computeForces();

        std::vector<double> forcesAndMoments(2*dim);

        // Now we have to extract the necessary dofs and postprocess them
        for (const auto& vdata : _vertices) {

            dealii::Tensor<1, dim> forces;

            const auto& vertex = vdata.second;

            for (unsigned int i=0; i<dim; i++) {
                forces[i] = _linearSOE->getRHS()[vertex.dofs[i]];
            }

            const auto& point = vertex.point;

            auto radius = point - _center;
            auto moments = dealii::cross_product_3d(radius, forces);

            // Unroll tensors to vectors
            dealii::Vector<double> vecForces(dim);
            dealii::Vector<double> vecMoments(dim);
            forces.unroll(vecForces);
            moments.unroll(vecMoments);

            std::vector<double> localForcesAndMoments(2*dim);
            for (unsigned int i=0; i<dim; i++) {
                localForcesAndMoments[i] = vecForces[i];
                localForcesAndMoments[dim+i] = vecMoments[i];
            }

            std::transform(forcesAndMoments.begin(), forcesAndMoments.end(), 
                localForcesAndMoments.begin(), forcesAndMoments.begin(), std::plus<double>());
        }

        return forcesAndMoments[_component_num];
    }

    virtual std::string name() const
    {
        return _name;
    }

private:

    void computeForces() const
    {
        const auto& approximation = _problem->getApproximation();

        // Compute RHS firstly

        // Compute velocities
        //_problem->getState().computeVelocity(_problem->getModel()->getCurrentTimeStep());

        // Handle ghost values
        _problem->getState().getSolution().update_ghost_values();

        // Nullify RHS
        _linearSOE->nullifyRHS();

        auto scratch_data = _model->buildScratchData(_problem);

        // allocate memory for element matrix, vector and index vector
        const unsigned int dofs_per_cell = approximation.getFESystem().dofs_per_cell;

        dealii::FullMatrix<double> cell_matrix;
        dealii::Vector<double> cell_rhs;
        cell_rhs.reinit(dofs_per_cell);

        //for (const auto &cell : _problem->getApproximation().getDoFHandler().active_cell_iterators()) {
        for (const auto &tria_cell : approximation.getDoFHandler().get_triangulation().active_cell_iterators()) {

            if (tria_cell->is_locally_owned()) {

                auto cell = typename dealii::DoFHandler<dim>::cell_iterator(*tria_cell, &approximation.getDoFHandler());

                // reset element vector
                cell_rhs = 0;

                _model->computeCell(tria_cell, cell_rhs, cell_matrix, scratch_data);

                // Work on boundaries
                for (unsigned int face = 0; face < tria_cell->n_faces(); ++face) {
                    if (tria_cell->at_boundary(face)) {
                        _model->computeBoundary(tria_cell, face, cell_rhs, cell_matrix, scratch_data);
                    }
                }

                cell->distribute_local_to_global(cell_rhs, _linearSOE->getRHS());
            }
        }

        _linearSOE->getRHS().compress(dealii::VectorOperation::add);
    }

    void buildVertices() const
    {
        _vertices.clear();

        double posTol = 1e-6;

        // Section is alogn x-axis - hardcoded so far
        unsigned int sectionOrientation = 0;

        for (const auto &cell : _problem->getApproximation().getDoFHandler().active_cell_iterators())
        {
            if (cell->is_locally_owned())
            {
                for (unsigned int v = 0; v < cell->n_vertices(); ++v)
                {
                    auto& point = cell->vertex(v);

                    
                    if (std::abs(point(sectionOrientation) - _center(sectionOrientation)) < posTol) {
                        std::vector<unsigned int> dofs;

                        for (unsigned int i = 0; i < dim; i++)  {
                            dofs.push_back(cell->vertex_dof_index(v, i));
                        }

                        VertexData vdata{point, dofs};

                        _vertices[cell->vertex_index(v)] = vdata;
                    }
                }
            }
        }
    }

};

#endif /* SCALAR_QUANTITY_SECTION_RESULTANT_H */
