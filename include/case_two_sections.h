#ifndef CASE_TWO_SECTIONS_H
#define CASE_TWO_SECTIONS_H

#include "study_case.h"
#include "scalar_quantity_section_resultant.h"

#include "core/linear_soe.h"
#include "core/problem_base.h"
#include "core/table_entry_scalar.h"
#include "core/table_output.h"

template <int dim>
struct Section {
    std::map<unsigned int, dealii::Point<dim>> vertices;
    dealii::Point<dim> center;
};

template <int dim, typename VectorType, typename MatrixType, typename ModelType>
class CaseTwoSections : public StudyCase<dim, VectorType>
{
private:
    std::shared_ptr<CORE::ProblemBase<dim, VectorType>> _problem;
    std::shared_ptr<const ModelType> _model;
    std::shared_ptr<CORE::LinearSOE<VectorType, MatrixType>> _linearSOE;

    std::vector<Section<dim>> _sections;

public:
    CaseTwoSections(std::shared_ptr<CORE::ProblemBase<dim, VectorType>> problem, 
    std::shared_ptr<const ModelType> model,
        std::shared_ptr<CORE::LinearSOE<VectorType, MatrixType>> linearSOE,
        const dealii::Point<dim>& O1, const dealii::Point<dim>& O2)
        : _problem(problem)
        , _model(model)
        , _linearSOE(linearSOE)
    {
        _sections.push_back(Section<dim>{{}, O1});
        _sections.push_back(Section<dim>{{}, O2});

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
                        _sections[0].vertices[cell->vertex_index(v)] = point;
                    } else if (std::abs(point(0) - O2(0)) < posTol) {
                        _sections[1].vertices[cell->vertex_index(v)] = point;
                    }
                }
            }
        }

/*
        for (auto ip : cLeft) {
            std::cout << ip->index() << std::endl;
        }
*/
        std::cout << "# of sections: " << _sections.size() << std::endl;
        std::cout << "# of left nodes:  " << _sections[0].vertices.size() << std::endl;
        std::cout << "# of right nodes: " << _sections[1].vertices.size() << std::endl;
    }

protected:

    const std::vector<Section<dim>>& getSections() const {
        return _sections;
    }

    std::shared_ptr<CORE::TableEntryScalar<double>> createTableEntry(const dealii::Point<dim>& center, unsigned int componentNum) {

        unsigned int precision = 10;

        auto scalarReaction = std::make_shared<ScalarQuantitySectionResultant<dim, VectorType, MatrixType, ModelType>>(
            _problem, _model, _linearSOE, center, componentNum);
            
        return std::make_shared<CORE::TableEntryScalar<double>>(scalarReaction, precision);
    }

};

#endif // CASE_TWO_SECTIONS_H