#ifndef STUDY_CASE_H
#define STUDY_CASE_H

#include "core/model.h"
#include "core/table_output.h"

#include <deal.II/base/point.h>

template <int dim, typename VectorType>
class StudyCase
{
public:
    virtual void imposeBoundaryConditions(std::shared_ptr<CORE::Model<dim, VectorType>> model) = 0;
    virtual void addReactions(std::shared_ptr<CORE::TableOutput> table) = 0;
    virtual double loadMagnitude() const = 0;
};

#endif // STUDY_CASE_H