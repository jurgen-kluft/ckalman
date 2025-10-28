#ifndef __C_KALMAN_MODELS_H__
#define __C_KALMAN_MODELS_H__
#include "ccore/c_target.h"
#ifdef USE_PRAGMA_ONCE
#    pragma once
#endif

#include "ckalman/c_math.h"

namespace ncore
{
    namespace nkalman
    {
        namespace nmodels
        {
            struct state_t
            {
                u64              m_Time;
                nmath::vector_t *m_State;
                nmath::matrix_t *m_Covariance;
            };

            struct measurement_t
            {
                nmath::matrix_t *m_Covariance;
                nmath::vector_t *m_Value;
                nmath::matrix_t *m_ObservationModel;
            };

            // model_t is used to initialize hidden states in the model and
            // provide transition matrices to the filter.
            // kalman/models provides commonly used models.
            class model_t
            {
            public:
                virtual void InitialState(state_t &state)                              = 0;
                virtual void Transition(u64 dt, nmath::matrix_t *&outMatrix)           = 0;
                virtual void CovarianceTransition(u64 dt, nmath::matrix_t *&outMatrix) = 0;
            };

            struct constantvelocitymodel_config_t
            {
                f64 m_InitialVariance;
                f64 m_ProcessVariance;
            };

            model_t *NewConstantVelocityModel(u64 model_t, nmath::vector_t *initialPosition, constantvelocitymodel_config_t cfg);

        }  // namespace nmodels
    }  // namespace nkalman
}  // namespace ncore
#endif  // __C_KALMAN_MODELS_H__
