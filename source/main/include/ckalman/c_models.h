#ifndef __C_KALMAN_MODELS_H__
#define __C_KALMAN_MODELS_H__
#include "ccore/c_target.h"
#ifdef USE_PRAGMA_ONCE
    #pragma once
#endif

#include "ckalman/c_math.h"
#include "ckalman/c_memory.h"

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

            // model is used to initialize hidden states in the model and
            // provide transition matrices to the filter.
            // kalman/models provides commonly used models.
            class model_t
            {
            public:
                virtual void             InitialState(state_t &state)                = 0;
                virtual nmath::matrix_t *Transition(memory_t *mem, u64 dt)           = 0;
                virtual nmath::matrix_t *CovarianceTransition(memory_t *mem, u64 dt) = 0;
            };

            // -----------------------------------------------------------------------------------------------------------------
            // -----------------------------------------------------------------------------------------------------------------
            // Constant Velocity Model

            struct constantvelocitymodel_config_t
            {
                f32 m_InitialVariance;
                f32 m_ProcessVariance;
            };

            class constantvelocitymodel_t : public model_t
            {
            public:
                virtual measurement_t    NewPositionMeasurement(memory_t *mem, nmath::vector_t *position, f32 measurementVariance) = 0;
                virtual nmath::vector_t *Position(memory_t *mem, nmath::vector_t *state)                                           = 0;
                virtual nmath::vector_t *Velocity(memory_t *mem, nmath::vector_t *state)                                           = 0;
            };

            constantvelocitymodel_t *NewConstantVelocityModel(memory_t *mem, u64 initialTime, nmath::vector_t *initialPosition, constantvelocitymodel_config_t cfg);

            // -----------------------------------------------------------------------------------------------------------------
            // -----------------------------------------------------------------------------------------------------------------
            // Simple Model

            struct simplemodel_config_t
            {
                f32 m_InitialVariance;
                f32 m_ProcessVariance;
                f32 m_ObservationVariance;
            };

            class simplemodel_t : public model_t
            {
            public:
                virtual measurement_t NewMeasurement(memory_t *mem, f32 value) = 0;
                virtual f32           Value(nmath::vector_t *state)            = 0;
            };

            simplemodel_t *NewSimpleModel(memory_t *mem, u64 initialTime, f32 initialValue, simplemodel_config_t cfg);

            // -----------------------------------------------------------------------------------------------------------------
            // -----------------------------------------------------------------------------------------------------------------
            // Brownian Motion Model

            struct brownianmodel_config_t
            {
                f32 m_InitialVariance;
                f32 m_ProcessVariance;
                f32 m_ObservationVariance;
            };

            class brownianmodel_t : public model_t
            {
            public:
                virtual measurement_t    NewMeasurement(memory_t *mem, nmath::vector_t *value) = 0;
                virtual nmath::vector_t *Value(nmath::vector_t *state)                         = 0;
            };

            brownianmodel_t *NewBrownianModel(memory_t *mem, u64 initialTime, nmath::vector_t *initialVector, brownianmodel_config_t cfg);

        }  // namespace nmodels
    }  // namespace nkalman
}  // namespace ncore
#endif  // __C_KALMAN_MODELS_H__
