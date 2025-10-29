#include "ccore/c_allocator.h"

#include "ckalman/c_kalman.h"
#include "ckalman/c_math.h"
#include "ckalman/c_memory.h"
#include "ckalman/c_models.h"

namespace ncore
{
    namespace nkalman
    {
        namespace nmodels
        {
            class constantvelocitymodel_imp_t : public constantvelocitymodel_t
            {
            public:
                virtual void             InitialState(state_t &state);
                virtual nmath::matrix_t *Transition(memory_t *mem, u64 dt);
                virtual nmath::matrix_t *CovarianceTransition(memory_t *mem, u64 dt);

                virtual measurement_t    NewPositionMeasurement(memory_t *mem, nmath::vector_t *position, f64 measurementVariance);
                virtual nmath::vector_t *Position(memory_t *mem, nmath::vector_t *state);
                virtual nmath::vector_t *Velocity(memory_t *mem, nmath::vector_t *state);

                DCORE_CLASS_PLACEMENT_NEW_DELETE

                state_t                       *m_initialState;
                memory_t                      *m_memory;
                s32                            m_dims;
                s32                            m_stateDims;
                constantvelocitymodel_config_t m_cfg;
            };

            constantvelocitymodel_t *NewConstantVelocityModel(memory_t *mem, u64 time, nmath::vector_t *initialPosition, constantvelocitymodel_config_t cfg)
            {
                const s32 dims      = initialPosition->Len();
                const s32 stateDims = dims * 2;

                nmath::matrix_t *initialCovariance = mem->AllocZeroMatrix(stateDims, stateDims);
                for (s32 i = 0; i < stateDims; i++)
                {
                    initialCovariance->Set(i, i, (f32)cfg.m_InitialVariance);
                }

                nmath::vector_t *initialStateVector = mem->AllocZeroVector(stateDims);
                for (s32 i = 0; i < dims; i++)
                {
                    initialStateVector->SetVec(i, initialPosition->AtVec(i));
                }

                state_t *initialState      = (state_t *)mem->AllocMemory(sizeof(state_t), alignof(state_t));
                initialState->m_Time       = time;
                initialState->m_State      = initialStateVector;
                initialState->m_Covariance = initialCovariance;

                void                        *model_mem = mem->AllocMemory(sizeof(constantvelocitymodel_imp_t), alignof(constantvelocitymodel_imp_t));
                constantvelocitymodel_imp_t *model     = new (model_mem) constantvelocitymodel_imp_t();
                model->m_initialState                  = initialState;
                model->m_dims                          = dims;
                model->m_stateDims                     = stateDims;
                model->m_cfg                           = cfg;

                return model;
            }

            void constantvelocitymodel_imp_t ::InitialState(state_t &state)
            {
                state.m_Time       = m_initialState->m_Time;
                state.m_State      = m_initialState->m_State;
                state.m_Covariance = m_initialState->m_Covariance;
            }

            nmath::matrix_t *constantvelocitymodel_imp_t::Transition(memory_t *mem, u64 dt)
            {
                nmath::matrix_t *outMatrix = mem->AllocZeroMatrix(m_stateDims, m_stateDims);
                for (s32 i = 0; i < m_stateDims; i++)
                {
                    outMatrix->Set(i, i, 1.0f);
                }
                f32 dts = static_cast<f32>(dt) * 0.001f;
                for (s32 i = 0; i < m_dims; i++)
                {
                    outMatrix->Set(i, m_dims + i, dts);
                }
            }

            nmath::matrix_t *constantvelocitymodel_imp_t::CovarianceTransition(memory_t *mem, u64 dt)
            {
                nmath::matrix_t *outMatrix = mem->AllocZeroMatrix(m_stateDims, m_stateDims);
                f32              v         = static_cast<f32>(dt) * 0.001f * static_cast<f32>(m_cfg.m_ProcessVariance);
                for (s32 i = 0; i < m_stateDims; i++)
                {
                    outMatrix->Set(i, i, v);
                }
            }

            measurement_t constantvelocitymodel_imp_t::NewPositionMeasurement(memory_t *mem, nmath::vector_t *position, f64 measurementVariance)
            {
                measurement_t measurement;
                measurement.m_Value            = nullptr;
                measurement.m_Covariance       = nullptr;
                measurement.m_ObservationModel = nullptr;

                if (position->Len() != m_dims)
                {
                    ASSERT(false);  // printf("position vector has incorrect number of entries: %d (expected %d)", position->Len(), m.dims))
                    return measurement;
                }
                nmath::matrix_t *covariance = mem->AllocZeroMatrix(m_dims, m_dims);
                for (s32 i = 0; i < m_dims; i++)
                {
                    covariance->Set(i, i, (f32)measurementVariance);
                }
                nmath::matrix_t *observationModel = mem->AllocZeroMatrix(m_dims, m_stateDims);
                for (s32 i = 0; i < m_dims; i++)
                {
                    observationModel->Set(i, i, 1.0f);
                }
                measurement.m_Value            = nmath::Duplicate(mem, position);
                measurement.m_Covariance       = covariance;
                measurement.m_ObservationModel = observationModel;
                return measurement;
            }

            nmath::vector_t *constantvelocitymodel_imp_t::Position(memory_t *mem, nmath::vector_t *state)
            {
                if (state->Len() != m_stateDims)
                {
                    ASSERT(false);  // printf("state vector has incorrect number of entries: %d (expected %d)", state->Len(), m_stateDims))
                    return nullptr;
                }
                nmath::vector_t *result = mem->AllocZeroVector(m_dims);
                for (s32 i = 0; i < m_dims; i++)
                {
                    result->SetVec(i, state->AtVec(i));
                }
                return result;
            }

            nmath::vector_t *constantvelocitymodel_imp_t::Velocity(memory_t *mem, nmath::vector_t *state)
            {
                if (state->Len() != m_stateDims)
                {
                    ASSERT(false);  // printf("state vector has incorrect number of entries: %d (expected %d)", state->Len(), m_stateDims))
                    return nullptr;
                }
                nmath::vector_t *result = mem->AllocZeroVector(m_dims);
                for (s32 i = 0; i < m_dims; i++)
                {
                    result->SetVec(i, state->AtVec(i + m_dims));
                }
                return result;
            }

        }  // namespace nmodels
    }  // namespace nkalman
}  // namespace ncore
