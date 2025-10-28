#include "ckalman/c_kalman.h"
#include "ckalman/c_models.h"

namespace ncore
{
    namespace nkalman
    {
        namespace nmodels
        {
            class ConstantVelocityModel : public model_t
            {
            public:
                virtual void InitialState(state_t &state) override
                {
                    state.m_Time = m_initialState->m_Time;
                    state.m_State = nmath::Copy(m_initialState->m_State);
                    state.m_Covariance = nmath::Copy(m_initialState->m_Covariance);
                }

                virtual void Transition(u64 dt, nmath::matrix_t *&outMatrix) override
                {
                }

                virtual void CovarianceTransition(u64 dt, nmath::matrix_t *&outMatrix) override
                {
                }

                state_t *m_initialState;
                s32 m_dims;
                s32 m_stateDims;
                constantvelocitymodel_config_t m_cfg;
            };

            model_t *NewConstantVelocityModel(u64 model_t, nmath::vector_t *initialPosition, constantvelocitymodel_config_t cfg)
            {
                return nullptr;
            }

        } // namespace nmodels
    } // namespace nkalman
} // namespace ncore
