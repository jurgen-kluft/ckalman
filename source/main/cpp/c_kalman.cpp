#include "ckalman/c_kalman.h"
#include "ckalman/c_models.h"
#include "ckalman/c_math.h"
#include "ckalman/c_memory.h"

#include "ccore/c_debug.h"
#include "ccore/c_memory.h"

namespace ncore
{
    namespace nkalman
    {
        // filter_t is responsible for prediction and filtering
        // of a given linear model. It is assumed that the process being modelled
        // is a time series, and that the time steps are non-uniform and specified
        // for each update and prediction operation.
        struct filter_t
        {
            nmodels::model_t *m_model;
            memory_t         *m_memory;
            s32               m_dims;
            u64               m_t;
            nmath::vector_t  *m_state;
            nmath::matrix_t  *m_covariance;
        };

        // NewFilter returns a new filter object for the given linear model.
        filter_t *NewFilter(nmodels::model_t *model, memory_t *mem)
        {
            nmodels::state_t initial;
            model->InitialState(initial);

            filter_t *kf = (filter_t *)mem->AllocZeroMemory(sizeof(filter_t), alignof(filter_t));

            kf->m_model      = model;
            kf->m_memory     = mem;
            kf->m_dims       = initial.m_State->Len();
            kf->m_t          = initial.m_Time;
            kf->m_state      = nmath::Duplicate(kf->m_memory, initial.m_State);
            kf->m_covariance = nmath::Duplicate(kf->m_memory, initial.m_Covariance);

            return kf;
        }

        nmath::vector_t *State(filter_t *kf) { return kf->m_state; }
        nmath::matrix_t *Covariance(filter_t *kf) { return kf->m_covariance; }
        void             SetCovariance(filter_t *kf, nmath::matrix_t *covariance) { kf->m_covariance->Copy(covariance); }
        void             SetState(filter_t *kf, nmath::vector_t *state) { kf->m_state->Copy(state); }
        u64              Time(filter_t *kf) { return kf->m_t; }

        bool Predict(filter_t *kf, u64 t)
        {
            if (t <= kf->m_t)
                return false;  // can't predict past

            if (t == kf->m_t)
                return true;

            const u64 dt = t - kf->m_t;
            kf->m_t      = t;

            nmath::matrix_t *T = kf->m_model->Transition(kf->m_memory, dt);
            nmath::matrix_t *Q = kf->m_model->CovarianceTransition(kf->m_memory, dt);
            nmath::matrix_t *P = kf->m_covariance;

            kf->m_state->MulVec(kf->m_memory, T, kf->m_state);

            kf->m_memory->PushScope();
            {
                nmath::matrix_t *newCovariance = nmath::NewMatrix(kf->m_memory, kf->m_dims, kf->m_dims, nullptr);
                nmath::matrix_t *transposedT   = nmath::NewMatrix(kf->m_memory, T->m_cols, T->m_rows, nullptr);
                transposedT->Transpose(kf->m_memory, T);
                newCovariance->Product(kf->m_memory, T, P, transposedT);
                kf->m_covariance->Add(kf->m_memory, newCovariance, Q);
            }
            kf->m_memory->PopScope();

            return true;
        }

        bool Update(filter_t *kf, u64 t, nmodels::measurement_t *m)
        {
            if (t <= kf->m_t)
                return false;  // can't predict past

            if (!Predict(kf, t))
                return false;

            nmath::vector_t *z = m->m_Value;
            nmath::matrix_t *R = m->m_Covariance;
            nmath::matrix_t *H = m->m_ObservationModel;
            nmath::matrix_t *P = kf->m_covariance;

            kf->m_memory->PushScope();
            {
                nmath::vector_t *preFitResidual = nmath::NewVector(kf->m_memory, z->Len(), nullptr);
                preFitResidual->MulVec(kf->m_memory, H, kf->m_state);
                preFitResidual->SubVec(kf->m_memory, z, preFitResidual);

                nmath::matrix_t *transposedH = nmath::NewMatrix(kf->m_memory, H->m_cols, H->m_rows, nullptr);
                transposedH->Transpose(kf->m_memory, H);

                nmath::matrix_t *gain = nmath::NewMatrix(kf->m_memory, kf->m_dims, z->Len(), nullptr);
                kf->m_memory->PushScope();
                {
                    nmath::matrix_t *preFitResidualCov = nmath::NewMatrix(kf->m_memory, z->Len(), z->Len(), nullptr);
                    preFitResidualCov->Product(kf->m_memory, H, P, transposedH);
                    preFitResidualCov->Add(kf->m_memory, preFitResidualCov, R);

                    nmath::matrix_t *preFitResidualCovInv = nmath::NewMatrix(kf->m_memory, z->Len(), z->Len(), nullptr);
                    preFitResidualCovInv->Inverse(kf->m_memory, preFitResidualCov);

                    gain->Product(kf->m_memory, P, transposedH, preFitResidualCovInv);
                }
                kf->m_memory->PopScope();

                nmath::vector_t *newState = nmath::NewVector(kf->m_memory, kf->m_dims, nullptr);
                newState->MulVec(kf->m_memory, gain, preFitResidual);
                newState->AddVec(kf->m_memory, kf->m_state, newState);

                nmath::matrix_t *newCovariance = nmath::NewMatrix(kf->m_memory, kf->m_dims, kf->m_dims, nullptr);
                newCovariance->Mul(kf->m_memory, gain, H);
                kf->m_memory->PushScope();
                {
                    nmath::matrix_t *eye = Eye(kf->m_memory, kf, kf->m_dims);
                    newCovariance->Sub(kf->m_memory, eye, newCovariance);
                }
                kf->m_memory->PopScope();
                newCovariance->Mul(kf->m_memory, newCovariance, P);

                kf->m_state->Copy(newState);
                kf->m_covariance->Copy(newCovariance);
                kf->m_t = t;
            }
            kf->m_memory->PopScope();

            return true;
        }

        nmath::matrix_t *Eye(memory_t *mem, filter_t *kf, s32 n)
        {
            nmath::matrix_t *result = nmath::NewMatrix(mem, n, n, nullptr);
            for (s32 i = 0; i < n; i++)
                result->Set(i, i, 1.0);
            return result;
        }

    }  // namespace nkalman
}  // namespace ncore
