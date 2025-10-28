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
        filter_t *NewFilter(nmodels::model_t *model)
        {
            nmodels::state_t initial;
            model->InitialState(initial);

            // TODO memory allocation
            filter_t *km = nullptr;

            km->m_model      = model;
            km->m_dims       = initial.m_State->Len();
            km->m_t          = initial.m_Time;
            km->m_state      = nmath::Copy(initial.m_State);
            km->m_covariance = nmath::Copy(initial.m_Covariance);

            return km;
        }

        nmath::vector_t *State(filter_t *kf) { return kf->m_state; }
        nmath::matrix_t *Covariance(filter_t *kf) { return kf->m_covariance; }
        void             SetCovariance(filter_t *kf, nmath::matrix_t *covariance) { nmath::CopyContent(kf->m_covariance, covariance); }
        void             SetState(filter_t *kf, nmath::vector_t *state) { nmath::CopyContent(kf->m_state, state); }
        u64              Time(filter_t *kf) { return kf->m_t; }

        bool Predict(filter_t *kf, u64 t)
        {
            if (t <= kf->m_t)
                return false;  // can't predict past

            if (t == kf->m_t)
                return true;

            const u64 dt = t - kf->m_t;
            kf->m_t      = t;

            nmath::matrix_t *T, *Q;
            kf->m_model->Transition(dt, T);
            kf->m_model->CovarianceTransition(dt, Q);
            nmath::matrix_t *P = kf->m_covariance;

            nmath::vector_t *currState = nmath::Copy(kf->m_state);
            kf->m_state->MulVec(T, currState);
            nmath::Free(currState);

            nmath::matrix_t *newCovariance = nmath::NewMatrix(kf->m_dims, kf->m_dims, nullptr);
            nmath::matrix_t *transposedT   = nmath::Transpose(T);
            newCovariance->Product(T, P, transposedT);
            kf->m_covariance->Add(newCovariance, Q);

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

            nmath::vector_t *preFitResidual = nmath::NewVector(z->Len(), nullptr);
            preFitResidual->MulVec(H, kf->m_state);
            preFitResidual->SubVec(z, preFitResidual);

            nmath::matrix_t *preFitResidualCov = nmath::NewMatrix(z->Len(), z->Len(), nullptr);
            nmath::matrix_t *transposedH       = nmath::Transpose(H);
            preFitResidualCov->Product(H, P, transposedH);
            preFitResidualCov->Add(preFitResidualCov, R);

            nmath::matrix_t *preFitResidualCovInv = nmath::NewMatrix(z->Len(), z->Len(), nullptr);
            preFitResidualCovInv->Inverse(preFitResidualCov);

            nmath::matrix_t *gain = nmath::NewMatrix(kf->m_dims, z->Len(), nullptr);
            gain->Product(P, transposedH, preFitResidualCovInv);

            nmath::vector_t *newState = nmath::NewVector(kf->m_dims, nullptr);
            newState->MulVec(gain, preFitResidual);
            newState->AddVec(kf->m_state, newState);

            nmath::matrix_t *newCovariance = nmath::NewMatrix(kf->m_dims, kf->m_dims, nullptr);
            newCovariance->Mul(gain, H);
            nmath::matrix_t *eye = Eye(kf, kf->m_dims);
            newCovariance->Sub(eye, newCovariance);
            newCovariance->Mul(newCovariance, P);
            nmath::Free(eye);

            nmath::Free(kf->m_covariance);
            nmath::Free(kf->m_state);

            kf->m_covariance = newCovariance;
            kf->m_state      = newState;
            kf->m_t          = t;

            return true;
        }

        nmath::matrix_t *Eye(filter_t *kf, s32 n)
        {
            nmath::matrix_t *result = nmath::NewMatrix(n, n, nullptr);
            for (s32 i = 0; i < n; i++)
                result->Set(i, i, 1.0);
            return result;
        }

    }  // namespace nkalman
}  // namespace ncore
