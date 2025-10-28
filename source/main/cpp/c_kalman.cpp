#include "ckalman/c_kalman.h"
#include "ckalman/c_models.h"
#include "ckalman/c_math.h"

#include "ccore/c_debug.h"
#include "ccore/c_memory.h"

namespace ncore
{
    namespace nkalman
    {
        struct memory_t
        {
            void *m_memory_base;
            u32   m_memory_size;
            u8   *m_memory_current;

            enum
            {
                MEMORY_SCOPE_MAX = 8
            };

            u8 *m_scopes[MEMORY_SCOPE_MAX];
            s32 m_scope_index;

            nmath::vector_t *AllocVector(s32 n)
            {
                u32 const size = sizeof(nmath::vector_t) + sizeof(float) * n;
                ASSERT((u32)(m_memory_current + size - (u8 *)m_memory_base) <= m_memory_size);
                nmath::vector_t *v = (nmath::vector_t *)m_memory_current;
                m_memory_current += sizeof(nmath::vector_t);
                v->N    = n;
                v->Inc  = 1;
                v->data = (float *)m_memory_current;
                m_memory_current += sizeof(float) * n;
                return v;
            }

            nmath::vector_t *AllocZeroVector(s32 n)
            {
                nmath::vector_t *v = AllocVector(n);
                nmem::memset(v->data, 0, sizeof(float) * n);
                return v;
            }

            nmath::matrix_t *AllocMatrix(s32 rows, s32 cols)
            {
                u32 const size = sizeof(nmath::matrix_t) + sizeof(float) * rows * cols;
                ASSERT((u32)(m_memory_current + size - (u8 *)m_memory_base) <= m_memory_size);
                nmath::matrix_t *m = (nmath::matrix_t *)m_memory_current;
                m_memory_current += size;
                m->Rows    = rows;
                m->Cols    = cols;
                m->stride  = cols;
                m->capRows = rows;
                m->capCols = cols;
                m->data    = (float *)m_memory_current;
                return m;
            }

            nmath::matrix_t *AllocZeroMatrix(s32 rows, s32 cols)
            {
                nmath::matrix_t *m = AllocMatrix(rows, cols);
                nmem::memset(m->data, 0, sizeof(float) * rows * cols);
                return m;
            }

            void PushScope()
            {
                ASSERT(m_scope_index + 1 < MEMORY_SCOPE_MAX);
                m_scopes[m_scope_index++] = m_memory_current;
            }

            void PopScope()
            {
                ASSERT(m_scope_index - 1 >= 0);
                m_memory_current = m_scopes[--m_scope_index];
#ifdef TARGET_DEBUG
                // In debug mode, clear memory to catch use-after-free bugs
                u32 const size = m_scopes[m_scope_index] - m_memory_current;
                nmem::memset(m_memory_current, 0xCD, (u64)size);
#endif
            }

            void Reset() { m_memory_current = (u8 *)m_memory_base; }
        };

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
            km->m_dims       = initial.State->Len();
            km->m_t          = initial.Time;
            km->m_state      = nmath::Copy(initial.State);
            km->m_covariance = nmath::Copy(initial.Covariance);

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

            nmath::vector_t *z = m->Value;
            nmath::matrix_t *R = m->Covariance;
            nmath::matrix_t *H = m->ObservationModel;
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
