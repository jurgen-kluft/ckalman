#include "ccore/c_allocator.h"
#include "ckalman/c_kalman.h"

namespace ncore
{
    namespace nmath
    {
        struct Vector
        {
            s32    N;
            s32    Inc;
            float* data;

            s32 Len() const { return N; }

            void MulVec(Matrix* m, Vector* v)
            {
                // TODO implement matrix-vector multiplication
                // Note: Incoming vector v can be the same as this!
            }

            void AddVec(Vector* m, Vector* v)
            {
                // TODO implement vector addition
                // Note: Incoming vector m or v can be the same as this!
            }

            void SubVec(Vector* m, Vector* v)
            {
                // TODO implement vector subtraction
                // Note: Incoming vector m or v can be the same as this!
            }
        };

        struct Matrix
        {
            s32    Rows;
            s32    Cols;
            float* data;
            s32    stride;
            s32    capRows;
            s32    capCols;

            void Set(s32 row, s32 col, float value) { data[row * stride + col] = value; }

            void Product(Matrix* T, Matrix* P, Matrix* transposedT)
            {
                // TODO implement matrix-matrix multiplication
            }

            void Add(Matrix* a, Matrix* b)
            {
                // TODO implement matrix addition
            }
        };

        Vector* NewVector(s32 n, float* data)
        {
            // TODO memory allocation
            Vector* v = nullptr;
            v->N      = n;
            v->Inc    = 1;
            if (data == nullptr)
            {
                // TODO memory allocation
                // v->data = (float*)malloc(sizeof(float) * n);
                // all elements initialized to zero
            }
            else
            {
                v->data = data;
            }
            return v;
        }

        Matrix* NewMatrix(s32 rows, s32 cols, float* data)
        {
            // TODO memory allocation
            Matrix* m = nullptr;
            m->Rows   = rows;
            m->Cols   = cols;
            m->stride = cols;
            if (data == nullptr)
            {
                // TODO memory allocation
                // m->data = (float*)malloc(sizeof(float) * rows * cols);
                // all elements initialized to zero
            }
            else
            {
                m->data = data;
            }
            return m;
        }

        Vector* Copy(Vector* v)
        {
            Vector* copy = NewVector(v->N, nullptr);
            for (s32 i = 0; i < v->N; i++)
            {
                copy->data[i] = v->data[i];
            }
            return copy;
        }

        void CopyContent(Vector* dest, Vector* src)
        {
            for (s32 i = 0; i < src->N; i++)
            {
                dest->data[i] = src->data[i];
            }
        }

        void CopyContent(Matrix* dest, Matrix* src)
        {
            for (s32 i = 0; i < src->Rows; i++)
            {
                for (s32 j = 0; j < src->Cols; j++)
                {
                    dest->data[i * dest->stride + j] = src->data[i * src->stride + j];
                }
            }
        }

        Matrix* Copy(Matrix* m)
        {
            Matrix* copy = NewMatrix(m->Rows, m->Cols, nullptr);
            for (s32 i = 0; i < m->Rows; i++)
            {
                for (s32 j = 0; j < m->Cols; j++)
                {
                    copy->data[i * copy->stride + j] = m->data[i * m->stride + j];
                }
            }
            return copy;
        }

        Matrix* Transpose(Matrix* m)
        {
            Matrix* transposed = NewMatrix(m->Cols, m->Rows, nullptr);
            for (s32 i = 0; i < m->Rows; i++)
            {
                for (s32 j = 0; j < m->Cols; j++)
                {
                    transposed->data[j * transposed->stride + i] = m->data[i * m->stride + j];
                }
            }
            return transposed;
        }

        void Free(Vector* v)
        {
            // TODO free memory
            // free(v->data);
            // free(v);
        }

        void Free(Matrix* m)
        {
            // TODO free memory
            // free(m->data);
            // free(m);
        }

    }  // namespace nmath

    namespace nkalman
    {
        // NewKalmanFilter returns a new KalmanFilter for the given linear model.
        KalmanFilter* NewKalmanFilter(nmodels::LinearModel* model)
        {
            nmodels::State initial;
            model->InitialState(initial);

            KalmanFilter* km = nullptr;

            km->model      = model;
            km->dims       = initial.State->Len();
            km->t          = initial.Time;
            km->state      = nmath::Copy(initial.State);
            km->covariance = nmath::Copy(initial.Covariance);

            return km;
        }

        nmath::Vector* State(KalmanFilter* kf) { return kf->state; }
        nmath::Matrix* Covariance(KalmanFilter* kf) { return kf->covariance; }

        void SetCovariance(KalmanFilter* kf, nmath::Matrix* covariance) { nmath::CopyContent(kf->covariance, covariance); }

        void SetState(KalmanFilter* kf, nmath::Vector* state) { nmath::CopyContent(kf->state, state); }

        u64 Time(KalmanFilter* kf) { return kf->t; }

        bool Predict(KalmanFilter* kf, u64 t)
        {
            if (t <= kf->t)
                return false;  // can't predict past

            if (t == kf->t)
                return true;

            const u64 dt = t - kf->t;
            kf->t        = t;

            nmath::Matrix *T, *Q;
            kf->model->Transition(dt, T);
            kf->model->CovarianceTransition(dt, Q);
            nmath::Matrix* P = kf->covariance;

            nmath::Vector* currState = nmath::Copy(kf->state);
            kf->state->MulVec(T, currState);
            nmath::Free(currState);

            nmath::Matrix* newCovariance = nmath::NewMatrix(kf->dims, kf->dims, nullptr);
            nmath::Matrix* transposedT   = nmath::Transpose(T);
            newCovariance->Product(T, P, transposedT);
            kf->covariance->Add(newCovariance, Q);

            return true;
        }

        // struct Measurement
        // {
        //     nmath::Matrix* Covariance;
        //     nmath::Vector* Value;
        //     nmath::Matrix* ObservationModel;
        // };

        bool Update(KalmanFilter* kf, u64 t, nmodels::Measurement* m)
        {
            if (t <= kf->t)
                return false;  // can't predict past

            if (!Predict(kf, t))
                return false;

            nmath::Vector* z = m->Value;
            nmath::Matrix* R = m->Covariance;
            nmath::Matrix* H = m->ObservationModel;
            nmath::Matrix* P = kf->covariance;

            nmath::Vector* preFitResidual = nmath::NewVector(z->Len(), nullptr);
            preFitResidual->MulVec(H, kf->state);
            preFitResidual->SubVec(z, preFitResidual);

            nmath::Matrix* preFitResidualCov = nmath::NewMatrix(z->Len(), z->Len(), nullptr);
            //     preFitResidualCov->Product(H, P, H->T());
            //     preFitResidualCov->Add(preFitResidualCov, R);

            nmath::Matrix* preFitResidualCovInv = nmath::NewMatrix(z->Len(), z->Len(), nullptr);
            //     preFitResidualCovInv->Inverse(preFitResidualCov);

            nmath::Matrix* gain = nmath::NewMatrix(kf->dims, z->Len(), nullptr);
            //     gain->Product(P, H->T(), preFitResidualCovInv);

            nmath::Vector* newState = nmath::NewVector(kf->dims, nullptr);
            newState->MulVec(gain, preFitResidual);
            newState->AddVec(kf->state, newState);

            nmath::Matrix* newCovariance = nmath::NewMatrix(kf->dims, kf->dims, nullptr);
            //     newCovariance->Mul(gain, H);
            //     newCovariance->Sub(eye(kf->dims), newCovariance);
            //     newCovariance->Mul(newCovariance, P);

            //     kf->covariance = newCovariance;
            //     kf->state = newState;
            //     kf->t = t;

            return true;
        }

        nmath::Matrix* Eye(KalmanFilter* kf, s32 n)
        {
            nmath::Matrix* result = nmath::NewMatrix(n, n, nullptr);
            for (s32 i = 0; i < n; i++)
                result->Set(i, i, 1.0);
            return result;
        }

    }  // namespace nkalman
}  // namespace ncore
