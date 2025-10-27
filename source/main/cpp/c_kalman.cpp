#include "ckalman/c_kalman.h"

namespace ncore
{
    namespace nkalman
    {
        namespace nmath
        {
            struct vector_t
            {
                s32    N;
                s32    Inc;
                float* data;

                s32 Len() const { return N; }

                void MulVec(matrix_t* m, vector_t* v)
                {
                    // TODO implement matrix-vector multiplication
                    // Note: Incoming vector v can be the same as this!
                }

                void AddVec(vector_t* m, vector_t* v)
                {
                    // TODO implement vector addition
                    // Note: Incoming vector m or v can be the same as this!
                }

                void SubVec(vector_t* m, vector_t* v)
                {
                    // TODO implement vector subtraction
                    // Note: Incoming vector m or v can be the same as this!
                }
            };

            struct matrix_t
            {
                s32    Rows;
                s32    Cols;
                float* data;
                s32    stride;
                s32    capRows;
                s32    capCols;

                void Set(s32 row, s32 col, float value) { data[row * stride + col] = value; }

                void Product(matrix_t* T, matrix_t* P, matrix_t* transposedT)
                {
                    // TODO implement matrix-matrix multiplication
                }

                void Add(matrix_t* a, matrix_t* b)
                {
                    // TODO implement matrix addition
                }

                void Sub(matrix_t* a, matrix_t* b)
                {
                    // TODO implement matrix subtraction
                }

                void Mul(matrix_t* a, matrix_t* b)
                {
                    // TODO implement matrix multiplication
                }

                void Inverse(matrix_t* m)
                {
                    // TODO implement matrix inversion
                }
            };

            vector_t* NewVector(s32 n, float* data)
            {
                // TODO memory allocation
                vector_t* v = nullptr;
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

            matrix_t* NewMatrix(s32 rows, s32 cols, float* data)
            {
                // TODO memory allocation
                matrix_t* m = nullptr;
                m->Rows     = rows;
                m->Cols     = cols;
                m->stride   = cols;
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

            vector_t* Copy(vector_t* v)
            {
                vector_t* copy = NewVector(v->N, nullptr);
                for (s32 i = 0; i < v->N; i++)
                {
                    copy->data[i] = v->data[i];
                }
                return copy;
            }

            void CopyContent(vector_t* dest, vector_t* src)
            {
                for (s32 i = 0; i < src->N; i++)
                {
                    dest->data[i] = src->data[i];
                }
            }

            void CopyContent(matrix_t* dest, matrix_t* src)
            {
                for (s32 i = 0; i < src->Rows; i++)
                {
                    for (s32 j = 0; j < src->Cols; j++)
                    {
                        dest->data[i * dest->stride + j] = src->data[i * src->stride + j];
                    }
                }
            }

            matrix_t* Copy(matrix_t* m)
            {
                matrix_t* copy = NewMatrix(m->Rows, m->Cols, nullptr);
                for (s32 i = 0; i < m->Rows; i++)
                {
                    for (s32 j = 0; j < m->Cols; j++)
                    {
                        copy->data[i * copy->stride + j] = m->data[i * m->stride + j];
                    }
                }
                return copy;
            }

            matrix_t* Transpose(matrix_t* m)
            {
                matrix_t* transposed = NewMatrix(m->Cols, m->Rows, nullptr);
                for (s32 i = 0; i < m->Rows; i++)
                {
                    for (s32 j = 0; j < m->Cols; j++)
                    {
                        transposed->data[j * transposed->stride + i] = m->data[i * m->stride + j];
                    }
                }
                return transposed;
            }

            void Free(vector_t* v)
            {
                // TODO free memory
                // free(v->data);
                // free(v);
            }

            void Free(matrix_t* m)
            {
                // TODO free memory
                // free(m->data);
                // free(m);
            }

        }  // namespace nmath

        // NewKalmanFilter returns a new filter_t for the given linear model.
        filter_t* NewFilter(nmodels::model_t* model)
        {
            nmodels::state_t initial;
            model->InitialState(initial);

            filter_t* km = nullptr;

            km->model      = model;
            km->dims       = initial.state_t->Len();
            km->t          = initial.Time;
            km->state      = nmath::Copy(initial.state_t);
            km->covariance = nmath::Copy(initial.Covariance);

            return km;
        }

        nmath::vector_t*   state_t(filter_t* kf) { return kf->state; }
        nmath::matrix_t* Covariance(filter_t* kf) { return kf->covariance; }

        void SetCovariance(filter_t* kf, nmath::matrix_t* covariance) { nmath::CopyContent(kf->covariance, covariance); }

        void SetState(filter_t* kf, nmath::vector_t* state) { nmath::CopyContent(kf->state, state); }

        u64 Time(filter_t* kf) { return kf->t; }

        bool Predict(filter_t* kf, u64 t)
        {
            if (t <= kf->t)
                return false;  // can't predict past

            if (t == kf->t)
                return true;

            const u64 dt = t - kf->t;
            kf->t        = t;

            nmath::matrix_t *T, *Q;
            kf->model->Transition(dt, T);
            kf->model->CovarianceTransition(dt, Q);
            nmath::matrix_t* P = kf->covariance;

            nmath::vector_t* currState = nmath::Copy(kf->state);
            kf->state->MulVec(T, currState);
            nmath::Free(currState);

            nmath::matrix_t* newCovariance = nmath::NewMatrix(kf->dims, kf->dims, nullptr);
            nmath::matrix_t* transposedT   = nmath::Transpose(T);
            newCovariance->Product(T, P, transposedT);
            kf->covariance->Add(newCovariance, Q);

            return true;
        }

        bool Update(filter_t* kf, u64 t, nmodels::measurement_t* m)
        {
            if (t <= kf->t)
                return false;  // can't predict past

            if (!Predict(kf, t))
                return false;

            nmath::vector_t*   z = m->Value;
            nmath::matrix_t* R = m->Covariance;
            nmath::matrix_t* H = m->ObservationModel;
            nmath::matrix_t* P = kf->covariance;

            nmath::vector_t* preFitResidual = nmath::NewVector(z->Len(), nullptr);
            preFitResidual->MulVec(H, kf->state);
            preFitResidual->SubVec(z, preFitResidual);

            nmath::matrix_t* preFitResidualCov = nmath::NewMatrix(z->Len(), z->Len(), nullptr);
            nmath::matrix_t* transposedH       = nmath::Transpose(H);
            preFitResidualCov->Product(H, P, transposedH);
            preFitResidualCov->Add(preFitResidualCov, R);

            nmath::matrix_t* preFitResidualCovInv = nmath::NewMatrix(z->Len(), z->Len(), nullptr);
            preFitResidualCovInv->Inverse(preFitResidualCov);

            nmath::matrix_t* gain = nmath::NewMatrix(kf->dims, z->Len(), nullptr);
            gain->Product(P, transposedH, preFitResidualCovInv);

            nmath::vector_t* newState = nmath::NewVector(kf->dims, nullptr);
            newState->MulVec(gain, preFitResidual);
            newState->AddVec(kf->state, newState);

            nmath::matrix_t* newCovariance = nmath::NewMatrix(kf->dims, kf->dims, nullptr);
            newCovariance->Mul(gain, H);
            nmath::matrix_t* eye = Eye(kf, kf->dims);
            newCovariance->Sub(eye, newCovariance);
            newCovariance->Mul(newCovariance, P);
            nmath::Free(eye);

            nmath::Free(kf->covariance);
            nmath::Free(kf->state);

            kf->covariance = newCovariance;
            kf->state      = newState;
            kf->t          = t;

            return true;
        }

        nmath::matrix_t* Eye(filter_t* kf, s32 n)
        {
            nmath::matrix_t* result = nmath::NewMatrix(n, n, nullptr);
            for (s32 i = 0; i < n; i++)
                result->Set(i, i, 1.0);
            return result;
        }

    }  // namespace nkalman
}  // namespace ncore
