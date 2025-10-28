#include "ckalman/c_math.h"

namespace ncore
{
    namespace nkalman
    {
        namespace nmath
        {
            void vector_t::MulVec(matrix_t *m, vector_t *v)
            {
                vector_t *result = NewVector(m->m_Rows, nullptr);
                for (s32 i = 0; i < m->m_Rows; i++)
                {
                    float f = 0.0f;
                    for (s32 j = 0; j < m->m_Cols; j++)
                    {
                        f += m->At(i, j) * v->AtVec(j);
                    }
                    result->SetVec(i, f);
                }
                for (s32 i = 0; i < result->Len(); i++)
                    SetVec(i, result->AtVec(i));
                Free(result);
            }

            void vector_t::AddVec(vector_t *a, vector_t *b)
            {
                const s32 ar     = a->Len();
                vector_t *result = NewVector(ar, nullptr);
                for (s32 i = 0; i < ar; i++)
                    result->SetVec(i, a->AtVec(i) + b->AtVec(i));
                for (s32 i = 0; i < ar; i++)
                    SetVec(i, result->AtVec(i));
                Free(result);
            }

            void vector_t::SubVec(vector_t *a, vector_t *b)
            {
                const s32 ar     = a->Len();
                vector_t *result = NewVector(ar, nullptr);
                for (s32 i = 0; i < ar; i++)
                    result->SetVec(i, a->AtVec(i) - b->AtVec(i));
                for (s32 i = 0; i < ar; i++)
                    SetVec(i, result->AtVec(i));
                Free(result);
            }

            void matrix_t::Product(matrix_t *T, matrix_t *P, matrix_t *transposedT)
            {
                matrix_t *result = NewMatrix(T->m_Rows, P->m_Cols, nullptr);
                result->Mul2(T, P);
                Mul2(result, transposedT);
                Free(result);
            }

            void matrix_t::Add(matrix_t *a, matrix_t *b)
            {
                matrix_t *result = NewMatrix(a->m_Rows, a->m_Cols, nullptr);
                for (s32 i = 0; i < a->m_Rows; i++)
                {
                    for (s32 j = 0; j < a->m_Cols; j++)
                        result->Set(i, j, a->At(i, j) + b->At(i, j));
                }
                CopyContent(this, result);
                Free(result);
            }

            void matrix_t::Sub(matrix_t *a, matrix_t *b)
            {
                matrix_t *result = NewMatrix(a->m_Rows, a->m_Cols, nullptr);
                for (s32 i = 0; i < a->m_Rows; i++)
                {
                    for (s32 j = 0; j < a->m_Cols; j++)
                        result->Set(i, j, a->At(i, j) - b->At(i, j));
                }
                CopyContent(this, result);
                Free(result);
            }

            void matrix_t::Mul(matrix_t *a, matrix_t *b)
            {
                const s32 ar = a->m_Rows;
                const s32 ac = a->m_Cols;

                matrix_t *result = NewMatrix(ar, ac, nullptr);
                float    *row    = NewBuffer(ac);

                for (s32 i = 0; i < ar; i++)
                {
                    for (s32 k = 0; k < ac; k++)
                        row[k] = a->At(i, k);

                    for (s32 j = 0; j < b->m_Cols; j++)
                    {
                        float f = 0.0f;
                        for (s32 k = 0; k < ac; k++)
                            f += row[k] * b->At(k, j);
                        result->Set(i, j, f);
                    }
                }
                CopyContent(this, result);
                FreeBuffer(row);
                Free(result);
            }

            // We know for sure that 'this' is not 'a' or 'b'
            void matrix_t::Mul2(matrix_t *a, matrix_t *b)
            {
                const s32 ar  = a->m_Rows;
                const s32 ac  = a->m_Cols;
                float    *row = NewBuffer(ac);
                for (s32 i = 0; i < ar; i++)
                {
                    for (s32 k = 0; k < ac; k++)
                        row[k] = a->At(i, k);

                    for (s32 j = 0; j < b->m_Cols; j++)
                    {
                        float f = 0.0f;
                        for (s32 k = 0; k < ac; k++)
                            f += row[k] * b->At(k, j);
                        Set(i, j, f);
                    }
                }
                FreeBuffer(row);
            }

            void matrix_t::Inverse(matrix_t *m)
            {
                // TODO implement matrix inversion
            }

            vector_t *NewVector(s32 n, float *data)
            {
                // TODO memory allocation
                vector_t *v = nullptr;
                v->m_N      = n;
                v->m_Inc    = 1;
                if (data == nullptr)
                {
                    // TODO memory allocation
                    // v->m_data = (float*)malloc(sizeof(float) * n);
                    // all elements initialized to zero
                }
                else
                {
                    v->m_data = data;
                }
                return v;
            }

            matrix_t *NewMatrix(s32 rows, s32 cols, float *data)
            {
                // TODO memory allocation
                matrix_t *m = nullptr;
                m->m_Rows   = rows;
                m->m_Cols   = cols;
                m->m_stride = cols;
                if (data == nullptr)
                {
                    // TODO memory allocation
                    // m->m_data = (float*)malloc(sizeof(float) * rows * cols);
                    // all elements initialized to zero
                }
                else
                {
                    m->m_data = data;
                }
                return m;
            }

            vector_t *Copy(vector_t *v)
            {
                vector_t *copy = NewVector(v->m_N, nullptr);
                for (s32 i = 0; i < v->m_N; i++)
                {
                    copy->m_data[i] = v->m_data[i];
                }
                return copy;
            }

            matrix_t *Copy(matrix_t *m)
            {
                matrix_t *copy = NewMatrix(m->m_Rows, m->m_Cols, nullptr);
                for (s32 i = 0; i < m->m_Rows; i++)
                {
                    for (s32 j = 0; j < m->m_Cols; j++)
                    {
                        copy->m_data[i * copy->m_stride + j] = m->m_data[i * m->m_stride + j];
                    }
                }
                return copy;
            }

            void CopyContent(vector_t *dest, vector_t *src)
            {
                for (s32 i = 0; i < src->m_N; i++)
                {
                    dest->m_data[i] = src->m_data[i];
                }
            }

            void CopyContent(matrix_t *dest, matrix_t *src)
            {
                for (s32 i = 0; i < src->m_Rows; i++)
                {
                    for (s32 j = 0; j < src->m_Cols; j++)
                    {
                        dest->m_data[i * dest->m_stride + j] = src->m_data[i * src->m_stride + j];
                    }
                }
            }

            matrix_t *Transpose(matrix_t *m)
            {
                matrix_t *transposed = NewMatrix(m->m_Cols, m->m_Rows, nullptr);
                for (s32 i = 0; i < m->m_Rows; i++)
                {
                    for (s32 j = 0; j < m->m_Cols; j++)
                    {
                        transposed->m_data[j * transposed->m_stride + i] = m->m_data[i * m->m_stride + j];
                    }
                }
                return transposed;
            }

            void Free(vector_t *v)
            {
                // TODO free memory
                // free(v->m_data);
                // free(v);
            }

            void Free(matrix_t *m)
            {
                // TODO free memory
                // free(m->m_data);
                // free(m);
            }

        }  // namespace nmath

    }  // namespace nkalman
}  // namespace ncore
