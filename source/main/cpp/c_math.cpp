#include "ckalman/c_math.h"
#include "ckalman/c_memory.h"

namespace ncore
{
    namespace nkalman
    {
        namespace nmath
        {
            void vector_t::MulVec(memory_t *mem, matrix_t *m, vector_t *v)
            {
                mem->PushScope();
                {
                    vector_t *result = NewVector(mem, m->m_Rows, nullptr);
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
                }
                mem->PopScope();
            }

            void vector_t::AddVec(memory_t *mem, vector_t *a, vector_t *b)
            {
                mem->PushScope();
                {
                    const s32 ar     = a->Len();
                    vector_t *result = NewVector(mem, ar, nullptr);
                    for (s32 i = 0; i < ar; i++)
                        result->SetVec(i, a->AtVec(i) + b->AtVec(i));
                    for (s32 i = 0; i < ar; i++)
                        SetVec(i, result->AtVec(i));
                }
                mem->PopScope();
            }

            void vector_t::SubVec(memory_t *mem, vector_t *a, vector_t *b)
            {
                mem->PushScope();
                {
                    const s32 ar     = a->Len();
                    vector_t *result = NewVector(mem, ar, nullptr);
                    for (s32 i = 0; i < ar; i++)
                        result->SetVec(i, a->AtVec(i) - b->AtVec(i));
                    for (s32 i = 0; i < ar; i++)
                        SetVec(i, result->AtVec(i));
                }
                mem->PopScope();
            }

            void matrix_t::Product(memory_t *mem, matrix_t *T, matrix_t *P, matrix_t *transposedT)
            {
                mem->PushScope();
                {
                    matrix_t *result = NewMatrix(mem, T->m_Rows, P->m_Cols, nullptr);
                    result->Mul2(mem, T, P);
                    Mul2(mem, result, transposedT);
                }
                mem->PopScope();
            }

            void matrix_t::Add(memory_t *mem, matrix_t *a, matrix_t *b)
            {
                mem->PushScope();
                {
                    matrix_t *result = NewMatrix(mem, a->m_Rows, a->m_Cols, nullptr);
                    for (s32 i = 0; i < a->m_Rows; i++)
                    {
                        for (s32 j = 0; j < a->m_Cols; j++)
                            result->Set(i, j, a->At(i, j) + b->At(i, j));
                    }
                    CopyContent(this, result);
                }
                mem->PopScope();
            }

            void matrix_t::Sub(memory_t *mem, matrix_t *a, matrix_t *b)
            {
                mem->PushScope();
                {
                    matrix_t *result = NewMatrix(mem, a->m_Rows, a->m_Cols, nullptr);
                    for (s32 i = 0; i < a->m_Rows; i++)
                    {
                        for (s32 j = 0; j < a->m_Cols; j++)
                            result->Set(i, j, a->At(i, j) - b->At(i, j));
                    }
                    CopyContent(this, result);
                }
                mem->PopScope();
            }

            void matrix_t::Mul(memory_t *mem, matrix_t *a, matrix_t *b)
            {
                mem->PushScope();
                {
                    const s32 ar = a->m_Rows;
                    const s32 ac = a->m_Cols;

                    matrix_t *result = NewMatrix(mem, ar, ac, nullptr);
                    float    *row    = NewBuffer(mem, ac);

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
                }
                mem->PopScope();
            }

            // We know for sure that 'this' is not 'a' or 'b'
            void matrix_t::Mul2(memory_t *mem, matrix_t *a, matrix_t *b)
            {
                mem->PushScope();
                {
                    const s32 ar  = a->m_Rows;
                    const s32 ac  = a->m_Cols;
                    float    *row = NewBuffer(mem, ac);
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
                }
                mem->PopScope();
            }

            void matrix_t::Inverse(memory_t *mem, matrix_t *m)
            {
                // TODO implement matrix inversion
            }

            void matrix_t::Transpose(memory_t *mem, matrix_t *m)
            {
                mem->PushScope();
                {
                    matrix_t *transposed = NewMatrix(mem, m->m_Cols, m->m_Rows, nullptr);
                    for (s32 i = 0; i < m->m_Rows; i++)
                    {
                        for (s32 j = 0; j < m->m_Cols; j++)
                        {
                            transposed->m_data[j * transposed->m_stride + i] = m->m_data[i * m->m_stride + j];
                        }
                    }
                    CopyContent(this, transposed);
                }
                mem->PopScope();
            }

            vector_t *NewVector(memory_t *mem, s32 n, float *data)
            {
                return mem->AllocVector(n);
            }

            matrix_t *NewMatrix(memory_t *mem, s32 rows, s32 cols, float *data)
            {
                return mem->AllocMatrix(rows, cols);
            }

            float* NewBuffer(memory_t *mem, s32 size)
            {
                return mem->AllocFloatArray((u32)size);
            }

            vector_t *Copy(memory_t *mem, vector_t *v)
            {
                vector_t *copy = NewVector(mem, v->m_N, nullptr);
                for (s32 i = 0; i < v->m_N; i++)
                {
                    copy->m_data[i] = v->m_data[i];
                }
                return copy;
            }

            matrix_t *Copy(memory_t *mem, matrix_t *m)
            {
                matrix_t *copy = NewMatrix(mem, m->m_Rows, m->m_Cols, nullptr);
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
        }  // namespace nmath

    }  // namespace nkalman
}  // namespace ncore
