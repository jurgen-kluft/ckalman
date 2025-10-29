#include "ckalman/c_memory.h"
#include "ckalman/c_math.h"
#include "ckalman/c_lud.h"

#include "ccore/c_debug.h"
#include "ccore/c_memory.h"

namespace ncore
{
    namespace nkalman
    {
        namespace nmath
        {
            struct lud_t
            {
                s32       m_N;
                s32       m_M;
                s32       m_PivotSign;
                s32*      m_Pivots;
                matrix_t* m_LU;
                memory_t* m_Memory;

                void Init(matrix_t* M);
            };

            lud_t* NewLUD(memory_t* mem)
            {
                lud_t* lud    = (lud_t*)mem->AllocZeroMemory(sizeof(lud_t), alignof(lud_t));
                lud->m_Memory = mem;
                return lud;
            }

            void lud_t::Init(matrix_t* M)
            {   
                m_LU        = M;
                m_M         = M->m_Rows;
                m_N         = M->m_Cols;
                m_PivotSign = 1;
                m_Pivots    = (s32*)m_Memory->AllocZeroMemory(sizeof(s32) * m_M, alignof(s32));
                for (s32 i = 0; i < m_M; i++)
                {
                    m_Pivots[i] = i;
                }

                float* LUcolj = m_Memory->AllocFloatArray((u32)m_M);

                // Outer loop.
                for (s32 j = 0; j < m_N; j++)
                {
                    // Make a copy of the j-th column to localize references.
                    for (s32 i = 0; i < m_M; i++)
                    {
                        LUcolj[i] = m_LU->At(i, j);
                    }

                    // Apply previous transformations.
                    for (s32 i = 0; i < m_M; i++)
                    {
                        float* LUrowi = &m_LU->m_data[i * m_LU->m_stride];

                        // Most of the time is spent in the following dot product.
                        s32   kmax = nmath::Min(i, j);
                        float s    = 0.0f;
                        for (s32 k = 0; k < kmax; k++)
                        {
                            s += LUrowi[k] * LUcolj[k];
                        }
                        LUrowi[j] = LUcolj[i] -= s;
                    }

                    // Find pivot and exchange if necessary.
                    s32 p = j;
                    for (s32 i = j + 1; i < m_M; i++)
                    {
                        if (nmath::Abs(LUcolj[i]) > nmath::Abs(LUcolj[p]))
                        {
                            p = i;
                        }
                    }

                    if (p != j)
                    {
                        for (s32 k = 0; k < m_N; k++)
                        {
                            float t = m_LU->At(p, k);
                            m_LU->Set(p, k, m_LU->At(j, k));
                            m_LU->Set(j, k, t);
                        }
                        s32 k       = m_Pivots[p];
                        m_Pivots[p] = m_Pivots[j];
                        m_Pivots[j] = k;
                        m_PivotSign = -m_PivotSign;
                    }

                    // Compute multipliers.
                    if (j < m_M && m_LU->At(j, j) != 0.0f)
                    {
                        for (s32 i = j + 1; i < m_M; i++)
                        {
                            m_LU->Set(i, j, m_LU->At(i, j) / m_LU->At(j, j));
                        }
                    }
                }
            }

            bool IsNonsingular(lud_t* lud)
            {
                for (s32 j = 0; j < lud->m_N; j++)
                {
                    if (lud->m_LU->At(j, j) == 0.0f)
                    {
                        return false;
                    }
                }
                return true;
            }

            void getPivotVector(lud_t* lud, s32*& pivots_out, s32& pivots_len)
            {
                pivots_out = (s32*)lud->m_Memory->AllocMemory(sizeof(s32) * lud->m_M, alignof(s32));
                for (s32 i = 0; i < lud->m_M; i++)
                {
                    pivots_out[i] = lud->m_Pivots[i];
                }
                pivots_len = lud->m_M;
            }

            float Determinant(lud_t* lud)
            {
                ASSERTS(lud->m_N == lud->m_M, "Matrix must be square to compute determinant");
                float det = (float)lud->m_PivotSign;
                for (s32 j = 0; j < lud->m_N; j++)
                {
                    det *= lud->m_LU->At(j, j);
                }
                return det;
            }

            matrix_t* IdentityMatrix(lud_t* lud, s32 size)
            {
                matrix_t* I = lud->m_Memory->AllocZeroMatrix(size, size);
                for (s32 i = 0; i < size; i++)
                {
                    I->Set(i, i, 1.0f);
                }
                return I;
            }

            matrix_t* SubMatrix(lud_t* lud, matrix_t* A, s32* rArray, s32 rLen, s32 j0, s32 j1)
            {
                matrix_t* B = lud->m_Memory->AllocZeroMatrix(rLen, j1 - j0 + 1);
                for (s32 i = 0; i < rLen; i++)
                {
                    for (s32 j = j0; j <= j1; j++)
                        B->Set(i, j - j0, A->At(rArray[i], j));
                }
                return B;
            }

            matrix_t* GetL(lud_t* lud)
            {
                matrix_t* L = lud->m_Memory->AllocMatrix(lud->m_M, lud->m_N);
                for (s32 i = 0; i < lud->m_M; i++)
                {
                    for (s32 j = 0; j < lud->m_N; j++)
                    {
                        if (i > j)
                        {
                            L->Set(i, j, lud->m_LU->At(i, j));
                        }
                        else if (i == j)
                        {
                            L->Set(i, j, 1.0f);
                        }
                        else
                        {
                            L->Set(i, j, 0.0f);
                        }
                    }
                }
                return L;
            }

            matrix_t* GetU(lud_t* lud)
            {
                matrix_t* U = lud->m_Memory->AllocMatrix(lud->m_M, lud->m_N);
                for (s32 i = 0; i < lud->m_M; i++)
                {
                    for (s32 j = 0; j < lud->m_N; j++)
                    {
                        if (i <= j)
                        {
                            U->Set(i, j, lud->m_LU->At(i, j));
                        }
                        else
                        {
                            U->Set(i, j, 0.0f);
                        }
                    }
                }
                return U;
            }

            // Solve A*X = B
            // @param  B   A Matrix with as many rows as A and any number of columns.
            // @return     X so that L*U*X = B(piv,:)
            // @required   Matrix row dimensions must agree.
            // @required   Matrix is not singular.
            matrix_t* Solve(lud_t* lud, matrix_t* B)
            {
                ASSERTS(B->m_Rows == lud->m_M, "Matrix row dimensions must agree.");
                ASSERTS(IsNonsingular(lud), "Matrix is singular.");

                // Copy right hand side with pivoting
                const s32 nx = B->m_Cols;
                matrix_t* X  = SubMatrix(lud, B, lud->m_Pivots, lud->m_N, 0, nx - 1);

                // Solve L*Y = B(piv,:)
                for (s32 k = 0; k < lud->m_N; k++)
                {
                    for (s32 i = k + 1; i < lud->m_N; i++)
                    {
                        for (s32 j = 0; j < nx; j++)
                        {
                            X->Set(i, j, X->At(i, j) - X->At(k, j) * lud->m_LU->At(i, k));
                        }
                    }
                }
                // Solve U*X = Y;
                for (s32 k = lud->m_N - 1; k >= 0; k--)
                {
                    for (s32 j = 0; j < nx; j++)
                    {
                        X->Set(k, j, X->At(k, j) / lud->m_LU->At(k, k));
                    }
                    for (s32 i = 0; i < k; i++)
                    {
                        for (s32 j = 0; j < nx; j++)
                        {
                            X->Set(i, j, X->At(i, j) - X->At(k, j) * lud->m_LU->At(i, k));
                        }
                    }
                }

                return X;
            }

            matrix_t* Inverse(lud_t* lud)
            {
                matrix_t* B      = IdentityMatrix(lud, lud->m_M);
                matrix_t* result = Solve(lud, B);
                return result;
            }

        }  // namespace nmath
    }  // namespace nkalman
}  // namespace ncore
