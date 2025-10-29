#include "ckalman/c_lud.h"
#include "ckalman/c_math.h"
#include "ckalman/c_memory.h"

#include "ccore/c_debug.h"
#include "ccore/c_memory.h"

namespace ncore
{
    namespace nkalman
    {
        namespace nmath
        {
            void begin(lud_t* lud, memory_t* mem, matrix_t* M)
            {
                lud->m_LU        = M;
                lud->m_M         = M->m_rows;
                lud->m_N         = M->m_cols;
                lud->m_PivotSign = 1;
                lud->m_Pivots    = (s8*)mem->AllocMemory(sizeof(s8) * lud->m_M, 4);
                for (s32 i = 0; i < lud->m_M; i++)
                    lud->m_Pivots[i] = i;

                mem->PushScope();
                float* LUcolj = mem->AllocFloatArray((u32)lud->m_M);

                // Outer loop.
                for (s32 j = 0; j < lud->m_N; j++)
                {
                    // Make a copy of the j-th column to localize references.
                    for (s32 i = 0; i < lud->m_M; i++)
                    {
                        LUcolj[i] = lud->m_LU->At(i, j);
                    }

                    // Apply previous transformations.
                    for (s32 i = 0; i < lud->m_M; i++)
                    {
                        float* LUrowi = &lud->m_LU->m_data[i * lud->m_LU->m_stride];

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
                    for (s32 i = j + 1; i < lud->m_M; i++)
                    {
                        if (nmath::Abs(LUcolj[i]) > nmath::Abs(LUcolj[p]))
                        {
                            p = i;
                        }
                    }

                    mem->PopScope();

                    if (p != j)
                    {
                        for (s32 k = 0; k < lud->m_N; k++)
                        {
                            float t = lud->m_LU->At(p, k);
                            lud->m_LU->Set(p, k, lud->m_LU->At(j, k));
                            lud->m_LU->Set(j, k, t);
                        }
                        s32 k       = lud->m_Pivots[p];
                        lud->m_Pivots[p] = lud->m_Pivots[j];
                        lud->m_Pivots[j] = k;
                        lud->m_PivotSign = -lud->m_PivotSign;
                    }

                    // Compute multipliers.
                    if (j < lud->m_M && lud->m_LU->At(j, j) != 0.0f)
                    {
                        for (s32 i = j + 1; i < lud->m_M; i++)
                        {
                            lud->m_LU->Set(i, j, lud->m_LU->At(i, j) / lud->m_LU->At(j, j));
                        }
                    }
                }
            }

            bool is_nonsingular(lud_t* lud)
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

            float determinant(lud_t* lud)
            {
                ASSERTS(lud->m_N == lud->m_M, "Matrix must be square to compute determinant");
                float det = (float)lud->m_PivotSign;
                for (s32 j = 0; j < lud->m_N; j++)
                {
                    det *= lud->m_LU->At(j, j);
                }
                return det;
            }

            matrix_t* identity_matrix(lud_t* lud, memory_t* mem, s32 size)
            {
                matrix_t* I = mem->AllocZeroMatrix(size, size);
                for (s32 i = 0; i < size; i++)
                {
                    I->Set(i, i, 1.0f);
                }
                return I;
            }

            matrix_t* sub_matrix(lud_t* lud, memory_t* mem, matrix_t* A, s8* rArray, s32 rLen, s32 j0, s32 j1)
            {
                matrix_t* B = mem->AllocZeroMatrix(rLen, j1 - j0 + 1);
                for (s32 i = 0; i < rLen; i++)
                {
                    for (s32 j = j0; j <= j1; j++)
                        B->Set(i, j - j0, A->At(rArray[i], j));
                }
                return B;
            }

            matrix_t* getL(lud_t* lud, memory_t* mem)
            {
                matrix_t* L = mem->AllocMatrix(lud->m_M, lud->m_N);
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

            matrix_t* getU(lud_t* lud, memory_t* mem)
            {
                matrix_t* U = mem->AllocMatrix(lud->m_M, lud->m_N);
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
            // @required   Matrix row dimensions must agree
            // @required   Matrix is not singular
            matrix_t* solve(lud_t* lud, memory_t* mem, matrix_t* B)
            {
                ASSERTS(B->m_rows == lud->m_M, "Matrix row dimensions must agree.");
                ASSERTS(is_nonsingular(lud), "Matrix is singular.");

                // Copy right hand side with pivoting
                const s32 nx = B->m_cols;
                matrix_t* X  = sub_matrix(lud, mem, B, lud->m_Pivots, lud->m_N, 0, nx - 1);

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

            matrix_t* inverse(lud_t* lud, memory_t* mem)
            {
                matrix_t* B      = identity_matrix(lud, mem, lud->m_M);
                matrix_t* result = solve(lud, mem, B);
                return result;
            }

        }  // namespace nmath
    }  // namespace nkalman
}  // namespace ncore
