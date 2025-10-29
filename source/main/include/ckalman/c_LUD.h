#ifndef __C_KALMAN_LUD_H__
#define __C_KALMAN_LUD_H__
#include "ccore/c_target.h"
#ifdef USE_PRAGMA_ONCE
    #pragma once
#endif

namespace ncore
{
    namespace nkalman
    {
        struct memory_t;

        namespace nmath
        {
            struct lud_t
            {
                inline lud_t()
                    : m_LU(nullptr)
                    , m_Pivots(nullptr)
                    , m_PivotSign(1)
                    , m_M(0)
                    , m_N(0)
                {
                }

                matrix_t* m_LU;         // User matrix M to decompose
                s8*       m_Pivots;     // Pivot permutation vector
                s8        m_PivotSign;  // Sign of the permutation
                s8        m_M;          // Number of rows
                s8        m_N;          // Number of columns
            };

            void      begin(lud_t* lud, memory_t* mem, matrix_t* M);
            bool      is_nonsingular(lud_t* lud);
            float     determinant(lud_t* lud);
            matrix_t* getL(lud_t* lud, memory_t* mem);
            matrix_t* getU(lud_t* lud, memory_t* mem);
            matrix_t* solve(lud_t* lud, memory_t* mem, matrix_t* B);
            matrix_t* inverse(lud_t* lud, memory_t* mem);

        }  // namespace nmath
    }  // namespace nkalman
}  // namespace ncore

#endif  // __C_KALMAN_LUD_H__
