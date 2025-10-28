#ifndef __C_KALMAN_MEMORY_H__
#define __C_KALMAN_MEMORY_H__
#include "ccore/c_target.h"
#ifdef USE_PRAGMA_ONCE
#    pragma once
#endif

namespace ncore
{
    namespace nkalman
    {
        namespace nmath
        {
            struct vector_t;
            struct matrix_t;
        }  // namespace nmath

        struct memory_t
        {
            enum
            {
                MEMORY_SCOPE_MAX = 8
            };

            void *m_memory_base;
            u32   m_memory_size;
            u8   *m_memory_current;
            u8   *m_scopes[MEMORY_SCOPE_MAX];
            s32   m_scope_index;

            void            *AllocMemory(u32 size, u32 alignment);
            void            *AllocZeroMemory(u32 size, u32 alignment);
            nmath::vector_t *AllocVector(s32 n);
            nmath::vector_t *AllocZeroVector(s32 n);
            nmath::matrix_t *AllocMatrix(s32 rows, s32 cols);
            nmath::matrix_t *AllocZeroMatrix(s32 rows, s32 cols);
            f32             *AllocFloatArray(u32 n);

            void PushScope();
            void PopScope();

            void Reset();
        };

    }  // namespace nkalman
}  // namespace ncore
#endif  // __C_KALMAN_MEMORY_H__
