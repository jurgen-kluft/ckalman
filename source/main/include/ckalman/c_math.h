#ifndef __C_KALMAN_MATH_H__
#define __C_KALMAN_MATH_H__
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
            struct matrix_t;

            inline f32 Abs(f32 f) { return (f < 0.0f) ? -f : f; }
            inline s32 Min(s32 a, s32 b) { return (a < b) ? a : b; }

            struct vector_t
            {
                s8   m_N;
                s8   m_Inc;
                f32 *m_data;

                inline s32  Len() const { return m_N; }
                inline f32  AtVec(s32 i) const { return m_data[i * m_Inc]; }
                inline void SetVec(s32 i, f32 value) { m_data[i * m_Inc] = value; }

                void MulVec(memory_t *mem, matrix_t *m, vector_t *v);
                void AddVec(memory_t *mem, vector_t *a, vector_t *b);
                void SubVec(memory_t *mem, vector_t *a, vector_t *b);

                void Copy(vector_t *src);
            };

            struct matrix_t
            {
                s8   m_rows;
                s8   m_cols;
                s8   m_stride;
                f32 *m_data;

                inline void Set(s32 row, s32 col, f32 value) { m_data[row * m_stride + col] = value; }
                inline f32  At(s32 row, s32 col) const { return m_data[row * m_stride + col]; }

                void Product(memory_t *mem, matrix_t *T, matrix_t *P, matrix_t *transposedT);
                void Add(memory_t *mem, matrix_t *a, matrix_t *b);
                void Sub(memory_t *mem, matrix_t *a, matrix_t *b);
                void Mul(memory_t *mem, matrix_t *a, matrix_t *b);
                void Mul2(memory_t *mem, matrix_t *a, matrix_t *b);
                void Inverse(memory_t *mem, matrix_t *m);
                void Transpose(memory_t *mem, matrix_t *m);

                void Copy(matrix_t *src);
            };

            vector_t *NewVector(memory_t *mem, s32 n, f32 *data);
            matrix_t *NewMatrix(memory_t *mem, s32 rows, s32 cols, f32 *data);
            f32      *NewBuffer(memory_t *mem, s32 size);
            vector_t *Duplicate(memory_t *mem, vector_t *v);
            matrix_t *Duplicate(memory_t *mem, matrix_t *m);

        }  // namespace nmath
    }  // namespace nkalman
}  // namespace ncore
#endif  // __C_KALMAN_MATH_H__
