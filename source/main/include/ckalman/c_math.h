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

            inline float Abs(float f) { return (f < 0.0f) ? -f : f; }
            inline s32   Min(s32 a, s32 b) { return (a < b) ? a : b; }

            struct vector_t
            {
                s32    m_N;
                s32    m_Inc;
                float *m_data;

                inline s32   Len() const { return m_N; }
                inline float AtVec(s32 i) const { return m_data[i * m_Inc]; }
                inline void  SetVec(s32 i, float value) { m_data[i * m_Inc] = value; }

                void MulVec(memory_t *mem, matrix_t *m, vector_t *v);
                void AddVec(memory_t *mem, vector_t *a, vector_t *b);
                void SubVec(memory_t *mem, vector_t *a, vector_t *b);
            };

            struct matrix_t
            {
                s32    m_Rows;
                s32    m_Cols;
                float *m_data;
                s32    m_stride;
                s32    m_capRows;
                s32    m_capCols;

                inline void  Set(s32 row, s32 col, float value) { m_data[row * m_stride + col] = value; }
                inline float At(s32 row, s32 col) const { return m_data[row * m_stride + col]; }

                void Product(memory_t *mem, matrix_t *T, matrix_t *P, matrix_t *transposedT);
                void Add(memory_t *mem, matrix_t *a, matrix_t *b);
                void Sub(memory_t *mem, matrix_t *a, matrix_t *b);
                void Mul(memory_t *mem, matrix_t *a, matrix_t *b);
                void Mul2(memory_t *mem, matrix_t *a, matrix_t *b);
                void Inverse(memory_t *mem, matrix_t *m);
                void Transpose(memory_t *mem, matrix_t *m);
            };

            vector_t *NewVector(memory_t *mem, s32 n, float *data);
            matrix_t *NewMatrix(memory_t *mem, s32 rows, s32 cols, float *data);
            float    *NewBuffer(memory_t *mem, s32 size);
            vector_t *Copy(memory_t *mem, vector_t *v);
            matrix_t *Copy(memory_t *mem, matrix_t *m);
            void      CopyContent(vector_t *dest, vector_t *src);
            void      CopyContent(matrix_t *dest, matrix_t *src);

        }  // namespace nmath
    }  // namespace nkalman
}  // namespace ncore
#endif  // __C_KALMAN_MATH_H__
