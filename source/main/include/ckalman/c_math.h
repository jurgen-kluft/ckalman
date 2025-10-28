#ifndef __C_KALMAN_MATH_H__
#define __C_KALMAN_MATH_H__
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
            struct matrix_t;

            struct vector_t
            {
                s32    m_N;
                s32    m_Inc;
                float *m_data;

                inline s32   Len() const { return m_N; }
                inline float AtVec(s32 i) const { return m_data[i * m_Inc]; }
                inline void  SetVec(s32 i, float value) { m_data[i * m_Inc] = value; }

                void MulVec(matrix_t *m, vector_t *v);
                void AddVec(vector_t *a, vector_t *b);
                void SubVec(vector_t *a, vector_t *b);
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

                void Product(matrix_t *T, matrix_t *P, matrix_t *transposedT);
                void Add(matrix_t *a, matrix_t *b);
                void Sub(matrix_t *a, matrix_t *b);
                void Mul(matrix_t *a, matrix_t *b);
                void Mul2(matrix_t *a, matrix_t *b);
                void Inverse(matrix_t *m);
            };

            vector_t *NewVector(s32 n, float *data);
            matrix_t *NewMatrix(s32 rows, s32 cols, float *data);
            float    *NewBuffer(s32 size);
            vector_t *Copy(vector_t *v);
            matrix_t *Copy(matrix_t *m);
            void      CopyContent(vector_t *dest, vector_t *src);
            void      CopyContent(matrix_t *dest, matrix_t *src);
            matrix_t *Transpose(matrix_t *m);
            void      Free(vector_t *v);
            void      Free(matrix_t *m);
            void      FreeBuffer(float *buf);

        }  // namespace nmath
    }  // namespace nkalman
}  // namespace ncore
#endif  // __C_KALMAN_MATH_H__
