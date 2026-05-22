#ifndef __C_KALMAN_FILTER_H__
#define __C_KALMAN_FILTER_H__
#include "ccore/c_target.h"
#ifdef USE_PRAGMA_ONCE
#    pragma once
#endif

namespace ncore
{
    namespace nkalman
    {
        template <i32 ROWS, i32 COLS>
        struct matrix_t
        {
            f32 data[ROWS][COLS];

            void setIdentity()
            {
                for (i32 i = 0; i < ROWS; ++i)
                {
                    for (i32 j = 0; j < COLS; ++j)
                    {
                        data[i][j] = (i == j) ? 1.0f : 0.0f;
                    }
                }
            }

            void clear()
            {
                for (i32 i = 0; i < ROWS; ++i)
                    for (i32 j = 0; j < COLS; ++j)
                        data[i][j] = 0.0f;
            }

            // matrix_t Addition: Out = This + Other
            void add(const matrix_t<ROWS, COLS>& other, matrix_t<ROWS, COLS>& out) const
            {
                for (i32 i = 0; i < ROWS; ++i)
                {
                    for (i32 j = 0; j < COLS; ++j)
                    {
                        out.data[i][j] = this->data[i][j] + other.data[i][j];
                    }
                }
            }

            // matrix_t Subtraction: Out = This - Other
            void subtract(const matrix_t<ROWS, COLS>& other, matrix_t<ROWS, COLS>& out) const
            {
                for (i32 i = 0; i < ROWS; ++i)
                {
                    for (i32 j = 0; j < COLS; ++j)
                    {
                        out.data[i][j] = this->data[i][j] - other.data[i][j];
                    }
                }
            }

            // matrix_t Multiplication: Out (ROWS x OTHERCOLS) = This (ROWS x COLS) * Other (COLS x OTHERCOLS)
            template <i32 OTHERCOLS>
            void multiply(const matrix_t<COLS, OTHERCOLS>& other, matrix_t<ROWS, OTHERCOLS>& out) const
            {
                out.clear();
                for (i32 i = 0; i < ROWS; ++i)
                {
                    for (i32 j = 0; j < OTHERCOLS; ++j)
                    {
                        for (i32 k = 0; k < COLS; ++k)
                        {
                            out.data[i][j] += this->data[i][k] * other.data[k][j];
                        }
                    }
                }
            }

            // Transpose: Out (COLS x ROWS) = This^T (ROWS x COLS)
            void transpose(matrix_t<COLS, ROWS>& out) const
            {
                for (i32 i = 0; i < ROWS; ++i)
                {
                    for (i32 j = 0; j < COLS; ++j)
                    {
                        out.data[j][i] = this->data[i][j];
                    }
                }
            }

            // Deterministic Stack Inversion (Supports up to 3x3 matrices)
            void invert(matrix_t<ROWS, COLS>& out) const
            {
                static_assert(ROWS == COLS, "Inversion requires a square matrix!");
                out.clear();

                if (ROWS == 1)
                {
                    out.data[0][0] = (data[0][0] != 0.0f) ? (1.0f / data[0][0]) : 1.0f;
                }
                else if (ROWS == 2)
                {
                    f32 det        = data[0][0] * data[1][1] - data[0][1] * data[1][0];
                    f32 invDet     = (det != 0.0f) ? (1.0f / det) : 1.0f;
                    out.data[0][0] = data[1][1] * invDet;
                    out.data[0][1] = -data[0][1] * invDet;
                    out.data[1][0] = -data[1][0] * invDet;
                    out.data[1][1] = data[0][0] * invDet;
                }
                else if (ROWS == 3)
                {
                    f32 det    = data[0][0] * (data[1][1] * data[2][2] - data[1][2] * data[2][1]) - data[0][1] * (data[1][0] * data[2][2] - data[1][2] * data[2][0]) + data[0][2] * (data[1][0] * data[2][1] - data[1][1] * data[2][0]);
                    f32 invDet = (det != 0.0f) ? (1.0f / det) : 1.0f;

                    out.data[0][0] = (data[1][1] * data[2][2] - data[1][2] * data[2][1]) * invDet;
                    out.data[0][1] = (data[0][2] * data[2][1] - data[0][1] * data[2][2]) * invDet;
                    out.data[0][2] = (data[0][1] * data[1][2] - data[0][2] * data[1][1]) * invDet;
                    out.data[1][0] = (data[1][2] * data[2][0] - data[1][0] * data[2][2]) * invDet;
                    out.data[1][1] = (data[0][0] * data[2][2] - data[0][2] * data[2][0]) * invDet;
                    out.data[1][2] = (data[0][2] * data[1][0] - data[0][0] * data[1][2]) * invDet;
                    out.data[2][0] = (data[1][0] * data[2][1] - data[1][1] * data[2][0]) * invDet;
                    out.data[2][1] = (data[0][1] * data[2][0] - data[0][0] * data[2][1]) * invDet;
                    out.data[2][2] = (data[0][0] * data[1][1] - data[0][1] * data[1][0]) * invDet;
                }
            }
        };

        // Simple 1D Kalman Filter Implementation
        struct kalman_1D_t
        {
            f32 q;  // process noise covariance
            f32 r;  // measurement noise covariance
            f32 p;  // estimation error covariance
            f32 x;  // state estimate
            f32 k;  // Kalman gain
        };

        static inline void initialize(kalman_1D_t& kf, f32 processNoise = 0.01f, f32 measurementNoise = 0.1f, f32 estimationError = 1.0f)
        {
            kf.q = (processNoise);
            kf.r = (measurementNoise);
            kf.p = (estimationError);
            kf.x = (0.0f);
            kf.k = (0.0f);
        }

        static inline void begin(kalman_1D_t& kf, f32 initialState, f32 initialError = 1.0f)
        {
            kf.x = initialState;
            kf.p = initialError;
        }

        static inline f32 update(kalman_1D_t& kf, f32 measurement)
        {
            kf.p = kf.p + kf.q;
            kf.k = kf.p / (kf.p + kf.r);
            kf.x = kf.x + kf.k * (measurement - kf.x);
            kf.p = (1.0f - kf.k) * kf.p;
            return kf.x;
        }

        static inline f32 getState(const kalman_1D_t& kf) { return kf.x; }

        // ============================================================================
        // MULTI-VARIABLE MATRIX KALMAN FILTER (Refactored to use matrix_t functions)
        // ============================================================================

        template <i32 N, i32 M>
        struct kalman_nd_t
        {
            matrix_t<N, N> F;  // State transition
            matrix_t<N, M> K;  // Kalman Gain
            matrix_t<N, N> P;  // Estimate error covariance
            matrix_t<N, N> Q;  // Process noise covariance
            matrix_t<M, M> R;  // Measurement noise covariance
            matrix_t<M, N> H;  // Measurement mapping
            matrix_t<N, 1> x;  // State vector (represented as Nx1 matrix_t)
        };

        template <i32 N, i32 M>
        static inline void initialize(kalman_nd_t<N, M>& kf)
        {
            kf.P.setIdentity();
            kf.Q.setIdentity();
            kf.R.setIdentity();
            kf.x.clear();
        }

        template <i32 N, i32 M>
        static inline void begin(kalman_nd_t<N, M>& kf, const f32 initial_states[N], f32 initial_uncertainty = 1.0f)
        {
            // 1. Set the initial tracking states
            for (int i = 0; i < N; ++i)
            {
                kf.x.data[i][0] = initial_states[i];
            }

            // 2. Reset the error covariance matrix to an identity layout scaled by uncertainty
            kf.P.setIdentity();
            for (int i = 0; i < N; ++i)
            {
                kf.P.data[i][i] = initial_uncertainty;
            }
        }

        template <i32 N, i32 M>
        static inline void update(kalman_nd_t<N, M>& kf, const f32 measurement[M])
        {
            // Wrap raw array into our reusable matrix_t object for computations
            matrix_t<M, 1> z;
            for (int i = 0; i < M; ++i)
                z.data[i][0] = measurement[i];

            // --- 1. PREDICT PHASE ---
            // x_pred = F * x
            matrix_t<N, 1> x_pred;
            kf.F.multiply(kf.x, x_pred);

            // P_pred = F * P * F^T + Q
            matrix_t<N, N> F_trans;
            kf.F.transpose(F_trans);
            matrix_t<N, N> FP;
            kf.F.multiply(kf.P, FP);
            matrix_t<N, N> FP_FT;
            FP.multiply(F_trans, FP_FT);
            matrix_t<N, N> P_pred;
            FP_FT.add(kf.Q, P_pred);

            // --- 2. MEASUREMENT UPDATE PHASE ---
            // y = z - H * x_pred
            matrix_t<M, 1> H_xpred;
            kf.H.multiply(x_pred, H_xpred);
            matrix_t<M, 1> y;
            z.subtract(H_xpred, y);

            // S = H * P_pred * H^T + R
            matrix_t<N, M> H_trans;
            kf.H.transpose(H_trans);
            matrix_t<M, N> HP;
            kf.H.multiply(P_pred, HP);
            matrix_t<M, M> HP_HT;
            HP.multiply(H_trans, HP_HT);
            matrix_t<M, M> S;
            HP_HT.add(kf.R, S);

            // K = P_pred * H^T * S^-1
            matrix_t<M, M> S_inv;
            S.invert(S_inv);
            matrix_t<N, M> P_HT;
            P_pred.multiply(H_trans, P_HT);
            P_HT.multiply(S_inv, kf.K);

            // x = x_pred + K * y
            matrix_t<N, 1> Ky;
            kf.K.multiply(y, Ky);
            x_pred.add(Ky, kf.x);

            // P = (I - K * H) * P_pred
            matrix_t<N, N> I;
            I.setIdentity();
            matrix_t<N, N> KH;
            kf.K.multiply(kf.H, KH);
            matrix_t<N, N> I_KH;
            I.subtract(KH, I_KH);
            I_KH.multiply(P_pred, kf.P);
        }

    }  // namespace nkalman
}  // namespace ncore
#endif  // __C_KALMAN_FILTER_H__
