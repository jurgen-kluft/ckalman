#ifndef __C_KALMAN_FILTER_RD03D_H__
#define __C_KALMAN_FILTER_RD03D_H__
#include "ccore/c_target.h"
#ifdef USE_PRAGMA_ONCE
#    pragma once
#endif

#include "ccore/c_math.h"
#include "ckalman/c_kalman.h"

namespace ncore
{
    namespace nkalman
    {
        enum config_t
        {
            MAX_TARGETS = 3,
            STATE_DIM   = 4,  // [X_pos, Y_pos, X_vel, Y_vel]
            MEASURE_DIM = 2   // [X_pos, Y_pos]
        };

        struct rd03d_t
        {
            kalman_nd_t<STATE_DIM, MEASURE_DIM> m_roomFilters[MAX_TARGETS];
            bool                                m_targetActive[MAX_TARGETS];
        };

        // dt = 0.05f -> 50ms intervals (20Hz)
        static inline void setup(rd03d_t& rd, f32 dt = 0.05f)
        {
            rd.m_targetActive[0] = false;
            rd.m_targetActive[1] = false;
            rd.m_targetActive[2] = false;

            for (i32 i = 0; i < MAX_TARGETS; i++)
            {
                // --- F MATRIX: State Transitions ---
                // New_X = X + (Vx * dt)
                // New_Y = Y + (Vy * dt)
                rd.m_roomFilters[i].F.data[0][0] = 1.0f;
                rd.m_roomFilters[i].F.data[0][1] = 0.0f;
                rd.m_roomFilters[i].F.data[0][2] = dt;
                rd.m_roomFilters[i].F.data[0][3] = 0.0f;
                rd.m_roomFilters[i].F.data[1][0] = 0.0f;
                rd.m_roomFilters[i].F.data[1][1] = 1.0f;
                rd.m_roomFilters[i].F.data[1][2] = 0.0f;
                rd.m_roomFilters[i].F.data[1][3] = dt;
                rd.m_roomFilters[i].F.data[2][0] = 0.0f;
                rd.m_roomFilters[i].F.data[2][1] = 0.0f;
                rd.m_roomFilters[i].F.data[2][2] = 1.0f;
                rd.m_roomFilters[i].F.data[2][3] = 0.0f;
                rd.m_roomFilters[i].F.data[3][0] = 0.0f;
                rd.m_roomFilters[i].F.data[3][1] = 0.0f;
                rd.m_roomFilters[i].F.data[3][2] = 0.0f;
                rd.m_roomFilters[i].F.data[3][3] = 1.0f;

                // --- H MATRIX: Map States to Measurements ---
                // We measure state index 0 (X position) and index 1 (Y position)
                rd.m_roomFilters[i].H.data[0][0] = 1.0f;
                rd.m_roomFilters[i].H.data[0][1] = 0.0f;
                rd.m_roomFilters[i].H.data[0][2] = 0.0f;
                rd.m_roomFilters[i].H.data[0][3] = 0.0f;
                rd.m_roomFilters[i].H.data[1][0] = 0.0f;
                rd.m_roomFilters[i].H.data[1][1] = 1.0f;
                rd.m_roomFilters[i].H.data[1][2] = 0.0f;
                rd.m_roomFilters[i].H.data[1][3] = 0.0f;

                // --- Q MATRIX: Process Noise (Target Acceleration dynamics) ---
                rd.m_roomFilters[i].Q.setIdentity();
                rd.m_roomFilters[i].Q.data[0][0] = 0.01f;  // Position X variance
                rd.m_roomFilters[i].Q.data[1][1] = 0.01f;  // Position Y variance
                rd.m_roomFilters[i].Q.data[2][2] = 0.10f;  // Velocity X variance
                rd.m_roomFilters[i].Q.data[3][3] = 0.10f;  // Velocity Y variance

                // --- R MATRIX: Sensor Measurement Noise Covariance ---
                // mmWave angle tracking is generally noisier than distance tracking
                rd.m_roomFilters[i].R.data[0][0] = 0.25f;  // X measurement variance (influenced heavily by angle)
                rd.m_roomFilters[i].R.data[1][1] = 0.15f;  // Y measurement variance
            }
        }

        struct target_t
        {
            u8   m_id;        // 1-based ID from sensor (1, 2, or 3)
            bool m_detected;  // Whether the target is currently detected in this frame
            f32  m_distance;  // in meters
            f32  m_angle;     // in degrees (Ensure your library converts raw bytes to degrees)
        };

        void processFrame(rd03d_t& rd, const target_t targets[], i32 count);

    }  // namespace nkalman
}  // namespace ncore

#endif  // __C_KALMAN_FILTER_RD03D_H__
