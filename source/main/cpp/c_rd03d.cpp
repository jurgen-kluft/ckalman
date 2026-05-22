#include "ckalman/c_kalman.h"
#include "ckalman/c_rd03d.h"

#include "ccore/c_debug.h"
#include "ccore/c_memory.h"
#include "ccore/c_math.h"

#include <cmath>

namespace ncore
{
    namespace nkalman
    {
        void processFrame(rd03d_t& rd, const target_t targets[], i32 count)
        {
            bool seenThisFrame[MAX_TARGETS] = {false, false, false};

            for (i32 i = 0; i < count; i++)
            {
                const i32 idx = targets[i].m_id - 1;
                if (idx < 0 || idx >= MAX_TARGETS)
                    continue;

                if (targets[i].m_detected)
                {
                    seenThisFrame[idx] = true;

                    // 1. Convert Polar (Distance, Angle) to Cartesian (X, Y)
                    const f32 rad  = targets[i].m_angle * (math::PI / 180.0f);  // Convert degrees to radians
                    const f32 posX = targets[i].m_distance * sin(rad);
                    const f32 posY = targets[i].m_distance * cos(rad);

                    // 2. Gating: New Target Arrival
                    if (!rd.m_targetActive[idx])
                    {
                        rd.m_targetActive[idx] = true;

                        // Initial State: Position X, Position Y, Velocity X (0), Velocity Y (0)
                        f32 initialStates[STATE_DIM] = {posX, posY, 0.0f, 0.0f};
                        begin(rd.m_roomFilters[idx], initialStates, 5.0f);  // Higher initial uncertainty for new targets

                        // Serial.print("🎯 Target ");
                        // Serial.print(targets[i].m_id);
                        // Serial.println(" spawned in room layout.");
                    }

                    // 3. Filter Execution
                    const f32 measurement[MEASURE_DIM] = {posX, posY};
                    update(rd.m_roomFilters[idx], measurement);

                    // 4. Extract Filtered 2D Trajectory Metrics
                    f32 cleanX   = rd.m_roomFilters[idx].x.data[0][0];
                    f32 cleanY   = rd.m_roomFilters[idx].x.data[1][0];
                    f32 speedX   = rd.m_roomFilters[idx].x.data[2][0];
                    f32 speedY   = rd.m_roomFilters[idx].x.data[3][0];
                    f32 netSpeed = sqrt(speedX * speedX + speedY * speedY);

                    // Output data stream for drawing paths or triggering rules
                    // Serial.print("[ID ");
                    // Serial.print(targets[i].m_id);
                    // Serial.print("] X: ");
                    // Serial.print(cleanX);
                    // Serial.print("m | Y: ");
                    // Serial.print(cleanY);
                    // Serial.print("m | Speed: ");
                    // Serial.print(netSpeed);
                    // Serial.println(" m/s");
                }
            }

            // 5. Gating: Target Dropout
            for (i32 i = 0; i < MAX_TARGETS; i++)
            {
                if (rd.m_targetActive[i] && !seenThisFrame[i])
                {
                    rd.m_targetActive[i] = false;
                    // Serial.print("❌ Target ");
                    // Serial.print(i + 1);
                    // Serial.println(" dropped from room.");
                }
            }
        }

    }  // namespace nkalman
}  // namespace ncore
