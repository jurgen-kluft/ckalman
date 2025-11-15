#include "ccore/c_allocator.h"
#include "ccore/c_math.h"
#include "ccore/c_printf.h"
#include "ccore/c_random.h"

#include "ckalman/c_kalman.h"
#include "ckalman/c_models.h"

#include "cunittest/cunittest.h"

#include <cmath>

using namespace ncore;

UNITTEST_SUITE_BEGIN(test_kalman)
{
    UNITTEST_FIXTURE(filter)
    {
        UNITTEST_ALLOCATOR;

        // Generate noisy observations around a linear trend
        void generateValues(int n, f32*& noisy, f32*& truth, f32 observationVariance)
        {
            f32 xNoisy = 0.0f;
            f32 xTruth = 0.0f;

            noisy = g_allocate_array_and_clear<f32>(Allocator, n);
            truth = g_allocate_array_and_clear<f32>(Allocator, n);

            rand_t rnd;
            for (int i = 0; i < n; i++)
            {
                f32 delta = static_cast<f32>(i) * 0.001f;
                // generate Gaussian noise using Box-Muller transform
                f32 u1    = g_random_f32(&rnd);
                f32 u2    = g_random_f32(&rnd);
                f32 noise = observationVariance * sqrtf(-2.0f * logf(u1)) * cosf(2.0f * 3.14159265f * u2);

                xTruth += delta;
                xNoisy = xTruth + noise;

                noisy[i] = xNoisy;
                truth[i] = xTruth;
            }

            // Same values as in Go to be able to debug
            noisy[0]  = -1.733117;
            noisy[1]  = 1.443965;
            noisy[2]  = 1.716114;
            noisy[3]  = -2.683945;
            noisy[4]  = 1.039937;
            noisy[5]  = 0.273394;
            noisy[6]  = -0.254378;
            noisy[7]  = 2.333305;
            noisy[8]  = -1.696381;
            noisy[9]  = 1.924037;
            noisy[10] = -0.735044;
            noisy[11] = -0.457046;
            noisy[12] = 1.273879;
            noisy[13] = 0.911561;
            noisy[14] = 1.333793;
            noisy[15] = 0.721973;
            noisy[16] = -0.328323;
            noisy[17] = 1.118762;
            noisy[18] = -1.903650;
            noisy[19] = 2.810987;

            truth[0]  = 0.000000;
            truth[1]  = 0.001000;
            truth[2]  = 0.003000;
            truth[3]  = 0.006000;
            truth[4]  = 0.010000;
            truth[5]  = 0.015000;
            truth[6]  = 0.021000;
            truth[7]  = 0.028000;
            truth[8]  = 0.036000;
            truth[9]  = 0.045000;
            truth[10] = 0.055000;
            truth[11] = 0.066000;
            truth[12] = 0.078000;
            truth[13] = 0.091000;
            truth[14] = 0.105000;
            truth[15] = 0.120000;
            truth[16] = 0.136000;
            truth[17] = 0.153000;
            truth[18] = 0.171000;
            truth[19] = 0.190000;
        }

        const f32 processVariance     = 0.01f;
        const f32 initialVariance     = 2.0f;
        const f32 observationVariance = 2.0f;

        UNITTEST_TEST(first)
        {
            f32* noisy;
            f32* truth;
            s32  n = 20;
            generateValues(n, noisy, truth, observationVariance);

            const u32         memory_size = 1024 * 1024;
            void*             memory_base = Allocator->allocate(memory_size);
            nkalman::memory_t memory(memory_base, memory_size);

            nkalman::nmodels::simplemodel_config_t simple_cfg = {initialVariance, processVariance, observationVariance};
            nkalman::nmodels::simplemodel_t*       model      = nkalman::nmodels::NewSimpleModel(&memory, 0, noisy[0], simple_cfg);

            nkalman::filter_t* kf       = nkalman::NewFilter(model, &memory);
            f32*               filtered = g_allocate_array_and_clear<f32>(Allocator, n);

            for (s32 i = 0; i < n; i++)
            {
                u64 t = i * 1000;  // time in milliseconds

                const bool pred_ok = nkalman::Predict(kf, t);
                CHECK_TRUE(pred_ok)  // Prediction failed

                memory.PushScope();
                {
                    nkalman::nmodels::measurement_t measurement = model->NewMeasurement(&memory, noisy[i]);

                    const bool upd_ok = nkalman::Update(kf, t, &measurement);
                    CHECK_TRUE(upd_ok);  // Update failed

                    nkalman::nmath::vector_t* state = nkalman::State(kf);
                    filtered[i]                     = state->AtVec(0);
                }
                memory.PopScope();
            }

            // Check that the filtered values are closer to the truth than the noisy values
            for (s32 i = 0; i < n; i++)
            {
                const f32 noisy_value    = noisy[i];
                const f32 filtered_value = filtered[i];
                const f32 truth_value    = truth[i];

                printf("t=%d noisy=%.6f filtered=%.6f truth=%.6f\n", va_t(i), va_t(noisy_value), va_t(filtered_value), va_t(truth_value));

                const f32 noisy_error    = fabsf(noisy_value - truth_value);
                const f32 filtered_error = fabsf(filtered_value - truth_value);
                CHECK_TRUE(filtered_error <= noisy_error); // Filtered value is not closer to truth than noisy value
            }
        }
    }
}
UNITTEST_SUITE_END
