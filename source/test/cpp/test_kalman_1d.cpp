#include "ccore/c_allocator.h"
#include "ccore/c_math.h"
#include "ccore/c_printf.h"
#include "ccore/c_random.h"

#include "ckalman/c_kalman.h"

#include "cunittest/cunittest.h"

#include <cmath>

using namespace ncore;

UNITTEST_SUITE_BEGIN(kalman_1D)
{
    UNITTEST_FIXTURE(tests)
    {
        UNITTEST_TEST(kalman_1d_initialize_begin_update)
        {
            nkalman::kalman_1D_t kf;
            nkalman::initialize(kf, 0.01f, 0.1f, 1.0f);

            CHECK_CLOSE(0.01f, kf.q, 0.00001f);
            CHECK_CLOSE(0.1f, kf.r, 0.00001f);
            CHECK_CLOSE(1.0f, kf.p, 0.00001f);
            CHECK_CLOSE(0.0f, kf.x, 0.00001f);
            CHECK_CLOSE(0.0f, kf.k, 0.00001f);

            nkalman::begin(kf, 0.0f, 1.0f);
            const f32 x = nkalman::update(kf, 1.0f);

            CHECK_CLOSE(0.9099099f, x, 0.0001f);
            CHECK_CLOSE(0.9099099f, nkalman::getState(kf), 0.0001f);
            CHECK_CLOSE(0.0909910f, kf.p, 0.0001f);
            CHECK_CLOSE(0.9099099f, kf.k, 0.0001f);
        }

        UNITTEST_TEST(kalman_1d_converges_to_constant_signal)
        {
            nkalman::kalman_1D_t kf;
            nkalman::initialize(kf, 0.01f, 0.1f, 1.0f);
            nkalman::begin(kf, 0.0f, 1.0f);

            for (s32 i = 0; i < 30; ++i)
            {
                nkalman::update(kf, 10.0f);
            }

            CHECK_CLOSE(10.0f, nkalman::getState(kf), 0.1f);
        }
    }
}
UNITTEST_SUITE_END
