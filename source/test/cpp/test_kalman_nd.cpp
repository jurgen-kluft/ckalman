#include "ccore/c_allocator.h"
#include "ccore/c_math.h"
#include "ccore/c_printf.h"
#include "ccore/c_random.h"

#include "ckalman/c_kalman.h"

#include "cunittest/cunittest.h"

#include <cmath>

using namespace ncore;

UNITTEST_SUITE_BEGIN(kalman_nd)
{
    UNITTEST_FIXTURE(tests)
    {
        UNITTEST_TEST(kalman_nd_initialize_begin_and_single_update)
        {
            nkalman::kalman_nd_t<2, 1> kf;
            nkalman::initialize(kf);

            CHECK_CLOSE(1.0f, kf.P.data[0][0], 0.00001f);
            CHECK_CLOSE(0.0f, kf.P.data[0][1], 0.00001f);
            CHECK_CLOSE(0.0f, kf.P.data[1][0], 0.00001f);
            CHECK_CLOSE(1.0f, kf.P.data[1][1], 0.00001f);
            CHECK_CLOSE(1.0f, kf.Q.data[0][0], 0.00001f);
            CHECK_CLOSE(1.0f, kf.R.data[0][0], 0.00001f);
            CHECK_CLOSE(0.0f, kf.x.data[0][0], 0.00001f);
            CHECK_CLOSE(0.0f, kf.x.data[1][0], 0.00001f);

            kf.F.setIdentity();
            kf.H.clear();
            kf.H.data[0][0] = 1.0f;
            kf.H.data[0][1] = 0.0f;
            kf.Q.clear();
            kf.R.setIdentity();

            const f32 initial[2] = {0.0f, 0.0f};
            nkalman::begin(kf, initial, 1.0f);

            const f32 measurement[1] = {1.0f};
            nkalman::update(kf, measurement);

            CHECK_CLOSE(0.5f, kf.x.data[0][0], 0.0001f);
            CHECK_CLOSE(0.0f, kf.x.data[1][0], 0.0001f);
            CHECK_CLOSE(0.5f, kf.P.data[0][0], 0.0001f);
            CHECK_CLOSE(1.0f, kf.P.data[1][1], 0.0001f);
            CHECK_CLOSE(0.5f, kf.K.data[0][0], 0.0001f);
            CHECK_CLOSE(0.0f, kf.K.data[1][0], 0.0001f);
        }




    }
}
UNITTEST_SUITE_END
