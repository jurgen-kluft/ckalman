#include "ccore/c_allocator.h"
#include "ccore/c_math.h"
#include "ccore/c_printf.h"
#include "ccore/c_random.h"

#include "ckalman/c_kalman.h"

#include "cunittest/cunittest.h"

#include <cmath>

using namespace ncore;

UNITTEST_SUITE_BEGIN(matrix)
{
    UNITTEST_FIXTURE(tests)
    {
        UNITTEST_TEST(matrix_identity_clear_add_subtract)
        {
            nkalman::matrix_t<2, 2> a;
            nkalman::matrix_t<2, 2> b;
            nkalman::matrix_t<2, 2> out;

            a.setIdentity();
            CHECK_CLOSE(1.0f, a.data[0][0], 0.00001f);
            CHECK_CLOSE(0.0f, a.data[0][1], 0.00001f);
            CHECK_CLOSE(0.0f, a.data[1][0], 0.00001f);
            CHECK_CLOSE(1.0f, a.data[1][1], 0.00001f);

            b.data[0][0] = 2.0f;
            b.data[0][1] = 3.0f;
            b.data[1][0] = 4.0f;
            b.data[1][1] = 5.0f;

            a.add(b, out);
            CHECK_CLOSE(3.0f, out.data[0][0], 0.00001f);
            CHECK_CLOSE(3.0f, out.data[0][1], 0.00001f);
            CHECK_CLOSE(4.0f, out.data[1][0], 0.00001f);
            CHECK_CLOSE(6.0f, out.data[1][1], 0.00001f);

            out.subtract(b, out);
            CHECK_CLOSE(1.0f, out.data[0][0], 0.00001f);
            CHECK_CLOSE(0.0f, out.data[0][1], 0.00001f);
            CHECK_CLOSE(0.0f, out.data[1][0], 0.00001f);
            CHECK_CLOSE(1.0f, out.data[1][1], 0.00001f);

            out.clear();
            CHECK_CLOSE(0.0f, out.data[0][0], 0.00001f);
            CHECK_CLOSE(0.0f, out.data[0][1], 0.00001f);
            CHECK_CLOSE(0.0f, out.data[1][0], 0.00001f);
            CHECK_CLOSE(0.0f, out.data[1][1], 0.00001f);
        }

        UNITTEST_TEST(matrix_multiply_transpose_invert)
        {
            nkalman::matrix_t<2, 2> a;
            nkalman::matrix_t<2, 2> b;
            nkalman::matrix_t<2, 2> mul;
            nkalman::matrix_t<2, 2> trans;
            nkalman::matrix_t<2, 2> inv;
            nkalman::matrix_t<2, 2> ident;

            a.data[0][0] = 4.0f;
            a.data[0][1] = 7.0f;
            a.data[1][0] = 2.0f;
            a.data[1][1] = 6.0f;

            b.data[0][0] = 1.0f;
            b.data[0][1] = 2.0f;
            b.data[1][0] = 3.0f;
            b.data[1][1] = 4.0f;

            a.multiply(b, mul);
            CHECK_CLOSE(25.0f, mul.data[0][0], 0.00001f);
            CHECK_CLOSE(36.0f, mul.data[0][1], 0.00001f);
            CHECK_CLOSE(20.0f, mul.data[1][0], 0.00001f);
            CHECK_CLOSE(28.0f, mul.data[1][1], 0.00001f);

            a.transpose(trans);
            CHECK_CLOSE(4.0f, trans.data[0][0], 0.00001f);
            CHECK_CLOSE(2.0f, trans.data[0][1], 0.00001f);
            CHECK_CLOSE(7.0f, trans.data[1][0], 0.00001f);
            CHECK_CLOSE(6.0f, trans.data[1][1], 0.00001f);

            a.invert(inv);
            CHECK_CLOSE(0.6f, inv.data[0][0], 0.00001f);
            CHECK_CLOSE(-0.7f, inv.data[0][1], 0.00001f);
            CHECK_CLOSE(-0.2f, inv.data[1][0], 0.00001f);
            CHECK_CLOSE(0.4f, inv.data[1][1], 0.00001f);

            a.multiply(inv, ident);
            CHECK_CLOSE(1.0f, ident.data[0][0], 0.0001f);
            CHECK_CLOSE(0.0f, ident.data[0][1], 0.0001f);
            CHECK_CLOSE(0.0f, ident.data[1][0], 0.0001f);
            CHECK_CLOSE(1.0f, ident.data[1][1], 0.0001f);
        }
    }
}
UNITTEST_SUITE_END
