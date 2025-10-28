#ifndef __C_KALMAN_FILTER_H__
#define __C_KALMAN_FILTER_H__
#include "ccore/c_target.h"
#ifdef USE_PRAGMA_ONCE
    #pragma once
#endif

namespace ncore
{
    namespace nkalman
    {
        namespace nmath
        {
            struct vector_t;
            struct matrix_t;
        }  // namespace nmath

        namespace nmodels
        {
            struct measurement_t;
            class model_t;
        }  // namespace nmodels

        struct filter_t;

        filter_t* NewFilter(nmodels::model_t* model);

        // state_t returns the current hidden state of the filter_t.
        // Example models provided with this package often provide functions
        // to extract meaningful information from the state vector, such as
        // .Velocity() for the provided constant velocity model.
        nmath::vector_t* State(filter_t* kf);

        // Covariance returns the current covaraince of the model.
        nmath::matrix_t* Covariance(filter_t* kf);

        // SetCovariance resets the covariance of the Kalman Filter to the given value.
        void SetCovariance(filter_t* kf, nmath::matrix_t* covariance);

        // SetState resets the state of the Kalman Filter to the given value.
        void SetState(filter_t* kf, nmath::vector_t* state);

        // Time returns the time for which the current hidden state is an estimate.
        // The time is monotone increasing.
        u64 Time(filter_t* kf);

        // Predict advances the filter_t from the internal current time
        // to the given time using the built-in linear model.
        // The state of the filter is updated and the current time is updated.
        // Each time can be no earlier than the current time of the filter.
        bool Predict(filter_t* kf, u64 t);

        // Update is used to take a new measurement from a sensor and fuse it to the model.
        // The time field must be no earlier than the current time of the filter.
        bool Update(filter_t* kf, u64 t, nmodels::measurement_t* m);

        nmath::matrix_t* Eye(filter_t* kf, s32 n);

    }  // namespace nkalman
}  // namespace ncore
#endif  // __C_KALMAN_FILTER_H__
