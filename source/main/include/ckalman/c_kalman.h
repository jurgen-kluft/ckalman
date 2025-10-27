#ifndef __C_KALMAN_FILTER_H__
#define __C_KALMAN_FILTER_H__
#include "ccore/c_target.h"
#ifdef USE_PRAGMA_ONCE
    #pragma once
#endif

namespace ncore
{
    namespace nmath
    {
        struct Vector;
        struct Matrix;
    }  // namespace nmath

    namespace nmodels
    {
        struct State
        {
            u64            Time;
            nmath::Vector* State;
            nmath::Matrix* Covariance;
        };

        struct Measurement
        {
            nmath::Matrix* Covariance;
            nmath::Vector* Value;
            nmath::Matrix* ObservationModel;
        };

        // LinearModel is used to initialize hidden states in the model and
        // provide transition matrices to the filter.
        // kalman/models provides commonly used models.
        class LinearModel
        {
        public:
            virtual void InitialState(State& state)                              = 0;
            virtual void Transition(u64 dt, nmath::Matrix*& outMatrix)           = 0;
            virtual void CovarianceTransition(u64 dt, nmath::Matrix*& outMatrix) = 0;
        };
    }  // namespace nmodels

    namespace nkalman
    {
        // KalmanFilter is responsible for prediction and filtering
        // of a given linear model. It is assumed that the process being modelled
        // is a time series, and that the time steps are non-uniform and specified
        // for each update and prediction operation.
        struct KalmanFilter
        {
            nmodels::LinearModel* model;
            s32                   dims;
            u64                   t;
            nmath::Vector*        state;
            nmath::Matrix*        covariance;
        };

        KalmanFilter* NewKalmanFilter(nmodels::LinearModel* model);

        // State returns the current hidden state of the KalmanFilter.
        // Example models provided with this package often provide functions
        // to extract meaningful information from the state vector, such as
        // .Velocity() for the provided constant velocity model.
        nmath::Vector* State(KalmanFilter* kf);

        // Covariance returns the current covaraince of the model.
        nmath::Matrix* Covariance(KalmanFilter* kf);

        // SetCovariance resets the covariance of the Kalman Filter to the given value.
        void SetCovariance(KalmanFilter* kf, nmath::Matrix* covariance);

        // SetState resets the state of the Kalman Filter to the given value.
        void SetState(KalmanFilter* kf, nmath::Vector* state);

        // Time returns the time for which the current hidden state is an estimate.
        // The time is monotone increasing.
        u64 Time(KalmanFilter* kf);

        // Predict advances the KalmanFilter from the internal current time
        // to the given time using the built-in linear model.
        // The state of the filter is updated and the current time is updated.
        // Each time can be no earlier than the current time of the filter.
        bool Predict(KalmanFilter* kf, u64 t);

        // Update is used to take a new measurement from a sensor and fuse it to the model.
        // The time field must be no earlier than the current time of the filter.
        bool Update(KalmanFilter* kf, u64 t, nmodels::Measurement* m);

        nmath::Matrix* Eye(KalmanFilter* kf, s32 n);

    }  // namespace nkalman
}  // namespace ncore
#endif  // __C_KALMAN_FILTER_H__
