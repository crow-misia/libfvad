/*
 *  Copyright (c) 2012 The WebRTC project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

/*
 * This header file includes all of the fix point signal processing library
 * (SPL) function descriptions and declarations. For specific function calls,
 * see bottom of file.
 */

#ifndef COMMON_AUDIO_SIGNAL_PROCESSING_INCLUDE_SIGNAL_PROCESSING_LIBRARY_H_
#define COMMON_AUDIO_SIGNAL_PROCESSING_INCLUDE_SIGNAL_PROCESSING_LIBRARY_H_

#include "../common.h"

// Macros specific for the fixed point implementation
#define WEBRTC_SPL_WORD16_MAX       32767
#define WEBRTC_SPL_WORD16_MIN       -32768
#define WEBRTC_SPL_WORD32_MAX       (int32_t)0x7fffffff
#define WEBRTC_SPL_WORD32_MIN       (int32_t)0x80000000

#define WEBRTC_SPL_MUL(a, b) \
    ((int32_t) ((int32_t)(a) * (int32_t)(b)))


// inline functions:
#include "spl_inl.h"

int16_t WebRtcSpl_GetScalingSquare(int16_t* in_vector,
                                   size_t in_vector_length,
                                   size_t times);

// Copy and set operations. Implementation in copy_set_operations.c.
// Descriptions at bottom of file.
void WebRtcSpl_MemSetW16(int16_t* vector,
                         int16_t set_value,
                         size_t vector_length);
void WebRtcSpl_MemSetW32(int32_t* vector,
                         int32_t set_value,
                         size_t vector_length);
void WebRtcSpl_MemCpyReversedOrder(int16_t* out_vector,
                                   int16_t* in_vector,
                                   size_t vector_length);
void WebRtcSpl_CopyFromEndW16(const int16_t* in_vector,
                              size_t in_vector_length,
                              size_t samples,
                              int16_t* out_vector);
void WebRtcSpl_ZerosArrayW16(int16_t* vector, size_t vector_length);
void WebRtcSpl_ZerosArrayW32(int32_t* vector, size_t vector_length);
// End: Copy and set operations.

// Minimum and maximum operation functions and their pointers.
// Implementation in min_max_operations.c.

// Returns the largest absolute value in a signed 16-bit vector.
//
// Input:
//      - vector : 16-bit input vector.
//      - length : Number of samples in vector.
//
// Return value  : Maximum absolute value in vector.
typedef int16_t (*MaxAbsValueW16)(const int16_t* vector, size_t length);
extern const MaxAbsValueW16 WebRtcSpl_MaxAbsValueW16;
int16_t WebRtcSpl_MaxAbsValueW16C(const int16_t* vector, size_t length);
#if defined(WEBRTC_HAS_NEON)
int16_t WebRtcSpl_MaxAbsValueW16Neon(const int16_t* vector, size_t length);
#endif
#if defined(MIPS32_LE)
int16_t WebRtcSpl_MaxAbsValueW16_mips(const int16_t* vector, size_t length);
#endif

// Returns the largest absolute value in a signed 32-bit vector.
//
// Input:
//      - vector : 32-bit input vector.
//      - length : Number of samples in vector.
//
// Return value  : Maximum absolute value in vector.
typedef int32_t (*MaxAbsValueW32)(const int32_t* vector, size_t length);
extern const MaxAbsValueW32 WebRtcSpl_MaxAbsValueW32;
int32_t WebRtcSpl_MaxAbsValueW32C(const int32_t* vector, size_t length);
#if defined(WEBRTC_HAS_NEON)
int32_t WebRtcSpl_MaxAbsValueW32Neon(const int32_t* vector, size_t length);
#endif
#if defined(MIPS_DSP_R1_LE)
int32_t WebRtcSpl_MaxAbsValueW32_mips(const int32_t* vector, size_t length);
#endif

// Returns the maximum value of a 16-bit vector.
//
// Input:
//      - vector : 16-bit input vector.
//      - length : Number of samples in vector.
//
// Return value  : Maximum sample value in |vector|.
typedef int16_t (*MaxValueW16)(const int16_t* vector, size_t length);
extern const MaxValueW16 WebRtcSpl_MaxValueW16;
int16_t WebRtcSpl_MaxValueW16C(const int16_t* vector, size_t length);
#if defined(WEBRTC_HAS_NEON)
int16_t WebRtcSpl_MaxValueW16Neon(const int16_t* vector, size_t length);
#endif
#if defined(MIPS32_LE)
int16_t WebRtcSpl_MaxValueW16_mips(const int16_t* vector, size_t length);
#endif

// Returns the maximum value of a 32-bit vector.
//
// Input:
//      - vector : 32-bit input vector.
//      - length : Number of samples in vector.
//
// Return value  : Maximum sample value in |vector|.
typedef int32_t (*MaxValueW32)(const int32_t* vector, size_t length);
extern const MaxValueW32 WebRtcSpl_MaxValueW32;
int32_t WebRtcSpl_MaxValueW32C(const int32_t* vector, size_t length);
#if defined(WEBRTC_HAS_NEON)
int32_t WebRtcSpl_MaxValueW32Neon(const int32_t* vector, size_t length);
#endif
#if defined(MIPS32_LE)
int32_t WebRtcSpl_MaxValueW32_mips(const int32_t* vector, size_t length);
#endif

// Returns the minimum value of a 16-bit vector.
//
// Input:
//      - vector : 16-bit input vector.
//      - length : Number of samples in vector.
//
// Return value  : Minimum sample value in |vector|.
typedef int16_t (*MinValueW16)(const int16_t* vector, size_t length);
extern const MinValueW16 WebRtcSpl_MinValueW16;
int16_t WebRtcSpl_MinValueW16C(const int16_t* vector, size_t length);
#if defined(WEBRTC_HAS_NEON)
int16_t WebRtcSpl_MinValueW16Neon(const int16_t* vector, size_t length);
#endif
#if defined(MIPS32_LE)
int16_t WebRtcSpl_MinValueW16_mips(const int16_t* vector, size_t length);
#endif

// Returns the minimum value of a 32-bit vector.
//
// Input:
//      - vector : 32-bit input vector.
//      - length : Number of samples in vector.
//
// Return value  : Minimum sample value in |vector|.
typedef int32_t (*MinValueW32)(const int32_t* vector, size_t length);
extern const MinValueW32 WebRtcSpl_MinValueW32;
int32_t WebRtcSpl_MinValueW32C(const int32_t* vector, size_t length);
#if defined(WEBRTC_HAS_NEON)
int32_t WebRtcSpl_MinValueW32Neon(const int32_t* vector, size_t length);
#endif
#if defined(MIPS32_LE)
int32_t WebRtcSpl_MinValueW32_mips(const int32_t* vector, size_t length);
#endif

// Returns both the minimum and maximum values of a 16-bit vector.
//
// Input:
//      - vector : 16-bit input vector.
//      - length : Number of samples in vector.
// Ouput:
//      - max_val : Maximum sample value in |vector|.
//      - min_val : Minimum sample value in |vector|.
void WebRtcSpl_MinMaxW16(const int16_t* vector,
                         size_t length,
                         int16_t* min_val,
                         int16_t* max_val);
#if defined(WEBRTC_HAS_NEON)
void WebRtcSpl_MinMaxW16Neon(const int16_t* vector,
                             size_t length,
                             int16_t* min_val,
                             int16_t* max_val);
#endif

// Returns the vector index to the largest absolute value of a 16-bit vector.
//
// Input:
//      - vector : 16-bit input vector.
//      - length : Number of samples in vector.
//
// Return value  : Index to the maximum absolute value in vector.
//                 If there are multiple equal maxima, return the index of the
//                 first. -32768 will always have precedence over 32767 (despite
//                 -32768 presenting an int16 absolute value of 32767).
size_t WebRtcSpl_MaxAbsIndexW16(const int16_t* vector, size_t length);

// Returns the element with the largest absolute value of a 16-bit vector. Note
// that this function can return a negative value.
//
// Input:
//      - vector : 16-bit input vector.
//      - length : Number of samples in vector.
//
// Return value  : The element with the largest absolute value. Note that this
//                 may be a negative value.
int16_t WebRtcSpl_MaxAbsElementW16(const int16_t* vector, size_t length);

// Returns the vector index to the maximum sample value of a 16-bit vector.
//
// Input:
//      - vector : 16-bit input vector.
//      - length : Number of samples in vector.
//
// Return value  : Index to the maximum value in vector (if multiple
//                 indexes have the maximum, return the first).
size_t WebRtcSpl_MaxIndexW16(const int16_t* vector, size_t length);

// Returns the vector index to the maximum sample value of a 32-bit vector.
//
// Input:
//      - vector : 32-bit input vector.
//      - length : Number of samples in vector.
//
// Return value  : Index to the maximum value in vector (if multiple
//                 indexes have the maximum, return the first).
size_t WebRtcSpl_MaxIndexW32(const int32_t* vector, size_t length);

// Returns the vector index to the minimum sample value of a 16-bit vector.
//
// Input:
//      - vector : 16-bit input vector.
//      - length : Number of samples in vector.
//
// Return value  : Index to the mimimum value in vector  (if multiple
//                 indexes have the minimum, return the first).
size_t WebRtcSpl_MinIndexW16(const int16_t* vector, size_t length);

// Returns the vector index to the minimum sample value of a 32-bit vector.
//
// Input:
//      - vector : 32-bit input vector.
//      - length : Number of samples in vector.
//
// Return value  : Index to the mimimum value in vector  (if multiple
//                 indexes have the minimum, return the first).
size_t WebRtcSpl_MinIndexW32(const int32_t* vector, size_t length);

// End: Minimum and maximum operations.

// Vector scaling operations. Implementation in vector_scaling_operations.c.
// Description at bottom of file.
void WebRtcSpl_VectorBitShiftW16(int16_t* out_vector,
                                 size_t vector_length,
                                 const int16_t* in_vector,
                                 int16_t right_shifts);
void WebRtcSpl_VectorBitShiftW32(int32_t* out_vector,
                                 size_t vector_length,
                                 const int32_t* in_vector,
                                 int16_t right_shifts);
void WebRtcSpl_VectorBitShiftW32ToW16(int16_t* out_vector,
                                      size_t vector_length,
                                      const int32_t* in_vector,
                                      int right_shifts);
void WebRtcSpl_ScaleVector(const int16_t* in_vector,
                           int16_t* out_vector,
                           int16_t gain,
                           size_t vector_length,
                           int16_t right_shifts);
void WebRtcSpl_ScaleVectorWithSat(const int16_t* in_vector,
                                  int16_t* out_vector,
                                  int16_t gain,
                                  size_t vector_length,
                                  int16_t right_shifts);
void WebRtcSpl_ScaleAndAddVectors(const int16_t* in_vector1,
                                  int16_t gain1,
                                  int right_shifts1,
                                  const int16_t* in_vector2,
                                  int16_t gain2,
                                  int right_shifts2,
                                  int16_t* out_vector,
                                  size_t vector_length);

// The functions (with related pointer) perform the vector operation:
//   out_vector[k] = ((scale1 * in_vector1[k]) + (scale2 * in_vector2[k])
//        + round_value) >> right_shifts,
//   where  round_value = (1 << right_shifts) >> 1.
//
// Input:
//      - in_vector1       : Input vector 1
//      - in_vector1_scale : Gain to be used for vector 1
//      - in_vector2       : Input vector 2
//      - in_vector2_scale : Gain to be used for vector 2
//      - right_shifts     : Number of right bit shifts to be applied
//      - length           : Number of elements in the input vectors
//
// Output:
//      - out_vector       : Output vector
// Return value            : 0 if OK, -1 if (in_vector1 == null
//                           || in_vector2 == null || out_vector == null
//                           || length <= 0 || right_shift < 0).
typedef int (*ScaleAndAddVectorsWithRound)(const int16_t* in_vector1,
                                           int16_t in_vector1_scale,
                                           const int16_t* in_vector2,
                                           int16_t in_vector2_scale,
                                           int right_shifts,
                                           int16_t* out_vector,
                                           size_t length);
extern const ScaleAndAddVectorsWithRound WebRtcSpl_ScaleAndAddVectorsWithRound;
int WebRtcSpl_ScaleAndAddVectorsWithRoundC(const int16_t* in_vector1,
                                           int16_t in_vector1_scale,
                                           const int16_t* in_vector2,
                                           int16_t in_vector2_scale,
                                           int right_shifts,
                                           int16_t* out_vector,
                                           size_t length);
#if defined(MIPS_DSP_R1_LE)
int WebRtcSpl_ScaleAndAddVectorsWithRound_mips(const int16_t* in_vector1,
                                               int16_t in_vector1_scale,
                                               const int16_t* in_vector2,
                                               int16_t in_vector2_scale,
                                               int right_shifts,
                                               int16_t* out_vector,
                                               size_t length);
#endif
// End: Vector scaling operations.

// iLBC specific functions. Implementations in ilbc_specific_functions.c.
// Description at bottom of file.
void WebRtcSpl_ReverseOrderMultArrayElements(int16_t* out_vector,
                                             const int16_t* in_vector,
                                             const int16_t* window,
                                             size_t vector_length,
                                             int16_t right_shifts);
void WebRtcSpl_ElementwiseVectorMult(int16_t* out_vector,
                                     const int16_t* in_vector,
                                     const int16_t* window,
                                     size_t vector_length,
                                     int16_t right_shifts);
void WebRtcSpl_AddVectorsAndShift(int16_t* out_vector,
                                  const int16_t* in_vector1,
                                  const int16_t* in_vector2,
                                  size_t vector_length,
                                  int16_t right_shifts);
void WebRtcSpl_AddAffineVectorToVector(int16_t* out_vector,
                                       const int16_t* in_vector,
                                       int16_t gain,
                                       int32_t add_constant,
                                       int16_t right_shifts,
                                       size_t vector_length);
void WebRtcSpl_AffineTransformVector(int16_t* out_vector,
                                     const int16_t* in_vector,
                                     int16_t gain,
                                     int32_t add_constant,
                                     int16_t right_shifts,
                                     size_t vector_length);
// End: iLBC specific functions.

// Signal processing operations.

// A 32-bit fix-point implementation of auto-correlation computation
//
// Input:
//      - in_vector        : Vector to calculate autocorrelation upon
//      - in_vector_length : Length (in samples) of |vector|
//      - order            : The order up to which the autocorrelation should be
//                           calculated
//
// Output:
//      - result           : auto-correlation values (values should be seen
//                           relative to each other since the absolute values
//                           might have been down shifted to avoid overflow)
//
//      - scale            : The number of left shifts required to obtain the
//                           auto-correlation in Q0
//
// Return value            : Number of samples in |result|, i.e. (order+1)
size_t WebRtcSpl_AutoCorrelation(const int16_t* in_vector,
                                 size_t in_vector_length,
                                 size_t order,
                                 int32_t* result,
                                 int* scale);

// A 32-bit fix-point implementation of the Levinson-Durbin algorithm that
// does NOT use the 64 bit class
//
// Input:
//      - auto_corr : Vector with autocorrelation values of length >= |order|+1
//      - order     : The LPC filter order (support up to order 20)
//
// Output:
//      - lpc_coef  : lpc_coef[0..order] LPC coefficients in Q12
//      - refl_coef : refl_coef[0...order-1]| Reflection coefficients in Q15
//
// Return value     : 1 for stable 0 for unstable
int16_t WebRtcSpl_LevinsonDurbin(const int32_t* auto_corr,
                                 int16_t* lpc_coef,
                                 int16_t* refl_coef,
                                 size_t order);

// Converts reflection coefficients |refl_coef| to LPC coefficients |lpc_coef|.
// This version is a 16 bit operation.
//
// NOTE: The 16 bit refl_coef -> lpc_coef conversion might result in a
// "slightly unstable" filter (i.e., a pole just outside the unit circle) in
// "rare" cases even if the reflection coefficients are stable.
//
// Input:
//      - refl_coef : Reflection coefficients in Q15 that should be converted
//                    to LPC coefficients
//      - use_order : Number of coefficients in |refl_coef|
//
// Output:
//      - lpc_coef  : LPC coefficients in Q12
void WebRtcSpl_ReflCoefToLpc(const int16_t* refl_coef,
                             int use_order,
                             int16_t* lpc_coef);

// Converts LPC coefficients |lpc_coef| to reflection coefficients |refl_coef|.
// This version is a 16 bit operation.
// The conversion is implemented by the step-down algorithm.
//
// Input:
//      - lpc_coef  : LPC coefficients in Q12, that should be converted to
//                    reflection coefficients
//      - use_order : Number of coefficients in |lpc_coef|
//
// Output:
//      - refl_coef : Reflection coefficients in Q15.
void WebRtcSpl_LpcToReflCoef(int16_t* lpc_coef,
                             int use_order,
                             int16_t* refl_coef);

// Calculates reflection coefficients (16 bit) from auto-correlation values
//
// Input:
//      - auto_corr : Auto-correlation values
//      - use_order : Number of coefficients wanted be calculated
//
// Output:
//      - refl_coef : Reflection coefficients in Q15.
void WebRtcSpl_AutoCorrToReflCoef(const int32_t* auto_corr,
                                  int use_order,
                                  int16_t* refl_coef);

// The functions (with related pointer) calculate the cross-correlation between
// two sequences |seq1| and |seq2|.
// |seq1| is fixed and |seq2| slides as the pointer is increased with the
// amount |step_seq2|. Note the arguments should obey the relationship:
// |dim_seq| - 1 + |step_seq2| * (|dim_cross_correlation| - 1) <
//      buffer size of |seq2|
//
// Input:
//      - seq1           : First sequence (fixed throughout the correlation)
//      - seq2           : Second sequence (slides |step_vector2| for each
//                            new correlation)
//      - dim_seq        : Number of samples to use in the cross-correlation
//      - dim_cross_correlation : Number of cross-correlations to calculate (the
//                            start position for |vector2| is updated for each
//                            new one)
//      - right_shifts   : Number of right bit shifts to use. This will
//                            become the output Q-domain.
//      - step_seq2      : How many (positive or negative) steps the
//                            |vector2| pointer should be updated for each new
//                            cross-correlation value.
//
// Output:
//      - cross_correlation : The cross-correlation in Q(-right_shifts)
typedef void (*CrossCorrelation)(int32_t* cross_correlation,
                                 const int16_t* seq1,
                                 const int16_t* seq2,
                                 size_t dim_seq,
                                 size_t dim_cross_correlation,
                                 int right_shifts,
                                 int step_seq2);
extern const CrossCorrelation WebRtcSpl_CrossCorrelation;
void WebRtcSpl_CrossCorrelationC(int32_t* cross_correlation,
                                 const int16_t* seq1,
                                 const int16_t* seq2,
                                 size_t dim_seq,
                                 size_t dim_cross_correlation,
                                 int right_shifts,
                                 int step_seq2);
#if defined(WEBRTC_HAS_NEON)
void WebRtcSpl_CrossCorrelationNeon(int32_t* cross_correlation,
                                    const int16_t* seq1,
                                    const int16_t* seq2,
                                    size_t dim_seq,
                                    size_t dim_cross_correlation,
                                    int right_shifts,
                                    int step_seq2);
#endif
#if defined(MIPS32_LE)
void WebRtcSpl_CrossCorrelation_mips(int32_t* cross_correlation,
                                     const int16_t* seq1,
                                     const int16_t* seq2,
                                     size_t dim_seq,
                                     size_t dim_cross_correlation,
                                     int right_shifts,
                                     int step_seq2);
#endif

// Creates (the first half of) a Hanning window. Size must be at least 1 and
// at most 512.
//
// Input:
//      - size      : Length of the requested Hanning window (1 to 512)
//
// Output:
//      - window    : Hanning vector in Q14.
void WebRtcSpl_GetHanningWindow(int16_t* window, size_t size);

// Calculates y[k] = sqrt(1 - x[k]^2) for each element of the input vector
// |in_vector|. Input and output values are in Q15.
//
// Inputs:
//      - in_vector     : Values to calculate sqrt(1 - x^2) of
//      - vector_length : Length of vector |in_vector|
//
// Output:
//      - out_vector    : Output values in Q15
void WebRtcSpl_SqrtOfOneMinusXSquared(int16_t* in_vector,
                                      size_t vector_length,
                                      int16_t* out_vector);
// End: Signal processing operations.

// Randomization functions. Implementations collected in
// randomization_functions.c and descriptions at bottom of this file.
int16_t WebRtcSpl_RandU(uint32_t* seed);
int16_t WebRtcSpl_RandN(uint32_t* seed);
int16_t WebRtcSpl_RandUArray(int16_t* vector,
                             int16_t vector_length,
                             uint32_t* seed);
// End: Randomization functions.

// Math functions
int32_t WebRtcSpl_Sqrt(int32_t value);

// Divisions. Implementations collected in division_operations.c and
// descriptions at bottom of this file.
int32_t WebRtcSpl_DivW32W16(int32_t num, int16_t den);
// End: Divisions.

int32_t WebRtcSpl_Energy(int16_t* vector,
                         size_t vector_length,
                         int* scale_factor);


/************************************************************
 *
 * RESAMPLING FUNCTIONS AND THEIR STRUCTS ARE DEFINED BELOW
 *
 ************************************************************/

/*******************************************************************
 * resample_fractional.c
 * Functions for internal use in the other resample functions
 *
 * Includes the following resampling combinations
 * 48 kHz -> 32 kHz
 *
 ******************************************************************/

void WebRtcSpl_Resample48khzTo32khz(const int32_t* In, int32_t* Out, size_t K);

/*******************************************************************
 * resample_48khz.c
 *
 * Includes the following resampling combinations
 * 48 kHz ->  8 kHz
 *
 ******************************************************************/

typedef struct {
  int32_t S_48_24[8];
  int32_t S_24_24[16];
  int32_t S_24_16[8];
  int32_t S_16_8[8];
} WebRtcSpl_State48khzTo8khz;

void WebRtcSpl_Resample48khzTo8khz(const int16_t* in,
                                   int16_t* out,
                                   WebRtcSpl_State48khzTo8khz* state,
                                   int32_t* tmpmem);

void WebRtcSpl_ResetResample48khzTo8khz(WebRtcSpl_State48khzTo8khz* state);

/************************************************************
 * END OF RESAMPLING FUNCTIONS
 ************************************************************/


#endif  // COMMON_AUDIO_SIGNAL_PROCESSING_INCLUDE_SIGNAL_PROCESSING_LIBRARY_H_
//
// WebRtcSpl_DivW32W16(...)
//
// Divides a int32_t |num| by a int16_t |den|.
//
// If |den|==0, (int32_t)0x7FFFFFFF is returned.
//
// Input:
//      - num       : Numerator
//      - den       : Denominator
//
// Return value     : Result of the division (as a int32_t), i.e., the
//                    integer part of num/den.
//

//
// WebRtcSpl_Energy(...)
//
// Calculates the energy of a vector
//
// Input:
//      - vector        : Vector which the energy should be calculated on
//      - vector_length : Number of samples in vector
//
// Output:
//      - scale_factor  : Number of left bit shifts needed to get the physical
//                        energy value, i.e, to get the Q0 value
//
// Return value         : Energy value in Q(-|scale_factor|)
//
