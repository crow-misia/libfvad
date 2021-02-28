/*
 *  Copyright (c) 2012 The WebRTC project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "../src/signal_processing/signal_processing_library.h"
#include "test_common.h"
#include <string.h>


#ifdef TEST_SPL_MACRO
void test_main() {
  // Macros with inputs.
  int B = 21;
  int a = -3;
  int b = WEBRTC_SPL_WORD32_MAX;

  EXPECT_EQ(-63, WEBRTC_SPL_MUL(a, B));
  EXPECT_EQ(-2147483645, WEBRTC_SPL_MUL(a, b)); // note: this check has been removed in upstream webrtc
}
#endif // TEST_SPL_MACRO


#ifdef TEST_SPL_INLINE
void test_main() {
  int32_t a32 = 111121;

  EXPECT_EQ(17, WebRtcSpl_GetSizeInBits(a32));

  EXPECT_EQ(0, WebRtcSpl_NormW32(0));
  EXPECT_EQ(31, WebRtcSpl_NormW32(-1));
  EXPECT_EQ(0, WebRtcSpl_NormW32(WEBRTC_SPL_WORD32_MIN));
  EXPECT_EQ(14, WebRtcSpl_NormW32(a32));

  EXPECT_EQ(0, WebRtcSpl_NormU32(0u));
  EXPECT_EQ(0, WebRtcSpl_NormU32(0xffffffff));
  EXPECT_EQ(15, WebRtcSpl_NormU32(a32));
}
#endif // TEST_SPL_INLINE

#ifdef TEST_SPL_LEADING_ZEROS
void test_main() {
  EXPECT_EQ(32, WebRtcSpl_CountLeadingZeros32(0));
  EXPECT_EQ(32, WebRtcSpl_CountLeadingZeros32_NotBuiltin(0));
  for (int i = 0; i < 32; ++i) {
    const uint32_t single_one = (uint32_t)1 << i;
    const uint32_t all_ones = 2 * single_one - 1;
    EXPECT_EQ(31 - i, WebRtcSpl_CountLeadingZeros32(single_one));
    EXPECT_EQ(31 - i, WebRtcSpl_CountLeadingZeros32_NotBuiltin(single_one));
    EXPECT_EQ(31 - i, WebRtcSpl_CountLeadingZeros32(all_ones));
    EXPECT_EQ(31 - i, WebRtcSpl_CountLeadingZeros32_NotBuiltin(all_ones));
  }
}
#endif // TEST_SPL_LEADING_ZEROS

#ifdef TEST_SPL_MATH_OPERATIONS
void test_main() {

  int32_t num = 117;
  int32_t den = -5;

  EXPECT_EQ(-23, WebRtcSpl_DivW32W16(num, (int16_t)den));
  EXPECT_EQ(23u, WebRtcSpl_DivU32U16(num, denU));
  EXPECT_EQ(0, WebRtcSpl_DivW32HiLow(128, 0, 256));
}

TEST(SplTest, BasicArrayOperationsTest) {
  const size_t kVectorSize = 4;
  int B[] = {4, 12, 133, 1100};
  int16_t b16[kVectorSize];
  int32_t b32[kVectorSize];

  int16_t bTmp16[kVectorSize];
  int32_t bTmp32[kVectorSize];

  WebRtcSpl_MemSetW16(b16, 3, kVectorSize);
  for (size_t kk = 0; kk < kVectorSize; ++kk) {
    EXPECT_EQ(3, b16[kk]);
  }
  WebRtcSpl_ZerosArrayW16(b16, kVectorSize);
  for (size_t kk = 0; kk < kVectorSize; ++kk) {
    EXPECT_EQ(0, b16[kk]);
  }
  WebRtcSpl_MemSetW32(b32, 3, kVectorSize);
  for (size_t kk = 0; kk < kVectorSize; ++kk) {
    EXPECT_EQ(3, b32[kk]);
  }
  WebRtcSpl_ZerosArrayW32(b32, kVectorSize);
  for (size_t kk = 0; kk < kVectorSize; ++kk) {
    EXPECT_EQ(0, b32[kk]);
  }
  for (size_t kk = 0; kk < kVectorSize; ++kk) {
    bTmp16[kk] = (int16_t)kk;
    bTmp32[kk] = (int32_t)kk;
  }
  WEBRTC_SPL_MEMCPY_W16(b16, bTmp16, kVectorSize);
  for (size_t kk = 0; kk < kVectorSize; ++kk) {
    EXPECT_EQ(b16[kk], bTmp16[kk]);
  }
  //    WEBRTC_SPL_MEMCPY_W32(b32, bTmp32, kVectorSize);
  //    for (int kk = 0; kk < kVectorSize; ++kk) {
  //        EXPECT_EQ(b32[kk], bTmp32[kk]);
  //    }
  WebRtcSpl_CopyFromEndW16(b16, kVectorSize, 2, bTmp16);
  for (size_t kk = 0; kk < 2; ++kk) {
    EXPECT_EQ(static_cast<int16_t>(kk + 2), bTmp16[kk]);
  }

  for (size_t kk = 0; kk < kVectorSize; ++kk) {
    b32[kk] = B[kk];
    b16[kk] = (int16_t)B[kk];
  }
  WebRtcSpl_VectorBitShiftW32ToW16(bTmp16, kVectorSize, b32, 1);
  for (size_t kk = 0; kk < kVectorSize; ++kk) {
    EXPECT_EQ((B[kk] >> 1), bTmp16[kk]);
  }
  WebRtcSpl_VectorBitShiftW16(bTmp16, kVectorSize, b16, 1);
  for (size_t kk = 0; kk < kVectorSize; ++kk) {
    EXPECT_EQ((B[kk] >> 1), bTmp16[kk]);
  }
  WebRtcSpl_VectorBitShiftW32(bTmp32, kVectorSize, b32, 1);
  for (size_t kk = 0; kk < kVectorSize; ++kk) {
    EXPECT_EQ((B[kk] >> 1), bTmp32[kk]);
  }

  WebRtcSpl_MemCpyReversedOrder(&bTmp16[3], b16, kVectorSize);
  for (size_t kk = 0; kk < kVectorSize; ++kk) {
    EXPECT_EQ(b16[3 - kk], bTmp16[kk]);
  }
}

TEST(SplTest, MinMaxOperationsTest) {
  const size_t kVectorSize = 17;

  // Vectors to test the cases where minimum values have to be caught
  // outside of the unrolled loops in ARM-Neon.
  int16_t vector16[kVectorSize] = {-1,
                                   7485,
                                   0,
                                   3333,
                                   -18283,
                                   0,
                                   12334,
                                   -29871,
                                   988,
                                   -3333,
                                   345,
                                   -456,
                                   222,
                                   999,
                                   888,
                                   8774,
                                   WEBRTC_SPL_WORD16_MIN};
  int32_t vector32[kVectorSize] = {-1,
                                   0,
                                   283211,
                                   3333,
                                   8712345,
                                   0,
                                   -3333,
                                   89345,
                                   -374585456,
                                   222,
                                   999,
                                   122345334,
                                   -12389756,
                                   -987329871,
                                   888,
                                   -2,
                                   WEBRTC_SPL_WORD32_MIN};

  EXPECT_EQ(WEBRTC_SPL_WORD16_MIN,
            WebRtcSpl_MinValueW16(vector16, kVectorSize));
  EXPECT_EQ(WEBRTC_SPL_WORD32_MIN,
            WebRtcSpl_MinValueW32(vector32, kVectorSize));
  EXPECT_EQ(kVectorSize - 1, WebRtcSpl_MinIndexW16(vector16, kVectorSize));
  EXPECT_EQ(kVectorSize - 1, WebRtcSpl_MinIndexW32(vector32, kVectorSize));
  EXPECT_EQ(WEBRTC_SPL_WORD16_MIN,
            WebRtcSpl_MaxAbsElementW16(vector16, kVectorSize));
  int16_t min_value, max_value;
  WebRtcSpl_MinMaxW16(vector16, kVectorSize, &min_value, &max_value);
  EXPECT_EQ(WEBRTC_SPL_WORD16_MIN, min_value);
  EXPECT_EQ(12334, max_value);

  // Test the cases where maximum values have to be caught
  // outside of the unrolled loops in ARM-Neon.
  vector16[kVectorSize - 1] = WEBRTC_SPL_WORD16_MAX;
  vector32[kVectorSize - 1] = WEBRTC_SPL_WORD32_MAX;

  EXPECT_EQ(WEBRTC_SPL_WORD16_MAX,
            WebRtcSpl_MaxAbsValueW16(vector16, kVectorSize));
  EXPECT_EQ(WEBRTC_SPL_WORD16_MAX,
            WebRtcSpl_MaxValueW16(vector16, kVectorSize));
  EXPECT_EQ(WEBRTC_SPL_WORD32_MAX,
            WebRtcSpl_MaxAbsValueW32(vector32, kVectorSize));
  EXPECT_EQ(WEBRTC_SPL_WORD32_MAX,
            WebRtcSpl_MaxValueW32(vector32, kVectorSize));
  EXPECT_EQ(kVectorSize - 1, WebRtcSpl_MaxAbsIndexW16(vector16, kVectorSize));
  EXPECT_EQ(kVectorSize - 1, WebRtcSpl_MaxIndexW16(vector16, kVectorSize));
  EXPECT_EQ(kVectorSize - 1, WebRtcSpl_MaxIndexW32(vector32, kVectorSize));
  EXPECT_EQ(WEBRTC_SPL_WORD16_MAX,
            WebRtcSpl_MaxAbsElementW16(vector16, kVectorSize));
  WebRtcSpl_MinMaxW16(vector16, kVectorSize, &min_value, &max_value);
  EXPECT_EQ(-29871, min_value);
  EXPECT_EQ(WEBRTC_SPL_WORD16_MAX, max_value);

  // Test the cases where multiple maximum and minimum values are present.
  vector16[1] = WEBRTC_SPL_WORD16_MAX;
  vector16[6] = WEBRTC_SPL_WORD16_MIN;
  vector16[11] = WEBRTC_SPL_WORD16_MIN;
  vector32[1] = WEBRTC_SPL_WORD32_MAX;
  vector32[6] = WEBRTC_SPL_WORD32_MIN;
  vector32[11] = WEBRTC_SPL_WORD32_MIN;

  EXPECT_EQ(WEBRTC_SPL_WORD16_MAX,
            WebRtcSpl_MaxAbsValueW16(vector16, kVectorSize));
  EXPECT_EQ(WEBRTC_SPL_WORD16_MAX,
            WebRtcSpl_MaxValueW16(vector16, kVectorSize));
  EXPECT_EQ(WEBRTC_SPL_WORD16_MIN,
            WebRtcSpl_MinValueW16(vector16, kVectorSize));
  EXPECT_EQ(WEBRTC_SPL_WORD32_MAX,
            WebRtcSpl_MaxAbsValueW32(vector32, kVectorSize));
  EXPECT_EQ(WEBRTC_SPL_WORD32_MAX,
            WebRtcSpl_MaxValueW32(vector32, kVectorSize));
  EXPECT_EQ(WEBRTC_SPL_WORD32_MIN,
            WebRtcSpl_MinValueW32(vector32, kVectorSize));
  EXPECT_EQ(6u, WebRtcSpl_MaxAbsIndexW16(vector16, kVectorSize));
  EXPECT_EQ(1u, WebRtcSpl_MaxIndexW16(vector16, kVectorSize));
  EXPECT_EQ(1u, WebRtcSpl_MaxIndexW32(vector32, kVectorSize));
  EXPECT_EQ(6u, WebRtcSpl_MinIndexW16(vector16, kVectorSize));
  EXPECT_EQ(6u, WebRtcSpl_MinIndexW32(vector32, kVectorSize));
  EXPECT_EQ(WEBRTC_SPL_WORD16_MIN,
            WebRtcSpl_MaxAbsElementW16(vector16, kVectorSize));
  WebRtcSpl_MinMaxW16(vector16, kVectorSize, &min_value, &max_value);
  EXPECT_EQ(WEBRTC_SPL_WORD16_MIN, min_value);
  EXPECT_EQ(WEBRTC_SPL_WORD16_MAX, max_value);

  // Test a one-element vector.
  int16_t single_element_vector = 0;
  EXPECT_EQ(0, WebRtcSpl_MaxAbsValueW16(&single_element_vector, 1));
  EXPECT_EQ(0, WebRtcSpl_MaxValueW16(&single_element_vector, 1));
  EXPECT_EQ(0, WebRtcSpl_MinValueW16(&single_element_vector, 1));
  EXPECT_EQ(0u, WebRtcSpl_MaxAbsIndexW16(&single_element_vector, 1));
  EXPECT_EQ(0u, WebRtcSpl_MaxIndexW16(&single_element_vector, 1));
  EXPECT_EQ(0u, WebRtcSpl_MinIndexW16(&single_element_vector, 1));
  EXPECT_EQ(0, WebRtcSpl_MaxAbsElementW16(&single_element_vector, 1));
  WebRtcSpl_MinMaxW16(&single_element_vector, 1, &min_value, &max_value);
  EXPECT_EQ(0, min_value);
  EXPECT_EQ(0, max_value);

  // Test a two-element vector with the values WEBRTC_SPL_WORD16_MIN and
  // WEBRTC_SPL_WORD16_MAX.
  int16_t two_element_vector[2] = {WEBRTC_SPL_WORD16_MIN,
                                   WEBRTC_SPL_WORD16_MAX};
  EXPECT_EQ(WEBRTC_SPL_WORD16_MAX,
            WebRtcSpl_MaxAbsValueW16(two_element_vector, 2));
  EXPECT_EQ(WEBRTC_SPL_WORD16_MAX,
            WebRtcSpl_MaxValueW16(two_element_vector, 2));
  EXPECT_EQ(WEBRTC_SPL_WORD16_MIN,
            WebRtcSpl_MinValueW16(two_element_vector, 2));
  EXPECT_EQ(0u, WebRtcSpl_MaxAbsIndexW16(two_element_vector, 2));
  EXPECT_EQ(1u, WebRtcSpl_MaxIndexW16(two_element_vector, 2));
  EXPECT_EQ(0u, WebRtcSpl_MinIndexW16(two_element_vector, 2));
  EXPECT_EQ(WEBRTC_SPL_WORD16_MIN,
            WebRtcSpl_MaxAbsElementW16(two_element_vector, 2));
  WebRtcSpl_MinMaxW16(two_element_vector, 2, &min_value, &max_value);
  EXPECT_EQ(WEBRTC_SPL_WORD16_MIN, min_value);
  EXPECT_EQ(WEBRTC_SPL_WORD16_MAX, max_value);
}

TEST(SplTest, VectorOperationsTest) {
  const size_t kVectorSize = 4;
  int B[] = {4, 12, 133, 1100};
  int16_t a16[kVectorSize];
  int16_t b16[kVectorSize];
  int16_t bTmp16[kVectorSize];

  for (size_t kk = 0; kk < kVectorSize; ++kk) {
    a16[kk] = B[kk];
    b16[kk] = B[kk];
  }

  WebRtcSpl_AffineTransformVector(bTmp16, b16, 3, 7, 2, kVectorSize);
  for (size_t kk = 0; kk < kVectorSize; ++kk) {
    EXPECT_EQ((B[kk] * 3 + 7) >> 2, bTmp16[kk]);
  }
  WebRtcSpl_ScaleAndAddVectorsWithRound(b16, 3, b16, 2, 2, bTmp16, kVectorSize);
  for (size_t kk = 0; kk < kVectorSize; ++kk) {
    EXPECT_EQ((B[kk] * 3 + B[kk] * 2 + 2) >> 2, bTmp16[kk]);
  }

  WebRtcSpl_AddAffineVectorToVector(bTmp16, b16, 3, 7, 2, kVectorSize);
  for (size_t kk = 0; kk < kVectorSize; ++kk) {
    EXPECT_EQ(((B[kk] * 3 + B[kk] * 2 + 2) >> 2) + ((b16[kk] * 3 + 7) >> 2),
              bTmp16[kk]);
  }

  WebRtcSpl_ScaleVector(b16, bTmp16, 13, kVectorSize, 2);
  for (size_t kk = 0; kk < kVectorSize; ++kk) {
    EXPECT_EQ((b16[kk] * 13) >> 2, bTmp16[kk]);
  }
  WebRtcSpl_ScaleVectorWithSat(b16, bTmp16, 13, kVectorSize, 2);
  for (size_t kk = 0; kk < kVectorSize; ++kk) {
    EXPECT_EQ((b16[kk] * 13) >> 2, bTmp16[kk]);
  }
  WebRtcSpl_ScaleAndAddVectors(a16, 13, 2, b16, 7, 2, bTmp16, kVectorSize);
  for (size_t kk = 0; kk < kVectorSize; ++kk) {
    EXPECT_EQ(((a16[kk] * 13) >> 2) + ((b16[kk] * 7) >> 2), bTmp16[kk]);
  }

  WebRtcSpl_AddVectorsAndShift(bTmp16, a16, b16, kVectorSize, 2);
  for (size_t kk = 0; kk < kVectorSize; ++kk) {
    EXPECT_EQ(B[kk] >> 1, bTmp16[kk]);
  }
  WebRtcSpl_ReverseOrderMultArrayElements(bTmp16, a16, &b16[3], kVectorSize, 2);
  for (size_t kk = 0; kk < kVectorSize; ++kk) {
    EXPECT_EQ((a16[kk] * b16[3 - kk]) >> 2, bTmp16[kk]);
  }
  WebRtcSpl_ElementwiseVectorMult(bTmp16, a16, b16, kVectorSize, 6);
  for (size_t kk = 0; kk < kVectorSize; ++kk) {
    EXPECT_EQ((a16[kk] * b16[kk]) >> 6, bTmp16[kk]);
  }

  WebRtcSpl_SqrtOfOneMinusXSquared(b16, kVectorSize, bTmp16);
  for (size_t kk = 0; kk < kVectorSize - 1; ++kk) {
    EXPECT_EQ(32767, bTmp16[kk]);
  }
  EXPECT_EQ(32749, bTmp16[kVectorSize - 1]);

  EXPECT_EQ(0, WebRtcSpl_GetScalingSquare(b16, kVectorSize, 1));
}

TEST(SplTest, EstimatorsTest) {
  const size_t kOrder = 2;
  const int32_t unstable_filter[] = {4, 12, 133, 1100};
  const int32_t stable_filter[] = {1100, 133, 12, 4};
  int16_t lpc[kOrder + 2] = {0};
  int16_t refl[kOrder + 2] = {0};
  int16_t lpc_result[] = {4096, -497, 15, 0};
  int16_t refl_result[] = {-3962, 123, 0, 0};

  EXPECT_EQ(0, WebRtcSpl_LevinsonDurbin(unstable_filter, lpc, refl, kOrder));
  EXPECT_EQ(1, WebRtcSpl_LevinsonDurbin(stable_filter, lpc, refl, kOrder));
  for (size_t i = 0; i < kOrder + 2; ++i) {
    EXPECT_EQ(lpc_result[i], lpc[i]);
    EXPECT_EQ(refl_result[i], refl[i]);
  }
}

TEST(SplTest, FilterTest) {
  const size_t kVectorSize = 4;
  const size_t kFilterOrder = 3;
  int16_t A[] = {1, 2, 33, 100};
  int16_t A5[] = {1, 2, 33, 100, -5};
  int16_t B[] = {4, 12, 133, 110};
  int16_t data_in[kVectorSize];
  int16_t data_out[kVectorSize];
  int16_t bTmp16Low[kVectorSize];
  int16_t bState[kVectorSize];
  int16_t bStateLow[kVectorSize];

  WebRtcSpl_ZerosArrayW16(bState, kVectorSize);
  WebRtcSpl_ZerosArrayW16(bStateLow, kVectorSize);

  for (size_t kk = 0; kk < kVectorSize; ++kk) {
    data_in[kk] = A[kk];
    data_out[kk] = 0;
  }

  // MA filters.
  // Note that the input data has |kFilterOrder| states before the actual
  // data (one sample).
  WebRtcSpl_FilterMAFastQ12(&data_in[kFilterOrder], data_out, B,
                            kFilterOrder + 1, 1);
  EXPECT_EQ(0, data_out[0]);
  // AR filters.
  // Note that the output data has |kFilterOrder| states before the actual
  // data (one sample).
  WebRtcSpl_FilterARFastQ12(data_in, &data_out[kFilterOrder], A,
                            kFilterOrder + 1, 1);
  EXPECT_EQ(0, data_out[kFilterOrder]);

  EXPECT_EQ(kVectorSize, WebRtcSpl_FilterAR(A5, 5, data_in, kVectorSize, bState,
                                            kVectorSize, bStateLow, kVectorSize,
                                            data_out, bTmp16Low, kVectorSize));
}

TEST(SplTest, RandTest) {
  const int kVectorSize = 4;
  int16_t BU[] = {3653, 12446, 8525, 30691};
  int16_t b16[kVectorSize];
  uint32_t bSeed = 100000;

  EXPECT_EQ(7086, WebRtcSpl_RandU(&bSeed));
  EXPECT_EQ(31565, WebRtcSpl_RandU(&bSeed));
  EXPECT_EQ(-9786, WebRtcSpl_RandN(&bSeed));
  EXPECT_EQ(kVectorSize, WebRtcSpl_RandUArray(b16, kVectorSize, &bSeed));
  for (int kk = 0; kk < kVectorSize; ++kk) {
    EXPECT_EQ(BU[kk], b16[kk]);
  }
}

TEST(SplTest, DotProductWithScaleTest) {
  EXPECT_EQ(605362796, WebRtcSpl_DotProductWithScale(vector16, vector16,
                                                     kVector16Size, 2));
}

TEST(SplTest, CrossCorrelationTest) {
  // Note the function arguments relation specificed by API.
  const size_t kCrossCorrelationDimension = 3;
  const int kShift = 2;
  const int kStep = 1;
  const size_t kSeqDimension = 6;

  const int16_t kVector16[kVector16Size] = {
      1,    4323, 1963, WEBRTC_SPL_WORD16_MAX, WEBRTC_SPL_WORD16_MIN + 5, -3333,
      -876, 8483, 142};
  int32_t vector32[kCrossCorrelationDimension] = {0};

  WebRtcSpl_CrossCorrelation(vector32, vector16, kVector16, kSeqDimension,
                             kCrossCorrelationDimension, kShift, kStep);

  // WebRtcSpl_CrossCorrelationC() and WebRtcSpl_CrossCorrelationNeon()
  // are not bit-exact.
  const int32_t kExpected[kCrossCorrelationDimension] = {-266947903, -15579555,
                                                         -171282001};
  const int32_t* expected = kExpected;
#if !defined(MIPS32_LE)
  const int32_t kExpectedNeon[kCrossCorrelationDimension] = {
      -266947901, -15579553, -171281999};
  if (WebRtcSpl_CrossCorrelation != WebRtcSpl_CrossCorrelationC) {
    expected = kExpectedNeon;
  }
#endif
  for (size_t i = 0; i < kCrossCorrelationDimension; ++i) {
    EXPECT_EQ(expected[i], vector32[i]);
  }
}
#endif // TEST_SPL_MATH_OPERATIONS

#ifdef TEST_SPL_SIGNAL_PROCESSING
void test_main() {
  const size_t kVectorSize = 4;
  int A[] = {1, 2, 33, 100};
  const int16_t kHanning[4] = {2399, 8192, 13985, 16384};
  int16_t b16[kVectorSize];

  int16_t bTmp16[kVectorSize];

  int bScale = 0;

  for (size_t kk = 0; kk < kVectorSize; ++kk) {
    b16[kk] = A[kk];
  }

  EXPECT_EQ(11094, WebRtcSpl_Energy(b16, kVectorSize, &bScale));
  EXPECT_EQ(0, bScale);
}
#endif // TEST_SPL_SIGNAL_PROCESSING

#ifdef TEST_SPL_RESAMPLE_48
void test_main() {
  // The test resamples 3*kBlockSize number of samples to 2*kBlockSize number
  // of samples.
  #define kBlockSize 16

  // Saturated input vector of 48 samples.
  const int32_t kVectorSaturated[3 * kBlockSize + 7] = {
      -32768, -32768, -32768, -32768, -32768, -32768, -32768, -32768,
      -32768, -32768, -32768, -32768, -32768, -32768, -32768, -32768,
      -32768, -32768, -32768, -32768, -32768, -32768, -32768, -32768,
      32767,  32767,  32767,  32767,  32767,  32767,  32767,  32767,
      32767,  32767,  32767,  32767,  32767,  32767,  32767,  32767,
      32767,  32767,  32767,  32767,  32767,  32767,  32767,  32767,
      32767,  32767,  32767,  32767,  32767,  32767,  32767};

  // All values in |out_vector| should be |kRefValue32kHz|.
  const int32_t kRefValue32kHz1 = -1077493760;
  const int32_t kRefValue32kHz2 = 1077493645;

  // After bit shift with saturation, |out_vector_w16| is saturated.

  //const int16_t kRefValue16kHz1 = -32768;
  //const int16_t kRefValue16kHz2 = 32767;
  // Vector for storing output.
  int32_t out_vector[2 * kBlockSize];
  // int16_t out_vector_w16[2 * kBlockSize];

  WebRtcSpl_Resample48khzTo32khz(kVectorSaturated, out_vector, kBlockSize);

  // Comparing output values against references. The values at position
  // 12-15 are skipped to account for the filter lag.
  for (size_t i = 0; i < 12; ++i) {
    EXPECT_EQ(kRefValue32kHz1, out_vector[i]);
    //EXPECT_EQ(kRefValue16kHz1, out_vector_w16[i]);
  }
  for (size_t i = 16; i < 2 * kBlockSize; ++i) {
    EXPECT_EQ(kRefValue32kHz2, out_vector[i]);
    //EXPECT_EQ(kRefValue16kHz2, out_vector_w16[i]);
  }
}
#endif // TEST_SPL_RESAMPLE_48
