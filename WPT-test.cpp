/*
 * Copyright 2026 L. Richard Moore Jr.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include "WPT.hpp"

using namespace dsp;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::Approx;

// Published Daubechies-2 (db2) detail/wavelet filter (Daubechies, 1992),
// shared by both db2 test cases below so the reference coefficients only
// appear once.
static const std::vector<float> db2Mother{-0.129409522550921f, -0.224143868041857f, 0.836516303737469f, -0.482962913144690f};

// This is an orthogonal wavelet family (Daubechies, of which Haar is the 2-tap
// special case). Given only the high-pass/detail ("mother") filter, the
// low-pass/scaling ("father") filter is derived via the standard quadrature
// mirror relationship g[n] = -(-1)^n * h[N-1-n]. Verify that derivation against
// hand-worked Haar coefficients and against the published Daubechies-2 (db2)
// filter pair.
TEST_CASE( "Wavelet derives the correct QMF low-pass filter - Haar" ) {
	Wavelet haar({1.0f, -1.0f});
	CHECK(haar.width == 2);
	CHECK_THAT(haar.hp, Approx(std::vector<float>{1.0f, -1.0f}).margin(1e-6));
	CHECK_THAT(haar.h, Approx(std::vector<float>{-1.0f, 1.0f}).margin(1e-6));
	CHECK_THAT(haar.lp, Approx(std::vector<float>{1.0f, 1.0f}).margin(1e-6));
	CHECK_THAT(haar.l, Approx(std::vector<float>{1.0f, 1.0f}).margin(1e-6));
}

TEST_CASE( "Wavelet derives the correct QMF low-pass filter - Daubechies-2" ) {
	Wavelet db2(db2Mother);
	// Published db2 scaling/low-pass filter - independent reference values,
	// not derived from this code.
	std::vector<float> expectedLp{0.482962913144690f, 0.836516303737469f, 0.224143868041857f, -0.129409522550921f};
	CHECK_THAT(db2.lp, Approx(expectedLp).margin(1e-6));
}

// Edge case: every other test uses an even-length filter (width 2 or 4).
// flipSigns negates every even index (0, 2, 4, ...), so an odd-length filter
// exercises a loop bound (index < size, stepping by 2) that never lands
// exactly on size-1 for even-length inputs but does for odd-length ones.
// Hand-derived: mother = {1, 2, 3}; h = mirror(mother) = {3, 2, 1};
// flipSigns negates indices 0 and 2 (both even) -> lp = {-3, 2, -1};
// l = mirror(lp) = {-1, 2, -3}.
TEST_CASE( "Wavelet derives the correct QMF low-pass filter - odd width" ) {
	Wavelet odd({1.0f, 2.0f, 3.0f});
	CHECK(odd.width == 3);
	CHECK_THAT(odd.h, Approx(std::vector<float>{3.0f, 2.0f, 1.0f}).margin(1e-6));
	CHECK_THAT(odd.lp, Approx(std::vector<float>{-3.0f, 2.0f, -1.0f}).margin(1e-6));
	CHECK_THAT(odd.l, Approx(std::vector<float>{-1.0f, 2.0f, -3.0f}).margin(1e-6));
}

// Edge case: envelope is allocated in WPT's constructor and zeroed at the
// top of every push(), but nothing zeroes it before the first push. Confirm
// it doesn't start out as uninitialized/garbage memory.
TEST_CASE( "WPT envelope is all zero before any push" ) {
	WPT wpt({1.0f, -1.0f}, 3);
	std::vector<float> zeros(8, 0.0f);
	CHECK_THAT(wpt.getEnvelope(), Approx(zeros).margin(1e-6));
}

// Edge case: width = 1 is the smallest legal (non-empty) filter. The
// circular buffer has exactly one slot, so oldestSignalIndex is always 0
// after every push ((0 + 1) % 1 == 0) -- each push fully replaces the
// buffer's only sample, unlike width >= 2 where a push blends old and new
// samples. Hand-derived: mother = {2}; h = mirror({2}) = {2};
// flipSigns negates index 0 (even) -> lp = {-2}; hp = mother = {2}.
// convolve(w) with a 1-slot buffer holding the just-pushed value v reduces
// to w[0] * v / 1, so low = -2v, high = 2v, with no memory of prior pushes.
TEST_CASE( "WPT single-tap (width=1) filter has no memory across pushes" ) {
	WPT wpt({2.0f}, 1);
	wpt.push(3.0f);
	CHECK_THAT(wpt.getEnvelope(), Approx(std::vector<float>{-6.0f, 6.0f}).margin(1e-6));
	wpt.push(-1.0f);
	CHECK_THAT(wpt.getEnvelope(), Approx(std::vector<float>{2.0f, -2.0f}).margin(1e-6));
}

// maxDepth = 1 collapses the tree to a single leaf Level, isolating the raw
// convolve()/circular-buffer mechanics from the packet routing. Expected
// values are hand-derived from the Haar filters above: for a two-sample
// window [old, new], low = (old + new) / 2 and high = (new - old) / 2.
TEST_CASE( "WPT single-level Haar impulse response matches hand-derived values" ) {
	WPT wpt({1.0f, -1.0f}, 1);
	struct Sample { float input; float low; float high; };
	std::vector<Sample> samples{
		{1.0f, 0.5f, -0.5f},
		{0.0f, 0.5f, 0.5f},
		{0.0f, 0.0f, 0.0f},
		{0.0f, 0.0f, 0.0f},
	};
	for(const auto& sample : samples) {
		wpt.push(sample.input);
		const auto& envelope = wpt.getEnvelope();
		CHECK_THAT(envelope[0], WithinAbs(sample.low, 1e-6));
		CHECK_THAT(envelope[1], WithinAbs(sample.high, 1e-6));
	}
}

// Same isolation as above but with the Daubechies-2 filters, so the expected
// values are a direct 4-tap convolution of an impulse against the published
// db2 lp/hp coefficients, computed by hand: for an impulse the only nonzero
// term is w[3] (the tap aligned with the newest, and only, sample), divided
// by the filter width (-0.129409522550921 / 4 and -0.482962913144690 / 4).
TEST_CASE( "WPT single-level Daubechies-2 impulse response matches hand-derived values" ) {
	WPT wpt(db2Mother, 1);
	wpt.push(1.0f);
	const auto& envelope = wpt.getEnvelope();
	CHECK_THAT(envelope[0], WithinAbs(-0.03235238063773025f, 1e-6));
	CHECK_THAT(envelope[1], WithinAbs(-0.12074072828617250f, 1e-6));
}

// Edge case: the test above only pushes once, so oldestSignalIndex only ever
// takes its initial (0 -> 1) step and the width-4 circular buffer never
// wraps back around to slot 0. Push 5 samples so oldestSignalIndex cycles
// through every slot once and back (0->1->2->3->0), meaning the 5th push
// overwrites the same slot the 1st push wrote to. Values 1-4 were confirmed
// against the single-push test above (a leading [1,0,0,0] impulse produces
// identical outputs either way); the 5th value was hand-derived separately:
// after 4 pushes the buffer holds [1,0,0,0] with oldestSignalIndex back at
// 0, so pushing 2.0 makes it [2,0,0,0] with oldestSignalIndex=1, and
// convolve only picks up the tap aligned with slot 0 (index 3 in the
// window starting at slot 1): low = lp[3]*2/4, high = hp[3]*2/4.
TEST_CASE( "WPT single-level Daubechies-2 circular buffer wraps correctly" ) {
	WPT wpt(db2Mother, 1);
	std::vector<float> input{1.0f, 0.0f, 0.0f, 0.0f, 2.0f};
	std::vector<std::vector<float>> expected{
		{-0.03235238063773025f, -0.12074072828617250f},
		{0.05603596701046425f, 0.20912907593436725f},
		{0.20912907593436725f, -0.05603596701046425f},
		{0.12074072828617250f, -0.03235238063773025f},
		{-0.06470476127546050f, -0.24148145657234500f},
	};
	REQUIRE(input.size() == expected.size());
	for(size_t index = 0; index < input.size(); ++index) {
		wpt.push(input[index]);
		CHECK_THAT(wpt.getEnvelope(), Approx(expected[index]).margin(1e-6));
	}
}

// maxDepth = 3 exercises the recursive even/odd packet routing (the part
// that's distinct from a plain DWT). Expected values were hand-traced
// leaf-by-leaf through the first push: the root splits [1] into
// low=0.5/high=-0.5, each of those is split again one level down, and again
// at the leaves, each split halving magnitude (1 -> 0.5 -> 0.25 -> 0.125) per
// the Haar low/high derivation used above.
TEST_CASE( "WPT three-level Haar packet routing matches hand-traced impulse response" ) {
	WPT wpt({1.0f, -1.0f}, 3);
	wpt.push(1.0f);
	std::vector<float> expected{0.125f, -0.125f, -0.125f, 0.125f, -0.125f, 0.125f, 0.125f, -0.125f};
	CHECK_THAT(wpt.getEnvelope(), Approx(expected).margin(1e-6));
}

// The single-push tests above only exercise each Level once, so pushToEven
// never flips and lowOdd/highOdd never fire. This test pushes a short
// sequence through maxDepth=2 (envelope size 4) so every depth-1 node
// alternates through both its Even and Odd children at least once, which
// is the part of the class the original header comment ("doesn't crash,
// AFAIK... not sure [if it's working]") was never confident about. The
// first two pushes' expected values (all 8 numbers) were hand-traced
// leaf-by-leaf the same way as the depth-3 test above; the remaining four
// pushes continue the same by-hand method and were cross-checked against
// an independent Python port of this exact algorithm.
TEST_CASE( "WPT depth-2 Haar packet routing matches hand-traced multi-push sequence" ) {
	WPT wpt({1.0f, -1.0f}, 2);
	std::vector<float> input{1.0f, 0.5f, -0.5f, 0.25f, -0.25f, 0.75f};
	std::vector<std::vector<float>> expected{
		{0.25f, -0.25f, -0.25f, 0.25f},
		{0.375f, -0.375f, 0.125f, -0.125f},
		{0.25f, 0.25f, 0.0f, -0.5f},
		{0.3125f, 0.4375f, -0.0625f, 0.3125f},
		{0.0f, 0.0f, 0.375f, 0.125f},
		{0.0625f, -0.1875f, -0.4375f, 0.0625f},
	};
	REQUIRE(input.size() == expected.size());
	for(size_t index = 0; index < input.size(); ++index) {
		wpt.push(input[index]);
		CHECK_THAT(wpt.getEnvelope(), Approx(expected[index]).margin(1e-6));
	}
}

// Throughput of steady-state push() at maxDepth=10, the depth used in this
// header's own doc-comment examples (wpt({1, -1}, 10)). Each push does
// O(maxDepth) recursive convolutions (4 children per non-leaf Level, each
// doing a width-sized convolve twice) plus a memset of the 2^maxDepth-entry
// envelope, so this is the number to watch if the recursion, the per-Level
// convolution, or the envelope-clearing strategy changes.
TEST_CASE( "WPT push() throughput" ) {
	WPT wpt({1.0f, -1.0f}, 10);
	BENCHMARK("push (maxDepth=10)") {
		wpt.push(1.0f);
	};
}