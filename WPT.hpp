#pragma once

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

#include <memory>
#include <vector>

// Ensure.hpp is an optional dependency: if it's available (either as part of this
// checkout or vendored alongside this header), use its ensure() for a formatted
// diagnostic on failure; otherwise fall back to plain assert() so this header still
// works standalone.
#if __has_include("commons/Ensure.hpp")
	#include "commons/Ensure.hpp"
#elif __has_include("Ensure.hpp")
	#include "Ensure.hpp"
#else
	#include <cassert>
	// Guard against Ensure.hpp having already been included under a path our
	// __has_include checks above don't know about (e.g. vendored elsewhere as
	// "3rdparty/Ensure.hpp") -- COMMONS_ENSURE_HPP is defined by Ensure.hpp
	// itself, so this still catches that case even under an unknown filename.
	#if !defined(COMMONS_ENSURE_HPP) && !defined(ensure)
		#define ensure(condition, ...) assert((condition))
	#endif
#endif

/**
 * @brief A full (not just single-level) Wavelet Packet Transform: recursively
 * splits both the low- and high-pass sub-bands at every level, unlike a
 * plain DWT which only recurses on the low/approximation branch.
 *
 * Give it a mother wavelet's decomposition (high-pass/detail) filter; the
 * low-pass (scaling) filter is derived automatically via the standard
 * quadrature mirror relationship, so any orthogonal wavelet family works --
 * Haar is just the 2-tap case:
 * 		wpt({1, -1}, 10)                              // Haar
 *		wpt({-0.1294, -0.2241, 0.8365, -0.4830}, 10)  // Daubechies-2
 *
 * push() is per-sample and online rather than block-based: each Level keeps
 * a small circular buffer and alternates which of two physical sub-levels
 * (Even/Odd) a given push descends into, computing both decimation phases
 * across time instead of discarding one. That's what makes this shift-
 * invariant: every leaf in getEnvelope() gets a fresh value on every push,
 * not just once every 2^depth samples the way a critically-sampled WPT
 * would -- at the cost of maintaining 2^depth Level objects and giving up
 * exact invertibility.
 *
 * Forward decomposition only, by design: this was built as a continuously-
 * updating streaming decomposition (an alternative to a STFT or a fixed
 * filter bank), not a compress/reconstruct codec. Wavelet still derives the
 * reconstruction (synthesis) filter pair (h/l) alongside the decomposition
 * pair (hp/lp) as a byproduct of deriving lp from hp, but nothing here
 * consumes h/l.
 */
namespace dsp {

	/**
	 * @brief Derives the full 4-filter analysis/synthesis set for an
	 * orthogonal wavelet from just its decomposition high-pass (detail)
	 * filter, via the standard quadrature mirror relationship.
	 */
	struct Wavelet {
		/**
		 * @brief Derives hp/h/lp/l/width from a mother wavelet's
		 * decomposition high-pass filter.
		 * @param mother The decomposition high-pass (detail) filter
		 *               coefficients, e.g. {1, -1} for Haar. Must be
		 *               non-empty.
		 * @note Terminates the process (via ensure()) if mother is empty.
		 */
		Wavelet(const std::vector<float>& mother)
			: hp(mother), h(mirror(mother)), lp(flipSigns(h)), l(mirror(lp)), width(mother.size()) {
			ensure(!mother.empty(), "Wavelet requires a non-empty mother filter.");
		}

		/// @brief Reverses a filter's tap order.
		static std::vector<float> mirror(const std::vector<float>& wavelet) {
			std::vector<float> mirror;
			mirror.reserve(wavelet.size());
			for(int index = wavelet.size() - 1; index >= 0; --index)
				mirror.push_back(wavelet[index]);
			return mirror;
		}

		/// @brief Negates every even-indexed tap (0, 2, 4, ...); half of the
		/// quadrature mirror relationship used to derive lp from h.
		static std::vector<float> flipSigns(const std::vector<float>& wavelet) {
			std::vector<float> mirror = wavelet;
			for(std::size_t index = 0; index < wavelet.size(); index += 2)
				mirror[index] = -mirror[index];
			return mirror;
		}

		// Declaration order is the initialization order (C++ ignores the
		// member-initializer-list order): each filter below is derived from
		// the one before it, so this order must track that dependency chain.
		// hp/lp are the decomposition pair Level actually convolves against;
		// h/l are the corresponding reconstruction pair for a future inverse
		// transform and are otherwise unused today.
		const std::vector<float> hp;	// mother ([1, -1 for Haar)
		const std::vector<float> h;	// mirror(hp)
		const std::vector<float> lp;	// flipSigns(h); father]
		const std::vector<float> l;	// mirror(lp)

		const unsigned int width;
	};

	/**
	 * @brief Interface for the type a Level tree delivers leaf outputs to: a
	 * single add(address, value) callback, invoked once per leaf address on
	 * every push().
	 */
	class FrameBuilder {
	public:
		virtual void add(std::uint32_t address, float value) = 0;
	};

	/**
	 * @brief One node of the recursive wavelet packet tree. Not intended to
	 * be constructed directly -- WPT builds and owns the whole tree via its
	 * topLevel member.
	 */
	class Level {
	public:
		/**
		 * @brief Constructs this node and, if it isn't a leaf, its four
		 * Even/Odd children.
		 * @param maxDepth     Total tree depth; this node is a leaf once
		 *                     depth == maxDepth - 1.
		 * @param depth        This node's depth from the root (0 at the root).
		 * @param wavelet      The shared filter set convolved against at
		 *                     every node in the tree.
		 * @param frameBuilder Receives leaf outputs via
		 *                     frameBuilder.add(address, value) once this
		 *                     node (or a descendant) reaches a leaf.
		 * @param address      This node's packet address, used to derive
		 *                     its children's addresses.
		 */
		Level(unsigned int maxDepth, unsigned int depth, const Wavelet& wavelet, FrameBuilder& frameBuilder, std::uint32_t address = 0)
			: depth(depth), wavelet(wavelet), lowAddress(address << 1), highAddress((address << 1) | 1),
			isLeaf(depth == maxDepth - 1), frameBuilder(frameBuilder), signal(wavelet.width) {
			if(!isLeaf) {
				lowEven = std::make_unique<Level>(maxDepth, depth + 1, wavelet, frameBuilder, lowAddress);
				highEven = std::make_unique<Level>(maxDepth, depth + 1, wavelet, frameBuilder, highAddress);
				lowOdd = std::make_unique<Level>(maxDepth, depth + 1, wavelet, frameBuilder, lowAddress);
				highOdd = std::make_unique<Level>(maxDepth, depth + 1, wavelet, frameBuilder, highAddress);
			}
		}

		/**
		 * @brief Pushes one sample through this node: convolves it into
		 * low/high outputs, then either delivers them to frameBuilder (if
		 * isLeaf) or recurses into one of this node's two Even/Odd child
		 * pairs, alternating each call to realize decimation-by-2 without
		 * buffering a block first.
		 */
		void push(float value) {
			// Store the value in our buffer and advance
			signal[oldestSignalIndex] = value;
			oldestSignalIndex = (oldestSignalIndex + 1) % wavelet.width;

			// Convolve with wavelets
			const float lowValue = convolve(wavelet.lp);
			const float highValue = convolve(wavelet.hp);

			if(isLeaf) {
				// Send to frame builders
				frameBuilder.add(lowAddress, lowValue);
				frameBuilder.add(highAddress, highValue);
			} else {
				// Send to lower levels
				if(pushToEven) {
					lowEven->push(lowValue);
					highEven->push(highValue);
				} else {
					lowOdd->push(lowValue);
					highOdd->push(highValue);
				}
				pushToEven = !pushToEven;
			}
		}

		/**
		 * @brief Convolves filter w against this node's circular sample
		 * buffer.
		 * @param w A filter of the same width as the shared Wavelet
		 *          (wavelet.lp or wavelet.hp).
		 * @return The convolution sum divided by the filter width.
		 */
		float convolve(const std::vector<float>& w) {
			float sum = 0;
			for(unsigned int index = 0; index < wavelet.width; ++index)
				sum += w[index] * signal[(oldestSignalIndex + index) % wavelet.width];
			return sum / wavelet.width;
		}

		const unsigned int depth;        ///< This node's depth from the root (0 at the root).
		const Wavelet& wavelet;          ///< The shared filter set convolved against at every node.
		const std::uint32_t lowAddress;  ///< Packet address this node's low-band output (leaf) or low subtree (non-leaf) writes to.
		const std::uint32_t highAddress; ///< Packet address this node's high-band output (leaf) or high subtree (non-leaf) writes to.
		const bool isLeaf;                ///< True once depth == maxDepth - 1; leaves call frameBuilder.add() instead of recursing.

	private:
		FrameBuilder& frameBuilder;

		bool pushToEven = true;
		std::unique_ptr<Level> lowEven;
		std::unique_ptr<Level> highEven;
		std::unique_ptr<Level> lowOdd;
		std::unique_ptr<Level> highOdd;

		std::vector<float> signal;
		unsigned int oldestSignalIndex = 0;
	};

	class WPT : public FrameBuilder {
	public:
		/**
		 * @brief Constructs a wavelet packet transform tree.
		 * @param mother   The mother wavelet's decomposition high-pass
		 *                 (detail) filter, e.g. {1, -1} for Haar. Must be
		 *                 non-empty.
		 * @param maxDepth Number of tree levels; getEnvelope() will hold
		 *                 2^maxDepth addresses. Must be > 0.
		 * @note Terminates the process (via ensure()) if mother is empty or
		 *       maxDepth is 0.
		 */
		WPT(const std::vector<float>& mother, unsigned int maxDepth = 10)
			: wavelet(mother), topLevel(checkedDepth(maxDepth), 0, wavelet, *this), envelope(1 << maxDepth) {
		}

		/**
		 * @brief Pushes one sample through the tree, refreshing every entry
		 * in getEnvelope() with this push's full packet decomposition.
		 */
		void push(float value) {
			// No need to zero envelope first: a single push() always calls
			// add() for every address in [0, envelope.size()) exactly once,
			// since Level's Even/Odd children of a given push always cover
			// the same full leaf-address range regardless of which phase
			// fires (verified empirically across maxDepth 1-6).
			topLevel.push(value);
		}

		/**
		 * @brief FrameBuilder callback invoked by the Level tree during
		 * push() -- not intended to be called directly; call push() instead.
		 */
		void add(std::uint32_t address, float value) override {
			envelope[address] = value;
		}

		/**
		 * @brief The most recent push()'s packet decomposition: one entry
		 * per leaf address in [0, 2^maxDepth).
		 */
		constexpr const std::vector<float>& getEnvelope() const {
			return envelope;
		}

	private:
		// maxDepth - 1 underflows (unsigned) if maxDepth is 0, which would make
		// Level::isLeaf never true and recurse without bound; called from the
		// member-initializer list so topLevel is never constructed with a bad depth.
		static unsigned int checkedDepth(unsigned int maxDepth) {
			ensure(maxDepth > 0, "WPT requires maxDepth > 0.");
			return maxDepth;
		}

		Wavelet wavelet;
		Level topLevel;

		std::vector<float> envelope;
	};

}
