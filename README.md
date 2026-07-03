# WPT.hpp

A full Wavelet Packet Transform: given a mother wavelet's decomposition filter, recursively
splits both the low- and high-pass sub-bands at every level -- unlike a plain DWT, which only
recurses on the low/approximation branch. `push()` is per-sample and online (not block-based):
each tree node keeps a small circular buffer and alternates which of two physical sub-levels
receives a given push, computing *both* decimation phases across time instead of discarding one.
That's what makes it shift-invariant: every leaf address in `getEnvelope()` gets a fresh value on
every single push, not just once every `2^depth` samples the way a critically-sampled WPT would --
at the cost of giving up exact invertibility (see Limitations below).

I originally wrote this while investigating whether a shift-invariant wavelet packet decomposition
might work better than a STFT or a fixed filter bank for some streaming analysis I was doing -- a
wavelet packet tree gives finer time resolution at the high-frequency leaves and finer frequency
resolution at the low ones "for free" from the recursive dyadic splitting, rather than needing to
hand-design each band the way a filter bank does. I ended up not adopting the approach for that
project, for reasons I don't fully remember, but the technique still seemed elegant enough to keep
around for the next time I have the same thought.

Furthermore, there really aren't a lot of simple C++ examples of WPT transforms floating around.
It's not a complicated concept, but having some working example to start with may be helpful --
particularly for the shift-invariant/streaming variant specifically, since most examples you'll
find elsewhere are the critically-sampled kind.

```cpp
#include "WPT.hpp"

dsp::WPT wpt({1.0f, -1.0f}, 4);   // Haar mother wavelet, 4 levels deep

for (float sample : signal)
    wpt.push(sample);

const std::vector<float>& envelope = wpt.getEnvelope();  // 2^4 = 16 packet coefficients,
                                                           // refreshed by every push()
```

Any orthogonal wavelet family works, not just Haar -- give it the mother wavelet's decomposition
high-pass (detail) filter and the low-pass (scaling) filter is derived automatically via the
standard quadrature mirror relationship:

```cpp
dsp::WPT wpt({-0.1294f, -0.2241f, 0.8365f, -0.4830f}, 10);  // Daubechies-2
```

## Requirements

- C++20 or later (uses concepts)
- Header-only -- copy `WPT.hpp` into your project and `#include` it
- Optional: [`Ensure.hpp`](https://github.com/lrmoorejr/ensure) for its `ensure()` helper; if
  it's not present, `WPT.hpp` falls back to plain `assert()` (see below)

## API

| Call | Behavior |
|---|---|
| `Wavelet(mother)` | Derives the full decomposition/reconstruction filter set from a mother wavelet's decomposition high-pass filter, e.g. `{1, -1}` for Haar. |
| `WPT(mother, maxDepth=10)` | Constructs a transform tree `maxDepth` levels deep from a mother wavelet filter. |
| `push(value)` | Pushes one sample through the tree, refreshing every entry in `getEnvelope()`. |
| `getEnvelope()` | The most recent `push()`'s packet decomposition: one entry per leaf address in `[0, 2^maxDepth)`. |

### Limitations

This is forward decomposition only -- there is no inverse/reconstruction transform. That's a
consequence of what it was built for (see above): a continuously-updating streaming decomposition
to read from, not a compress/reconstruct codec. If you need the original signal back, this isn't
the transform for that -- you'd want a standard critically-sampled WPT (one child per level
instead of two) paired with its synthesis filters. `Wavelet` still derives that reconstruction
(synthesis) filter pair (`h`/`l`) alongside the decomposition pair (`hp`/`lp`) as a byproduct of
deriving `lp` from `hp`, but nothing in this header consumes `h`/`l`.

### Validation

The constructors validate their input via `ensure()`, not by throwing: `Wavelet` requires a
non-empty mother filter, and `WPT` requires `maxDepth > 0`. Either violation terminates the
process with a diagnostic (or, without `Ensure.hpp` present, calls `assert()`) rather than
raising a catchable exception -- these are treated as programmer errors, not recoverable runtime
conditions.

### Ensure.hpp fallback

`WPT.hpp` uses one helper from [`Ensure.hpp`](https://github.com/lrmoorejr/ensure): `ensure()`,
a runtime assertion that prints a diagnostic and terminates on failure. If `Ensure.hpp` is
available -- either checked out alongside `WPT.hpp`, or reachable as `commons/Ensure.hpp` --
that's what gets used. Otherwise `WPT.hpp` falls back to plain `assert()`, so it works standalone
either way.

## License

Apache License 2.0 -- see [LICENSE](LICENSE).
