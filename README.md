# WPT.hpp

A full Wavelet Packet Transform: given a mother wavelet's decomposition filter, recursively
splits both the low- and high-pass sub-bands at every level -- unlike a plain DWT, which only
recurses on the low/approximation branch. `push()` is per-sample and online (not block-based):
each tree node keeps a small circular buffer and alternates which of two physical sub-levels
receives a given push, realizing decimation-by-2 without buffering a block first.

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

This is forward decomposition only -- there is no inverse/reconstruction transform. `Wavelet`
derives the reconstruction (synthesis) filter pair alongside the decomposition pair, but nothing
in this header currently consumes it.

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
