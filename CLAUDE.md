# Distributed GHZ Project

## Project Context

See [context.md](./context.md) for the full scope and design of this project (distributed quantum SE via Bell and GHZ-based measurements using stim/sinter).

## Reference Project

The project `~/Projects/long-range-parity` has the same general structure we're following here:

- `src/` — library code (noise models, circuit utilities, decoders)
- `tools/` — standalone scripts (decoding, circuit generation, scanning)
- `notebooks/` — analysis and visualization
- `in/` — input circuits
- `out/` — output data
- `pyproject.toml` — packaging with optional dependency groups

Use it as a reference for conventions, project layout, and the decode/sinter pipeline.
