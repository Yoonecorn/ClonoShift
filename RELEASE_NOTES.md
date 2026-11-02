# ClonoShift 1.0.0 source release

This repository provides the ClonoShift 1.0.0 source release, public documentation, manuscript PDF, and scientific figures. It does not include a prebuilt ClonoShift executable.

Build from source on the qualified platform, macOS on Apple Silicon. ClonoShift is distributed under GPL-3.0-or-later.

The `analyze` command runs ClonoShift Extended and ClonoShift Core separately, with descriptive comparison and explicit model statuses. It preserves a valid result when the other model refuses or fails, and publishes outputs using an atomic no-replace operation. Extended uses a qualified Core fit in its conditional formulation; Core is a NoisET-compatible Rust implementation. See PROVENANCE.md for software and scientific attribution.

Neither model is a calibrated biological caller. Synthetic fixtures are illustrative, and the repository contains no real biological data. Research interfaces retain their prior scope and are not promoted.
