# Scientific and software provenance

ClonoShift Extended is the ClonoShift-developed corrected workflow. ClonoShift Core is a NoisET-compatible Rust reimplementation/port of the bounded Poisson fit and response formulation. The two models remain separate results, and neither is a calibrated biological directional caller. ClonoShift is written in Rust and is distributed under GPL-3.0-or-later; see LICENSE.

NoisET is the upstream scientific and software reference for Core. The compatibility target is NoisET commit `2b6a45fcd26e7f9eefb3004d81725355cf6d66c9`, available at <https://github.com/statbiophys/NoisET>. The upstream source header identifies copyright © 2021 Meriem Bensouda Koraichi and permits redistribution and modification under GNU GPL version 3 or, at the recipient's option, any later version. The upstream repository's license file contains GPL version 3. This notice records attribution and license information; it does not make a legal determination about historical derivation.

Please cite both NoisET and the scientific method:

- Meriem Bensouda Koraichi, Maximilian Puelma Touzel, Andrea Mazzolini, Thierry Mora and Aleksandra M. Walczak, “NoisET: Noise learning and Expansion detection of T-cell receptors,” arXiv:2102.03568, <https://doi.org/10.48550/arXiv.2102.03568>.
- Maximilian Puelma Touzel, Aleksandra M. Walczak and Thierry Mora, “Inferring the immune response from repertoire sequencing,” PLOS Computational Biology (2020), <https://doi.org/10.1371/journal.pcbi.1007873>.

This source package contains no NoisET Python source files or upstream-invoking oracle scripts. It contains no real TCR or BCR datasets or precompiled binaries. Dataset acquisition guidance and third-party dependency notices are in docs/DATA_SOURCES.md and THIRD_PARTY_NOTICES.md. ClonoShift 1.0.0 is the finalized source release. The release includes the qualified source package, public documentation, manuscript PDF, and scientific figures. Release status does not imply journal publication, peer review, manuscript acceptance, broader platform qualification, or a calibrated biological consensus call.
