//! rustyomestats — genome statistics, codon density, and U50-family
//! assembly metrics.
//!
//! Modules:
//! * [`io_utils`] — FASTA discovery / loading
//! * [`stats`]    — length / GC / N-L stats (polars output)
//! * [`codon`]    — 6-frame translation + absolute/predicted codon density
//! * [`fgs`]      — FragGeneScanRs subprocess wrapper for predicted ORFs
//! * [`u50`]      — Castro et al. (2016) N50/U50 assembly metrics from a
//!                  reference FASTA + mapped-contigs BED

pub mod codon;
pub mod fgs;
pub mod io_utils;
pub mod stats;
pub mod u50;
pub mod visualization;
