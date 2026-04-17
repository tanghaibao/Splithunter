// Shannon-like trinucleotide entropy, matching the C++ implementation in
// Splithunter.cpp.  The result is scaled to [0, 100] so `< min_entropy` in
// callers retains the same semantics as the legacy tool (MINENT = 50).

use std::collections::HashMap;

pub fn trinucleotide_entropy(seq: &[u8]) -> f64 {
    if seq.len() < 3 {
        return 0.0;
    }
    let l = seq.len() - 2;
    let k = l.min(64);
    if k <= 1 {
        return 0.0;
    }
    let log_k = (k as f64).ln();
    if log_k == 0.0 {
        return 0.0;
    }

    let mut counts: HashMap<[u8; 3], u32> = HashMap::with_capacity(l);
    for i in 0..l {
        let tri = [seq[i], seq[i + 1], seq[i + 2]];
        *counts.entry(tri).or_insert(0) += 1;
    }

    let l_f = l as f64;
    let mut h = 0.0;
    for &c in counts.values() {
        let f = c as f64 / l_f;
        h += f * f.ln() / log_k;
    }
    -h * 100.0
}
