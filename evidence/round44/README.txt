Round 44: GPT-6.1 Sol (Codex), dna-firehorse, 2026-11-09

Main 3935eeac uses the existing ten-term Bernoulli correction for GAMMALN
above 8, with stored PC64 operations. The lower rational tables remain open.
Bank: 3,887,877 -> 3,887,881 / 3,910,773 (+4 net). Corpus: unchanged
75,050/76,957, inclusive failures 1,907. measurement-bindings.json names the
numeric epochs; subsequent UI/font changes are separate.

files.json lists immutable payloads with original and compressed hashes.
Original questions, frozen predictions and original replies are retained.
The early direct client did not request input echoes: its echo slots are
null and original socket frames were not saved. Do not describe them as
having complete echo/frame provenance. The later audited-cache envelopes
retain original frames, exact echoes and full environment definitions.
fingerprint.accept validates all six envelopes. The 2,048 audited random
arguments are new; 62 audited controls (2 new, 60 retained) and four further
order controls (2 new, 2 retained) match. The last set includes the input
that separates both nine and eleven terms from the supported ten terms.

The earlier second frozen set has 9,312 new arguments and 900 retained
boundary controls: 9,875 -> 10,212/10,212 total, 8,980 -> 9,312 new-only.
The registered application route matches all 10,212. Individual Horner
stores and higher coefficient last bits remain undetermined; 24 selected
questions support the leading 1/12 load. No fitted coefficients were added.

Replay from the Research checkout, using new writable output paths:
  python3 racer/x87_programs.py build /path/native.so
  python3 racer/gamma_stages.py ../DnaCalc /path/proof.json \
    --cache evidence/round44/retained/audited-cache --library /path/native.so \
    --gammaln-stirling-proof --count 2048 --seed 44005 --offline
The current helper supports explicit extra controls and refuses uncached
offline questions. Acquisition scripts in retained/ preserve their original
scratch paths; frozen numeric payloads and producer snapshots are primary.
normal_quantile.py subtraction-lattice replays the normal constraint using
independent rational inverse rounding. Beta and YIELD changes are diagnostic
only; no corresponding application kernel was changed. See design/.

Checks: main dependent 2,265/0; final light 705/0. Research Rust 154/0,
Python 162/0; current unit and discrete-disabled compilation pass. UBS has
zero critical findings; Python SQL-flow analysis remains partial because
the scanner does not support a pre-existing list-aliasing pattern.

Large corpus outputs are outside Git: dnacalc-evidence/research/round44/
20261109T083900Z-dna-firehorse, 2 files, 17,880,250 bytes, checksum copy and
one-way check with 0 differences. r2-evidence.json gives object hashes.
Excel waits above three seconds: none (waits.tsv is empty). Excel was used
only through the question service; no Windows or browser sweep was run.
