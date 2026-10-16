# Version 4 verification materials

See the root README.txt and `revision_audit_v4.md` for the current status.
The eight v3 numerical/algebraic JSON files and their engines are unchanged.
New files are `critical_norm_exact.json`, `critical_norm_primary.json`, and
`critical_norm_audit.json`. Full mathematical proofs are in `proofs/`.

The standard-library runner adds the independent critical-norm reconstruction
to the former nine jobs. It regenerates nine JSON files in an empty temporary
directory and compares their complete contents with the saved records.
The primary runner executes 15 jobs using SymPy/mpmath; its new two records
are compared exactly as well. It does not need Bessel quadrature.

The proof of the modular bridge is in `proofs/critical_norm_bridge.md`,
with the complete rational eta vectors in `proofs/critical_norm_certificate.md`.
No Bessel or critical-period value initializes the new proof checks.

The logs shipped in this version are current v4 replays. Earlier revision
narratives are clearly segregated under `history/`; an earlier statement
that (3.5) is open is superseded only by the new proof, not by numerical data.
