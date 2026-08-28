# GM12878 chromosome 1 N=400 regression fixture

This directory contains a processed experimental Hi-C contact map for the
10 Mb `chr1:31,000,000-41,000,000` region at 25 kb resolution and a converged
noise-aware COV reference Gram matrix. The contact data originate from
ENCODE accession `ENCFF065LSP`; `metadata.json` records the extraction,
preprocessing, hashes, model parameters, and reference result.

The ordinary CPU test suite loads the full N=400 data, reconstructs the target
squared-distance observations, and independently verifies both the COV
objective and the dual-eliminated KKT residual of the reference solution. This
is intentionally a scientific golden-result regression rather than a smaller
synthetic proxy.

On a machine with CuPy and an accessible CUDA GPU, the test marked `real_data`
additionally runs the complete public contact-map workflow with the default
hybrid optimizer: Rouse initialization, PDHG to a relative KKT residual of
`1e-2`, then direct-Gram monotone FISTA to the final `1e-5` certificate. The
PDHG phase is the canonical variance-whitened, inverse-free implementation;
FISTA uses its physical Gram matrix directly without recalibration.

```bash
pytest -q -m real_data tests/test_covariance_pdhg.py
```

On an NVIDIA GeForce RTX 4070, the direct-Gram `1e-2` hybrid reached an
independent KKT residual of `9.972e-6` with a 27.046-second median wall time;
continued PDHG required 34.364 seconds in the matched experiment. Relative
errors against the stored solution were `9.230e-10` for the objective,
`2.382e-5` for the Gram matrix, and `1.609e-5` for connectivity. Both paths
reserved 52.50 MiB in the CuPy memory pool. On CPU-only systems the end-to-end
test is skipped, while the fast full-size objective/KKT regression still runs.

Historical provenance: the previous public hybrid used a `1e-3` PDHG handoff
and Newton-CG refinement. Those measurements describe the removed optimizer
and are not the current workflow.
