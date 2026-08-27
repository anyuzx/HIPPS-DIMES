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
`1e-3`, then Newton-CG to the final `1e-5` certificate.

```bash
pytest -q -m real_data tests/test_covariance_pdhg.py
```

The standalone PDHG reference took about 115 seconds on an NVIDIA GeForce RTX
4070; the hybrid regression took about 70 seconds on the same GPU. On CPU-only
systems the end-to-end test is skipped, while the fast full-size objective/KKT
regression still runs.
