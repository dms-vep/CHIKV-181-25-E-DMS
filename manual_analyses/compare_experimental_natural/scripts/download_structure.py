"""Download a file (the 6NK6 biological assembly) from a URL, failing fast on error.

Run via the Snakemake `script:` directive; reads `snakemake.params.url` and writes
`snakemake.output[0]`.
"""

import sys

import requests

snakemake = snakemake  # noqa: F821  (injected by the Snakemake script directive)
sys.stdout = sys.stderr = open(snakemake.log[0], "w")  # capture all output to the log

url = snakemake.params.url
out = snakemake.output[0]

response = requests.get(url, timeout=120)
response.raise_for_status()
if not response.content:
    raise ValueError(f"empty download from {url}")
with open(out, "wb") as f:
    f.write(response.content)
print(f"downloaded {len(response.content)} bytes from {url}")
