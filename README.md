# LAIMethodComparison

Comparison of metrics to estimate LAI from UAV- and ground-based sampling methods.

# Project internal usage notes
To compute UAV-based metrics, log into one of the GPU nodes, pull the latest version of this repo, and run
```bash
./git/LAIMethodComparison/shell/run_ms.sh &
```
or the LiDAR version, respectively. Note that the LiDAR version needs to finish first, since the multispectral metrics partly rely on the LiDAR-derived tree cover data.
