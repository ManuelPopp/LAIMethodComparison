# LAIMethodComparison

Comparison of metrics to estimate LAI from UAV- and ground-based sampling methods.

# Hemispherical camera lens calibration
Hemisfer requires some calibration parameters to estimate LAI from hemispherical photos.
The lens parameters can be estimated through a polynomial of the form:
$y = \sum_{k=1}^{p} \beta_k \theta^k$
where $\theta$ is the angle (in radians) between the camera's optical axis and a point in the captured scene, and $y$ is the radial distance of the corresponding image point from the image centre, expressed as a proportion of the radial distance at a specified maximum angle.

# Project internal usage notes
To compute UAV-based metrics, log into one of the GPU nodes, pull the latest version of this repo, and run
```bash
./git/LAIMethodComparison/shell/run_ms.sh &
```
or the LiDAR version, respectively. Note that the LiDAR version needs to finish first, since the multispectral metrics partly rely on the LiDAR-derived tree cover data.
