# LAIMethodComparison

Comparison of metrics to estimate LAI from UAV- and ground-based sampling methods.

# Hemispherical camera lens calibration
To derive measures such as LAI from hemispherical photos, programs such as Hemispher require some calibration to set how image space relates to object space (Schwalbe, 2005).
The relationship may be approximated through a polynomial of the form:
$y = \sum_{k=1}^{p} \beta_k \theta^k$
where $\theta$ is the angle (in radians) between the camera's optical axis and a point in the captured scene, and $y$ is the radial distance of the corresponding image point from the image centre, expressed as a proportion of the radial distance at a specified maximum angle. [./data/camera_calibration](https://github.com/ManuelPopp/LAIMethodComparison/tree/main/data/camera_calibration) contains files and a script used to estimate the parameters for the lenses we used in this study. The coefficients of the polynomial model can be entered into software such as [Hemisfer](https://www.schleppi.ch/patrick/hemisfer/) (Schleppi, 2007).

# References
Schwalbe, E., 2005, February. Geometric modelling and calibration of fisheye lens camera systems. *Proc. 2nd Panoramic Photogrammetry Workshop, Int. Archives of Photogrammetry and Remote Sensing.* 36, 5. [url](https://www.isprs.org/proceedings/xxxvi/5-w8/Paper/PanoWS_Berlin2005_Schwalbe.pdf)
Schleppi, P., Conedera, M., Sedivy, I., Thimonier, A., 2007. Correcting non-linearity and slope effects in the estimation of the leaf area index of forests from hemispherical photographs. *Agricultural and Forest Meteorology* 144, 236–242. doi: [10.1016/j.agrformet.2007.02.004](https://doi.org/10.1016/j.agrformet.2007.02.004)

# Project internal usage notes
To compute UAV-based metrics, log into one of the GPU nodes, pull the latest version of this repo, and run
```bash
./git/LAIMethodComparison/shell/run_ms.sh &
```
or the LiDAR version, respectively. Note that the LiDAR version needs to finish first, since the multispectral metrics partly rely on the LiDAR-derived tree cover data.
