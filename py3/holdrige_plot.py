#!/usr/bin/env python3
"""
Holdridge life zones (HLZ) for the LAI method comparison plots.
Usage
    python holdridge_plots.py                    # Classical HLZ, table + figure
    python holdridge_plots.py --variant Penman-Monteith
    python holdridge_plots.py --no-plot          # table only
    python holdridge_plots.py --root /path/to/LAIMethodComparison
"""

from __future__ import annotations

import argparse
import io
from pathlib import Path

import numpy as np
import pandas as pd
import pygmt
import pooch
import rasterio
from rasterio.warp import transform as warp_transform

# ----------------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------------
PROJECT_DIRNAME = "LAIMethodComparison"
PLOT_CSV_REL = Path("data") / "plot_info.csv"
HLZ_DIR_REL = Path("data") / "holdridge"
FIG_DIR_REL = Path("figures")
FIG_STEM = "holdridge_life_zones"

# figshare files (article 28236242, version 1) with their published MD5 sums
FIGSHARE = {
    "Classical": {
        "tif": ("HLZ_Classical.tif",
                "https://ndownloader.figshare.com/files/52378187",
                "md5:14aa251e651730fd20d20f8350831623"),
        "codes": ("HLZ_Codes_Classical.csv",
                  "https://ndownloader.figshare.com/files/52378184",
                  "md5:4cc5c156a51091637a5235b0c336711b"),
    },
    "Penman-Monteith": {
        "tif": ("HLZ_Penman-Monteith.tif",
                "https://ndownloader.figshare.com/files/52378190",
                "md5:af5c4173abe37a9dd68b0678822645e7"),
        "codes": ("HLZ_Codes_Penman-Monteith.csv",
                  "https://ndownloader.figshare.com/files/52378181",
                  "md5:89de47d36b1cdb1abb791df2d42a0b1a"),
    },
}

# Holdridge hexagon definitions (hlz_defs.csv from the author's repository).
# Row index (1-based) == veg_class in the HLZ code. Columns: lower bound of
# annual biotemperature (degC), annual precipitation (mm), PET ratio, name.
HLZ_DEFS_CSV = """abt,tap,per,veg_class
0,0,0,none
0.75,62.5,0.5,desert
0.75,125,0.25,desert
0.75,250,0.125,desert
1.5,62.5,1,dry tundra
1.5,125,0.5,moist tundra
1.5,250,0.25,wet tundra
1.5,500,0.125,rain tundra
3,62.5,2,desert
3,125,1,dry scrub
3,250,0.5,moist forest
3,500,0.25,wet forest
3,1000,0.125,rain forest
6,62.5,4,desert
6,125,2,desert scrub
6,250,1,steppe
6,500,0.5,moist forest
6,1000,0.25,wet forest
6,2000,0.125,rain forest
12,62.5,8,desert
12,125,4,desert scrub
12,250,2,thorn steppe
12,500,1,dry forest
12,1000,0.5,moist forest
12,2000,0.25,wet forest
12,4000,0.125,rain forest
24,62.5,16,desert
24,125,8,desert scrub
24,250,4,thorn woodland
24,500,2,very dry forest
24,1000,1,dry forest
24,2000,0.5,moist forest
24,4000,0.25,wet forest
24,8000,0.125,rain forest
12,62.5,8,desert
12,125,4,desert scrub
12,250,2,thorn woodland
12,500,1,dry forest
12,1000,0.5,moist forest
12,2000,0.25,wet forest
12,4000,0.125,rain forest
"""
HLZ_DEFS = pd.read_csv(io.StringIO(HLZ_DEFS_CSV))
HLZ_DEFS.index = np.arange(1, len(HLZ_DEFS) + 1)  # veg_class is 1-based

# Plot styling (same look as the Biome.jl GMT figures)
WATER_RGB = "100/149/190"
LAND_OUTSIDE = "232/232/228"
OUTSIDE_HLZ_RGB = (0.86, 0.86, 0.86)     # code 1: ice / no vegetation
HLZ_TRANSPARENCY = 25                    # percent, HLZ over hillshade
GLOBAL_STEP_DEG = 0.1                    # decimated global grid spacing
FIG_WIDTH_CM = 18.0
INSET_GAP_CM = 1.2                       # room for the latitude annotations
NEIGHBOR_DEG = 5.0                       # country polygons kept if this close to plots
PAD_FRAC = 0.04                          # inset padding around the country
TRIANGLE_WIDTH_CM = FIG_WIDTH_CM         # Holdridge chart incl. side labels
SITE_DOT = "c0.2c"                       # plot dots on the insets
SITE_FONT = "6p,Helvetica-Bold,black"
LABEL_SIDE = {}                          # optional override, e.g. {"Labi": "left"}

# Holdridge chart color scheme
# Hue: perceptually ordered 'turbo' colormap (matplotlib), red at the highest
# PET ratio through yellow and green to blue at the lowest PET ratio.
# Opacity over white grows with precipitation, and the wettest hexagons are
# also darkened, so neighbors along every axis of the chart differ clearly.
COLORMAP = "turbo"
CMAP_RANGE = (0.96, 0.10)                # colormap position at PET ratio 16 .. 0.125
ALPHA_RANGE = (0.14, 1.00)               # opacity at 62.5 mm .. 8000 mm precipitation
ALPHA_GAMMA = 0.85
DARKEN_WET = 0.22                        # extra darkening of the wettest hexagons

# Holdridge chart guide lines
PRECIP_PEN = "0.4p,blue3,1_1.5"
PET_PEN = "0.4p,red3,3_1.5"
BIOTEMP_PEN = "0.35p,darkgreen,1_2"


# ----------------------------------------------------------------------------
# Paths and data access
# ----------------------------------------------------------------------------
def find_project_root(explicit: str | None) -> Path:
    """Return the LAIMethodComparison directory (explicit > script parents > cwd)."""
    if explicit:
        return Path(explicit).expanduser().resolve()
    candidates = [Path(__file__).resolve(), Path.cwd().resolve()]
    for start in candidates:
        for p in [start, *start.parents]:
            if p.name.lower() == PROJECT_DIRNAME.lower():
                return p
    for start in candidates:
        for p in [start, *start.parents]:
            if (p / PLOT_CSV_REL).is_file():
                return p
    raise FileNotFoundError(
        f"Could not locate {PROJECT_DIRNAME}/{PLOT_CSV_REL}; pass --root explicitly.")


def fetch_hlz(root: Path, variant: str) -> tuple[Path, Path]:
    """Download (or reuse) the HLZ GeoTIFF and code table, MD5-checked by pooch."""
    out_dir = root / HLZ_DIR_REL
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    for key in ("tif", "codes"):
        fname, url, md5 = FIGSHARE[variant][key]
        paths.append(Path(pooch.retrieve(url=url, known_hash=md5, fname=fname,
                                         path=out_dir, progressbar=True)))
    return paths[0], paths[1]


def read_plot_table(csv_path: Path) -> pd.DataFrame:
    """Read plot_info.csv, which pads fields with spaces around the delimiters."""
    df = pd.read_csv(csv_path, skipinitialspace=True, encoding="utf-8")
    df.columns = df.columns.str.strip()
    for col in df.columns:
        if df[col].dtype == object or pd.api.types.is_string_dtype(df[col]):
            df[col] = df[col].astype("string").str.strip()
    df["Latitude"] = pd.to_numeric(df["Latitude"])
    df["Longitude"] = pd.to_numeric(df["Longitude"])
    return df


def read_code_table(codes_csv: Path) -> pd.Series:
    t = pd.read_csv(codes_csv)
    t.columns = t.columns.str.strip()
    return pd.Series(t["Name"].astype(str).values, index=t["Code"].astype(int).values)


# ----------------------------------------------------------------------------
# Step 2: sample HLZ codes at plot coordinates
# ----------------------------------------------------------------------------
def sample_hlz(tif: Path, lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
    with rasterio.open(tif) as src:
        xs, ys = lon, lat
        if src.crs is not None and src.crs.to_epsg() != 4326:
            xs, ys = warp_transform("EPSG:4326", src.crs, lon, lat)
        vals = np.array([v[0] for v in src.sample(zip(xs, ys), indexes=1)])
    return vals.astype(np.int64)


def add_holdridge_columns(csv_path: Path, tif: Path, names: pd.Series) -> pd.DataFrame:
    df = read_plot_table(csv_path)
    codes = sample_hlz(tif, df["Longitude"].to_numpy(), df["Latitude"].to_numpy())
    df["Holdridge_code"] = codes
    df["Holdridge_life_zone"] = [names.get(int(c), f"Unknown code {c}") for c in codes]
    df.to_csv(csv_path, index=False, encoding="utf-8")
    return df


# ----------------------------------------------------------------------------
# Holdridge hexagons: geometry, colors, CPT
# ----------------------------------------------------------------------------
def decode(code: int) -> tuple[int, int, int, int]:
    return code // 1000, (code % 1000) // 100, (code % 100) // 10, code % 10


def build_hexagons() -> tuple[pd.DataFrame, np.ndarray]:
    """One row per hexagon of the Holdridge chart (unique tap/per pair), and a
    lookup veg_class -> hexagon id (0 = outside HLZ / none).

    Chart coordinates (unit hexagon spacing):
        u = log2(tap / 62.5)     precipitation step, 0..7
        v = log2(per / 0.125)    PET-ratio step, 0..7
        row = u + v - 2          0 polar .. 5 tropical (biotemperature belt)
        x = (u - v) / 2,  y = -row * sqrt(3) / 2
    """
    defs = HLZ_DEFS.iloc[1:].copy()
    defs["u"] = np.round(np.log2(defs["tap"] / 62.5)).astype(int)
    defs["v"] = np.round(np.log2(defs["per"] / 0.125)).astype(int)
    hexes = []
    veg_to_hex = np.zeros(len(HLZ_DEFS) + 1, dtype=np.int64)
    for (u, v), grp in defs.groupby(["u", "v"], sort=False):
        names = list(dict.fromkeys(grp["veg_class"]))
        name = names[0] if len(names) == 1 else " / ".join(
            [names[0]] + [n.split()[-1] for n in names[1:]])      # "thorn steppe / woodland"
        hid = len(hexes) + 1
        hexes.append(dict(hex_id=hid, u=u, v=v, tap=grp["tap"].iloc[0],
                          per=grp["per"].iloc[0], name=name, row=u + v - 2,
                          x=(u - v) / 2.0, y=-(u + v - 2) * np.sqrt(3) / 2.0))
        veg_to_hex[grp.index.to_numpy()] = hid
    return pd.DataFrame(hexes).set_index("hex_id"), veg_to_hex


HEXES, VEG_TO_HEX = build_hexagons()


def hex_rgb(hex_id: int) -> tuple[float, float, float]:
    """Hue from PET ratio (red -> yellow -> green -> blue, 'turbo'), opacity and
    depth from precipitation (dry hexagons near transparent, wet ones saturated)."""
    from matplotlib import colormaps

    if hex_id == 0:
        return OUTSIDE_HLZ_RGB
    h = HEXES.loc[hex_id]
    humid = np.clip((4.0 - np.log2(h["per"])) / 7.0, 0.0, 1.0)  # 0 superarid .. 1 superhumid
    wet = np.clip(h["u"] / 7.0, 0.0, 1.0)                        # 0 = 62.5 mm .. 1 = 8000 mm
    pos = CMAP_RANGE[0] + (CMAP_RANGE[1] - CMAP_RANGE[0]) * humid
    pure = np.array(colormaps[COLORMAP](pos)[:3])
    alpha = ALPHA_RANGE[0] + (ALPHA_RANGE[1] - ALPHA_RANGE[0]) * wet ** ALPHA_GAMMA
    rgb = (1.0 - alpha * (1.0 - pure)) * (1.0 - DARKEN_WET * wet ** 2)  # "transparency" over white
    return tuple(np.clip(rgb, 0.0, 1.0))


def text_color(rgb) -> str:
    """Black or white label, whichever reads better on the fill."""
    lum = 0.2126 * rgb[0] + 0.7152 * rgb[1] + 0.0722 * rgb[2]
    return "white" if lum < 0.42 else "gray5"


def gmt_rgb(rgb) -> str:
    return "/".join(str(int(round(255 * c))) for c in rgb)


def codes_to_hex(codes: np.ndarray) -> np.ndarray:
    """Raster codes -> CPT key: hexagon id (1..33), 0 for code 1, NaN for no data."""
    codes = codes.astype(np.int64)
    veg = np.clip(codes // 1000, 0, len(VEG_TO_HEX) - 1)
    keys = VEG_TO_HEX[veg].astype(np.float32)
    keys[codes == 1] = 0
    keys[codes == 0] = np.nan
    return keys


def write_hex_cpt(path: Path) -> Path:
    lines = ["# COLOR_MODEL = RGB", f"0\t{gmt_rgb(OUTSIDE_HLZ_RGB)}\t;outside HLZ parameters"]
    for hid, h in HEXES.iterrows():
        lines.append(f"{hid}\t{gmt_rgb(hex_rgb(hid))}\t;{h['name']}")
    lines += ["B\t255/255/255", "F\t255/255/255", "N\t255/255/255"]
    path.write_text("\n".join(lines) + "\n")
    return path


# ----------------------------------------------------------------------------
# Rasters as xarray grids for PyGMT
# ----------------------------------------------------------------------------
def to_grid(values: np.ndarray, west: float, north: float, dx: float, dy: float):
    import xarray as xr

    ny, nx = values.shape
    lon = west + (np.arange(nx) + 0.5) * dx
    lat = north - (np.arange(ny) + 0.5) * dy
    da = xr.DataArray(values[::-1, :], coords={"lat": lat[::-1], "lon": lon},
                      dims=("lat", "lon"))
    da.gmt.registration = 1   # pixel
    da.gmt.gtype = 1          # geographic
    return da


def read_window(tif: Path, region):
    import rasterio
    from rasterio.windows import from_bounds

    w, e, s, n = region
    with rasterio.open(tif) as src:
        win = from_bounds(w, s, e, n, src.transform).round_offsets().round_lengths()
        arr = src.read(1, window=win, boundless=True, fill_value=0)
        t = src.window_transform(win)
    return arr, t.c, t.f, t.a, -t.e


def read_global(tif: Path, cache_dir: Path, step: float = GLOBAL_STEP_DEG):
    """Mode-decimated global grid (categorical-safe), cached as .npz."""
    import rasterio
    from rasterio.enums import Resampling

    cache = cache_dir / f"{tif.stem}_mode_{step:g}deg.npz"
    if cache.is_file():
        z = np.load(cache)
        return z["arr"], *z["meta"]
    with rasterio.open(tif) as src:
        b = src.bounds
        nx = int(round((b.right - b.left) / step))
        ny = int(round((b.top - b.bottom) / step))
        arr = src.read(1, out_shape=(ny, nx), resampling=Resampling.mode)
        meta = np.array([b.left, b.top, (b.right - b.left) / nx, (b.top - b.bottom) / ny])
    np.savez_compressed(cache, arr=arr, meta=meta)
    return arr, *meta


# ----------------------------------------------------------------------------
# Country geometry from GMT's Digital Chart of the World
# ----------------------------------------------------------------------------
def iso2(country: str) -> str:
    import pycountry

    try:
        return pycountry.countries.lookup(country).alpha_2
    except LookupError:
        return pycountry.countries.search_fuzzy(country)[0].alpha_2


def country_region(code: str, lon: np.ndarray, lat: np.ndarray, tmp: Path):
    """Bounding box of the DCW polygons of `code` lying near the plots
    (drops remote islands such as the Galapagos)."""
    import pygmt

    out = tmp / f"dcw_{code}.txt"
    with pygmt.clib.Session() as lib:
        lib.call_module("coast", f"-E{code} -M ->{out}")
    segs, cur = [], []
    for line in out.read_text().splitlines():
        if line.startswith(">") or not line.strip():
            if cur:
                segs.append(np.array(cur))
            cur = []
        elif not line.startswith("#"):
            cur.append([float(v) for v in line.split()[:2]])
    if cur:
        segs.append(np.array(cur))

    pw, pe = lon.min() - NEIGHBOR_DEG, lon.max() + NEIGHBOR_DEG
    ps, pn = lat.min() - NEIGHBOR_DEG, lat.max() + NEIGHBOR_DEG
    keep = [s for s in segs
            if s[:, 0].max() >= pw and s[:, 0].min() <= pe
            and s[:, 1].max() >= ps and s[:, 1].min() <= pn]
    pts = np.vstack(keep) if keep else np.column_stack([lon, lat])
    w, e = min(pts[:, 0].min(), lon.min()), max(pts[:, 0].max(), lon.max())
    s, n = min(pts[:, 1].min(), lat.min()), max(pts[:, 1].max(), lat.max())
    pad = PAD_FRAC * max(e - w, n - s) + 0.05
    return [w - pad, e + pad, s - pad, n + pad]


def inset_aspect(region) -> float:
    """Height/width of a -JQ frame whose standard parallel is the region center."""
    w, e, s, n = region
    return (n - s) / ((e - w) * np.cos(np.deg2rad(0.5 * (s + n))))


def match_aspect(region, target: float):
    """Grow a region in lon or lat (about its center) to the target aspect."""
    w, e, s, n = region
    coslat = np.cos(np.deg2rad(0.5 * (s + n)))
    a = inset_aspect(region)
    if a < target:
        half = 0.5 * target * (e - w) * coslat
        c = 0.5 * (s + n)
        s, n = c - half, c + half
    else:
        half = 0.5 * (n - s) / (target * coslat)
        c = 0.5 * (w + e)
        w, e = c - half, c + half
    return [w, e, max(s, -89.9), min(n, 89.9)]


def nice_interval(span: float, n_target: int = 3) -> float:
    raw = span / n_target
    steps = np.array([0.1, 0.25, 0.5, 1, 2, 2.5, 5, 10, 20])
    return float(steps[np.argmin(np.abs(np.log(steps / raw)))])


# ----------------------------------------------------------------------------
# Plotting
# ----------------------------------------------------------------------------
def hillshade(region, resolution):
    import pygmt

    relief = pygmt.datasets.load_earth_relief(resolution=resolution, region=region,
                                              registration="pixel")
    shade = pygmt.grdgradient(grid=relief, azimuth=315, normalize="e0.6")
    return shade


def relief_resolution(region) -> str:
    span = max(region[1] - region[0], region[3] - region[2])
    return "15s" if span < 3 else "30s" if span < 15 else "01m"


def draw_hlz_panel(fig, region, projection, hlz_grid, cpt, shade, clip_dcw=None,
                   frame=("WSne", "af")):
    fig.basemap(region=region, projection=projection, frame=list(frame))
    fig.coast(land=LAND_OUTSIDE, water=WATER_RGB, resolution="h" if clip_dcw else "l")
    if clip_dcw:
        fig.coast(dcw=f"{clip_dcw}+c")
    fig.grdimage(grid=shade, cmap="gray")
    fig.grdimage(grid=hlz_grid, cmap=str(cpt), nan_transparent=True,
                 transparency=HLZ_TRANSPARENCY)
    if clip_dcw:
        fig.coast(Q=True)
    fig.coast(water=WATER_RGB, shorelines="0.2p,gray35",
              borders="1/0.35p,gray35", resolution="h" if clip_dcw else "l")
    if clip_dcw:
        fig.coast(dcw=f"{clip_dcw}+p0.6p,black")
    fig.basemap(frame=list(frame))      # redraw frame on top of the raster


def site_hex(df) -> pd.Series:
    """Hexagon id of every plot (0 = outside HLZ, -1 = no data)."""
    return pd.Series(np.nan_to_num(codes_to_hex(df["Holdridge_code"].to_numpy()), nan=-1)
                     .astype(int), index=df.index)


def plot_points(fig, df, size=SITE_DOT):
    """Plots as dots filled with the color of their Holdridge life zone."""
    hexes = site_hex(df)
    for hid, sub in df.groupby(hexes, sort=False):
        fill = gmt_rgb(hex_rgb(hid)) if hid >= 0 else "white"
        fig.plot(x=sub["Longitude"], y=sub["Latitude"], style=size, fill=fill,
                 pen="0.5p,black")


def label_sites(fig, df):
    """One name per field site, placed beside its cluster of plots. Sites west of
    the country's median site longitude are labeled on the left, the others on
    the right (override with LABEL_SIDE)."""
    centers = df.groupby("Location", sort=False).agg(
        lon=("Longitude", "mean"), lat=("Latitude", "mean"),
        lon_min=("Longitude", "min"), lon_max=("Longitude", "max"))
    median_lon = centers["lon"].median()
    for loc, c in centers.iterrows():
        side = LABEL_SIDE.get(loc, "left" if c["lon"] < median_lon else "right")
        x, just, off = ((c["lon_min"], "MR", "-0.18c/0c") if side == "left"
                        else (c["lon_max"], "ML", "0.18c/0c"))
        fig.text(x=x, y=c["lat"], text=loc, font=SITE_FONT, justify=just, offset=off,
                 fill="white@25", clearance="0.04c/0.03c")


def geo_box(region, n=60):
    """Region outline densified along meridians/parallels (curves in Robinson)."""
    w, e, s, n_ = region
    xs = np.r_[np.linspace(w, e, n), np.full(n, e), np.linspace(e, w, n), np.full(n, w)]
    ys = np.r_[np.full(n, s), np.linspace(s, n_, n), np.full(n, n_), np.linspace(n_, s, n)]
    return xs, ys


def hexagon_xy(x0, y0, r=1.0 / np.sqrt(3)):
    ang = np.deg2rad(90 + 60 * np.arange(7))
    return x0 + r * np.cos(ang), y0 + r * np.sin(ang)


def uv_to_xy(u, v):
    """Chart coordinates (see build_hexagons) to plot units."""
    u, v = np.asarray(u, float), np.asarray(v, float)
    return (u - v) / 2.0, -(u + v - 2.0) * np.sqrt(3) / 2.0


def clip_polygon(xs, ys, clip_xy):
    """Sutherland-Hodgman clipping of a polygon by a convex polygon (CCW)."""
    pts = list(zip(xs, ys))
    cx = list(clip_xy)
    for (ax, ay), (bx, by) in zip(cx, cx[1:] + cx[:1]):
        inside = lambda p: (bx - ax) * (p[1] - ay) - (by - ay) * (p[0] - ax) >= 0
        out, prev = [], pts[-1]
        for cur in pts:
            if inside(cur):
                if not inside(prev):
                    out.append(_intersect(prev, cur, (ax, ay), (bx, by)))
                out.append(cur)
            elif inside(prev):
                out.append(_intersect(prev, cur, (ax, ay), (bx, by)))
            prev = cur
        pts = out
        if not pts:
            break
    pts = pts + pts[:1]
    return np.array([p[0] for p in pts]), np.array([p[1] for p in pts])


def _intersect(p1, p2, a, b):
    x1, y1 = p1
    x2, y2 = p2
    x3, y3 = a
    x4, y4 = b
    den = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
    t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / den
    return x1 + t * (x2 - x1), y1 + t * (y2 - y1)


def draw_holdridge_chart(fig, df, width_cm):
    """Classical Holdridge triangle with the map colors.

    Guide lines follow the classifier thresholds (hlz_classify.m): each hexagon
    spans one doubling of precipitation, PET ratio and biotemperature, so the
    threshold lines lie half a step from the hexagon centers:
        precipitation  u = k - 1/2  -> 62.5 * 2**k mm   (blue, right edge)
        PET ratio      v = k - 1/2  -> 0.125 * 2**k     (red, left edge)
        biotemperature rows r + 1/2 -> 1.5 * 2**r degC  (green, horizontal)
    Hexagons containing field plots are outlined and list the site names.
    """
    s32 = np.sqrt(3) / 2
    region = [-8.4, 9.9, -5.5 * s32 - 1.75, 3 * s32 + 1.05]
    x_lat, x_alt, x_bt = region[0] + 0.1, region[1] - 1.5, region[1] - 1.4
    scale = width_cm / (region[1] - region[0])
    fig.basemap(region=region, projection=f"x{scale}c", frame="+n")

    # outer triangle (apex: 62.5 mm / PET ratio 0.125; base: 24 degC belt bottom)
    apex, b_left, b_right = uv_to_xy(-0.5, -0.5), uv_to_xy(-0.5, 8.0), uv_to_xy(8.0, -0.5)
    tri = [tuple(b_left), tuple(b_right), tuple(apex)]                      # CCW

    # sites per hexagon (each site assigned to its most frequent hexagon)
    hexes = site_hex(df)
    occ = (pd.DataFrame({"Location": df["Location"], "hex": hexes}).query("hex > 0")
           .groupby("Location", sort=False)["hex"].agg(lambda s: s.mode().iloc[0]))
    locs_in_hex = {}
    for loc, hid in occ.items():
        locs_in_hex.setdefault(int(hid), []).append(loc)

    # hexagons, clipped to the triangle
    for hid, h in HEXES.iterrows():
        xs, ys = clip_polygon(*hexagon_xy(h["x"], h["y"]), tri)
        fig.plot(x=xs, y=ys, fill=gmt_rgb(hex_rgb(hid)), close=True, pen="0.3p,gray45")

    # biotemperature lines (horizontal) with labels on the far right
    alt_belts = ["alvar", "alpine", "subalpine", "montane",
                 "lower montane / premontane", "basal"]
    lat_regions = ["polar", "subpolar", "boreal", "cool temperate",
                   "warm temperate / subtropical", "tropical"]
    for r in range(5):
        y = -(r + 0.5) * s32
        fig.plot(x=[x_lat, x_bt - 0.05], y=[y, y], pen=BIOTEMP_PEN)
        fig.text(x=x_bt, y=y, text=f"{1.5 * 2 ** r:g} @.C", font="5p,Helvetica,darkgreen",
                 justify="LM")
    for r in range(6):
        y = -r * s32
        fig.text(x=x_lat, y=y, text=lat_regions[r], font="5.5p,Helvetica,gray10",
                 justify="LM")
        fig.text(x=x_alt, y=y, text=alt_belts[r], font="5.5p,Helvetica,gray10", justify="RM")

    # precipitation (blue) and PET ratio (red) lines, labels outside the edges
    ext = 0.45
    for k in range(9):
        c = k - 0.5
        x, y = uv_to_xy([c, c], [-0.5 - ext, 7.5 - c])                      # precipitation
        fig.plot(x=x, y=y, pen=PRECIP_PEN)
        fig.text(x=x[0], y=y[0], text=f"{62.5 * 2 ** k:g}", font="5p,Helvetica,blue3",
                 angle=60, justify="LM", offset="0.03c/0.03c")
        x, y = uv_to_xy([-0.5 - ext, 7.5 - c], [c, c])                      # PET ratio
        fig.plot(x=x, y=y, pen=PET_PEN)
        fig.text(x=x[0], y=y[0], text=f"{0.125 * 2 ** k:g}", font="5p,Helvetica,red3",
                 angle=-60, justify="RM", offset="-0.03c/0.03c")
    fig.plot(x=[p[0] for p in tri] + [tri[0][0]], y=[p[1] for p in tri] + [tri[0][1]],
             pen="0.8p,black")

    # edge titles
    mid_l = uv_to_xy(-0.5 - 1.75, 3.75)
    mid_r = uv_to_xy(3.75, -0.5 - 1.75)
    fig.text(x=mid_l[0], y=mid_l[1], text="Potential evapotranspiration ratio",
             font="6.5p,Helvetica-Bold,red3", angle=60, justify="CM")
    fig.text(x=mid_r[0], y=mid_r[1], text="Annual precipitation (mm)",
             font="6.5p,Helvetica-Bold,blue3", angle=-60, justify="CM")
    fig.text(x=region[1] - 0.2, y=-2.5 * s32, text="Biotemperature",
             font="6.5p,Helvetica-Bold,darkgreen", angle=90, justify="CM")

    # hexagon names, and site names in occupied hexagons
    for hid, h in HEXES.iterrows():
        fill_rgb = hex_rgb(hid)
        col = text_color(fill_rgb)
        words = h["name"].split(" ")
        lines = [" ".join(words[:-1]), words[-1]] if len(words) > 1 else [words[0]]
        sites = locs_in_hex.get(hid, [])
        y0 = h["y"] + (0.2 if sites else 0.0)
        for k, txt in enumerate(lines):
            dy = (0.075 if k == 0 else -0.075) if len(lines) == 2 else 0.0
            fig.text(x=h["x"], y=y0 + dy, text=txt, font=f"4p,Helvetica,{col}", justify="CM")
        if sites:
            xs, ys = clip_polygon(*hexagon_xy(h["x"], h["y"]), tri)
            fig.plot(x=xs, y=ys, pen="1p,black")
            for j, loc in enumerate(sites[:3]):
                fig.text(x=h["x"], y=h["y"] - 0.06 - 0.12 * j, text=loc,
                         font=f"3.6p,Helvetica-BoldOblique,{col}", justify="CM")

    # humidity provinces under the base
    provinces = [("super-", "arid"), ("per-", "arid"), ("arid",), ("semi-", "arid"),
                 ("sub-", "humid"), ("humid",), ("per-", "humid"), ("super-", "humid")]
    y_bot = -5.5 * s32 - 0.3
    for k, p in enumerate(provinces):
        for j, txt in enumerate(p):
            fig.text(x=-3.5 + k, y=y_bot - 0.17 * j, text=txt, font="5.5p,Helvetica,gray10",
                     justify="CM")
    fig.text(x=0, y=y_bot - 0.72, text="Humidity province", font="6.5p,Helvetica-Bold",
             justify="CM")
    fig.text(x=x_lat, y=1.0, text="Latitudinal region", font="6.5p,Helvetica-Bold",
             justify="LM")
    fig.text(x=x_alt, y=1.0, text="Altitudinal belt", font="6.5p,Helvetica-Bold",
             justify="RM")
    return (region[3] - region[2]) * scale


def make_figure(df, tif, names, root):
    fig_dir = root / FIG_DIR_REL
    fig_dir.mkdir(parents=True, exist_ok=True)
    work = root / HLZ_DIR_REL / "cache"
    work.mkdir(parents=True, exist_ok=True)

    countries = list(dict.fromkeys(df["Country"]))

    # Global grid
    g_arr, g_w, g_n, g_dx, g_dy = read_global(tif, work)

    # Country regions, all grown to one common frame aspect
    insets = []
    for c in countries:
        sub = df[df["Country"] == c]
        code = iso2(c)
        reg = country_region(code, sub["Longitude"].to_numpy(), sub["Latitude"].to_numpy(), work)
        insets.append(dict(country=c, iso=code, region=reg, sub=sub))
    target = float(np.median([inset_aspect(i["region"]) for i in insets]))
    for i in insets:
        i["region"] = match_aspect(i["region"], target)
        arr, w, n, dx, dy = read_window(tif, i["region"])
        i["keys"] = codes_to_hex(arr)
        i["geo"] = (w, n, dx, dy)

    cpt = write_hex_cpt(work / "hlz_hexagons.cpt")

    letters = "abcdefghij"
    n_in = len(insets)
    inset_w = (FIG_WIDTH_CM - (n_in - 1) * INSET_GAP_CM) / n_in
    inset_h = inset_w * target

    fig = pygmt.Figure()
    with pygmt.config(MAP_FRAME_TYPE="plain", MAP_FRAME_PEN="0.6p,black",
                      FONT_ANNOT_PRIMARY="6p,Helvetica,gray20", FONT_LABEL="7p,Helvetica",
                      MAP_TICK_LENGTH_PRIMARY="2p", MAP_ANNOT_OFFSET_PRIMARY="1.5p",
                      MAP_GRID_PEN_PRIMARY="0.2p,gray60", FORMAT_GEO_MAP="ddd:mmF",
                      PS_PAGE_COLOR="white"):
        # -- a: global map, Robinson
        g_region = "d"
        g_proj = f"N0/{FIG_WIDTH_CM}c"
        g_grid = to_grid(codes_to_hex(g_arr), g_w, g_n, g_dx, g_dy)
        with pygmt.config(FORMAT_GEO_MAP="dddF"):
            draw_hlz_panel(fig, g_region, g_proj, g_grid, cpt,
                           hillshade([-180, 180, -90, 90], "10m"),
                           frame=("WSne", "xg60", "yg30"))
        for k, i in enumerate(insets):
            bx, by = geo_box(i["region"])
            fig.plot(x=bx, y=by, pen="0.9p,black")
            fig.text(x=i["region"][1], y=i["region"][3], text=letters[k + 1],
                     font="8p,Helvetica-Bold,black", justify="BL", offset="0.06c/0.04c",
                     fill="white@20", clearance="0.03c/0.03c")
        plot_points(fig, df, size="c0.11c")
        fig.text(position="TL", text="a", font="9p,Helvetica-Bold", offset="0.05c/-0.05c")
        g_h = FIG_WIDTH_CM / 1.97      # Robinson height/width

        # -- b, c, ...: country insets, identical frames in one row
        fig.shift_origin(yshift=f"-{inset_h + 1.0}c")
        for k, i in enumerate(insets):
            if k > 0:
                fig.shift_origin(xshift=f"{inset_w + INSET_GAP_CM}c")
            w, e, s, n = i["region"]
            step = nice_interval(min(e - w, n - s))
            fmt = "dddF" if step >= 1 else "ddd:mmF"
            w0, n0, dx, dy = i["geo"]
            grid = to_grid(i["keys"], w0, n0, dx, dy)
            proj = f"Q{0.5 * (w + e):.4f}/{0.5 * (s + n):.4f}/{inset_w}c"
            with pygmt.config(FORMAT_GEO_MAP=fmt):
                draw_hlz_panel(fig, i["region"], proj, grid, cpt,
                               hillshade(i["region"], relief_resolution(i["region"])),
                               clip_dcw=i["iso"],
                               frame=("WSne", f"xa{step}f{step / 2}",
                                      f"ya{step}f{step / 2}"))
            plot_points(fig, i["sub"])
            label_sites(fig, i["sub"])
            fig.text(position="TL", text=f"{letters[k + 1]}  {i['country']}",
                     font="8p,Helvetica-Bold", offset="0.1c/-0.1c", fill="white@20",
                     clearance="0.05c/0.05c")
        fig.shift_origin(xshift=f"-{(n_in - 1) * (inset_w + INSET_GAP_CM)}c")

        # -- Holdridge chart as legend
        s32 = np.sqrt(3) / 2
        tri_h = ((3 + 5.5) * s32 + 2.8) * TRIANGLE_WIDTH_CM / 18.3
        fig.shift_origin(yshift=f"-{tri_h + 0.8}c")
        draw_holdridge_chart(fig, df, TRIANGLE_WIDTH_CM)

    out = fig_dir / FIG_STEM
    fig.savefig(f"{out}.png", dpi=400)
    fig.savefig(f"{out}.pdf")
    print(f"Figure written to {out}.png and .pdf (global map height {g_h:.1f} cm)")


# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--root", default=None, help="LAIMethodComparison directory")
    ap.add_argument("--variant", default="Classical", choices=list(FIGSHARE))
    ap.add_argument("--no-plot", action="store_true", help="only update the CSV")
    args = ap.parse_args()

    root = find_project_root(args.root)
    csv_path = root / PLOT_CSV_REL
    tif, codes_csv = fetch_hlz(root, args.variant)
    names = read_code_table(codes_csv)

    df = add_holdridge_columns(csv_path, tif, names)
    print(df[["PlotID", "Country", "Location", "Holdridge_code", "Holdridge_life_zone"]]
          .to_string(index=False))
    print(f"\nUpdated {csv_path}")

    if not args.no_plot:
        make_figure(df, tif, names, root)


if __name__ == "__main__":
    main()