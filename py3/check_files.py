from pathlib import Path
import pandas as pd

# ------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------
BASE = Path("/dme/eunomia/plotdata/measurments")
OUTPUT = Path("/lud11/poppman/plot_measurement_inventory.csv")

AIRBORNE = BASE / "airborne"
MULTISPECTRAL = AIRBORNE / "multispectral"
THERMAL = AIRBORNE / "thermal"
LIDAR = AIRBORNE / "lidar"
LICOR = BASE / "ground" / "par" / "licor"

# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------
def has_file(directory, extensions):
    """True if directory contains at least one non-empty file
    with one of the specified extensions."""
    if not directory.is_dir():
        return False

    return any(
        p.is_file()
        and p.stat().st_size > 0
        and p.suffix.lower() in extensions
        for p in directory.iterdir()
    )


def has_nonempty_folder(directory):
    """True if directory contains at least one subdirectory
    that contains a file or another subdirectory."""
    if not directory.is_dir():
        return False

    for folder in directory.iterdir():
        if folder.is_dir() and any(folder.iterdir()):
            return True

    return False


def plot_directories(directory):
    """Return plot IDs from subdirectories containing plot data."""
    if not directory.is_dir():
        return set()

    return {
        p.name
        for p in directory.iterdir()
        if p.is_dir()
    }


# ------------------------------------------------------------------
# Define measurements and their expected file types
# ------------------------------------------------------------------
measurements = {}

def add_raster_measurements(group, root):
    if not root.is_dir():
        return

    for method_dir in root.iterdir():
        if not method_dir.is_dir() or method_dir.name == "rawdata":
            continue

        method = method_dir.name

        def check(plot_id, method_dir=method_dir, method=method):
            plot_dir = method_dir / plot_id
            return has_file(plot_dir, {".tif", ".tiff"})

        measurements[f"{group}_{method}"] = {
            "plots": plot_directories(root / method),
            "check": check,
        }


def add_rawdata_measurements(group, root):
    raw_root = root / "rawdata"

    if not raw_root.is_dir():
        return

    def check(plot_id):
        plot_dir = raw_root / plot_id

        if not plot_dir.is_dir():
            return False

        return any(
            item.is_file() and item.stat().st_size > 0
            or item.is_dir() and any(item.iterdir())
            for item in plot_dir.iterdir()
        )

    measurements[f"{group}_rawdata"] = {
        "plots": set(),
        "check": check,
    }


# Multispectral and thermal rasters
add_raster_measurements("multispectral", MULTISPECTRAL)
add_rawdata_measurements("multispectral", MULTISPECTRAL)

add_raster_measurements("thermal", THERMAL)
add_rawdata_measurements("thermal", THERMAL)

# LiDAR rasters (DEM, DOM, DSM, etc.)
add_raster_measurements("lidar", LIDAR)

# LiDAR point clouds
pointcloud_root = LIDAR / "pointclouds"

measurements["lidar_pointclouds"] = {
    "plots": plot_directories(pointcloud_root),
    "check": lambda plot_id: has_file(
        pointcloud_root / plot_id, {".las", ".laz"}
    ),
}

add_rawdata_measurements("lidar", LIDAR)

# Ground PAR / Licor CSV files
measurements["ground_par_licor"] = {
    "plots": plot_directories(LICOR),
    "check": lambda plot_id: has_file(
        LICOR / plot_id, {".csv"}
    ),
}

# ------------------------------------------------------------------
# Discover all plot IDs
# ------------------------------------------------------------------
all_plots = set()

for measurement in measurements.values():
    all_plots.update(measurement["plots"])

all_plots = sorted(all_plots)

# ------------------------------------------------------------------
# Build Boolean inventory
# ------------------------------------------------------------------
rows = []

for plot_id in all_plots:
    row = {"plot_id": plot_id}

    for name, measurement in measurements.items():
        row[name] = bool(measurement["check"](plot_id))

    rows.append(row)

df = pd.DataFrame(rows)

# Ensure Boolean columns are written as True / False.
if not df.empty:
    df = df.astype({
        column: bool
        for column in df.columns
        if column != "plot_id"
    })

df.to_csv(OUTPUT, index=False)

print(f"Found {len(all_plots)} plots.")
print(f"Inventory saved to: {OUTPUT.resolve()}")
print(f"Number of measurement columns: {len(measurements)}")
print(df.head())