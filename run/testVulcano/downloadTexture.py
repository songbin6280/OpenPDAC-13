#
# This script downloads a satellite image that geographically matches an
# ESRI ASCII raster file (.asc) and reprojects it to have the exact same
# coordinate system and bounds as the source raster.
# It includes an upscale factor to generate a final texture with a higher
# resolution than the source DEM.
#

import rasterio
from rasterio.io import MemoryFile
from rasterio.warp import reproject, Resampling
from rasterio.transform import from_bounds
import contextily as cx
import numpy as np
from PIL import Image
from pyproj import Transformer, CRS

# ----- USER SETTINGS -----

# Path to your ESRI ASCII DEM file.
RASTER_FILE_PATH = "./constant/DEM/DEMcropped.asc"

# Filename for the final texture output.
OUTPUT_TEXTURE_PATH = "./constant/DEM/texture.jpg"

# Fallback Coordinate Reference System (CRS) if it's not found in the raster file.
# 'EPSG:32633' is for WGS 84 / UTM zone 33N.
SOURCE_CRS_FALLBACK = "EPSG:32633"

# Tile provider for the satellite imagery. Esri.WorldImagery is a great choice.
MAP_PROVIDER = cx.providers.Esri.WorldImagery

# --- SETTINGS FOR HIGHER RESOLUTION ---

# UPSCALE_FACTOR: Determines the final output size.
# 1 = Same size as DEM (e.g., 800x800)
# 2 = 2x the DEM size (e.g., 1600x1600)
# 4 = 4x the DEM size (e.g., 3200x3200)
UPSCALE_FACTOR = 2

# FORCED_ZOOM: Determines the quality of the downloaded source image.
# For higher upscale factors, you should use a higher zoom level to provide
# enough detail. A good rule of thumb is zoom 16 or 17 for a 2x-4x upscale.
FORCED_ZOOM = 16

# ----- END OF USER SETTINGS -----


# The target CRS for web map tile providers is always Web Mercator.
TARGET_CRS_MERCATOR = "EPSG:32633"

try:
    # --- STEP 1: Read geographic metadata from the source .asc file ---
    print(f"1. Reading geographic metadata from source raster: '{RASTER_FILE_PATH}'")
    with rasterio.open(RASTER_FILE_PATH) as dst_raster:
        dst_crs = dst_raster.crs if dst_raster.crs else CRS(SOURCE_CRS_FALLBACK)
        dst_transform = dst_raster.transform
        dst_width = dst_raster.width
        dst_height = dst_raster.height
        dst_bounds = dst_raster.bounds
        
        print(f"   - Base DEM dimensions: {dst_width}x{dst_height} pixels")
        print(f"   - Target CRS: {dst_crs.name}")

    # --- STEP 2: Create an upscaled destination profile ---
    # Apply the upscale factor to the dimensions and transform
    final_width = dst_width * UPSCALE_FACTOR
    final_height = dst_height * UPSCALE_FACTOR

    # The new transform scales the pixel size down by the upscale factor
    scaled_transform = rasterio.Affine(
        dst_transform.a / UPSCALE_FACTOR, dst_transform.b, dst_transform.c,
        dst_transform.d, dst_transform.e / UPSCALE_FACTOR, dst_transform.f
    )

    dst_profile = {
        'driver': 'GTiff', 'crs': dst_crs, 'transform': scaled_transform,
        'width': final_width, 'height': final_height,
        'count': 3, 'dtype': 'uint8'
    }
    print(f"   - Final texture dimensions will be: {final_width}x{final_height} pixels")

    # --- STEP 3: Download the source image ---
    print(f"\n2. Downloading source image from '{MAP_PROVIDER.name}' (Zoom: {FORCED_ZOOM})...")
    transformer = Transformer.from_crs(dst_crs, TARGET_CRS_MERCATOR, always_xy=True)
    west_merc, south_merc, east_merc, north_merc = transformer.transform_bounds(*dst_bounds)
    
    src_array, src_extent = cx.bounds2img(
        west_merc, south_merc, east_merc, north_merc,
        source=MAP_PROVIDER, ll=False, zoom=FORCED_ZOOM
    )
    print(f"   - Download complete. Source image dimensions: {src_array.shape[1]}x{src_array.shape[0]} pixels")

    if src_array.shape[2] == 4:
        print("   - Alpha channel (RGBA) detected, converting to RGB.")
        src_array = src_array[:, :, :3]

    # --- STEP 4: Perform the reprojection ---
    print("\n3. Reprojecting the texture onto the upscaled DEM grid...")
    src_west, src_east, src_south, src_north = src_extent
    
    src_profile = dict(
        driver='GTiff', crs=CRS(TARGET_CRS_MERCATOR),
        transform=from_bounds(src_west, src_south, src_east, src_north, width=src_array.shape[1], height=src_array.shape[0]),
        count=3, width=src_array.shape[1], height=src_array.shape[0], dtype='uint8'
    )
    
    with MemoryFile() as dst_memfile:
        with dst_memfile.open(**dst_profile) as dst:
            with MemoryFile() as src_memfile:
                with src_memfile.open(**src_profile) as src:
                    src.write(np.moveaxis(src_array, -1, 0))
                    
                    reproject(
                        source=rasterio.band(src, [1, 2, 3]),
                        destination=rasterio.band(dst, [1, 2, 3]),
                        resampling=Resampling.cubic_spline
                    )
            final_image_array = dst.read()

    # --- STEP 5: Save the final image ---
    print("\n4. Saving the final image...")
    final_image_array = np.moveaxis(final_image_array, 0, -1)
    img = Image.fromarray(final_image_array)
    img.save(OUTPUT_TEXTURE_PATH, quality=95)
    
    print("\nOperation completed successfully!")
    print(f"Texture saved to '{OUTPUT_TEXTURE_PATH}' with final dimensions {final_width}x{final_height}.")

except FileNotFoundError:
    print(f"ERROR: The file '{RASTER_FILE_PATH}' was not found.")
except Exception as e:
    print(f"An unexpected error occurred: {e}")
