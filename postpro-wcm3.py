import logging  # Use logging for traceable runtime output
import sys
import numpy as np
from netCDF4 import Dataset
from util.Interpolator import Interp2D, depths
from util.Distributor import Distrib3D
from util.Wacomm import Wacomm

# Configure a simple logger for the script
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)


def compute_sfconc(conc, depth_limit, mask2d, depths, fill_value=1e37):
    """
    Compute the vertical sum of concentration up to a given depth,
    ignoring fill values and masked areas.

    Parameters:
        conc (np.ndarray): 3D concentration array with shape (depth, lat, lon).
        depth_limit (float): Maximum depth (in meters) for the integration.
        mask2d (np.ndarray): 2D mask of shape (lat, lon) with 1 for water, 0 for land.
        depths (list or np.ndarray): List of depth levels corresponding to the first axis of `conc`.
        fill_value (float): Value used to indicate missing or invalid data.

    Returns:
        np.ndarray: 2D array (lat, lon) containing the summed concentration
                    from the surface down to `depth_limit`.
    """
    # Ensure depths is a NumPy array
    depth_arr = np.array(depths)

    # Create boolean mask of levels <= depth_limit
    depth_mask = depth_arr <= depth_limit  # shape (n_depths,)

    # Expand spatial mask to 3D and combine with depth mask
    mask3d = depth_mask[:, None, None] * mask2d[None, :, :]

    # Replace fill values with zero so they don’t inflate the sum
    conc_clean = np.where(conc == fill_value, 0.0, conc)

    # Sum along the depth axis
    sfconc = np.sum(conc_clean * mask3d, axis=0)  # shape (lat, lon)

    # Restore fill_value on land points
    sfconc[mask2d == 0] = fill_value

    return sfconc


if __name__ == '__main__':
    # Ensure the user supplied all required arguments
    if len(sys.argv) != 5:
        logger.error("Usage: python %s initialization_date source_file history_dir destination_file", sys.argv[0])
        sys.exit(-1)

    # Parse the CLI arguments for clarity
    iDate = sys.argv[1]  # Initialization date string
    src = sys.argv[2]  # Input NetCDF path
    history_dir = sys.argv[3]  # Historical folder (unused placeholder)
    dst = sys.argv[4]  # Output archive path

    # Log the inputs so users can trace what the script is processing
    logger.info("iDate:%s src:%s history:%s dst:%s", iDate, src, history_dir, dst)

    # Open the NetCDF file containing wave chemistry model data
    ncsrcfile = Dataset(src)

    # Read time and grid variables from the source dataset
    time = ncsrcfile.variables["ocean_time"][:]  # Time steps
    Xlat = ncsrcfile["lat_rho"][:]  # Latitude grid on Rho points
    Xlon = ncsrcfile["lon_rho"][:]  # Longitude grid on Rho points
    s_rho = ncsrcfile["s_rho"][:]  # Sigma levels for the vertical coordinate
    mask_rho = ncsrcfile["mask_rho"][:]  # Land/sea mask
    H = ncsrcfile["h"][:]  # Bathymetry

    # Build destination longitude and latitude arrays with even spacing
    dstLon = np.linspace(Xlon.min(), Xlon.max(), len(Xlon[0]))
    dstLat = np.linspace(Xlat.min(), Xlat.max(), len(Xlat))

    # Instantiate a Wacomm archive writer for the destination grid
    wacomm = Wacomm(dst, time, depths, dstLon, dstLat)

    # Create 2D bilinear interpolator for surface fields
    interpolator2DRho = Interp2D(Xlon, Xlat, dstLon, dstLat)

    # Create a 3D distributor that handles vertical stretching and land masking
    distributor3DRho = Distrib3D(Xlon, Xlat, dstLon, dstLat, s_rho, mask_rho, H)

    # Read and distribute the 3D concentration field onto the destination grid
    logger.info("Distributing concentration field")
    conc = ncsrcfile.variables["conc"][:]
    conc = distributor3DRho.distrib(conc)

    # Compute surface and depth-integrated concentration metrics
    logger.info("Calculating surface and integrated concentration totals")
    sfconc = conc[0, 0]
    sfconc_10m = compute_sfconc(conc[0], 10.0, distributor3DRho.mask, depths)
    sfconc_30m = compute_sfconc(conc[0], 30.0, distributor3DRho.mask, depths)

    # Persist outputs to the archive file
    logger.info("Saving processed concentration fields")
    wacomm.mask = distributor3DRho.mask
    wacomm.conc = conc
    wacomm.sfconc = sfconc
    wacomm.sfconc_10m = sfconc_10m
    wacomm.sfconc_30m = sfconc_30m
    wacomm.write()

    # Close the open dataset handles to flush to disk
    ncsrcfile.close()
    wacomm.close()
