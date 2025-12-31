import logging  # Provides configurable logging instead of ad-hoc prints
import sys
import numpy as np
from netCDF4 import Dataset
from util.Interpolator import Interp2D, Interp3D, depths
from util.ROMS import ROMS

# Configure a basic logger for the script runtime
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)


if __name__ == '__main__':
    # Validate that the script received all required arguments
    if len(sys.argv) != 5:
        logger.error("Usage: python %s initialization_date source_file history_dir destination_file", sys.argv[0])
        sys.exit(-1)

    # Extract command-line arguments for clarity
    iDate = sys.argv[1]  # Initialization date (string)
    src = sys.argv[2]  # Path to the ROMS source NetCDF file
    history_dir = sys.argv[3]  # Historical directory (currently unused but kept for parity)
    dst = sys.argv[4]  # Destination path for the processed archive

    # Log the received arguments to help trace inputs
    logger.info("iDate:%s src:%s history:%s dst:%s", iDate, src, history_dir, dst)

    # Open the NetCDF source file to access simulation outputs
    ncsrcfile = Dataset(src)

    # Read temporal and grid variables from the source file
    time = ncsrcfile.variables["ocean_time"][:]  # Simulation time steps
    Xlat = ncsrcfile["lat_rho"][:]  # Latitude coordinates on Rho points
    Xlon = ncsrcfile["lon_rho"][:]  # Longitude coordinates on Rho points
    Ulat = ncsrcfile["lat_u"][:]  # Latitude coordinates on U points
    Ulon = ncsrcfile["lon_u"][:]  # Longitude coordinates on U points
    Vlat = ncsrcfile["lat_v"][:]  # Latitude coordinates on V points
    Vlon = ncsrcfile["lon_v"][:]  # Longitude coordinates on V points
    s_rho = ncsrcfile["s_rho"][:]  # Vertical sigma levels
    mask_rho = ncsrcfile["mask_rho"][:]  # Land/sea mask for Rho points
    mask_u = ncsrcfile["mask_u"][:]  # Land/sea mask for U points
    mask_v = ncsrcfile["mask_v"][:]  # Land/sea mask for V points
    H = ncsrcfile["h"][:]  # Bathymetry depth

    # Build evenly spaced destination longitude and latitude arrays
    dstLon = np.linspace(Xlon.min(), Xlon.max(), len(Xlon[0]))
    dstLat = np.linspace(Xlat.min(), Xlat.max(), len(Xlat))

    # Instantiate a ROMS archive file writer with the destination grid
    roms = ROMS(dst, time, depths, dstLon, dstLat)

    # Create 2D bilinear interpolators for Rho, U, and V grids
    interpolator2DRho = Interp2D(Xlon, Xlat, dstLon, dstLat)
    interpolator2DU = Interp2D(Ulon, Ulat, dstLon, dstLat)
    interpolator2DV = Interp2D(Vlon, Vlat, dstLon, dstLat)

    # Create 3D interpolators that also account for vertical sigma levels and masks
    interpolator3DRho = Interp3D(Xlon, Xlat, dstLon, dstLat, s_rho, mask_rho, H)
    interpolator3DU = Interp3D(Ulon, Ulat, dstLon, dstLat, s_rho, mask_u, H)
    interpolator3DV = Interp3D(Vlon, Vlat, dstLon, dstLat, s_rho, mask_v, H)

    # Interpolate bathymetry depth onto the destination grid
    logger.info("Interpolating bathymetry depth (h)")
    H = interpolator2DRho.interp(H)

    # Interpolate sea surface height
    logger.info("Interpolating sea surface height (zeta)")
    zeta = ncsrcfile["zeta"][:]
    zeta = interpolator2DRho.interp(zeta)

    # Interpolate temperature on all vertical levels
    logger.info("Interpolating 3D temperature field")
    temp = ncsrcfile["temp"][:]
    temp = interpolator3DRho.interp(temp)

    # Extract temperature at the ocean bottom for diagnostics
    logger.info("Extracting bottom temperature")
    temp_at_bottom = interpolator3DRho.bottomValues(temp)

    # Extract temperature at the ocean surface for diagnostics
    logger.info("Extracting surface temperature")
    temp_at_surface = interpolator3DRho.surfaceValues(temp)

    # Interpolate salinity on all vertical levels
    logger.info("Interpolating 3D salinity field")
    salt = ncsrcfile["salt"][:]
    salt = interpolator3DRho.interp(salt)

    # Extract bottom salinity
    logger.info("Extracting bottom salinity")
    salt_at_bottom = interpolator3DRho.bottomValues(salt)

    # Extract surface salinity
    logger.info("Extracting surface salinity")
    salt_at_surface = interpolator3DRho.surfaceValues(salt)

    # Interpolate 3D U-velocity component
    logger.info("Interpolating 3D U current")
    u = ncsrcfile["u"][:]
    u = interpolator3DU.interp(u)

    # Extract bottom U current
    logger.info("Extracting bottom U current")
    U_at_bottom = interpolator3DRho.bottomValues(u)

    # Extract surface U current with slight amplification near surface
    logger.info("Extracting surface U current")
    U_at_surface = interpolator3DRho.surfaceValues(u, factor=1.2)

    # Interpolate depth-averaged U current
    logger.info("Interpolating depth-averaged U current (ubar)")
    ubar = ncsrcfile["ubar"][:]
    ubar = interpolator2DU.interp(ubar)

    # Interpolate 3D V-velocity component
    logger.info("Interpolating 3D V current")
    v = ncsrcfile["v"][:]
    v = interpolator3DV.interp(v)

    # Extract bottom V current
    logger.info("Extracting bottom V current")
    V_at_bottom = interpolator3DRho.bottomValues(v)

    # Extract surface V current with slight amplification near surface
    logger.info("Extracting surface V current")
    V_at_surface = interpolator3DRho.surfaceValues(v, factor=1.2)

    # Interpolate depth-averaged V current
    logger.info("Interpolating depth-averaged V current (vbar)")
    vbar = ncsrcfile["vbar"][:]
    vbar = interpolator2DV.interp(vbar)

    # Save all processed variables into the ROMS archive
    logger.info("Saving interpolated fields to ROMS archive")
    roms.h = H
    roms.temp = temp
    roms.tempBottom = temp_at_bottom
    roms.tempSurface = temp_at_surface
    roms.salt = salt
    roms.saltBottom = salt_at_bottom
    roms.saltSurface = salt_at_surface
    roms.zeta = zeta
    roms.U = u
    roms.uBottom = U_at_bottom
    roms.uSurface = U_at_surface
    roms.V = v
    roms.vBottom = V_at_bottom
    roms.vSurface = V_at_surface
    roms.ubar = ubar
    roms.vbar = vbar
    roms.write()

    # Close the NetCDF and archive handles to flush output
    ncsrcfile.close()
    roms.close()
