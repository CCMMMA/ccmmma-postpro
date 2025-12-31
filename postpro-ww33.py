import logging  # Structured logging for traceability
import sys
import numpy as np
from netCDF4 import Dataset
from util.WW33 import WW33
from util.Interpolator import Interp2D

# Initialize a logger with a simple format
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)


if __name__ == '__main__':
    # Verify the expected number of arguments
    if len(sys.argv) != 5:
        logger.error("Usage: python %s initialization_date source_file history_dir destination_file", sys.argv[0])
        sys.exit(-1)

    # Parse the command-line inputs
    iDate = sys.argv[1]  # Initialization date string
    src = sys.argv[2]  # Path to the WW3 NetCDF source
    history_dir = sys.argv[3]  # Historical path (placeholder)
    dst = sys.argv[4]  # Output archive path

    # Log the inputs so users know what is being processed
    logger.info("iDate:%s src:%s history:%s dst:%s", iDate, src, history_dir, dst)

    # Open the NetCDF source file containing wave model outputs
    ncsrcfile = Dataset(src)

    # Read time and coordinate variables
    time = ncsrcfile.variables["time"][:]  # Forecast lead times
    srcLats = ncsrcfile["latitude"][:]  # Source latitude array
    srcLons = ncsrcfile["longitude"][:]  # Source longitude array

    # Build 2D meshgrid for interpolation
    Xlon, Xlat = np.meshgrid(srcLons, srcLats)
    dLon = (srcLons[1] - srcLons[0]) * 0.75  # Destination longitudinal spacing
    dLat = (srcLats[1] - srcLats[0]) * 0.75  # Destination latitudinal spacing
    dstLat = np.arange(Xlat.min(), Xlat.max(), dLat)  # Target latitude grid
    dstLon = np.arange(Xlon.min(), Xlon.max(), dLon)  # Target longitude grid

    # Instantiate a WW33 archive writer on the destination grid
    ww33 = WW33(dst, time, dstLon, dstLat)

    # Set up 2D bilinear interpolator for all scalar fields
    interpolator2D = Interp2D(Xlon, Xlat, dstLon, dstLat)

    # Interpolate bathymetry depth
    logger.info("Interpolating bathymetry depth (dpt)")
    dpt = ncsrcfile.variables["dpt"][:]
    dpt = interpolator2D.interp(dpt)

    # Interpolate significant wave height
    logger.info("Interpolating significant wave height (hs)")
    hs = ncsrcfile.variables["hs"][:]
    hs = interpolator2D.interp(hs)

    # Interpolate mean wave length
    logger.info("Interpolating mean wave length (lm)")
    lm = ncsrcfile.variables["lm"][:]
    lm = interpolator2D.interp(lm)

    # Interpolate peak wave frequency
    logger.info("Interpolating peak frequency (fp)")
    fp = ncsrcfile.variables["fp"][:]
    fp = interpolator2D.interp(fp)

    # Interpolate mean wave direction
    logger.info("Interpolating mean wave direction (dir)")
    dir = ncsrcfile.variables["dir"][:]
    dir = interpolator2D.interp(dir)

    # Interpolate mean wave period
    logger.info("Interpolating mean wave period (t0m1)")
    t0m1 = ncsrcfile.variables["t0m1"][:]
    t0m1 = interpolator2D.interp(t0m1)

    # Persist all processed wave diagnostics
    logger.info("Saving processed wave fields to archive")
    ww33.dpt = dpt
    ww33.hs = hs
    ww33.lm = lm
    ww33.fp = fp
    ww33.dir = dir
    ww33.period = t0m1
    ww33.write()

    # Close file handles to flush changes
    ncsrcfile.close()
    ww33.close()
