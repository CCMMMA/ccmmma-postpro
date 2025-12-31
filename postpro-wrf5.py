import logging  # Use structured logging for traceable runtime feedback
import sys
import numpy as np
import math
from os.path import basename
from netCDF4 import Dataset
from datetime import timedelta, datetime
from wrf import getvar, interplevel
from util.WRF import WRF
from util.Interpolator import Interp2D

# Configure module-level logging with timestamps
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)


def get_date_time(date):
    # Build a compact YYYYMMDDZHHMM timestamp for filenames
    dateTime = format(date.year,'04') + format(date.month,'02') + format(date.day,'02') + "Z"+format(date.hour,'02') + format(date.minute,'02')
    return dateTime

def get_date_time_path(date):
    # Build YYYY/MM/DD directory structure used by history storage
    dateTimePath = format(date.year,'04') + "/" + format(date.month,'02') + "/"+format(date.day,'02')
    return dateTimePath

def getBoundaries(Xlon, Xlat):
    # Compute the geographic bounding box of the input grid
    row_lat = len(Xlat) - 1
    col_lat = len(Xlat[0]) - 1

    row_long = len(Xlon) - 1
    col_long = len(Xlon[0]) - 1

    A = [Xlat[0][0], Xlon[0][0]]
    B = [Xlat[0][col_lat], Xlon[0][col_long]]
    C = [Xlat[row_lat][col_lat], Xlon[row_long][col_long]]
    D = [Xlat[row_lat][0], Xlon[row_long][0]]

    min_lat = Xlat[0][0]
    minI = 0

    ''' from A to B '''
    for i in range(col_lat,-1,-1):
        np1 = [Xlat[0][i], Xlon[0][i]]
        if np1[0] > min_lat:
            minI = i
            min_lat = np1[0]

    max_lat = Xlat[row_lat][col_lat]
    maxI = col_lat

    ''' from C to D '''
    for i in range(col_lat,-1,-1):
        np1 = [Xlat[row_lat][i], Xlon[row_long][i]]
        if np1[0] < max_lat:
            maxI = i
            max_lat = np1[0]

    min_long = Xlon[0][0]
    minJ = 0

    ''' from A to D '''
    for i in range(row_lat,-1,-1):
        np1 = [Xlat[i][0], Xlon[i][0]]
        if np1[1] > min_long:
            minJ = i
            min_long = np1[1]

    max_long = Xlon[0][col_long]
    maxJ = row_lat

    ''' from B to C '''
    for i in range(row_lat,-1,-1):
        np1 = [Xlat[i][col_lat], Xlon[i][col_lat]]
        if np1[1] < max_long:
            maxJ = i
            max_long = np1[1]


    #minLat=np.asscalar(min_lat)
    #maxLat=np.asscalar(max_lat)
    #minLon=np.asscalar(min_long)
    #maxLon=np.asscalar(max_long)
    minLat = min_lat
    maxLat = max_lat
    minLon = min_long
    maxLon = max_long

    return minLon, minLat, maxLon, maxLat


if __name__ == '__main__':
    # Confirm the required arguments are supplied
    if len(sys.argv) != 6:
        logger.error("Usage: python %s initialization_date source_file source_file_1hago source_file_00 destination_file", sys.argv[0])
        sys.exit(-1)

    # Parse CLI arguments for readability and future validation
    iDate = sys.argv[1]  # Initialization date string
    src = sys.argv[2]  # Path to current WRF NetCDF output
    src_1hago = sys.argv[3]  # Path to previous hour WRF output
    src_00 = sys.argv[4]  # Path to 00Z WRF output
    dst = sys.argv[5]  # Destination archive path

    # Log incoming paths to aid debugging
    logger.info("iDate:%s src:%s src_1hago:%s src_00:%s dst:%s", iDate, src, src_1hago, src_00, dst)

    # Open the NetCDF file for the current forecast time
    ncsrcfile = Dataset(src)

    # Extract and decode the model timestamp into datetime components
    timeVariable = [el.decode('UTF-8') for el in ncsrcfile.variables["Times"][:][0]]
    datetimeStr = ''.join(timeVariable).split("_")
    dateStr = datetimeStr[0].split("-")
    timeStr = datetimeStr[1].split(":")

    # Build datetime objects for the current hour, previous hour, and midnight
    datetime_current = datetime(int(dateStr[0]), int(dateStr[1]), int(dateStr[2]), int(timeStr[0]), int(timeStr[1]), int(timeStr[2]))
    datetime_1h_ago = datetime_current - timedelta(hours=1)
    datetime_00 = datetime(int(dateStr[0]), int(dateStr[1]), int(dateStr[2]), 0, 0, 0)

    # Store time as list for compatibility with WRF archive writer
    time = [datetime_current]
    logger.info("Dates -- current:%s 1h ago:%s today:%s", datetime_current, datetime_1h_ago, datetime_00)

    # Log reference files to clarify which inputs are used for deltas
    logger.info("Prev hour file:%s", src_1hago)
    logger.info("Current day file:%s", src_00)

    # Read latitude and longitude grids from the WRF dataset
    Xlat = np.array(getvar(ncsrcfile, "XLAT", meta=False))
    Xlon = np.array(getvar(ncsrcfile, "XLONG", meta=False))

    # Compute mean location to derive grid spacing in degrees
    lon = np.average(Xlon)
    lat = np.average(Xlat)

    # Earth's radius (meters) used for geographic calculations
    R = 6378137

    # Extract native grid spacing in meters from the dataset
    dn = ncsrcfile.DY
    de = ncsrcfile.DX

    # Convert grid spacing to degrees for the destination grid
    dLat = 0.5 * (dn/R) * 180/math.pi
    dLon = 0.5 * (de/(R * math.cos(math.pi * lat/180))) * 180/math.pi

    # Determine min/max bounds to build target latitude and longitude arrays
    minLon, minLat, maxLon, maxLat = getBoundaries(Xlon, Xlat)
    dstLat = np.arange(minLat, maxLat, dLat)
    dstLon = np.arange(minLon, maxLon, dLon)

    # Initialize the output WRF archive writer with the target grid
    wrf = WRF(dst, time, dstLon, dstLat)

    # Create a bilinear 2D interpolator for horizontal fields
    interpolator2DRho = Interp2D(Xlon, Xlat, dstLon, dstLat)

    # Extract atmospheric state variables at multiple pressure levels
    p = getvar(ncsrcfile, "pressure")
    z = getvar(ncsrcfile, "z", units="dm")
    tc = getvar(ncsrcfile, "temp", units="degC")
    td = getvar(ncsrcfile, "td", units="degC")
    ter = getvar(ncsrcfile,"ter", units="m")
    theta_e = getvar(ncsrcfile, "theta_e", units="degC")
    theta_w = getvar(ncsrcfile, "twb", units="degC")
    rh = getvar(ncsrcfile, "rh")
    uvmet = getvar(ncsrcfile, "uvmet")

    # Interpolate thermodynamic and wind fields to standard pressure levels
    tc1000  = interplevel(tc, p, 1000)
    rh1000  = interplevel(rh, p, 1000)
    u1000  = interplevel(uvmet[0], p, 1000)
    v1000  = interplevel(uvmet[1], p, 1000)
    tc975  = interplevel(tc, p, 975)
    rh975  = interplevel(rh, p, 975)
    u975  = interplevel(uvmet[0], p, 975)
    v975  = interplevel(uvmet[1], p, 975)
    tc950 = interplevel(tc, p, 950)
    rh950  = interplevel(rh, p, 950)
    u950  = interplevel(uvmet[0], p, 950)
    v950  = interplevel(uvmet[1], p, 950)
    tc925  = interplevel(tc, p, 925)
    rh925  = interplevel(rh, p, 925)
    u925  = interplevel(uvmet[0], p, 925)
    v925  = interplevel(uvmet[1], p, 925)
    gph850 = interplevel(z, p, 850)
    theta_e850 = interplevel(theta_e, p, 850)
    theta_w850 = interplevel(theta_w, p, 850)
    td850 = interplevel(td, p, 850)
    tc850 = interplevel(tc, p, 850)
    rh850  = interplevel(rh, p, 850)
    u850  = interplevel(uvmet[0], p, 850)
    v850  = interplevel(uvmet[1], p, 850)
    td700 = interplevel(td, p, 700)
    tc700 = interplevel(tc, p, 700)
    rh700  = interplevel(rh, p, 700)
    u700  = interplevel(uvmet[0], p, 700)
    v700  = interplevel(uvmet[1], p, 700)
    theta_e500 = interplevel(theta_e, p, 500)
    gph500 = interplevel(z, p, 500)
    tc500 = interplevel(tc, p, 500)
    rh500  = interplevel(rh, p, 500)
    u500  = interplevel(uvmet[0], p, 500)
    v500  = interplevel(uvmet[1], p, 500)
    tc300 = interplevel(tc, p, 300)
    rh300  = interplevel(rh, p, 300)
    u300  = interplevel(uvmet[0], p, 300)
    v300  = interplevel(uvmet[1], p, 300)

    # Compute stability indices
    tt = tc850 + td850 - 2 * tc500  
    ki = (tc850 - tc500) + td850 - (tc700 - td700)
    delta_theta = theta_e500 - theta_e850

    # Gather severe weather diagnostics
    cape_2d = getvar(ncsrcfile, "cape_2d", meta=False)
    updraft_helicity = getvar(ncsrcfile, "updraft_helicity", meta=False)
    helicity = getvar(ncsrcfile, "helicity", meta=False)

    # Get the sea level pressure
    slp = getvar(ncsrcfile, "slp", meta=False)
    rh2 = getvar(ncsrcfile, "rh2", meta=False)
    pw = getvar(ncsrcfile, "pw", meta=False)

    # Cloud fraction as the maximum of the low and mid layers
    cloudfrac = getvar(ncsrcfile, "cloudfrac", meta=False)
    clf = np.maximum(cloudfrac[0], cloudfrac[1])

    # Read the temperature at 2m in celsius
    t2c = ncsrcfile["T2"][:] - 273.15

    # Read the snow rate
    sr = ncsrcfile["SR"][:]

    # Get the wind at 10m u and v components (meteo oriented)
    uvmet10 = getvar(ncsrcfile, "uvmet10", meta=False)

    # Get the wind speed and wind dir at 10m (meteo oriented)
    uvmet10_wspd_wdir = getvar(ncsrcfile, "uvmet10_wspd_wdir", meta=False)

    # Read the simulation cumulated rain from the current file
    rain = ncsrcfile["RAINC"][:] + ncsrcfile["RAINNC"][:] + ncsrcfile["RAINSH"][:]

    # Set placeholders for the 00Z comparison fields
    raini_00 = None
    wspd10i_00 = None
    wdir10i_00 = None

    # Interpolate current accumulated rain and snow rate
    logger.info("Interpolating accumulated rain and snow rate")
    raini = interpolator2DRho.interp(rain[0])
    sri = interpolator2DRho.interp(sr[0])

    # Interpolate 10 m wind speed and direction
    logger.info("Interpolating 10 m wind speed and direction")
    wspd10i = interpolator2DRho.interp(uvmet10_wspd_wdir[0])
    wdir10i = interpolator2DRho.interp(uvmet10_wspd_wdir[1])

    try:
        # Try to open the midnight dataset to compute daily deltas
        ncsrc_00 = Dataset(src_00)
        logger.info("Calculating daily deltas against %s", src_00)

        # Read and interpolate accumulated rain at 00Z
        rain_00 = ncsrc_00["RAINC"][:] + ncsrc_00["RAINNC"][:] + ncsrc_00["RAINSH"][:]
        raini_00 = interpolator2DRho.interp(rain_00[0])

        # Close the dataset after use
        ncsrc_00.close()
        logger.info("Daily dataset processed successfully")
    except Exception:
        # If the midnight dataset is missing, fall back to current totals
        logger.exception("Unable to process daily dataset %s; using current rain totals", src_00)
        raini_00 = raini

    # Calculate daily accumulated rain by differencing against 00Z
    logger.info("Calculating daily accumulated rain")
    draini = raini - raini_00

    # Initialize placeholders for previous hour comparisons
    raini_1hago = None
    wspd10i_1hago = None
    wdir10i_1hago = None

    try:
        # Try to open the previous hour dataset to compute hourly deltas
        ncsrc_1hago = Dataset(src_1hago)
        logger.info("Calculating 1 hour deltas against %s", src_1hago)

        # Read and interpolate accumulated rain one hour ago
        rain_1hago = ncsrc_1hago["RAINC"][:] + ncsrc_1hago["RAINNC"][:] + ncsrc_1hago["RAINSH"][:]
        raini_1hago = interpolator2DRho.interp(rain_1hago[0])

        # Read and interpolate previous hour 10 m wind speed/direction
        uvmet10_wspd_wdir_1hago = getvar(ncsrc_1hago, "uvmet10_wspd_wdir", meta=False)
        wspd10i_1hago = interpolator2DRho.interp(uvmet10_wspd_wdir_1hago[0])
        wdir10i_1hago = interpolator2DRho.interp(uvmet10_wspd_wdir_1hago[1])

        # Close the previous hour dataset
        ncsrc_1hago.close()
        logger.info("Previous hour dataset processed successfully")
    except Exception:
        # If the prior dataset is unavailable, reuse current interpolations
        logger.exception("Unable to process previous hour dataset %s; using current fields", src_1hago)
        raini_1hago = raini
        wspd10i_1hago = wspd10i
        wdir10i_1hago = wdir10i

    # Calculate hourly rainfall accumulation
    logger.info("Calculating hourly accumulated rain")
    hraini = raini - raini_1hago

    # Calculate wind shift over the last hour
    logger.info("Calculating 1-hour wind shift")
    dwspd10i = wspd10i - wspd10i_1hago
    dwdir10i = wdir10i - wdir10i_1hago

    # Calculate snow water equivalent from hourly rain and snow rate
    hswei = np.array(hraini * (sri - 0.75) * 5)
    hswei[hswei < 0] = 0

    # Interpolate precipitation and thermodynamic diagnostics
    logger.info("Interpolating atmospheric diagnostics")
    pw = interpolator2DRho.interp(pw)
    rh2 = interpolator2DRho.interp(rh2)
    t2c = interpolator2DRho.interp(t2c[0])
    uh = interpolator2DRho.interp(updraft_helicity)
    srh = interpolator2DRho.interp(helicity)
    mcape = interpolator2DRho.interp(cape_2d[0])
    mcin = interpolator2DRho.interp(cape_2d[1])

    # Interpolate wind and thermodynamic fields at multiple pressure levels
    logger.info("Interpolating wind and temperature at mandatory levels")
    u1000 = interpolator2DRho.interp(u1000)
    v1000 = interpolator2DRho.interp(v1000)
    tc1000 = interpolator2DRho.interp(tc1000)
    rh1000 = interpolator2DRho.interp(rh1000)
    u975 = interpolator2DRho.interp(u975)
    v975 = interpolator2DRho.interp(v975)
    tc975 = interpolator2DRho.interp(tc975)
    rh975 = interpolator2DRho.interp(rh975)
    u950 = interpolator2DRho.interp(u950)
    v950 = interpolator2DRho.interp(v950)
    tc950 = interpolator2DRho.interp(tc950)
    rh950 = interpolator2DRho.interp(rh950)
    u925 = interpolator2DRho.interp(u925)
    v925 = interpolator2DRho.interp(v925)
    tc925 = interpolator2DRho.interp(tc925)
    rh925 = interpolator2DRho.interp(rh925)
    u850 = interpolator2DRho.interp(u850)
    v850 = interpolator2DRho.interp(v850)
    tc850 = interpolator2DRho.interp(tc850)
    rh850 = interpolator2DRho.interp(rh850)
    u700 = interpolator2DRho.interp(u700)
    v700 = interpolator2DRho.interp(v700)
    tc700 = interpolator2DRho.interp(tc700)
    rh700 = interpolator2DRho.interp(rh700)
    u500 = interpolator2DRho.interp(u500)
    v500 = interpolator2DRho.interp(v500)
    tc500 = interpolator2DRho.interp(tc500)
    rh500 = interpolator2DRho.interp(rh500)
    u300 = interpolator2DRho.interp(u300)
    v300 = interpolator2DRho.interp(v300)
    tc300 = interpolator2DRho.interp(tc300)
    rh300 = interpolator2DRho.interp(rh300)

    # Interpolate geopotential heights and surface diagnostics
    logger.info("Interpolating geopotential heights and surface diagnostics")
    gph500 = interpolator2DRho.interp(gph500)
    gph850 = interpolator2DRho.interp(gph850)
    slp = interpolator2DRho.interp(slp)
    clf = interpolator2DRho.interp(clf)
    u10m = interpolator2DRho.interp(uvmet10[0])
    v10m = interpolator2DRho.interp(uvmet10[1])

    # Interpolate stability indices
    logger.info("Interpolating stability indices")
    tt = interpolator2DRho.interp(tt)
    ki = interpolator2DRho.interp(ki)
    theta_e850 = interpolator2DRho.interp(theta_e850)
    theta_w850 = interpolator2DRho.interp(theta_w850)
    delta_theta = interpolator2DRho.interp(delta_theta)
    
    # Persist all processed diagnostics into the WRF archive
    logger.info("Saving processed WRF diagnostics to archive")
    wrf.dwspd10 = dwspd10i
    wrf.dwdir10 = dwdir10i
    wrf.drain = draini
    wrf.hrain = hraini
    wrf.hswe = hswei
    wrf.pw = pw
    wrf.rh2 = rh2
    wrf.t2c = t2c
    wrf.uh = uh
    wrf.srh = srh
    wrf.mcape = mcape
    wrf.mcin = mcin
    wrf.u1000 = u1000
    wrf.v1000 = v1000
    wrf.tc1000 = tc1000
    wrf.rh1000 = rh1000
    wrf.u975 = u975
    wrf.v975 = v975
    wrf.tc975 = tc975
    wrf.rh975 = rh975
    wrf.u950 = u950
    wrf.v950 = v950
    wrf.tc950 = tc950
    wrf.rh950 = rh950
    wrf.u925 = u925
    wrf.v925 = v925
    wrf.tc925 = tc925
    wrf.rh925 = rh925
    wrf.u850 = u850
    wrf.v850 = v850
    wrf.tc850 = tc850
    wrf.rh850 = rh850
    wrf.u700 = u700
    wrf.v700 = v700
    wrf.tc700 = tc700
    wrf.rh700 = rh700
    wrf.u500 = u500
    wrf.v500 = v500
    wrf.tc500 = tc500
    wrf.rh500 = rh500
    wrf.u300 = u300
    wrf.v300 = v300
    wrf.tc300 = tc300
    wrf.rh300 = rh300
    wrf.tt = tt
    wrf.ki = ki
    wrf.theta_e850 = theta_e850
    wrf.theta_w850 = theta_w850
    wrf.delta_theta = delta_theta
    wrf.gph500 = gph500
    wrf.gph850 = gph850
    wrf.slp = slp
    wrf.clf = clf
    wrf.u10m = u10m
    wrf.v10m = v10m
    wrf.wspd10 = wspd10i
    wrf.wdir10 = wdir10i

    wrf.write()

    # Close the NetCDF and archive objects to ensure data is flushed
    ncsrcfile.close()
    wrf.close()
