#include <stdarg.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

#define LE_FLG_JPLEPH 2

#define LE_FLG_SWIEPH 4

#define LE_FLG_MOSEPH 8

#define LE_FLG_VSOP2013 1

#define LE_FLG_HELIO 16

#define LE_FLG_BARYHEL 32

#define LE_FLG_TOPOCTR 64

#define LE_FLG_XYZ 128

#define LE_FLG_SPEED 256

#define LE_FLG_NOABERR 512

#define LE_FLG_NOGDEFL 1024

#define LE_FLG_NOBIRR 2048

#define LE_FLG_J2000 4096

#define LE_FLG_JPLHOR 8192

#define LE_FLG_SIDEREAL 16384

#define LE_FLG_ICRS 32768

#define LE_FLG_EQUATORIAL 65536

#define LE_FLG_ECLIPTIC 131072

#define LE_FLG_TRUE 262144

#define LE_FLG_USER 524288

#define LE_FLG_NONUT 1048576

#define LE_FLG_SPEED3 2097152

#define LE_FLG_RADIANT 4194304

#define LE_FLG_CENTER_BODY 8388608

#define LE_FLG_TRUE_VELOCITY 16777216

#define LE_SUN 0

#define LE_MOON 1

#define LE_MERCURY 2

#define LE_VENUS 3

#define LE_MARS 4

#define LE_JUPITER 5

#define LE_SATURN 6

#define LE_URANUS 7

#define LE_NEPTUNE 8

#define LE_PLUTO 9

#define LE_CHIRON 10

#define LE_MEAN_BARY 11

#define LE_TRUE_NODE 12

#define LE_MEAN_NODE 13

#define LE_MEAN_APOG 14

#define LE_OSC_APOG 15

#define LE_INT_APOG 16

#define LE_EARTH 17

#define LE_NPLANETS 12

#define LE_NSIDEREAL_PLANETS 19

#define LE_NALL_PLANETS 19

#define LE_ECL_GEO 256

#define LE_ECL_HEL 512

#define LE_ECL_MAG 1024

#define LE_ECL_BARY 2048

#define LE_EQ_GEO 4096

#define LE_EQ_HEL 8192

#define LE_HOR 16384

#define LE_SIDM_FAGAN_BRADLEY 0

#define LE_SIDM_LAHIRI 1

#define LE_SIDM_DELUCE 2

#define LE_SIDM_RAMAN 3

#define LE_SIDM_USHASHASHI 4

#define LE_SIDM_KRISHNAMURTI 5

#define LE_SIDM_DJWHAL_KHUL 6

#define LE_SIDM_YUKTESHWAR 7

#define LE_SIDM_JN_BHASIN 8

#define LE_SIDM_BABYLONIAN_KUGLER1 9

#define LE_SIDM_BABYLONIAN_KUGLER2 10

#define LE_SIDM_BABYLONIAN_KUGLER3 11

#define LE_SIDM_BABYLONIAN_HUBER 12

#define LE_SIDM_BABYLONIAN_ETPSC 13

#define LE_SIDM_ALDEBARAN_15TAU 14

#define LE_SIDM_HIPPARCHUS 15

#define LE_SIDM_SASSANIAN 16

#define LE_SIDM_GALACTIC_CENTER_0SAG 17

#define LE_SIDM_J2000 18

#define LE_SIDM_J1900 19

#define LE_SIDM_B1950 20

#define LE_SIDM_SURYASIDDHANTA 21

#define LE_SIDM_ARYABHATA 22

#define LE_SIDM_SS_CHITRAPAKSHA 23

#define LE_SIDM_SS_REVATI 24

#define LE_SIDM_TRUE_CITRA 25

#define LE_SIDM_TRUE_REVATI 26

#define LE_SIDM_TRUE_PUSHYA 27

#define LE_SIDM_TRUE_ASHVINI 28

#define LE_SIDM_TRUE_MAGHA 29

#define LE_SIDM_TRUE_MULA 30

#define LE_SIDM_GALCENT_CENTER_0SAG 31

#define LE_SIDM_GALCENT_MULA 32

#define LE_SIDM_GALCENT_RADIO_0SAG 33

#define LE_SIDM_GALCENT_0CAP 34

#define LE_SIDM_TRUE_ANTARES 35

#define LE_SIDM_TRUE_FOMALHAUT 36

#define LE_SIDM_VALENTINE 37

#define LE_SIDM_USER 38

#define LE_SIDM_LARGE 39

#define LE_SIDM_TRUE_ALDEBARAN 40

#define LE_SIDM_STEINER_SET 41

#define LE_SIDM_RAMAN_A 42

#define LE_SIDM_SURYAPAKSHA_A 43

#define LE_SIDM_TRUE_VEGA 44

#define LE_SIDM_TRUE_ZETA_PSC 45

#define LE_SIDM_SURYAPAKSHA_B 46

#define LE_NMODES_AYANAMSA 47

#define LE_PREC_IAU_1976 0

#define LE_PREC_LASKAR_1986 1

#define LE_PREC_WILLIAMS_1994 2

#define LE_PREC_SIMON_1994 3

#define LE_PREC_IAU_2000 4

#define LE_PREC_BRETAGNON_2003 5

#define LE_PREC_IAU_2006 6

#define LE_PREC_VONDRAK_2011 7

#define LE_PREC_OWEN_1990 8

#define LE_PREC_NEWCOMB 9

#define LE_PREC_IAU_2006_REDUCED 10

#define LE_NMODES_PREC 11

#define LE_NUT_IAU_1980 0

#define LE_NUT_IAU_1980_HERRING 1

#define LE_NUT_IAU_2000A 2

#define LE_NUT_IAU_2000B 3

#define LE_NUT_WOOLARD_1953 4

#define LE_NMODES_NUT 5

#define LE_DT_STEPHENSON_1984 0

#define LE_DT_STEPHENSON_1997 1

#define LE_DT_ESPENAK_MEEUS_2006 2

#define LE_DT_STEPHENSON_2016 3

#define LE_DT_SCHOCH 4

#define LE_DT_USER 5

#define LE_NMODES_DT 6

#define LE_ST_IAU_1976 0

#define LE_ST_IAU_2006 1

#define LE_ST_IERS_2010 2

#define LE_ST_LONG_TERM 3

#define LE_NMODES_ST 4

#define LE_BIAS_NONE 0

#define LE_BIAS_IAU_2000 1

#define LE_BIAS_IAU_2006 2

#define LE_NMODES_BIAS 3

#define LE_HOR_APPROX_NONE 0

#define LE_HOR_APPROX_STANDARD 1

#define LE_HOR_APPROX_REFINED 2

#define LE_OK 0

#define LE_ERR -1

#define LE_ERR_INVALID_PARAMS -2

#define ERR_FILE_NOT_FOUND -3

#define ERR_OUT_OF_RANGE -4

#define ERR_NO_EPHEMERIS -5

#define ERR_IO -6

#define ERR_ENGINE -7

#define ERR_INVALID_PLANET -8

#define ERR_INVALID_FLAG -9

#define ERR_MEMORY -10

#define ERR_NOT_IMPLEMENTED -11

#define ERR_UNKNOWN_BODY -12

#define LE_X 0

#define LE_Y 1

#define LE_Z 2

#define LE_VX 3

#define LE_VY 4

#define LE_VZ 5

#define LE_RA 0

#define LE_DEC 1

#define LE_DIST 2

#define LE_RA_DOT 3

#define LE_DEC_DOT 4

#define LE_DIST_DOT 5

#define LE_LON 0

#define LE_LAT 1

#define LE_DISTANCE 2

#define LE_LON_DOT 3

#define LE_LAT_DOT 4

#define LE_DISTANCE_DOT 5

#define LE_MAX_ERR_LEN 256

#define LE_GAUSS_G 0.01720209895

#define LE_J1970 2440587.5

#define LE_J2000 2451545.0

#define LE_B1950 2433282.4235

#define LE_J1900 2415020.0

#define LE_DAY_PER_YEAR 365.25

#define LE_DAY_PER_CENTURY 36525.0

#define LE_DAY_PER_MILLENNIUM 365250.0

#define LE_DEG 0.017453292519943295

#define LE_RAD 57.29577951308232

#define LE_ARCSEC 4.84813681109536e-6

#define LE_ARCMIN 2.908882086657216e-4

#define LE_CLIGHT 173.1446326846693

/**
 * Valid range flags for Julian day conversions.
 */
#define LE_GREG_CAL 1

#define LE_JUL_CAL 0

#define HS_PLACIDUS 80

#define HS_KOCH 75

#define HS_EQUAL_ASC 69

#define HS_EQUAL_MC 65

#define HS_CAMPANUS 67

#define HS_REGIOMONTANUS 82

#define HS_MORINUS 77

#define HS_TOPOCENTRIC 84

#define HS_ALCABITIUS 66

#define HS_PORPHYRIUS 80

#define HS_GAUQUELIN 71

#define HS_HORIZONTAL 72

#define HS_WHOLE_SIGN 87

#define HS_AXIAL_ROT 88

#define HS_APC 89

#define HS_VERTEX 86

#define HS_SUNSHINE 83

#define HS_KRUSINSKI 85

#define HS_MERIDIAN 78

#define HS_POLAR 76

#define HS_AZIMUTHAL 90

#define HS_DEFAULT 80

/**
 * Number of stars in the embedded catalog.
 */
#define LE_STAR_COUNT 200

/**
 * Rise/set/transit event codes (matching `rsmi` parameter conventions).
 */
#define RISE 0

#define SET 1

#define TRANSIT_UPPER 2

#define TRANSIT_LOWER 3

/**
 * A fixed star entry in the catalog.
 */
typedef struct LeStar {
  /**
   * Bayer/Flamsteed name
   */
  int8_t name[32];
  /**
   * RA J2000 in radians
   */
  double ra;
  /**
   * Dec J2000 in radians
   */
  double dec;
  /**
   * Parallax in arcseconds
   */
  double parallax;
  /**
   * Proper motion in RA (radians/year)
   */
  double pm_ra;
  /**
   * Proper motion in Dec (radians/year)
   */
  double pm_dec;
  /**
   * Apparent visual magnitude
   */
  double mag;
} LeStar;

/**
 * C ABI: set ephemeris search paths (colon/semicolon-separated).
 */
void le_set_ephe_path(const int8_t *path);

/**
 * C ABI: set JPL DE ephemeris file path.
 */
void le_set_jpl_file(const int8_t *fname);

/**
 * C ABI: set observer position (longitude deg, latitude deg, altitude m).
 */
void le_set_topo(double geolon, double geolat, double geoalt);

/**
 * C ABI: set sidereal (ayanamsa) mode.
 */
void le_set_sid_mode(int32_t sid_mode, double t0, double ayan_t0);

/**
 * C ABI: override Delta T (TT - UT) in seconds.
 */
void le_set_delta_t_user(double dt);

/**
 * C ABI: set tidal acceleration (arcsec/cy²).
 */
void le_set_tid_acc(double acc);

/**
 * C ABI: load a VSOP2013 data file for a given planet (1 = Mercury .. 9 = Pluto).
 */
int32_t le_set_vsop2013_file(int32_t planet, const int8_t *path);

/**
 * C ABI: load ELP-MPP02 data files from a directory.
 */
int32_t le_set_elpmpp02_dir(const int8_t *dir);

/**
 * C ABI: close ephemeris files and reset context.
 */
void le_close(void);

/**
 * C ABI: set ELP-MPP02 adjustable coefficients.
 *
 * Pass a pointer to a 7-element `[f64; 7]` array:
 * `[fa, fb1, fb2, fb3, fb4, fb5, fb6]`.
 * Takes effect on the next call to le_set_elpmpp02_dir().
 */
int32_t le_set_elpmpp02_paras(const double *paras);

/**
 * C ABI: return library version string.
 */
const int8_t *le_version(void);

/**
 * Convert date to Julian day number, Gregorian calendar assumed.
 */
double le_julday(int32_t year, int32_t month, double day, int32_t gregflag);

/**
 * Convert Julian day number to calendar date.
 */
void le_revjul(double jd, int32_t gregflag, int32_t *year, int32_t *month, double *day);

/**
 * Compute day of week from Julian day.
 * Returns 0=Monday ... 6=Sunday.
 */
int32_t le_day_of_week(double jd);

/**
 * C ABI: coordinate transformation between systems.
 * iflag bits select source and target.
 * 0: equatorial J2000 -> ecliptic J2000
 * 1: ecliptic J2000 -> equatorial J2000
 * 2: equatorial J2000 -> equatorial of date
 * etc.
 */
void le_cotrans(const double *x,
                const double *y,
                const double *z,
                double eps,
                double *xout,
                double *yout,
                double *zout);

/**
 * C ABI: coordinate transformation for speed vectors (6-vector).
 */
void le_cotrans_sp(const double *x,
                   const double *y,
                   const double *z,
                   double eps,
                   double *xout,
                   double *yout,
                   double *zout);

/**
 * C ABI: split decimal degrees.
 */
void le_split_deg(double ddeg,
                  int32_t _roundflag,
                  int32_t *ideg,
                  int32_t *imin,
                  double *isec,
                  double *dsecfrac);

/**
 * C ABI: normalize angle to [0, 360).
 */
double le_csnorm(double d);

/**
 * C ABI: difference of two angles.
 */
double le_difdeg2n(double d1, double d2);

/**
 * C ABI: main ephemeris calculation (ET input).
 */
int32_t le_calc(double tjd, int32_t ipl, int32_t iflag, double *xx, int8_t *serr);

/**
 * C ABI: main ephemeris calculation (UT input).
 *
 * Computes the position and velocity of a celestial body at a given
 * Universal Time. Internally converts UT to ET using Delta T.
 *
 * # Arguments
 * * `tjd_ut` - Julian day in Universal Time
 * * `ipl` - Planet index (0=Sun, 1=Moon, 2=Mercury, ..., 9=Pluto, 10=Chiron, 17=Earth)
 * * `iflag` - Bitmask of LE_FLG_* flags controlling output frame, coordinate system, corrections
 * * `xx` - Output array of 24 doubles: [x, y, z, vx, vy, vz, ...]
 * * `serr` - Error string buffer (256 bytes), zeroed on success
 *
 * # Returns
 * 0 on success, negative on error.
 *
 * # Flags
 * - `LE_FLG_XYZ`: Cartesian output (default: polar)
 * - `LE_FLG_J2000`: J2000 frame (default: date)
 * - `LE_FLG_ECLIPTIC`: Ecliptic coordinates (default: equatorial)
 * - `LE_FLG_HELIO`: Heliocentric frame (default: geocentric)
 * - `LE_FLG_BARYHEL`: Barycentric frame
 * - `LE_FLG_NOABERR`: Skip annual aberration correction
 * - `LE_FLG_NOGDEFL`: Skip gravitational deflection
 * - `LE_FLG_NOBIRR`: Skip frame bias
 * - `LE_FLG_NONUT`: Skip nutation
 * - `LE_FLG_SPEED`: Include velocity in output
 * - `LE_FLG_SIDEREAL`: Apply sidereal (ayanamsa) correction
 * - `LE_FLG_TOPOCTR`: Topocentric position
 * - `LE_FLG_VSOP2013`: Use VSOP2013/ELP-MPP02 file-based engine
 * - `LE_FLG_JPLEPH`: Use JPL DE ephemeris
 */
int32_t le_calc_ut(double tjd_ut, int32_t ipl, int32_t iflag, double *xx, int8_t *serr);

/**
 * C ABI: solar eclipse calculation.
 */
int32_t le_sol_eclipse_how(double tjd, int32_t _iflag, double *attr, int8_t *serr);

/**
 * C ABI: lunar eclipse calculation.
 */
int32_t le_lun_eclipse_how(double tjd, int32_t _iflag, double *attr, int8_t *serr);

/**
 * C ABI: compute Delta T (TT - UT) in seconds for a given UT Julian day.
 */
double le_deltat(double tjd);

/**
 * C ABI: compute Delta T with extended interface (iflag, error string).
 */
double le_deltat_ex(double tjd, int32_t _iflag, int8_t *serr);

/**
 * C ABI: compute ayanamsa.
 */
double le_get_ayanamsa(double tjd_et);

/**
 * C ABI: compute ayanamsa with error string.
 */
double le_get_ayanamsa_ex(double tjd_et, int32_t _iflag, int8_t *serr);

/**
 * C ABI: get ayanamsa name.
 */
const int8_t *le_get_ayanamsa_name(int32_t isidmode);

/**
 * C ABI: compute houses.
 */
int32_t le_houses(double tjd_ut,
                  double geolat,
                  double geolon,
                  int32_t hsys,
                  double *cusps,
                  double *ascmc);

/**
 * C ABI: compute houses with ephemeris flag and error string.
 */
int32_t le_houses_ex(double tjd_ut,
                     int32_t _iflag,
                     double geolat,
                     double geolon,
                     int32_t hsys,
                     double *cusps,
                     double *ascmc,
                     int8_t *serr);

/**
 * C ABI: house position of an ecliptic longitude point.
 */
double le_house_pos(double _armc,
                    double _geolat,
                    double _eps,
                    int32_t _hsys,
                    double *xpin,
                    int8_t *serr);

/**
 * C ABI: get house system name.
 *
 * Returns a pointer to a static null-terminated string; no free needed.
 * Uses a static array indexed by (hsys & 0xFF) as a simple cache.
 */
const int8_t *le_house_name(int32_t hsys);

/**
 * C ABI: get star position at a given Julian day.
 */
int32_t le_fixstar(const int8_t *star_name, double jd_et, int32_t iflag, double *xx, int8_t *serr);

/**
 * C ABI: get star count in catalog.
 */
int32_t le_star_count(void);

/**
 * C ABI: get star data by index.
 */
int32_t le_star_data(int32_t index, struct LeStar *star_out);

/**
 * C ABI: compute rise, set, or transit time.
 *
 * For fast-moving bodies (Moon), iterates up to 3 times,
 * recomputing the body position at each estimated event time.
 */
double le_rise_trans(double tjd_ut,
                     int32_t ipl,
                     int32_t _iflag,
                     int32_t rsmi,
                     double geo_lon,
                     double geo_lat,
                     double *_attr,
                     int8_t *serr);
