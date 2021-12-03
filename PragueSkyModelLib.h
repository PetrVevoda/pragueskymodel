#define PSM_SUN_RADIUS                          0.2667 * 0.01745329251994329576923690768488612713443
#define PSM_PLANET_RADIUS                       6378000.0
#define PSM_PLANET_RADIUS_SQR                   PSM_PLANET_RADIUS * PSM_PLANET_RADIUS

#define PSM_MIN_ALTITUDE                              0.0
#define PSM_MAX_ALTITUDE                          15000.0
#define PSM_LIGHTCOLLECTION_VERTICAL_STEPSIZE       250.0
#define PSM_ARRAYSIZE                                61 // = PSM_MAX_ALTITUDE / PSM_LIGHTCOLLECTION_VERTICAL_STEPSIZE + 1

//#include "ART_Foundation_Geometry.h"
//#include "ART_Foundation_ColourAndSpectra.h"

#include <cmath>

class Vector3 {
public:
    double x, y, z;

    Vector3() {
        this->x = 0.0;
        this->y = 0.0;
        this->z = 0.0;
    }

    Vector3(double x, double y, double z) {
        this->x = x;
        this->y = y;
        this->z = z;
    }

	Vector3 operator+(const Vector3& other) const {
		return Vector3(x + other.x, y + other.y, z + other.z);
	}
	Vector3 operator-(const Vector3& other) const {
		return Vector3(x - other.x, y - other.y, z - other.z);
	}

	Vector3 operator*(const double factor) const {
		return Vector3(x * factor, y * factor, z * factor);
	}

	Vector3 operator/(const double factor) const {
		return Vector3(x / factor, y / factor, z / factor);
	}
};

float dot(const Vector3& v1, const Vector3& v2);

double magnitude(const Vector3& vector);

Vector3 normalize(const Vector3& vector);

//   This computes the canonical angles of the model from
//   a normalised view vector and solar elevation.

void arpragueskymodel_compute_angles(
        const Vector3 & viewpoint,
        const Vector3 & viewDirection,
        const double    groundLevelSolarElevationAtOrigin,
        const double    groundLevelSolarAzimuthAtOrigin,
              double  * solarElevationAtViewpoint,
              double  * altitudeOfViewpoint,
              double  * theta,
              double  * gamma,
              double  * shadow,
              double  * zero
        );

//   One blob of floats for each wavelength and task

typedef struct ArPragueSkyModelState
{
	// Radiance metadata

	int turbidities;
	double * turbidity_vals;

	int albedos;
	double * albedo_vals;

	int altitudes;
	double * altitude_vals;

	int elevations;
	double * elevation_vals;

	int channels;
	double channel_start;
	double channel_width;

	int tensor_components;

	int sun_nbreaks;
	int sun_offset;
	int sun_stride;
	double * sun_breaks;

	int zenith_nbreaks;
	int zenith_offset;
	int zenith_stride;
	double * zenith_breaks;

	int emph_nbreaks;
	int emph_offset;
	double * emph_breaks;

	int total_coefs_single_config; // this is for one specific configuration
	int total_coefs_all_configs;
	int total_configs;

	// Radiance data

	double * radiance_dataset;



    // Tranmittance metadata

	int     trans_n_a;
	int     trans_n_d;
	int     trans_turbidities;
	int     trans_altitudes;
	int     trans_rank;
	float * transmission_altitudes;
	float * transmission_turbities;

    // Tranmittance data

	float * transmission_dataset_U;
	float * transmission_dataset_V;



    // Polarisation metadata

	int tensor_components_pol;
     
    int sun_nbreaks_pol;
	int sun_offset_pol;
	int sun_stride_pol;
	double * sun_breaks_pol;

	int zenith_nbreaks_pol;
	int zenith_offset_pol;
	int zenith_stride_pol;
	double * zenith_breaks_pol;

	int total_coefs_single_config_pol; // this is for one specific configuration
	int total_coefs_all_configs_pol;

	// Polarisation data

	double * polarisation_dataset;
}
ArPragueSkyModelState;

ArPragueSkyModelState  * arpragueskymodelstate_alloc_init(
        const char  * library_path
        );

void arpragueskymodelstate_free(
        ArPragueSkyModelState  * state
        );

//   theta  - zenith angle
//   gamma  - sun angle
//   shadow - angle from the shadow point, which is further 90 degrees above the sun
//   zero   - angle from the zero point, which lies at the horizon below the sun
//   altitude
//   wavelength

double arpragueskymodel_radiance(
        const ArPragueSkyModelState  * state,
        const double                   theta,
        const double                   gamma,
        const double                   shadow,
        const double                   zero,
        const double                   elevation,
        const double                   altitude,
        const double                   turbidity,
        const double                   albedo,
        const double                   wavelength
        );

/* ----------------------------------------------------------------------------

    arpragueskymodel_solar_radiance
    ---------------------------

    This computes transmittance between a point at 'altitude' and infinity in
    the direction 'theta' at a wavelength 'wavelength'.

---------------------------------------------------------------------------- */


double arpragueskymodel_solar_radiance(
        const ArPragueSkyModelState  * state,
        const double                   theta,
        const double                   gamma,
        const double                   shadow,
        const double                   zero,
        const double                   elevation,
        const double                   altitude,
        const double                   turbidity,
        const double                   albedo,
        const double                   wavelength
        );

double arpragueskymodel_polarisation(
        const ArPragueSkyModelState  * state,
        const double                   theta,
        const double                   gamma,
        const double                   elevation,
        const double                   altitude,
        const double                   turbidity,
        const double                   albedo,
        const double                   wavelength
        );

/* ----------------------------------------------------------------------------

    arpragueskymodel_tau
    ------------------------------

    This computes transmittance between a point at 'altitude' and infinity in
    the direction 'theta' at a wavelength 'wavelength'.

---------------------------------------------------------------------------- */

void arpragueskymodel_toAD(
	double theta,
	double distance,
	double altitude,
	double *a,
	double *d
);

double arpragueskymodel_tau(
        const ArPragueSkyModelState  * state,
        const double                   theta,
        const double                   altitude,
        const double                   turbidity,
        const double                   wavelength,
        const double                   distance
        );