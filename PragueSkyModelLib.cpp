#include "PragueSkyModelLib.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cstring>
//#include <stdlib.h>
#include <io.h>

//   Some macro definitions that occur elsewhere in ART, and that have to be
//   replicated to make this a stand-alone module.

#ifndef MATH_PI
#define MATH_PI                     3.141592653589793
#endif

#ifndef MATH_RAD_TO_DEG
#define MATH_RAD_TO_DEG             ( 180.0 / MATH_PI )
#endif

#ifndef ALLOC
#define ALLOC(_struct)              ((_struct *)malloc(sizeof(_struct)))
#endif

#ifndef ALLOC_ARRAY
#define ALLOC_ARRAY(_struct,_number) \
            ((_struct *) malloc(sizeof(_struct) * (_number)))
#endif

#ifndef M_MAX
#define M_MAX(_a,_b)            ((_a) > (_b) ? (_a) : (_b))
#endif

#define MATH_HUGE_DOUBLE        5.78960446186580977117855E+76

#ifndef FREE
//#define FREE(_pointer)                ((void *)(free(_pointer), NULL))
#define FREE(_pointer) \
do { \
    void *_ptr=(void *)(_pointer); \
    free(_ptr); \
    _ptr=NULL; \
    _pointer=NULL; \
} while (0)
#endif

float dot(const Vector3 & v1, const Vector3 & v2) {
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

double magnitude(const Vector3& vector) {
	return std::sqrt(dot(vector, vector));
}

Vector3 normalize(const Vector3& vector) {
	return vector / magnitude(vector);
}




const int transsvdrank = 12;
const double safety_altitude = 50.0;

typedef unsigned short     uint16;
typedef unsigned int       uint32;
typedef unsigned long long uint64;

double double_from_half(const uint16 value)
{
	uint32 hi = uint32(value & 0x8000) << 16;
	uint16 abs = value & 0x7FFF;
	if (abs)
	{
		hi |= 0x3F000000 << uint16(abs >= 0x7C00);
		for (; abs < 0x400; abs <<= 1, hi -= 0x100000);
		hi += uint32(abs) << 10;
	}
	uint64 dbits = uint64(hi) << 32;
	double out;
	std::memcpy(&out, &dbits, sizeof(double));
	return out;
}

/*double double_from_half(const unsigned short value)
{
	unsigned long hi = (unsigned long)(value&0x8000) << 16;
	unsigned int abs = value & 0x7FFF;
	if(abs)
	{
		hi |= 0x3F000000 << (unsigned)(abs>=0x7C00);
		for(; abs<0x400; abs<<=1,hi-=0x100000) ;
		hi += (unsigned long)(abs) << 10;
	}
	unsigned long dbits = (unsigned long)(hi) << 32;
	double out;
    std::memcpy(&out, &dbits, sizeof(double));
	return out;
}*/

/*typedef unsigned short     uint16;
typedef unsigned int       uint32;
typedef unsigned long long uint64;

double doubleFromHalf(const uint16 halfBits) {
	uint32 hi = uint32(halfBits & 0x8000) << 16;
	uint16 abs = halfBits & 0x7FFF;
	if (abs) {
		hi |= 0x3F000000 << uint16(abs >= 0x7C00);
		for (; abs < 0x400; abs <<= 1, hi -= 0x100000) {
		}
		hi += uint32(abs) << 10;
	}
	const uint64 doubleBits = uint64(hi) << 32;
	return bitCast<double>(doubleBits);
}*/

int compute_pp_coefs_from_half(const int nbreaks, const double * breaks, const unsigned short * values, double * coefs, const int offset, const double scale)
{
	for (int i = 0; i < nbreaks - 1; ++i)
	{
		const double val1 = double_from_half(values[i+1]) / scale;
		const double val2 = double_from_half(values[i]) / scale;
		const double diff = val1 - val2;

		coefs[offset + 2 * i] = diff / (breaks[i+1] - breaks[i]);
		coefs[offset + 2 * i + 1]  = val2;
	}
	return 2 * nbreaks - 2;
}

int compute_pp_coefs_from_float(const int nbreaks, const double * breaks, const float * values, double * coefs, const int offset)
{
	for (int i = 0; i < nbreaks - 1; ++i)
	{
		coefs[offset + 2 * i] = ((double)values[i+1] - (double)values[i]) / (breaks[i+1] - breaks[i]);
		coefs[offset + 2 * i + 1] = (double)values[i];
	}
	return 2 * nbreaks - 2;
}

void printErrorAndExit(const char * message) {
	fprintf(stderr, "%s", message);
	fprintf(stderr, "\n");
	fflush(stderr);
	exit(-1);
}

void read_radiance(ArPragueSkyModelState * state, FILE * handle)
{
	// Read metadata

	// Structure of the metadata part of the data file:
	// turbidities       (1 * int),  turbidity_vals (turbidities * double),
	// albedos           (1 * int),  albedo_vals    (albedos * double),
	// altitudes         (1 * int),  altitude_vals  (altitudes * double),
	// elevations        (1 * int),  elevation_vals (elevations * double),
	// channels          (1 * int),  channel_start  (1 * double), channel_width (1 * double),
	// tensor_components (1 * int),
        // sun_nbreaks       (1 * int),  sun_breaks     (sun_nbreaks * double),
	// zenith_nbreaks    (1 * int),  zenith_breaks  (zenith_nbreaks * double),
	// emph_nbreaks      (1 * int),  emph_breaks    (emph_nbreaks * double)

	int valsRead;

	valsRead = fread(&state->turbidities, sizeof(int), 1, handle);
	if (valsRead != 1 || state->turbidities < 1) printErrorAndExit("Error reading sky model data: turbidities");

	state->turbidity_vals = ALLOC_ARRAY(double, state->turbidities);
	valsRead = fread(state->turbidity_vals, sizeof(double), state->turbidities, handle);
	if (valsRead != state->turbidities) printErrorAndExit("Error reading sky model data: turbidity_vals");

	valsRead = fread(&state->albedos, sizeof(int), 1, handle);
	if (valsRead != 1 || state->albedos < 1) printErrorAndExit("Error reading sky model data: albedos");

	state->albedo_vals = ALLOC_ARRAY(double, state->albedos);
	valsRead = fread(state->albedo_vals, sizeof(double), state->albedos, handle);
	if (valsRead != state->albedos) printErrorAndExit("Error reading sky model data: albedo_vals");

	valsRead = fread(&state->altitudes, sizeof(int), 1, handle);
	if (valsRead != 1 || state->altitudes < 1) printErrorAndExit("Error reading sky model data: altitudes");

	state->altitude_vals = ALLOC_ARRAY(double, state->altitudes);
	valsRead = fread(state->altitude_vals, sizeof(double), state->altitudes, handle);
	if (valsRead != state->altitudes) printErrorAndExit("Error reading sky model data: altitude_vals");

	valsRead = fread(&state->elevations, sizeof(int), 1, handle);
	if (valsRead != 1 || state->elevations < 1) printErrorAndExit("Error reading sky model data: elevations");

	state->elevation_vals = ALLOC_ARRAY(double, state->elevations);
	valsRead = fread(state->elevation_vals, sizeof(double), state->elevations, handle);
	if (valsRead != state->elevations) printErrorAndExit("Error reading sky model data: elevation_vals");

	valsRead = fread(&state->channels, sizeof(int), 1, handle);
	if (valsRead != 1 || state->channels < 1) printErrorAndExit("Error reading sky model data: channels");

	valsRead = fread(&state->channel_start, sizeof(double), 1, handle);
	if (valsRead != 1 || state->channel_start < 0) printErrorAndExit("Error reading sky model data: channel_start");

	valsRead = fread(&state->channel_width, sizeof(double), 1, handle);
	if (valsRead != 1 || state->channel_width <= 0) printErrorAndExit("Error reading sky model data: channel_width");

	valsRead = fread(&state->tensor_components, sizeof(int), 1, handle);
	if (valsRead != 1 || state->tensor_components < 1) printErrorAndExit("Error reading sky model data: tensor_components");

	valsRead = fread(&state->sun_nbreaks, sizeof(int), 1, handle);
	if (valsRead != 1 || state->sun_nbreaks < 2) printErrorAndExit("Error reading sky model data: sun_nbreaks");

	state->sun_breaks = ALLOC_ARRAY(double, state->sun_nbreaks);
	valsRead = fread(state->sun_breaks, sizeof(double), state->sun_nbreaks, handle);
	if (valsRead != state->sun_nbreaks) printErrorAndExit("Error reading sky model data: sun_breaks");

	valsRead = fread(&state->zenith_nbreaks, sizeof(int), 1, handle);
	if (valsRead != 1 || state->zenith_nbreaks < 2) printErrorAndExit("Error reading sky model data: zenith_nbreaks");

	state->zenith_breaks = ALLOC_ARRAY(double, state->zenith_nbreaks);
	valsRead = fread(state->zenith_breaks, sizeof(double), state->zenith_nbreaks, handle);
	if (valsRead != state->zenith_nbreaks) printErrorAndExit("Error reading sky model data: zenith_breaks");

	valsRead = fread(&state->emph_nbreaks, sizeof(int), 1, handle);
	if (valsRead != 1 || state->emph_nbreaks < 2) printErrorAndExit("Error reading sky model data: emph_nbreaks");

	state->emph_breaks = ALLOC_ARRAY(double, state->emph_nbreaks);
	valsRead = fread(state->emph_breaks, sizeof(double), state->emph_nbreaks, handle);
	if (valsRead != state->emph_nbreaks) printErrorAndExit("Error reading sky model data: emph_breaks");

	// Calculate offsets and strides

	state->sun_offset = 0;
	state->sun_stride = 2 * state->sun_nbreaks - 2 + 2 * state->zenith_nbreaks - 2;

	state->zenith_offset = state->sun_offset + 2 * state->sun_nbreaks - 2;
	state->zenith_stride = state->sun_stride;

	state->emph_offset = state->sun_offset + state->tensor_components * state->sun_stride;

	state->total_coefs_single_config = state->emph_offset + 2 * state->emph_nbreaks - 2; // this is for one specific configuration
	state->total_configs = state->channels * state->elevations * state->altitudes * state->albedos * state->turbidities;
	state->total_coefs_all_configs = state->total_coefs_single_config * state->total_configs;

	// Read data

	// Structure of the data part of the data file:
	// [[[[[[ sun_coefs (sun_nbreaks * half), zenith_scale (1 * double), zenith_coefs (zenith_nbreaks * half) ] * tensor_components, emph_coefs (emph_nbreaks * half) ]
	//   * channels ] * elevations ] * altitudes ] * albedos ] * turbidities

	int offset = 0;
	state->radiance_dataset = ALLOC_ARRAY(double, state->total_coefs_all_configs);

	if (1) {
		unsigned short * radiance_temp = ALLOC_ARRAY(unsigned short, M_MAX(state->sun_nbreaks, M_MAX(state->zenith_nbreaks, state->emph_nbreaks)));

		for (int con = 0; con < state->total_configs; ++con)
		{
			for (int tc = 0; tc < state->tensor_components; ++tc)
			{
				const double sun_scale = 1.0;
				valsRead = fread(radiance_temp, sizeof(unsigned short), state->sun_nbreaks, handle);
				if (valsRead != state->sun_nbreaks) printErrorAndExit("Error reading sky model data: sun_coefs");
				offset += compute_pp_coefs_from_half(state->sun_nbreaks, state->sun_breaks, radiance_temp, state->radiance_dataset, offset, sun_scale);

				double zenith_scale;
				valsRead = fread(&zenith_scale, sizeof(double), 1, handle);
				if (valsRead != 1) printErrorAndExit("Error reading sky model data: zenith_scale");

				valsRead = fread(radiance_temp, sizeof(unsigned short), state->zenith_nbreaks, handle);
				if (valsRead != state->zenith_nbreaks) printErrorAndExit("Error reading sky model data: zenith_coefs");
				offset += compute_pp_coefs_from_half(state->zenith_nbreaks, state->zenith_breaks, radiance_temp, state->radiance_dataset, offset, zenith_scale);
			}

			const double emph_scale = 1.0;
			valsRead = fread(radiance_temp, sizeof(unsigned short), state->emph_nbreaks, handle);
			if (valsRead != state->emph_nbreaks) printErrorAndExit("Error reading sky model data: emph_coefs");
			offset += compute_pp_coefs_from_half(state->emph_nbreaks, state->emph_breaks, radiance_temp, state->radiance_dataset, offset, emph_scale);
		}

		free(radiance_temp);
	}
	else {
		float * radiance_temp = ALLOC_ARRAY(float, M_MAX(state->sun_nbreaks, M_MAX(state->zenith_nbreaks, state->emph_nbreaks)));

		for (int con = 0; con < state->total_configs; ++con)
		{
			for (int tc = 0; tc < state->tensor_components; ++tc)
			{
				fread(radiance_temp, sizeof(float), state->sun_nbreaks, handle);
				offset += compute_pp_coefs_from_float(state->sun_nbreaks, state->sun_breaks, radiance_temp, state->radiance_dataset, offset);

				fread(radiance_temp, sizeof(float), state->zenith_nbreaks, handle);
				offset += compute_pp_coefs_from_float(state->zenith_nbreaks, state->zenith_breaks, radiance_temp, state->radiance_dataset, offset);
			}

			fread(radiance_temp, sizeof(float), state->emph_nbreaks, handle);
			offset += compute_pp_coefs_from_float(state->emph_nbreaks, state->emph_breaks, radiance_temp, state->radiance_dataset, offset);
		}

		free(radiance_temp);
	}
}

void read_transmittance(ArPragueSkyModelState * state, FILE * handle)
{
	// Read metadata

	int valsRead;

	valsRead = fread(&state->trans_n_d, sizeof(int), 1, handle);
	if (valsRead != 1 || state->trans_n_d < 1) printErrorAndExit("Error reading sky model data: trans_n_d");

	valsRead = fread(&state->trans_n_a, sizeof(int), 1, handle);
	if (valsRead != 1 || state->trans_n_a < 1) printErrorAndExit("Error reading sky model data: trans_n_a");

	valsRead = fread(&state->trans_turbidities, sizeof(int), 1, handle);
	if (valsRead != 1 || state->trans_turbidities < 1) printErrorAndExit("Error reading sky model data: trans_turbidities");

	valsRead = fread(&state->trans_altitudes, sizeof(int), 1, handle);
	if (valsRead != 1 || state->trans_altitudes < 1) printErrorAndExit("Error reading sky model data: trans_altitudes");

	valsRead = fread(&state->trans_rank, sizeof(int), 1, handle);
	if (valsRead != 1 || state->trans_rank < 1) printErrorAndExit("Error reading sky model data: trans_rank");

	state->transmission_altitudes = ALLOC_ARRAY(float, state->trans_altitudes);
	valsRead = fread(state->transmission_altitudes, sizeof(float), state->trans_altitudes, handle);
	if (valsRead != state->trans_altitudes) printErrorAndExit("Error reading sky model data: transmission_altitudes");

	state->transmission_turbities = ALLOC_ARRAY(float, state->trans_turbidities);
	valsRead = fread(state->transmission_turbities, sizeof(float), state->trans_turbidities, handle);
	if (valsRead != state->trans_turbidities) printErrorAndExit("Error reading sky model data: transmission_turbities");

	const int total_coefs_U = state->trans_n_d * state->trans_n_a * state->trans_rank * state->trans_altitudes;
	const int total_coefs_V = state->trans_turbidities * state->trans_rank * 11 * state->trans_altitudes;

	// Read data

	state->transmission_dataset_U = ALLOC_ARRAY(float, total_coefs_U);
	valsRead = fread(state->transmission_dataset_U, sizeof(float), total_coefs_U, handle);
	if (valsRead != total_coefs_U) printErrorAndExit("Error reading sky model data: transmission_dataset_U");

	state->transmission_dataset_V = ALLOC_ARRAY(float, total_coefs_V);
	valsRead = fread(state->transmission_dataset_V, sizeof(float), total_coefs_V, handle);
	if (valsRead != total_coefs_V) printErrorAndExit("Error reading sky model data: transmission_dataset_V");
}

void read_polarisation(ArPragueSkyModelState * state, FILE * handle)
{
	// Read metadata

	// Structure of the metadata part of the data file:
	// tensor_components_pol (1 * int),
        // sun_nbreaks_pol       (1 * int),  sun_breaks_pol     (sun_nbreaks_pol * double),
	// zenith_nbreaks_pol    (1 * int),  zenith_breaks_pol  (zenith_nbreaks_pol * double),
	// emph_nbreaks_pol      (1 * int),  emph_breaks_pol    (emph_nbreaks_pol * double)

	int valsRead;

	valsRead = fread(&state->tensor_components_pol, sizeof(int), 1, handle);
	if (valsRead != 1)
	{
		// Polarisation dataset not present
		state->tensor_components_pol = 0;
        printErrorAndExit("No polarisation dataset available!\n");
		return;
	}

	valsRead = fread(&state->sun_nbreaks_pol, sizeof(int), 1, handle);
	if (valsRead != 1 || state->sun_nbreaks_pol < 1) printErrorAndExit("Error reading sky model data: sun_nbreaks_pol");

	state->sun_breaks_pol = ALLOC_ARRAY(double, state->sun_nbreaks_pol);
	valsRead = fread(state->sun_breaks_pol, sizeof(double), state->sun_nbreaks_pol, handle);
	if (valsRead != state->sun_nbreaks_pol) printErrorAndExit("Error reading sky model data: sun_breaks_pol");

	valsRead = fread(&state->zenith_nbreaks_pol, sizeof(int), 1, handle);
	if (valsRead != 1 || state->zenith_nbreaks_pol < 1) printErrorAndExit("Error reading sky model data: zenith_nbreaks_pol");

	state->zenith_breaks_pol = ALLOC_ARRAY(double, state->zenith_nbreaks_pol);
	valsRead = fread(state->zenith_breaks_pol, sizeof(double), state->zenith_nbreaks_pol, handle);
	if (valsRead != state->zenith_nbreaks_pol) printErrorAndExit("Error reading sky model data: zenith_breaks_pol");

	// Calculate offsets and strides

	state->sun_offset_pol = 0;
	state->sun_stride_pol = 2 * state->sun_nbreaks_pol - 2 + 2 * state->zenith_nbreaks_pol - 2;

	state->zenith_offset_pol = state->sun_offset_pol + 2 * state->sun_nbreaks_pol - 2;
	state->zenith_stride_pol = state->sun_stride_pol;

	state->total_coefs_single_config_pol = state->sun_offset_pol + state->tensor_components_pol * state->sun_stride_pol; // this is for one specific configuration
	state->total_coefs_all_configs_pol = state->total_coefs_single_config_pol * state->total_configs;

	// Read data

	// Structure of the data part of the data file:
	// [[[[[[ sun_coefs_pol (sun_nbreaks_pol * float), zenith_coefs_pol (zenith_nbreaks_pol * float) ] * tensor_components_pol]
	//   * channels ] * elevations ] * altitudes ] * albedos ] * turbidities

	int offset = 0;
	state->polarisation_dataset = ALLOC_ARRAY(double, state->total_coefs_all_configs_pol);
	float * polarisation_temp = ALLOC_ARRAY(float, M_MAX(state->sun_nbreaks_pol, state->zenith_nbreaks_pol));

	for (int con = 0; con < state->total_configs; ++con)
	{
		for (int tc = 0; tc < state->tensor_components_pol; ++tc)
		{
			valsRead = fread(polarisation_temp, sizeof(float), state->sun_nbreaks_pol, handle);
			if (valsRead != state->sun_nbreaks_pol) printErrorAndExit("Error reading sky model data: sun_coefs_pol");
			offset += compute_pp_coefs_from_float(state->sun_nbreaks_pol, state->sun_breaks_pol, polarisation_temp, state->polarisation_dataset, offset);

			valsRead = fread(polarisation_temp, sizeof(float), state->zenith_nbreaks_pol, handle);
			if (valsRead != state->zenith_nbreaks_pol) printErrorAndExit("Error reading sky model data: zenith_coefs_pol");
			offset += compute_pp_coefs_from_float(state->zenith_nbreaks_pol, state->zenith_breaks_pol, polarisation_temp, state->polarisation_dataset, offset);
		}
	}

	free(polarisation_temp);
}

ArPragueSkyModelState  * arpragueskymodelstate_alloc_init(
	const char  * library_path
	)
{
	ArPragueSkyModelState * state = ALLOC(ArPragueSkyModelState);

	char filename[1024];
    
	sprintf(filename, "%s", library_path);
    
    if ( _access(filename, 0 | 4) != 0 )
    {
        printErrorAndExit("Sky model dataset not found");
    }
    
	FILE * handle = fopen(filename, "rb");

	// Read data
	read_radiance(state, handle);
	read_transmittance(state, handle);
	read_polarisation(state, handle);

	fclose(handle);

	return state;
}

void arpragueskymodelstate_free(
	ArPragueSkyModelState  * state
	)
{
	free(state->turbidity_vals);
	free(state->albedo_vals);
	free(state->altitude_vals);
	free(state->elevation_vals);

	free(state->sun_breaks);
	free(state->zenith_breaks);
	free(state->emph_breaks);
	free(state->radiance_dataset);

	free(state->transmission_dataset_U);
	free(state->transmission_dataset_V);
	free(state->transmission_altitudes);
	free(state->transmission_turbities);

	if (state->tensor_components_pol > 0)
	{
		free(state->sun_breaks_pol);
		free(state->zenith_breaks_pol);
		free(state->polarisation_dataset);
	}

	FREE(state);
}

void arpragueskymodel_compute_altitude_and_elevation(
        const Vector3 & viewpoint,
        const double    groundLevelSolarElevationAtOrigin,
        const double    groundLevelSolarAzimuthAtOrigin,
              double  * solarElevationAtViewpoint,
              double  * altitudeOfViewpoint,
              double  * distanceToView,
			  Vector3 * directionToZenithN,
	          Vector3 * directionToSunN
        )
{
    // Direction to zenith

    Vector3 centerOfTheEarth = Vector3(0.0, 0.0, -PSM_PLANET_RADIUS);
	Vector3 directionToZenith = viewpoint - centerOfTheEarth;
	*directionToZenithN = normalize(directionToZenith);

    // Altitude of viewpoint

    *distanceToView = magnitude(directionToZenith);
    //ASSERT_DOUBLE_LARGER_THAN(*distanceToView, -0.0001);
    *distanceToView = M_MAX(*distanceToView, 0.0);

    *altitudeOfViewpoint = *distanceToView - PSM_PLANET_RADIUS;
    *altitudeOfViewpoint = M_MAX( *altitudeOfViewpoint, 0.0 );


    // Direction to sun

	directionToSunN->x = cos( groundLevelSolarAzimuthAtOrigin )
                       * cos( groundLevelSolarElevationAtOrigin );
	directionToSunN->y = sin( groundLevelSolarAzimuthAtOrigin )
                       * cos( groundLevelSolarElevationAtOrigin );
	directionToSunN->z = sin( groundLevelSolarElevationAtOrigin );


    // Solar elevation at viewpoint (more precisely, solar elevation at the point on the ground directly below viewpoint)

    const double dotZenithSun =
        dot(
            *directionToZenithN,
            *directionToSunN
            );

    *solarElevationAtViewpoint = 0.5 * MATH_PI - acos(dotZenithSun);
}

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
        )
{
    //ASSERT_VALID_DOUBLE(groundLevelSolarElevationAtOrigin);
    //ASSERT_VALID_DOUBLE(groundLevelSolarAzimuthAtOrigin);

    // Shift viewpoint about safety altitude up

	const Vector3 centerOfTheEarth = Vector3(0.0, 0.0, -PSM_PLANET_RADIUS);
	Vector3 toViewpoint = viewpoint - centerOfTheEarth;
	Vector3 toViewpointN = normalize(toViewpoint);
    const double distanceToViewTmp = magnitude(toViewpoint) + safety_altitude;
	Vector3 toShiftedViewpoint = toViewpointN * distanceToViewTmp;
	Vector3 shiftedViewpoint = centerOfTheEarth + toShiftedViewpoint;

	Vector3 viewDirectionN = normalize(viewDirection);

    double distanceToView;
	Vector3 directionToZenithN;
	Vector3 directionToSunN;

    arpragueskymodel_compute_altitude_and_elevation(
          shiftedViewpoint,
          groundLevelSolarElevationAtOrigin,
          groundLevelSolarAzimuthAtOrigin,
          solarElevationAtViewpoint,
          altitudeOfViewpoint,
        & distanceToView,
        & directionToZenithN,
        & directionToSunN
        );


    // Altitude-corrected view direction

	Vector3 correctViewN;
    if (distanceToView > PSM_PLANET_RADIUS)
    {
		Vector3 lookAt = shiftedViewpoint + viewDirectionN;

        const double correction = sqrt(distanceToView * distanceToView - PSM_PLANET_RADIUS_SQR) / distanceToView;

		Vector3 toNewOrigin = directionToZenithN * (distanceToView - correction);
		Vector3 newOrigin = centerOfTheEarth + toNewOrigin;
		Vector3 correctView = lookAt - newOrigin;

		correctViewN = normalize(correctView);
    }
    else
    {
        correctViewN = viewDirectionN;
    }


    // Sun angle (gamma) - no correction

    double  dotProductSun =
        dot(
            viewDirectionN,
            directionToSunN
            );

    *gamma = acos(dotProductSun);


    // Shadow angle - requires correction

    const double effectiveElevation = groundLevelSolarElevationAtOrigin;
    const double effectiveAzimuth = groundLevelSolarAzimuthAtOrigin;
    const double shadow_angle = effectiveElevation + MATH_PI * 0.5;

	Vector3 shadowDirectionN = Vector3(cos(shadow_angle) * cos(effectiveAzimuth), cos(shadow_angle) * sin(effectiveAzimuth), sin(shadow_angle));

    const double  dotProductShadow  =
        dot(
            correctViewN,
            shadowDirectionN
            );

    *shadow = acos(dotProductShadow);


    // Zenith angle (theta) - corrected version stored in otherwise unused zero angle

    double  cosThetaCor =
        dot(
            correctViewN,
            directionToZenithN
            );

    *zero  = acos(cosThetaCor);


    // Zenith angle (theta) - uncorrected version goes outside

    double  cosTheta =
        dot(
            viewDirectionN,
            directionToZenithN
            );

    *theta  = acos(cosTheta);
}

double lerp(const double from, const double to, const double factor)
{
	return (1.0 - factor) * from + factor * to;
}

int find_segment(const double x, const int nbreaks, const double* breaks)
{
	int segment = 0;
	for (segment = 0; segment < nbreaks; ++segment)
	{
		if (breaks[segment+1] >= x)
		break;
	}
	return segment;
}

double eval_pp(const double x, const int segment, const double * breaks, const double * coefs)
{
	const double x0 = x - breaks[segment];
	const double * sc = coefs + 2 * segment; // segment coefs
	return sc[0] * x0 + sc[1];
}

const double * control_params_single_config(
	const ArPragueSkyModelState * state,
	const double *                dataset,
	const int                     total_coefs_single_config,
	const int                     elevation,
	const int                     altitude,
	const int                     turbidity,
	const int                     albedo,
	const int                     wavelength
)
{
	return dataset + (total_coefs_single_config * (
		wavelength +
		state->channels*elevation +
		state->channels*state->elevations*altitude +
		state->channels*state->elevations*state->altitudes*albedo +
		state->channels*state->elevations*state->altitudes*state->albedos*turbidity
	));
}

double reconstruct(
  const ArPragueSkyModelState  * state,
  const double                   gamma,
  const double                   alpha,
  const double                   zero,
  const int                      gamma_segment,
  const int                      alpha_segment,
  const int                      zero_segment,
  const double                 * control_params
)
{
  double res = 0.0;
  for (int t = 0; t < state->tensor_components; ++t) {
	const double sun_val_t = eval_pp(gamma, gamma_segment, state->sun_breaks, control_params + state->sun_offset + t * state->sun_stride);
	const double zenith_val_t = eval_pp(alpha, alpha_segment, state->zenith_breaks, control_params + state->zenith_offset + t * state->zenith_stride);
	res += sun_val_t * zenith_val_t;
  }
  const double emph_val_t = eval_pp(zero, zero_segment, state->emph_breaks, control_params + state->emph_offset);
  res *= emph_val_t;

  return M_MAX(res, 0.0);
}

double map_parameter(const double param, const int value_count, const double * values)
{
	double mapped = 0.0;
    
	if (param < values[0])
	{
		mapped = 0.0;
	}
	else if (param > values[value_count - 1])
	{
		mapped = (double)value_count - 1.0;
	}
	else
	{
		for (int v = 0; v < value_count; ++v)
		{
			const double val = values[v];
			if (fabs(val - param) < 1e-6)
			{
				mapped = v;
				break;
			}
			else if (param < val)
			{
				mapped = v - ((val - param) / (val - values[v - 1]));
				break;
			}
		}
	}

	return mapped;
}


///////////////////////////////////////////////
// Radiance mono version
///////////////////////////////////////////////


double interpolate_elevation(
  const ArPragueSkyModelState  * state,
  double                  elevation,
  int                     altitude,
  int                     turbidity,
  int                     albedo,
  int                     wavelength,
  double                  gamma,
  double                  alpha,
  double                  zero,
  int                     gamma_segment,
  int                     alpha_segment,
  int                     zero_segment
)
{
  const int elevation_low = (int)elevation;
  const double factor = elevation - (double)elevation_low;

  const double * control_params_low = control_params_single_config(
    state,
    state->radiance_dataset,
    state->total_coefs_single_config,
    elevation_low,
    altitude,
    turbidity,
    albedo,
    wavelength);

  double res_low = reconstruct(
    state,
    gamma,
    alpha,
    zero,
    gamma_segment,
    alpha_segment,
    zero_segment,
    control_params_low);

  if (factor < 1e-6 || elevation_low >= (state->elevations - 1))
  {
    return res_low;
  }

  const double * control_params_high = control_params_single_config(
    state,
    state->radiance_dataset,
    state->total_coefs_single_config,
    elevation_low+1,
    altitude,
    turbidity,
    albedo,
    wavelength);

  double res_high = reconstruct(
    state,
    gamma,
    alpha,
    zero,
    gamma_segment,
    alpha_segment,
    zero_segment,
    control_params_high);

  return lerp(res_low, res_high, factor);
}

double interpolate_altitude(
  const ArPragueSkyModelState  * state,
  double                  elevation,
  double                  altitude,
  int                     turbidity,
  int                     albedo,
  int                     wavelength,
  double                  gamma,
  double                  alpha,
  double                  zero,
  int                     gamma_segment,
  int                     alpha_segment,
  int                     zero_segment
)
{
  const int altitude_low = (int)altitude;
  const double factor = altitude - (double)altitude_low;

  double res_low = interpolate_elevation(
    state,
    elevation,
    altitude_low,
    turbidity,
    albedo,
    wavelength,
    gamma,
    alpha,
    zero,
    gamma_segment,
    alpha_segment,
    zero_segment);

  if (factor < 1e-6 || altitude_low >= (state->altitudes - 1))
  {
    return res_low;
  }

  double res_high = interpolate_elevation(
    state,
    elevation,
    altitude_low + 1,
    turbidity,
    albedo,
    wavelength,
    gamma,
    alpha,
    zero,
    gamma_segment,
    alpha_segment,
    zero_segment);

  return lerp(res_low, res_high, factor);
}

double interpolate_turbidity(
  const ArPragueSkyModelState  * state,
  double                  elevation,
  double                  altitude,
  double                  turbidity,
  int                     albedo,
  int                     wavelength,
  double                  gamma,
  double                  alpha,
  double                  zero,
  int                     gamma_segment,
  int                     alpha_segment,
  int                     zero_segment
)
{
  const int turbidity_low = (int)turbidity;
  const double factor = turbidity - (double)turbidity_low;

  double res_low = interpolate_altitude(
    state,
    elevation,
    altitude,
    turbidity_low,
    albedo,
    wavelength,
    gamma,
    alpha,
    zero,
    gamma_segment,
    alpha_segment,
    zero_segment);

  if (factor < 1e-6 || turbidity_low >= (state->turbidities - 1))
  {
    return res_low;
  }

  double res_high = interpolate_altitude(
    state,
    elevation,
    altitude,
    turbidity_low + 1,
    albedo,
    wavelength,
    gamma,
    alpha,
    zero,
    gamma_segment,
    alpha_segment,
    zero_segment);

  return lerp(res_low, res_high, factor);
}

double interpolate_albedo(
  const ArPragueSkyModelState  * state,
  double                  elevation,
  double                  altitude,
  double                  turbidity,
  double                  albedo,
  int                     wavelength,
  double                  gamma,
  double                  alpha,
  double                  zero,
  int                     gamma_segment,
  int                     alpha_segment,
  int                     zero_segment
)
{
  const int albedo_low = (int)albedo;
  const double factor = albedo - (double)albedo_low;

  double res_low = interpolate_turbidity(
    state,
    elevation,
    altitude,
    turbidity,
    albedo_low,
    wavelength,
    gamma,
    alpha,
    zero,
    gamma_segment,
    alpha_segment,
    zero_segment);

  if (factor < 1e-6 || albedo_low >= (state->albedos - 1))
  {
    return res_low;
  }

  double res_high = interpolate_turbidity(
    state,
    elevation,
    altitude,
    turbidity,
    albedo_low + 1,
    wavelength,
    gamma,
    alpha,
    zero,
    gamma_segment,
    alpha_segment,
    zero_segment);

  return lerp(res_low, res_high, factor);
}

double interpolate_wavelength(
  const ArPragueSkyModelState  * state,
  double                  elevation,
  double                  altitude,
  double                  turbidity,
  double                  albedo,
  double                  wavelength,
  double                  gamma,
  double                  alpha,
  double                  zero,
  int                     gamma_segment,
  int                     alpha_segment,
  int                     zero_segment
)
{
  // Don't interpolate, use the bin it belongs to

  return interpolate_albedo(
    state,
    elevation,
    altitude,
    turbidity,
    albedo,
    (int)wavelength,
    gamma,
    alpha,
    zero,
    gamma_segment,
    alpha_segment,
    zero_segment);
}

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
        )
{
  (void) theta;
    
  // Translate parameter values to indices

  const double turbidity_control = map_parameter(turbidity, state->turbidities, state->turbidity_vals);
  const double albedo_control    = map_parameter(albedo, state->albedos, state->albedo_vals);
  const double altitude_control  = map_parameter(altitude, state->altitudes, state->altitude_vals);
  const double elevation_control = map_parameter(elevation * MATH_RAD_TO_DEG, state->elevations, state->elevation_vals);

  const double channel_control = (wavelength - state->channel_start) / state->channel_width;

  if ( channel_control >= state->channels || channel_control < 0.) return 0.;

  // Get params corresponding to the indices, reconstruct result and interpolate

  const double alpha = elevation < 0 ? shadow : zero;

  const int gamma_segment = find_segment(gamma, state->sun_nbreaks, state->sun_breaks);
  const int alpha_segment = find_segment(alpha, state->zenith_nbreaks, state->zenith_breaks);
  const int zero_segment = find_segment(zero, state->emph_nbreaks, state->emph_breaks);

  const double res = interpolate_wavelength(
    state,
    elevation_control,
    altitude_control,
    turbidity_control,
    albedo_control,
    channel_control,
    gamma,
    alpha,
    zero,
    gamma_segment,
    alpha_segment,
    zero_segment);

  //ASSERT_NONNEGATIVE_DOUBLE(res);

  return res;
}


const double SunRadStartWL = 310;
const double SunRadIncrementWL = 1;
const double SunRadTable[] =
{
9829.41, 10184., 10262.6, 10375.7, 10276., 10179.3, 10156.6, 10750.7, 11134., 11463.6, 11860.4, 12246.2, 12524.4, 12780., 13187.4, 13632.4, 13985.9, 13658.3, 13377.4, 13358.3, 13239., 13119.8, 13096.2, 13184., 13243.5, 13018.4, 12990.4, 13159.1, 13230.8, 13258.6, 13209.9, 13343.2, 13404.8, 13305.4, 13496.3, 13979.1, 14153.8, 14188.4, 14122.7, 13825.4, 14033.3, 13914.1, 13837.4, 14117.2, 13982.3, 13864.5, 14118.4, 14545.7, 15029.3, 15615.3, 15923.5, 16134.8, 16574.5, 16509., 16336.5, 16146.6, 15965.1, 15798.6, 15899.8, 16125.4, 15854.3, 15986.7, 15739.7, 15319.1, 15121.5, 15220.2, 15041.2, 14917.7, 14487.8, 14011., 14165.7, 14189.5, 14540.7, 14797.5, 14641.5, 14761.6, 15153.7, 14791.8, 14907.6, 15667.4, 16313.5, 16917., 17570.5, 18758.1, 20250.6, 21048.1, 21626.1, 22811.6, 23577.2, 23982.6, 24062.1, 23917.9, 23914.1, 23923.2, 24052.6, 24228.6, 24360.8, 24629.6, 24774.8, 24648.3, 24666.5, 24938.6, 24926.3, 24693.1, 24613.5, 24631.7, 24569.8, 24391.5, 24245.7, 24084.4, 23713.7, 22985.4, 22766.6, 22818.9, 22834.3, 22737.9, 22791.6, 23086.3, 23377.7, 23461., 23935.5, 24661.7, 25086.9, 25520.1, 25824.3, 26198., 26350.2, 26375.4, 26731.2, 27250.4, 27616., 28145.3, 28405.9, 28406.8, 28466.2, 28521.5, 28783.8, 29025.1, 29082.6, 29081.3, 29043.1, 28918.9, 28871.6, 29049., 29152.5, 29163.2, 29143.4, 28962.7, 28847.9, 28854., 28808.7, 28624.1, 28544.2, 28461.4, 28411.1, 28478., 28469.8, 28513.3, 28586.5, 28628.6, 28751.5, 28948.9, 29051., 29049.6, 29061.7, 28945.7, 28672.8, 28241.5, 27903.2, 27737., 27590.9, 27505.6, 27270.2, 27076.2, 26929.1, 27018.2, 27206.8, 27677.2, 27939.9, 27923.9, 27899.2, 27725.4, 27608.4, 27599.4, 27614.6, 27432.4, 27460.4, 27392.4, 27272., 27299.1, 27266.8, 27386.5, 27595.9, 27586.9, 27504.8, 27480.6, 27329.8, 26968.4, 26676.3, 26344.7, 26182.5, 26026.3, 25900.3, 25842.9, 25885.4, 25986.5, 26034.5, 26063.5, 26216.9, 26511.4, 26672.7, 26828.5, 26901.8, 26861.5, 26865.4, 26774.2, 26855.8, 27087.1, 27181.3, 27183.1, 27059.8, 26834.9, 26724.3, 26759.6, 26725.9, 26724.6, 26634.5, 26618.5, 26560.1, 26518.7, 26595.3, 26703.2, 26712.7, 26733.9, 26744.3, 26764.4, 26753.2, 26692.7, 26682.7, 26588.1, 26478., 26433.7, 26380.7, 26372.9, 26343.3, 26274.7, 26162.3, 26160.5, 26210., 26251.2, 26297.9, 26228.9, 26222.3, 26269.7, 26295.6, 26317.9, 26357.5, 26376.1, 26342.4, 26303.5, 26276.7, 26349.2, 26390., 26371.6, 26346.7, 26327.6, 26274.2, 26247.3, 26228.7, 26152.1, 25910.3, 25833.2, 25746.5, 25654.3, 25562., 25458.8, 25438., 25399.1, 25324.3, 25350., 25514., 25464.9, 25398.5, 25295.2, 25270.2, 25268.4, 25240.6, 25184.9, 25149.6, 25123.9, 25080.3, 25027.9, 25012.3, 24977.9, 24852.6, 24756.4, 24663.5, 24483.6, 24398.6, 24362.6, 24325.1, 24341.7, 24288.7, 24284.2, 24257.3, 24178.8, 24097.6, 24175.6, 24175.7, 24139.7, 24088.1, 23983.2, 23902.7, 23822.4, 23796.2, 23796.9, 23814.5, 23765.5, 23703., 23642., 23592.6, 23552., 23514.6, 23473.5, 23431., 23389.3, 23340., 23275.1, 23187.3, 23069.5, 22967., 22925.3, 22908.9, 22882.5, 22825., 22715.4, 22535.5, 22267.1, 22029.4, 21941.6, 21919.5, 21878.8, 21825.6, 21766., 21728.9, 21743.2, 21827.1, 21998.7, 22159.4, 22210., 22187.2, 22127.2, 22056.2, 22000.2, 21945.9, 21880.2, 21817.1, 21770.3, 21724.3, 21663.2, 21603.3, 21560.4, 21519.8, 21466.2, 21401.6, 21327.7, 21254.2, 21190.7, 21133.6, 21079.3, 21024., 20963.7, 20905.5, 20856.6, 20816.6, 20785.2, 20746.7, 20685.3, 20617.8, 20561.1, 20500.4, 20421.2, 20333.4, 20247., 20175.3, 20131.4, 20103.2, 20078.5, 20046.8, 19997.2, 19952.9, 19937.2, 19930.8, 19914.4, 19880.8, 19823., 19753.8, 19685.9, 19615.3, 19537.5, 19456.8, 19377.6, 19309.4, 19261.9, 19228., 19200.5, 19179.5, 19164.8, 19153.1, 19140.6, 19129.2, 19120.6, 19104.5, 19070.6, 19023.9, 18969.3, 18911.4, 18855., 18798.6, 18740.8, 18672.7, 18585.2, 18501., 18442.4, 18397.5, 18353.9, 18313.2, 18276.8, 18248.3, 18231.2, 18224., 18225.4, 18220.1, 18192.6, 18155.1, 18119.8, 18081.6, 18035.6, 17987.4, 17942.8, 17901.7, 17864.2, 17831.1, 17802.9, 17771.5, 17728.6, 17669.7, 17590.1, 17509.5, 17447.4, 17396., 17347.4, 17300.3, 17253.2, 17206.1, 17159., 17127.6, 17127.6, 17133.6, 17120.4, 17097.2, 17073.3, 17043.7, 17003.4, 16966.3, 16946.3, 16930.9, 16907.7, 16882.7, 16862., 16837.8, 16802.1, 16759.2, 16713.6, 16661.8, 16600.8, 16542.6, 16499.4, 16458.7, 16408., 16360.6, 16329.5, 16307.4, 16286.7, 16264.9, 16239.6, 16207.8, 16166.8, 16118.2, 16064., 16011.2, 15966.9, 15931.9, 15906.9, 15889.1, 15875.5, 15861.2, 15841.3, 15813.1, 15774.2, 15728.8, 15681.4, 15630., 15572.9, 15516.5, 15467.2, 15423., 15381.6, 15354.4, 15353., 15357.3, 15347.3, 15320.2, 15273.1, 15222., 15183.1, 15149.6, 15114.6, 15076.8, 15034.6, 14992.9
};

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
        )
{
    (void) gamma;
    (void) shadow;
    (void) zero;
    (void) elevation;
    (void) albedo;
        double idx = (wavelength - SunRadStartWL) / SunRadIncrementWL;
        double sunRadiance = 0.0;

        if ( idx >= 0.0 )
        {
            int lowIdx = floor(idx);
            double idxFloat = idx - floor(idx);

            sunRadiance = SunRadTable[lowIdx] * (1.0 - idxFloat) + SunRadTable[lowIdx+1] * idxFloat;
            //ASSERT_POSITIVE_DOUBLE(sunRadiance);
        }

	double tau = arpragueskymodel_tau(
	      state,
	      theta,
	      altitude,
	      turbidity,
	      wavelength,
	      MATH_HUGE_DOUBLE
	    );
        //ASSERT_UNIT_RANGE_DOUBLE(tau);

	return sunRadiance * tau;
}



///////////////////////////////////////////////
// Transmittance
///////////////////////////////////////////////


int arpragueskymodel_circleBounds2D(
	double x_v,
	double y_v,
	double y_c,
	double radius,
	double *d
)
{
	double qa = (x_v * x_v) + (y_v * y_v);
	double qb = 2.0 * y_c * y_v;
	double qc = (y_c * y_c) - (radius * radius);
	double n = (qb * qb) - (4.0 * qa * qc);
	if (n <= 0)
	{
		return 0;
	}
	float d1;
	float d2;
	n = sqrt(n);
	d1 = (-qb + n) / (2.0 * qa);
	d2 = (-qb - n) / (2.0 * qa);
	*d = (d1 > 0 && d2 > 0) ? (d1 < d2 ? d1 : d2) : (d1 > d2 ? d1 : d2); // It fits in one line.
	if (*d <= 0)
	{
		return 0;
	}
	return 1;
}

void arpragueskymodel_scaleAD(
	double x_p,
	double y_p,
	double *a,
	double *d
)
{
	double n;
	n = sqrt((x_p * x_p) + (y_p * y_p));
	*a = n - PSM_PLANET_RADIUS;
	*a = *a > 0 ? *a : 0;
	*a = pow(*a / 100000.0, 1.0 / 3.0);
	*d = acos(y_p / n) * PSM_PLANET_RADIUS;
	*d = *d / 1571524.413613; // Maximum distance to the edge of the atmosphere in the transmittance model
	*d = pow(*d, 0.25);
	*d = *d > 1.0 ? 1.0 : *d;
}

void arpragueskymodel_toAD(
	double theta,
	double distance,
	double altitude,
	double *a,
	double *d
)
{
	// Ray circle intersection
	double x_v = sin(theta);
	double y_v = cos(theta);
	double y_c = PSM_PLANET_RADIUS + altitude;
	double atmo_edge = PSM_PLANET_RADIUS + 90000;
	double n;
        if (altitude < 0.001) // Handle altitudes close to 0 separately to avoid reporting intersections on the other side of the planet
	{
		if (theta <= 0.5 * MATH_PI)
		{
			if (arpragueskymodel_circleBounds2D(x_v, y_v, y_c, atmo_edge, &n) == 0)
			{
				// Then we have a problem!
				// Return something, but this should never happen so long as the camera is inside the atmosphere
				// Which it should be in this work
				*a = 0;
				*d = 0;
				return;
			}
		}
		else
		{
			n = 0;
		}
	}
	else
	{
		if (arpragueskymodel_circleBounds2D(x_v, y_v, y_c, PSM_PLANET_RADIUS, &n) == 1) // Check for planet intersection
		{
			if (n <= distance) // We do intersect the planet so return a and d at the surface
			{
				double x_p = x_v * n;
				double y_p = (y_v * n) + PSM_PLANET_RADIUS + altitude;
				arpragueskymodel_scaleAD(x_p, y_p, a, d);
				return;
			}
		}
		if (arpragueskymodel_circleBounds2D(x_v, y_v, y_c, atmo_edge, &n) == 0)
		{
			// Then we have a problem!
			// Return something, but this should never happen so long as the camera is inside the atmosphere
			// Which it should be in this work
			*a = 0;
			*d = 0;
			return;
		}
	}
	double distance_corrected = n;
	// Use the smaller of the distances
	distance_corrected = distance < distance_corrected ? distance : distance_corrected;
	// Points in world space
	double x_p = x_v * distance_corrected;
	double y_p = (y_v * distance_corrected) + PSM_PLANET_RADIUS + altitude;
	arpragueskymodel_scaleAD(x_p, y_p, a, d);
}

float *arpragueskymodel_transmittanceCoefsIndex(
	const ArPragueSkyModelState  * state,
	const int turbidity,
	const int altitude,
	const int wavelength
)
{
	int transmittance_values_per_turbidity = state->trans_rank * 11 * state->trans_altitudes;
	return &state->transmission_dataset_V[(turbidity * transmittance_values_per_turbidity) + (((altitude * 11) + wavelength) * state->trans_rank)];
}

double clamp0_1(const double x)
{
	return (x < 0 ? 0 : (x > 1.0 ? 1.0 : x));
}

void arpragueskymodel_transmittanceInterpolateWaveLength(
	const ArPragueSkyModelState  * state,
	const int turbidity,
	const int altitude,
	const int wavelength_low,
	const int wavelength_inc,
	const double wavelength_w,
	double *coefficients
)
{
	float *wll = arpragueskymodel_transmittanceCoefsIndex(state, turbidity, altitude, wavelength_low);
	float *wlu = arpragueskymodel_transmittanceCoefsIndex(state, turbidity, altitude, wavelength_low + wavelength_inc);
	for (int i = 0; i < state->trans_rank; i++)
	{
		coefficients[i] = lerp(wll[i], wlu[i], wavelength_w);
	}
}

double arpragueskymodel_calc_transmittance_svd_altitude(
	const ArPragueSkyModelState *state,
	const int turbidity,
	const int altitude,
	const int wavelength_low,
	const int wavelength_inc,
	const double wavelength_factor,
	const int a_int,
	const int d_int,
	const int a_inc,
	const int d_inc,
	const double wa,
	const double wd)
{
	float t[4] = { 0.0, 0.0, 0.0, 0.0 };
	double interpolated_coefficients[transsvdrank];
	arpragueskymodel_transmittanceInterpolateWaveLength(state, turbidity, altitude, wavelength_low, wavelength_inc, wavelength_factor, interpolated_coefficients);
	int index = 0;
	// Calculate pow space values
	for (int al = a_int; al <= a_int + a_inc; al++)
	{
		for (int dl = d_int; dl <= d_int + d_inc; dl++)
		{
			for (int i = 0; i < state->trans_rank; i++)
			{
				t[index] = t[index] + (state->transmission_dataset_U[(altitude * state->trans_n_a * state->trans_n_d * state->trans_rank) + (((dl * state->trans_n_a) + al) * state->trans_rank) + i] * interpolated_coefficients[i]);
			}
			index++;
		}
	}
	if (d_inc == 1)
	{
		t[0] = lerp(t[0], t[1], wd);
		t[1] = lerp(t[2], t[3], wd);
	}
	if (a_inc == 1)
	{
		t[0] = lerp(t[0], t[1], wa);
	}
	return t[0];
}

double nonlinlerp(const double a, const double b, const double w, const double p)
{
	double c1 = pow(a, p);
	double c2 = pow(b, p);
	return ((pow(w, p) - c1) / (c2 - c1));
}

double arpragueskymodel_calc_transmittance_svd(
	const ArPragueSkyModelState *state,
	const double a,
	const double d,
	const int turbidity,
	const int wavelength_low,
	const int wavelength_inc,
	const double wavelength_factor,
	const int altitude_low,
	const int altitude_inc,
	const double altitude_factor)
{
	int a_int = (int)floor(a * (double)state->trans_n_a);
	int d_int = (int)floor(d * (double)state->trans_n_d);
	int a_inc = 0;
	int d_inc = 0;
	double wa = (a * (double)state->trans_n_a) - (double)a_int;
	double wd = (d * (double)state->trans_n_d) - (double)d_int;
	if (a_int < (state->trans_n_a - 1))
	{
		a_inc = 1;
		wa = nonlinlerp((double)a_int / (double)state->trans_n_a, (double)(a_int + a_inc) / (double)state->trans_n_a, a, 3.0);
	} else
	{
		a_int = state->trans_n_a - 1;
		wa = 0;
	}
	if (d_int < (state->trans_n_d - 1))
	{
		d_inc = 1;
		wd = nonlinlerp((double)d_int / (double)state->trans_n_d, (double)(d_int + d_inc) / (double)state->trans_n_d, d, 4.0);
	} else
	{
		d_int = state->trans_n_d - 1;
		wd = 0;
	}
	wa = wa < 0 ? 0 : wa;
	wa = wa > 1.0 ? 1.0 : wa;
	wd = wd < 0 ? 0 : wd;
	wd = wd > 1.0 ? 1.0 : wd;
	double trans[2];
	trans[0] = arpragueskymodel_calc_transmittance_svd_altitude(state, turbidity, altitude_low, wavelength_low, wavelength_inc, wavelength_factor, a_int, d_int, a_inc, d_inc, wa, wd);
	if (altitude_inc == 1)
	{
		trans[1] = arpragueskymodel_calc_transmittance_svd_altitude(state, turbidity, altitude_low + altitude_inc, wavelength_low, wavelength_inc, wavelength_factor, a_int, d_int, a_inc, d_inc, wa, wd);
		trans[0] = lerp(trans[0], trans[1], altitude_factor);
	}
	return trans[0];
}

void arpragueskymodel_findInArray(const float *arr, const int arrLength, const double value, int *index, int *inc, double *w)
{
	*inc = 0;
	if (value <= arr[0])
	{
		*index = 0;
		*w = 1.0;
		return;
	}
	if (value >= arr[arrLength - 1])
	{
		*index = arrLength - 1;
		*w = 0;
		return;
	}
	for (int i = 1; i < arrLength; i++)
	{
		if (value < arr[i])
		{
			*index = i - 1;
			*inc = 1;
			*w = (value - arr[i - 1]) / (arr[i] - arr[i - 1]); // Assume linear
			return;
		}
	}
}

double arpragueskymodel_tau(
	const ArPragueSkyModelState  * state,
	const double                   theta,
	const double                   altitude,
        const double                   turbidity,
	const double                   wavelength,
	const double                   distance
)
{
	/*ASSERT_DOUBLE_WITHIN_RANGE(theta, 0.0, MATH_PI);
	ASSERT_DOUBLE_WITHIN_RANGE(altitude, 0.0, 15000.0);
	ASSERT_NONNEGATIVE_DOUBLE(turbidity);
	ASSERT_POSITIVE_DOUBLE(wavelength);
	ASSERT_POSITIVE_DOUBLE(distance);*/

	const double wavelength_norm = (wavelength - state->channel_start) / state->channel_width;
	if (wavelength_norm >= state->channels || wavelength_norm < 0.)
		return 0.;
	const int wavelength_low = (int)wavelength_norm;
	const double wavelength_factor = 0.0;
	const int wavelength_inc = wavelength_low < 10 ? 1 : 0;
    /*ASSERT_INTEGER_WITHIN_RANGE(wavelength_low, 0, 10);
    ASSERT_DOUBLE_WITHIN_RANGE(wavelength_factor, 0.0, 1.0);
    ASSERT_INTEGER_WITHIN_RANGE(wavelength_inc, 0, 1);*/

	int altitude_low;
	double altitude_factor;
	int altitude_inc;
	arpragueskymodel_findInArray(state->transmission_altitudes, state->trans_altitudes, altitude, &altitude_low, &altitude_inc, &altitude_factor);
	/*ASSERT_INTEGER_WITHIN_RANGE(altitude_low, 0, 21);
    ASSERT_DOUBLE_WITHIN_RANGE(altitude_factor, 0.0, 1.0);
    ASSERT_INTEGER_WITHIN_RANGE(altitude_inc, 0, 1);*/

	int turb_low;
	double turb_w;
	int turb_inc;
	arpragueskymodel_findInArray(state->transmission_turbities, state->trans_turbidities, turbidity, &turb_low, &turb_inc, &turb_w);
	/*ASSERT_INTEGER_WITHIN_RANGE(turb_low, 0, 2);
    ASSERT_DOUBLE_WITHIN_RANGE(turb_w, 0.0, 1.0);
    ASSERT_INTEGER_WITHIN_RANGE(turb_inc, 0, 1);*/

	// Calculate normalized and non-linearly scaled position in the atmosphere
	double a;
	double d;
	arpragueskymodel_toAD(theta, distance, altitude, &a, &d);
	/*ASSERT_NONNEGATIVE_DOUBLE(a);
	ASSERT_NONNEGATIVE_DOUBLE(d);*/

    // Evaluate basis at low turbidity
	double trans_low = arpragueskymodel_calc_transmittance_svd(state, a, d, turb_low, wavelength_low, wavelength_inc, wavelength_factor, altitude_low, altitude_inc, altitude_factor);

    // Evaluate basis at high turbidity
	double trans_high = arpragueskymodel_calc_transmittance_svd(state, a, d, turb_low + turb_inc, wavelength_low, wavelength_inc, wavelength_factor, altitude_low, altitude_inc, altitude_factor);

	// Return interpolated transmittance values
	double trans = lerp(trans_low, trans_high, turb_w);
    //ASSERT_VALID_DOUBLE(trans);

	trans = clamp0_1(trans);
	trans = trans * trans;
    //ASSERT_UNIT_RANGE_DOUBLE(trans);

	return trans;
}


///////////////////////////////////////////////
// Polarisation mono version
///////////////////////////////////////////////

double reconstruct_pol(
  const ArPragueSkyModelState * state,
  const double                  gamma,
  const double                  alpha,
  const int                     gamma_segment,
  const int                     alpha_segment,
  const double *                control_params
)
{
  double res = 0;
  for (int t = 0; t < state->tensor_components_pol; ++t) {
	const double sun_val_t = eval_pp(gamma, gamma_segment, state->sun_breaks_pol, control_params + state->sun_offset_pol + t * state->sun_stride_pol);
	const double zenith_val_t = eval_pp(alpha, alpha_segment, state->zenith_breaks_pol, control_params + state->zenith_offset_pol + t * state->zenith_stride_pol);
	res += sun_val_t * zenith_val_t;
  }
  return res;
}

double interpolate_elevation_pol(
  const ArPragueSkyModelState  * state,
  double                  elevation,
  int                     altitude,
  int                     turbidity,
  int                     albedo,
  int                     wavelength,
  double                  gamma,
  double                  alpha,
  int                     gamma_segment,
  int                     alpha_segment
)
{
  const int elevation_low = (int)elevation;
  const double factor = elevation - (double)elevation_low;

  const double * control_params_low = control_params_single_config(
    state,
    state->polarisation_dataset,
    state->total_coefs_single_config_pol,
    elevation_low,
    altitude,
    turbidity,
    albedo,
    wavelength);

  double res_low = reconstruct_pol(
    state,
    gamma,
    alpha,
    gamma_segment,
    alpha_segment,
    control_params_low);

  if (factor < 1e-6 || elevation_low >= (state->elevations - 1))
  {
    return res_low;
  }

  const double * control_params_high = control_params_single_config(
    state,
    state->polarisation_dataset,
    state->total_coefs_single_config_pol,
    elevation_low+1,
    altitude,
    turbidity,
    albedo,
    wavelength);

  double res_high = reconstruct_pol(
    state,
    gamma,
    alpha,
    gamma_segment,
    alpha_segment,
    control_params_high);

  return lerp(res_low, res_high, factor);
}

double interpolate_altitude_pol(
  const ArPragueSkyModelState  * state,
  double                  elevation,
  double                  altitude,
  int                     turbidity,
  int                     albedo,
  int                     wavelength,
  double                  gamma,
  double                  alpha,
  int                     gamma_segment,
  int                     alpha_segment
)
{
  const int altitude_low = (int)altitude;
  const double factor = altitude - (double)altitude_low;

  double res_low = interpolate_elevation_pol(
    state,
    elevation,
    altitude_low,
    turbidity,
    albedo,
    wavelength,
    gamma,
    alpha,
    gamma_segment,
    alpha_segment);

  if (factor < 1e-6 || altitude_low >= (state->altitudes - 1))
  {
    return res_low;
  }

  double res_high = interpolate_elevation_pol(
    state,
    elevation,
    altitude_low + 1,
    turbidity,
    albedo,
    wavelength,
    gamma,
    alpha,
    gamma_segment,
    alpha_segment);

  return lerp(res_low, res_high, factor);
}

double interpolate_turbidity_pol(
  const ArPragueSkyModelState  * state,
  double                  elevation,
  double                  altitude,
  double                  turbidity,
  int                     albedo,
  int                     wavelength,
  double                  gamma,
  double                  alpha,
  int                     gamma_segment,
  int                     alpha_segment
)
{
  // Ignore turbidity

  return interpolate_altitude_pol(
    state,
    elevation,
    altitude,
    (int)turbidity,
    albedo,
    wavelength,
    gamma,
    alpha,
    gamma_segment,
    alpha_segment);
}

double interpolate_albedo_pol(
  const ArPragueSkyModelState  * state,
  double                  elevation,
  double                  altitude,
  double                  turbidity,
  double                  albedo,
  int                     wavelength,
  double                  gamma,
  double                  alpha,
  int                     gamma_segment,
  int                     alpha_segment
)
{
  const int albedo_low = (int)albedo;
  const double factor = albedo - (double)albedo_low;

  double res_low = interpolate_turbidity_pol(
    state,
    elevation,
    altitude,
    turbidity,
    albedo_low,
    wavelength,
    gamma,
    alpha,
    gamma_segment,
    alpha_segment);

  if (factor < 1e-6 || albedo_low >= (state->albedos - 1))
  {
    return res_low;
  }

  double res_high = interpolate_turbidity_pol(
    state,
    elevation,
    altitude,
    turbidity,
    albedo_low + 1,
    wavelength,
    gamma,
    alpha,
    gamma_segment,
    alpha_segment);

  return lerp(res_low, res_high, factor);
}

double interpolate_wavelength_pol(
  const ArPragueSkyModelState  * state,
  double                  elevation,
  double                  altitude,
  double                  turbidity,
  double                  albedo,
  double                  wavelength,
  double                  gamma,
  double                  alpha,
  int                     gamma_segment,
  int                     alpha_segment
)
{
  // Don't interpolate, use the bin it belongs to

  return interpolate_albedo_pol(
    state,
    elevation,
    altitude,
    turbidity,
    albedo,
    (int)wavelength,
    gamma,
    alpha,
    gamma_segment,
    alpha_segment);
}

double arpragueskymodel_polarisation(
        const ArPragueSkyModelState  * state,
        const double                   theta,
        const double                   gamma,
        const double                   elevation,
        const double                   altitude,
        const double                   turbidity,
        const double                   albedo,
        const double                   wavelength
        )
{
  // If no polarisation data available
  if (state->tensor_components_pol == 0)
  {
    return 0.0;
  }

  // Translate parameter values to indices

  const double turbidity_control = map_parameter(turbidity, state->turbidities, state->turbidity_vals);
  const double albedo_control    = map_parameter(albedo, state->albedos, state->albedo_vals);
  const double altitude_control  = map_parameter(altitude, state->altitudes, state->altitude_vals);
  const double elevation_control = map_parameter(elevation * MATH_RAD_TO_DEG, state->elevations, state->elevation_vals);

  const double channel_control = (wavelength - state->channel_start) / state->channel_width;
  if ( channel_control >= state->channels || channel_control < 0.) return 0.;

  // Get params corresponding to the indices, reconstruct result and interpolate

  const int gamma_segment = find_segment(gamma, state->sun_nbreaks_pol, state->sun_breaks_pol);
  const int theta_segment = find_segment(theta, state->zenith_nbreaks_pol, state->zenith_breaks_pol);

  return -interpolate_wavelength_pol(
    state,
    elevation_control,
    altitude_control,
    turbidity_control,
    albedo_control,
    channel_control,
    gamma,
    theta,
    gamma_segment,
    theta_segment);
}
