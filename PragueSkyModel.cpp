#include <algorithm>
#include <cmath>
#include <io.h>
#include <limits>
#include "PragueSkyModel.h"

#ifndef ALLOC
#    define ALLOC(_struct) ((_struct*)malloc(sizeof(_struct)))
#endif

#ifndef ALLOC_ARRAY
#    define ALLOC_ARRAY(_struct, _number) ((_struct*)malloc(sizeof(_struct) * (_number)))
#endif

///////////////////////////////////////////////
// Constants
///////////////////////////////////////////////

constexpr double PI              = 3.141592653589793;
constexpr double PLANET_RADIUS   = 6378000.0;
constexpr double SAFETY_ALTITUDE = 50.0;
constexpr double SUN_RADIUS      = 0.004654793; // = 0.2667 degrees
constexpr double SUN_RAD_START   = 310;
constexpr double SUN_RAD_STEP    = 1;
constexpr double SUN_RAD_TABLE[] = {
    9829.41, 10184.,  10262.6, 10375.7, 10276.,  10179.3, 10156.6, 10750.7, 11134.,  11463.6, 11860.4,
    12246.2, 12524.4, 12780.,  13187.4, 13632.4, 13985.9, 13658.3, 13377.4, 13358.3, 13239.,  13119.8,
    13096.2, 13184.,  13243.5, 13018.4, 12990.4, 13159.1, 13230.8, 13258.6, 13209.9, 13343.2, 13404.8,
    13305.4, 13496.3, 13979.1, 14153.8, 14188.4, 14122.7, 13825.4, 14033.3, 13914.1, 13837.4, 14117.2,
    13982.3, 13864.5, 14118.4, 14545.7, 15029.3, 15615.3, 15923.5, 16134.8, 16574.5, 16509.,  16336.5,
    16146.6, 15965.1, 15798.6, 15899.8, 16125.4, 15854.3, 15986.7, 15739.7, 15319.1, 15121.5, 15220.2,
    15041.2, 14917.7, 14487.8, 14011.,  14165.7, 14189.5, 14540.7, 14797.5, 14641.5, 14761.6, 15153.7,
    14791.8, 14907.6, 15667.4, 16313.5, 16917.,  17570.5, 18758.1, 20250.6, 21048.1, 21626.1, 22811.6,
    23577.2, 23982.6, 24062.1, 23917.9, 23914.1, 23923.2, 24052.6, 24228.6, 24360.8, 24629.6, 24774.8,
    24648.3, 24666.5, 24938.6, 24926.3, 24693.1, 24613.5, 24631.7, 24569.8, 24391.5, 24245.7, 24084.4,
    23713.7, 22985.4, 22766.6, 22818.9, 22834.3, 22737.9, 22791.6, 23086.3, 23377.7, 23461.,  23935.5,
    24661.7, 25086.9, 25520.1, 25824.3, 26198.,  26350.2, 26375.4, 26731.2, 27250.4, 27616.,  28145.3,
    28405.9, 28406.8, 28466.2, 28521.5, 28783.8, 29025.1, 29082.6, 29081.3, 29043.1, 28918.9, 28871.6,
    29049.,  29152.5, 29163.2, 29143.4, 28962.7, 28847.9, 28854.,  28808.7, 28624.1, 28544.2, 28461.4,
    28411.1, 28478.,  28469.8, 28513.3, 28586.5, 28628.6, 28751.5, 28948.9, 29051.,  29049.6, 29061.7,
    28945.7, 28672.8, 28241.5, 27903.2, 27737.,  27590.9, 27505.6, 27270.2, 27076.2, 26929.1, 27018.2,
    27206.8, 27677.2, 27939.9, 27923.9, 27899.2, 27725.4, 27608.4, 27599.4, 27614.6, 27432.4, 27460.4,
    27392.4, 27272.,  27299.1, 27266.8, 27386.5, 27595.9, 27586.9, 27504.8, 27480.6, 27329.8, 26968.4,
    26676.3, 26344.7, 26182.5, 26026.3, 25900.3, 25842.9, 25885.4, 25986.5, 26034.5, 26063.5, 26216.9,
    26511.4, 26672.7, 26828.5, 26901.8, 26861.5, 26865.4, 26774.2, 26855.8, 27087.1, 27181.3, 27183.1,
    27059.8, 26834.9, 26724.3, 26759.6, 26725.9, 26724.6, 26634.5, 26618.5, 26560.1, 26518.7, 26595.3,
    26703.2, 26712.7, 26733.9, 26744.3, 26764.4, 26753.2, 26692.7, 26682.7, 26588.1, 26478.,  26433.7,
    26380.7, 26372.9, 26343.3, 26274.7, 26162.3, 26160.5, 26210.,  26251.2, 26297.9, 26228.9, 26222.3,
    26269.7, 26295.6, 26317.9, 26357.5, 26376.1, 26342.4, 26303.5, 26276.7, 26349.2, 26390.,  26371.6,
    26346.7, 26327.6, 26274.2, 26247.3, 26228.7, 26152.1, 25910.3, 25833.2, 25746.5, 25654.3, 25562.,
    25458.8, 25438.,  25399.1, 25324.3, 25350.,  25514.,  25464.9, 25398.5, 25295.2, 25270.2, 25268.4,
    25240.6, 25184.9, 25149.6, 25123.9, 25080.3, 25027.9, 25012.3, 24977.9, 24852.6, 24756.4, 24663.5,
    24483.6, 24398.6, 24362.6, 24325.1, 24341.7, 24288.7, 24284.2, 24257.3, 24178.8, 24097.6, 24175.6,
    24175.7, 24139.7, 24088.1, 23983.2, 23902.7, 23822.4, 23796.2, 23796.9, 23814.5, 23765.5, 23703.,
    23642.,  23592.6, 23552.,  23514.6, 23473.5, 23431.,  23389.3, 23340.,  23275.1, 23187.3, 23069.5,
    22967.,  22925.3, 22908.9, 22882.5, 22825.,  22715.4, 22535.5, 22267.1, 22029.4, 21941.6, 21919.5,
    21878.8, 21825.6, 21766.,  21728.9, 21743.2, 21827.1, 21998.7, 22159.4, 22210.,  22187.2, 22127.2,
    22056.2, 22000.2, 21945.9, 21880.2, 21817.1, 21770.3, 21724.3, 21663.2, 21603.3, 21560.4, 21519.8,
    21466.2, 21401.6, 21327.7, 21254.2, 21190.7, 21133.6, 21079.3, 21024.,  20963.7, 20905.5, 20856.6,
    20816.6, 20785.2, 20746.7, 20685.3, 20617.8, 20561.1, 20500.4, 20421.2, 20333.4, 20247.,  20175.3,
    20131.4, 20103.2, 20078.5, 20046.8, 19997.2, 19952.9, 19937.2, 19930.8, 19914.4, 19880.8, 19823.,
    19753.8, 19685.9, 19615.3, 19537.5, 19456.8, 19377.6, 19309.4, 19261.9, 19228.,  19200.5, 19179.5,
    19164.8, 19153.1, 19140.6, 19129.2, 19120.6, 19104.5, 19070.6, 19023.9, 18969.3, 18911.4, 18855.,
    18798.6, 18740.8, 18672.7, 18585.2, 18501.,  18442.4, 18397.5, 18353.9, 18313.2, 18276.8, 18248.3,
    18231.2, 18224.,  18225.4, 18220.1, 18192.6, 18155.1, 18119.8, 18081.6, 18035.6, 17987.4, 17942.8,
    17901.7, 17864.2, 17831.1, 17802.9, 17771.5, 17728.6, 17669.7, 17590.1, 17509.5, 17447.4, 17396.,
    17347.4, 17300.3, 17253.2, 17206.1, 17159.,  17127.6, 17127.6, 17133.6, 17120.4, 17097.2, 17073.3,
    17043.7, 17003.4, 16966.3, 16946.3, 16930.9, 16907.7, 16882.7, 16862.,  16837.8, 16802.1, 16759.2,
    16713.6, 16661.8, 16600.8, 16542.6, 16499.4, 16458.7, 16408.,  16360.6, 16329.5, 16307.4, 16286.7,
    16264.9, 16239.6, 16207.8, 16166.8, 16118.2, 16064.,  16011.2, 15966.9, 15931.9, 15906.9, 15889.1,
    15875.5, 15861.2, 15841.3, 15813.1, 15774.2, 15728.8, 15681.4, 15630.,  15572.9, 15516.5, 15467.2,
    15423.,  15381.6, 15354.4, 15353.,  15357.3, 15347.3, 15320.2, 15273.1, 15222.,  15183.1, 15149.6,
    15114.6, 15076.8, 15034.6, 14992.9
};
constexpr int SVD_RANK = 12;

///////////////////////////////////////////////
// Conversion functions
///////////////////////////////////////////////

double radiansToDegrees(const double radians) {
    return radians * 180.0 / PI;
}

typedef unsigned short     uint16;
typedef unsigned int       uint32;
typedef unsigned long long uint64;

double doubleFromHalf(const uint16 half) {
    uint32 hi  = uint32(half & 0x8000) << 16;
    uint16 abs = half & 0x7FFF;
    if (abs) {
        hi |= 0x3F000000 << uint16(abs >= 0x7C00);
        for (; abs < 0x400; abs <<= 1, hi -= 0x100000)
            ;
        hi += uint32(abs) << 10;
    }
    uint64 dbits = uint64(hi) << 32;
    double out;
    std::memcpy(&out, &dbits, sizeof(double));
    return out;
}

///////////////////////////////////////////////
// Vector3 operations
///////////////////////////////////////////////

float dot(const PragueSkyModel::Vector3& v1, const PragueSkyModel::Vector3& v2) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

double magnitude(const PragueSkyModel::Vector3& vector) {
    return std::sqrt(dot(vector, vector));
}

PragueSkyModel::Vector3 normalize(const PragueSkyModel::Vector3& vector) {
    return vector / magnitude(vector);
}

///////////////////////////////////////////////
// Helper functions
///////////////////////////////////////////////

double lerp(const double from, const double to, const double factor) {
    return (1.0 - factor) * from + factor * to;
}

double nonlinlerp(const double a, const double b, const double w, const double p) {
    double c1 = pow(a, p);
    double c2 = pow(b, p);
    return ((pow(w, p) - c1) / (c2 - c1));
}

double clamp01(const double x) {
    return (x < 0 ? 0 : (x > 1.0 ? 1.0 : x));
}

///////////////////////////////////////////////
// Data reading
///////////////////////////////////////////////

int unpackCoefsFromHalf(const int             nbreaks,
                        const double*         breaks,
                        const unsigned short* values,
                        double*               coefs,
                        const int             offset,
                        const double          scale) {
    for (int i = 0; i < nbreaks - 1; ++i) {
        const double val1 = doubleFromHalf(values[i + 1]) / scale;
        const double val2 = doubleFromHalf(values[i]) / scale;
        const double diff = val1 - val2;

        coefs[offset + 2 * i]     = diff / (breaks[i + 1] - breaks[i]);
        coefs[offset + 2 * i + 1] = val2;
    }
    return 2 * nbreaks - 2;
}

int unpackCoefsFromFloat(const int     nbreaks,
                         const double* breaks,
                         const float*  values,
                         double*       coefs,
                         const int     offset) {
    for (int i = 0; i < nbreaks - 1; ++i) {
        coefs[offset + 2 * i]     = ((double)values[i + 1] - (double)values[i]) / (breaks[i + 1] - breaks[i]);
        coefs[offset + 2 * i + 1] = (double)values[i];
    }
    return 2 * nbreaks - 2;
}

void printErrorAndExit(const char* message) {
    fprintf(stderr, "%s", message);
    fprintf(stderr, "\n");
    fflush(stderr);
    exit(-1);
}

void PragueSkyModel::readRadiance(FILE* handle) {
    // Read metadata

    // Structure of the metadata part of the data file:
    // turbidities       (1 * int),  turbidity_vals (turbidities * double),
    // albedos           (1 * int),  albedo_vals    (albedos * double),
    // altitudes         (1 * int),  altitude_vals  (altitudes * double),
    // elevations        (1 * int),  elevation_vals (elevations * double),
    // channels          (1 * int),  channel_start  (1 * double), channel_width (1
    // * double), tensor_components (1 * int), sun_nbreaks       (1 * int),
    // sun_breaks     (sun_nbreaks * double), zenith_nbreaks    (1 * int),
    // zenith_breaks  (zenith_nbreaks * double), emph_nbreaks      (1 * int),
    // emph_breaks    (emph_nbreaks * double)

    int valsRead;

    valsRead = fread(&turbidities, sizeof(int), 1, handle);
    if (valsRead != 1 || turbidities < 1)
        printErrorAndExit("Error reading sky model data: turbidities");

    turbidity_vals = ALLOC_ARRAY(double, turbidities);
    valsRead       = fread(turbidity_vals, sizeof(double), turbidities, handle);
    if (valsRead != turbidities)
        printErrorAndExit("Error reading sky model data: turbidity_vals");

    valsRead = fread(&albedos, sizeof(int), 1, handle);
    if (valsRead != 1 || albedos < 1)
        printErrorAndExit("Error reading sky model data: albedos");

    albedo_vals = ALLOC_ARRAY(double, albedos);
    valsRead    = fread(albedo_vals, sizeof(double), albedos, handle);
    if (valsRead != albedos)
        printErrorAndExit("Error reading sky model data: albedo_vals");

    valsRead = fread(&altitudes, sizeof(int), 1, handle);
    if (valsRead != 1 || altitudes < 1)
        printErrorAndExit("Error reading sky model data: altitudes");

    altitude_vals = ALLOC_ARRAY(double, altitudes);
    valsRead      = fread(altitude_vals, sizeof(double), altitudes, handle);
    if (valsRead != altitudes)
        printErrorAndExit("Error reading sky model data: altitude_vals");

    valsRead = fread(&elevations, sizeof(int), 1, handle);
    if (valsRead != 1 || elevations < 1)
        printErrorAndExit("Error reading sky model data: elevations");

    elevation_vals = ALLOC_ARRAY(double, elevations);
    valsRead       = fread(elevation_vals, sizeof(double), elevations, handle);
    if (valsRead != elevations)
        printErrorAndExit("Error reading sky model data: elevation_vals");

    valsRead = fread(&channels, sizeof(int), 1, handle);
    if (valsRead != 1 || channels < 1)
        printErrorAndExit("Error reading sky model data: channels");

    valsRead = fread(&channel_start, sizeof(double), 1, handle);
    if (valsRead != 1 || channel_start < 0)
        printErrorAndExit("Error reading sky model data: channel_start");

    valsRead = fread(&channel_width, sizeof(double), 1, handle);
    if (valsRead != 1 || channel_width <= 0)
        printErrorAndExit("Error reading sky model data: channel_width");

    valsRead = fread(&tensor_components, sizeof(int), 1, handle);
    if (valsRead != 1 || tensor_components < 1)
        printErrorAndExit("Error reading sky model data: tensor_components");

    valsRead = fread(&sun_nbreaks, sizeof(int), 1, handle);
    if (valsRead != 1 || sun_nbreaks < 2)
        printErrorAndExit("Error reading sky model data: sun_nbreaks");

    sun_breaks = ALLOC_ARRAY(double, sun_nbreaks);
    valsRead   = fread(sun_breaks, sizeof(double), sun_nbreaks, handle);
    if (valsRead != sun_nbreaks)
        printErrorAndExit("Error reading sky model data: sun_breaks");

    valsRead = fread(&zenith_nbreaks, sizeof(int), 1, handle);
    if (valsRead != 1 || zenith_nbreaks < 2)
        printErrorAndExit("Error reading sky model data: zenith_nbreaks");

    zenith_breaks = ALLOC_ARRAY(double, zenith_nbreaks);
    valsRead      = fread(zenith_breaks, sizeof(double), zenith_nbreaks, handle);
    if (valsRead != zenith_nbreaks)
        printErrorAndExit("Error reading sky model data: zenith_breaks");

    valsRead = fread(&emph_nbreaks, sizeof(int), 1, handle);
    if (valsRead != 1 || emph_nbreaks < 2)
        printErrorAndExit("Error reading sky model data: emph_nbreaks");

    emph_breaks = ALLOC_ARRAY(double, emph_nbreaks);
    valsRead    = fread(emph_breaks, sizeof(double), emph_nbreaks, handle);
    if (valsRead != emph_nbreaks)
        printErrorAndExit("Error reading sky model data: emph_breaks");

    // Calculate offsets and strides

    sun_offset = 0;
    sun_stride = 2 * sun_nbreaks - 2 + 2 * zenith_nbreaks - 2;

    zenith_offset = sun_offset + 2 * sun_nbreaks - 2;
    zenith_stride = sun_stride;

    emph_offset = sun_offset + tensor_components * sun_stride;

    total_coefs_single_config = emph_offset + 2 * emph_nbreaks - 2; // this is for one specific configuration
    total_configs             = channels * elevations * altitudes * albedos * turbidities;
    total_coefs_all_configs   = total_coefs_single_config * total_configs;

    // Read data

    // Structure of the data part of the data file:
    // [[[[[[ sun_coefs (sun_nbreaks * half), zenith_scale (1 * double),
    // zenith_coefs (zenith_nbreaks * half) ] * tensor_components, emph_coefs
    // (emph_nbreaks * half) ]
    //   * channels ] * elevations ] * altitudes ] * albedos ] * turbidities

    int offset       = 0;
    radiance_dataset = ALLOC_ARRAY(double, total_coefs_all_configs);

    unsigned short* radiance_temp =
        ALLOC_ARRAY(unsigned short, std::max(sun_nbreaks, std::max(zenith_nbreaks, emph_nbreaks)));

    for (int con = 0; con < total_configs; ++con) {
        for (int tc = 0; tc < tensor_components; ++tc) {
            const double sun_scale = 1.0;
            valsRead               = fread(radiance_temp, sizeof(unsigned short), sun_nbreaks, handle);
            if (valsRead != sun_nbreaks)
                printErrorAndExit("Error reading sky model data: sun_coefs");
            offset += unpackCoefsFromHalf(sun_nbreaks,
                                          sun_breaks,
                                          radiance_temp,
                                          radiance_dataset,
                                          offset,
                                          sun_scale);

            double zenith_scale;
            valsRead = fread(&zenith_scale, sizeof(double), 1, handle);
            if (valsRead != 1)
                printErrorAndExit("Error reading sky model data: zenith_scale");

            valsRead = fread(radiance_temp, sizeof(unsigned short), zenith_nbreaks, handle);
            if (valsRead != zenith_nbreaks)
                printErrorAndExit("Error reading sky model data: zenith_coefs");
            offset += unpackCoefsFromHalf(zenith_nbreaks,
                                          zenith_breaks,
                                          radiance_temp,
                                          radiance_dataset,
                                          offset,
                                          zenith_scale);
        }

        const double emph_scale = 1.0;
        valsRead                = fread(radiance_temp, sizeof(unsigned short), emph_nbreaks, handle);
        if (valsRead != emph_nbreaks)
            printErrorAndExit("Error reading sky model data: emph_coefs");
        offset += unpackCoefsFromHalf(emph_nbreaks,
                                      emph_breaks,
                                      radiance_temp,
                                      radiance_dataset,
                                      offset,
                                      emph_scale);
    }

    free(radiance_temp);
}

void PragueSkyModel::readTransmittance(FILE* handle) {
    // Read metadata

    int valsRead;

    valsRead = fread(&trans_n_d, sizeof(int), 1, handle);
    if (valsRead != 1 || trans_n_d < 1)
        printErrorAndExit("Error reading sky model data: trans_n_d");

    valsRead = fread(&trans_n_a, sizeof(int), 1, handle);
    if (valsRead != 1 || trans_n_a < 1)
        printErrorAndExit("Error reading sky model data: trans_n_a");

    valsRead = fread(&trans_turbidities, sizeof(int), 1, handle);
    if (valsRead != 1 || trans_turbidities < 1)
        printErrorAndExit("Error reading sky model data: trans_turbidities");

    valsRead = fread(&trans_altitudes, sizeof(int), 1, handle);
    if (valsRead != 1 || trans_altitudes < 1)
        printErrorAndExit("Error reading sky model data: trans_altitudes");

    valsRead = fread(&trans_rank, sizeof(int), 1, handle);
    if (valsRead != 1 || trans_rank < 1)
        printErrorAndExit("Error reading sky model data: trans_rank");

    transmission_altitudes = ALLOC_ARRAY(float, trans_altitudes);
    valsRead               = fread(transmission_altitudes, sizeof(float), trans_altitudes, handle);
    if (valsRead != trans_altitudes)
        printErrorAndExit("Error reading sky model data: transmission_altitudes");

    transmission_turbities = ALLOC_ARRAY(float, trans_turbidities);
    valsRead               = fread(transmission_turbities, sizeof(float), trans_turbidities, handle);
    if (valsRead != trans_turbidities)
        printErrorAndExit("Error reading sky model data: transmission_turbities");

    const int total_coefs_U = trans_n_d * trans_n_a * trans_rank * trans_altitudes;
    const int total_coefs_V = trans_turbidities * trans_rank * 11 * trans_altitudes;

    // Read data

    transmission_dataset_U = ALLOC_ARRAY(float, total_coefs_U);
    valsRead               = fread(transmission_dataset_U, sizeof(float), total_coefs_U, handle);
    if (valsRead != total_coefs_U)
        printErrorAndExit("Error reading sky model data: transmission_dataset_U");

    transmission_dataset_V = ALLOC_ARRAY(float, total_coefs_V);
    valsRead               = fread(transmission_dataset_V, sizeof(float), total_coefs_V, handle);
    if (valsRead != total_coefs_V)
        printErrorAndExit("Error reading sky model data: transmission_dataset_V");
}

void PragueSkyModel::readPolarisation(FILE* handle) {
    // Read metadata

    // Structure of the metadata part of the data file:
    // tensor_components_pol (1 * int),
    // sun_nbreaks_pol       (1 * int),  sun_breaks_pol     (sun_nbreaks_pol *
    // double), zenith_nbreaks_pol    (1 * int),  zenith_breaks_pol
    // (zenith_nbreaks_pol * double), emph_nbreaks_pol      (1 * int),
    // emph_breaks_pol    (emph_nbreaks_pol * double)

    int valsRead;

    valsRead = fread(&tensor_components_pol, sizeof(int), 1, handle);
    if (valsRead != 1) {
        // Polarisation dataset not present
        tensor_components_pol = 0;
        printErrorAndExit("No polarisation dataset available!\n");
        return;
    }

    valsRead = fread(&sun_nbreaks_pol, sizeof(int), 1, handle);
    if (valsRead != 1 || sun_nbreaks_pol < 1)
        printErrorAndExit("Error reading sky model data: sun_nbreaks_pol");

    sun_breaks_pol = ALLOC_ARRAY(double, sun_nbreaks_pol);
    valsRead       = fread(sun_breaks_pol, sizeof(double), sun_nbreaks_pol, handle);
    if (valsRead != sun_nbreaks_pol)
        printErrorAndExit("Error reading sky model data: sun_breaks_pol");

    valsRead = fread(&zenith_nbreaks_pol, sizeof(int), 1, handle);
    if (valsRead != 1 || zenith_nbreaks_pol < 1)
        printErrorAndExit("Error reading sky model data: zenith_nbreaks_pol");

    zenith_breaks_pol = ALLOC_ARRAY(double, zenith_nbreaks_pol);
    valsRead          = fread(zenith_breaks_pol, sizeof(double), zenith_nbreaks_pol, handle);
    if (valsRead != zenith_nbreaks_pol)
        printErrorAndExit("Error reading sky model data: zenith_breaks_pol");

    // Calculate offsets and strides

    sun_offset_pol = 0;
    sun_stride_pol = 2 * sun_nbreaks_pol - 2 + 2 * zenith_nbreaks_pol - 2;

    zenith_offset_pol = sun_offset_pol + 2 * sun_nbreaks_pol - 2;
    zenith_stride_pol = sun_stride_pol;

    total_coefs_single_config_pol =
        sun_offset_pol + tensor_components_pol * sun_stride_pol; // this is for one specific configuration
    total_coefs_all_configs_pol = total_coefs_single_config_pol * total_configs;

    // Read data

    // Structure of the data part of the data file:
    // [[[[[[ sun_coefs_pol (sun_nbreaks_pol * float), zenith_coefs_pol
    // (zenith_nbreaks_pol * float) ] * tensor_components_pol]
    //   * channels ] * elevations ] * altitudes ] * albedos ] * turbidities

    int offset               = 0;
    polarisation_dataset     = ALLOC_ARRAY(double, total_coefs_all_configs_pol);
    float* polarisation_temp = ALLOC_ARRAY(float, std::max(sun_nbreaks_pol, zenith_nbreaks_pol));

    for (int con = 0; con < total_configs; ++con) {
        for (int tc = 0; tc < tensor_components_pol; ++tc) {
            valsRead = fread(polarisation_temp, sizeof(float), sun_nbreaks_pol, handle);
            if (valsRead != sun_nbreaks_pol)
                printErrorAndExit("Error reading sky model data: sun_coefs_pol");
            offset += unpackCoefsFromFloat(sun_nbreaks_pol,
                                           sun_breaks_pol,
                                           polarisation_temp,
                                           polarisation_dataset,
                                           offset);

            valsRead = fread(polarisation_temp, sizeof(float), zenith_nbreaks_pol, handle);
            if (valsRead != zenith_nbreaks_pol)
                printErrorAndExit("Error reading sky model data: zenith_coefs_pol");
            offset += unpackCoefsFromFloat(zenith_nbreaks_pol,
                                           zenith_breaks_pol,
                                           polarisation_temp,
                                           polarisation_dataset,
                                           offset);
        }
    }

    free(polarisation_temp);
}

///////////////////////////////////////////////
// Constructor & destructor
///////////////////////////////////////////////

PragueSkyModel::PragueSkyModel(const char* library_path) {
    char filename[1024];

    sprintf(filename, "%s", library_path);

    if (_access(filename, 0 | 4) != 0) {
        printErrorAndExit("Sky model dataset not found");
    }

    FILE* handle = fopen(filename, "rb");

    // Read data
    readRadiance(handle);
    readTransmittance(handle);
    readPolarisation(handle);

    fclose(handle);
}

PragueSkyModel::~PragueSkyModel() {
    free(turbidity_vals);
    free(albedo_vals);
    free(altitude_vals);
    free(elevation_vals);

    free(sun_breaks);
    free(zenith_breaks);
    free(emph_breaks);
    free(radiance_dataset);

    free(transmission_dataset_U);
    free(transmission_dataset_V);
    free(transmission_altitudes);
    free(transmission_turbities);

    if (tensor_components_pol > 0) {
        free(sun_breaks_pol);
        free(zenith_breaks_pol);
        free(polarisation_dataset);
    }
}

///////////////////////////////////////////////
// Angles
///////////////////////////////////////////////

void computeAltitudeAndElevation(const PragueSkyModel::Vector3& viewpoint,
                                 const double                   groundLevelSolarElevationAtOrigin,
                                 const double                   groundLevelSolarAzimuthAtOrigin,
                                 double*                        solarElevationAtViewpoint,
                                 double*                        altitudeOfViewpoint,
                                 double*                        distanceToView,
                                 PragueSkyModel::Vector3*       directionToZenithN,
                                 PragueSkyModel::Vector3*       directionToSunN) {
    // Direction to zenith

    PragueSkyModel::Vector3 centerOfTheEarth  = PragueSkyModel::Vector3(0.0, 0.0, -PLANET_RADIUS);
    PragueSkyModel::Vector3 directionToZenith = viewpoint - centerOfTheEarth;
    *directionToZenithN                       = normalize(directionToZenith);

    // Altitude of viewpoint

    *distanceToView = magnitude(directionToZenith);
    // ASSERT_DOUBLE_LARGER_THAN(*distanceToView, -0.0001);
    *distanceToView = std::max(*distanceToView, 0.0);

    *altitudeOfViewpoint = *distanceToView - PLANET_RADIUS;
    *altitudeOfViewpoint = std::max(*altitudeOfViewpoint, 0.0);

    // Direction to sun

    directionToSunN->x = cos(groundLevelSolarAzimuthAtOrigin) * cos(groundLevelSolarElevationAtOrigin);
    directionToSunN->y = sin(groundLevelSolarAzimuthAtOrigin) * cos(groundLevelSolarElevationAtOrigin);
    directionToSunN->z = sin(groundLevelSolarElevationAtOrigin);

    // Solar elevation at viewpoint (more precisely, solar elevation at the point
    // on the ground directly below viewpoint)

    const double dotZenithSun = dot(*directionToZenithN, *directionToSunN);

    *solarElevationAtViewpoint = 0.5 * PI - acos(dotZenithSun);
}

void PragueSkyModel::computeAngles(const Vector3& viewpoint,
                                   const Vector3& viewDirection,
                                   const double   groundLevelSolarElevationAtOrigin,
                                   const double   groundLevelSolarAzimuthAtOrigin,
                                   double*        solarElevationAtViewpoint,
                                   double*        altitudeOfViewpoint,
                                   double*        theta,
                                   double*        gamma,
                                   double*        shadow,
                                   double*        zero) const {
    // ASSERT_VALID_DOUBLE(groundLevelSolarElevationAtOrigin);
    // ASSERT_VALID_DOUBLE(groundLevelSolarAzimuthAtOrigin);

    // Shift viewpoint about safety altitude up

    const Vector3 centerOfTheEarth   = Vector3(0.0, 0.0, -PLANET_RADIUS);
    Vector3       toViewpoint        = viewpoint - centerOfTheEarth;
    Vector3       toViewpointN       = normalize(toViewpoint);
    const double  distanceToViewTmp  = magnitude(toViewpoint) + SAFETY_ALTITUDE;
    Vector3       toShiftedViewpoint = toViewpointN * distanceToViewTmp;
    Vector3       shiftedViewpoint   = centerOfTheEarth + toShiftedViewpoint;

    Vector3 viewDirectionN = normalize(viewDirection);

    double  distanceToView;
    Vector3 directionToZenithN;
    Vector3 directionToSunN;

    computeAltitudeAndElevation(shiftedViewpoint,
                                groundLevelSolarElevationAtOrigin,
                                groundLevelSolarAzimuthAtOrigin,
                                solarElevationAtViewpoint,
                                altitudeOfViewpoint,
                                &distanceToView,
                                &directionToZenithN,
                                &directionToSunN);

    // Altitude-corrected view direction

    Vector3 correctViewN;
    if (distanceToView > PLANET_RADIUS) {
        Vector3 lookAt = shiftedViewpoint + viewDirectionN;

        const double correction =
            sqrt(distanceToView * distanceToView - PLANET_RADIUS * PLANET_RADIUS) / distanceToView;

        Vector3 toNewOrigin = directionToZenithN * (distanceToView - correction);
        Vector3 newOrigin   = centerOfTheEarth + toNewOrigin;
        Vector3 correctView = lookAt - newOrigin;

        correctViewN = normalize(correctView);
    } else {
        correctViewN = viewDirectionN;
    }

    // Sun angle (gamma) - no correction

    double dotProductSun = dot(viewDirectionN, directionToSunN);

    *gamma = acos(dotProductSun);

    // Shadow angle - requires correction

    const double effectiveElevation = groundLevelSolarElevationAtOrigin;
    const double effectiveAzimuth   = groundLevelSolarAzimuthAtOrigin;
    const double shadow_angle       = effectiveElevation + PI * 0.5;

    Vector3 shadowDirectionN = Vector3(cos(shadow_angle) * cos(effectiveAzimuth),
                                       cos(shadow_angle) * sin(effectiveAzimuth),
                                       sin(shadow_angle));

    const double dotProductShadow = dot(correctViewN, shadowDirectionN);

    *shadow = acos(dotProductShadow);

    // Zenith angle (theta) - corrected version stored in otherwise unused zero
    // angle

    double cosThetaCor = dot(correctViewN, directionToZenithN);

    *zero = acos(cosThetaCor);

    // Zenith angle (theta) - uncorrected version goes outside

    double cosTheta = dot(viewDirectionN, directionToZenithN);

    *theta = acos(cosTheta);
}

///////////////////////////////////////////////
// Parameters
///////////////////////////////////////////////

void findInArray(const float* arr, const int arrLength, const double value, int* index, int* inc, double* w) {
    *inc = 0;
    if (value <= arr[0]) {
        *index = 0;
        *w     = 1.0;
        return;
    }
    if (value >= arr[arrLength - 1]) {
        *index = arrLength - 1;
        *w     = 0;
        return;
    }
    for (int i = 1; i < arrLength; i++) {
        if (value < arr[i]) {
            *index = i - 1;
            *inc   = 1;
            *w     = (value - arr[i - 1]) / (arr[i] - arr[i - 1]); // Assume linear
            return;
        }
    }
}

double mapParameter(const double param, const int value_count, const double* values) {
    double mapped = 0.0;

    if (param < values[0]) {
        mapped = 0.0;
    } else if (param > values[value_count - 1]) {
        mapped = (double)value_count - 1.0;
    } else {
        for (int v = 0; v < value_count; ++v) {
            const double val = values[v];
            if (fabs(val - param) < 1e-6) {
                mapped = v;
                break;
            } else if (param < val) {
                mapped = v - ((val - param) / (val - values[v - 1]));
                break;
            }
        }
    }

    return mapped;
}

int findSegment(const double x, const int nbreaks, const double* breaks) {
    int segment = 0;
    for (segment = 0; segment < nbreaks; ++segment) {
        if (breaks[segment + 1] >= x)
            break;
    }
    return segment;
}

double evalPP(const double x, const int segment, const double* breaks, const double* coefs) {
    const double  x0 = x - breaks[segment];
    const double* sc = coefs + 2 * segment; // segment coefs
    return sc[0] * x0 + sc[1];
}

const double* PragueSkyModel::controlParams(const double* dataset,
                                            const int     total_coefs_single_config,
                                            const int     elevation,
                                            const int     altitude,
                                            const int     turbidity,
                                            const int     albedo,
                                            const int     wavelength) const {
    return dataset + (total_coefs_single_config *
                      (wavelength + channels * elevation + channels * elevations * altitude +
                       channels * elevations * altitudes * albedo +
                       channels * elevations * altitudes * albedos * turbidity));
}

double PragueSkyModel::reconstruct(const double  gamma,
                                   const double  alpha,
                                   const double  zero,
                                   const int     gamma_segment,
                                   const int     alpha_segment,
                                   const int     zero_segment,
                                   const double* control_params) const {
    double res = 0.0;
    for (int t = 0; t < tensor_components; ++t) {
        const double sun_val_t =
            evalPP(gamma, gamma_segment, sun_breaks, control_params + sun_offset + t * sun_stride);
        const double zenith_val_t =
            evalPP(alpha, alpha_segment, zenith_breaks, control_params + zenith_offset + t * zenith_stride);
        res += sun_val_t * zenith_val_t;
    }
    const double emph_val_t = evalPP(zero, zero_segment, emph_breaks, control_params + emph_offset);
    res *= emph_val_t;

    return std::max(res, 0.0);
}

double PragueSkyModel::reconstructPol(const double  gamma,
                                      const double  alpha,
                                      const int     gamma_segment,
                                      const int     alpha_segment,
                                      const double* control_params) const {
    double res = 0;
    for (int t = 0; t < tensor_components_pol; ++t) {
        const double sun_val_t    = evalPP(gamma,
                                        gamma_segment,
                                        sun_breaks_pol,
                                        control_params + sun_offset_pol + t * sun_stride_pol);
        const double zenith_val_t = evalPP(alpha,
                                           alpha_segment,
                                           zenith_breaks_pol,
                                           control_params + zenith_offset_pol + t * zenith_stride_pol);
        res += sun_val_t * zenith_val_t;
    }
    return res;
}

///////////////////////////////////////////////
// Sky radiance
///////////////////////////////////////////////

double PragueSkyModel::interpolateElevation(double elevation,
                                            int    altitude,
                                            int    turbidity,
                                            int    albedo,
                                            int    wavelength,
                                            double gamma,
                                            double alpha,
                                            double zero,
                                            int    gamma_segment,
                                            int    alpha_segment,
                                            int    zero_segment) const {
    const int    elevation_low = (int)elevation;
    const double factor        = elevation - (double)elevation_low;

    const double* control_params_low = controlParams(radiance_dataset,
                                                     total_coefs_single_config,
                                                     elevation_low,
                                                     altitude,
                                                     turbidity,
                                                     albedo,
                                                     wavelength);

    double res_low =
        reconstruct(gamma, alpha, zero, gamma_segment, alpha_segment, zero_segment, control_params_low);

    if (factor < 1e-6 || elevation_low >= (elevations - 1)) {
        return res_low;
    }

    const double* control_params_high = controlParams(radiance_dataset,
                                                      total_coefs_single_config,
                                                      elevation_low + 1,
                                                      altitude,
                                                      turbidity,
                                                      albedo,
                                                      wavelength);

    double res_high =
        reconstruct(gamma, alpha, zero, gamma_segment, alpha_segment, zero_segment, control_params_high);

    return lerp(res_low, res_high, factor);
}

double PragueSkyModel::interpolateAltitude(double elevation,
                                           double altitude,
                                           int    turbidity,
                                           int    albedo,
                                           int    wavelength,
                                           double gamma,
                                           double alpha,
                                           double zero,
                                           int    gamma_segment,
                                           int    alpha_segment,
                                           int    zero_segment) const {
    const int    altitude_low = (int)altitude;
    const double factor       = altitude - (double)altitude_low;

    double res_low = interpolateElevation(elevation,
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

    if (factor < 1e-6 || altitude_low >= (altitudes - 1)) {
        return res_low;
    }

    double res_high = interpolateElevation(elevation,
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

double PragueSkyModel::interpolateVisibility(double elevation,
                                             double altitude,
                                             double turbidity,
                                             int    albedo,
                                             int    wavelength,
                                             double gamma,
                                             double alpha,
                                             double zero,
                                             int    gamma_segment,
                                             int    alpha_segment,
                                             int    zero_segment) const {
    const int    turbidity_low = (int)turbidity;
    const double factor        = turbidity - (double)turbidity_low;

    double res_low = interpolateAltitude(elevation,
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

    if (factor < 1e-6 || turbidity_low >= (turbidities - 1)) {
        return res_low;
    }

    double res_high = interpolateAltitude(elevation,
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

double PragueSkyModel::interpolateAlbedo(double elevation,
                                         double altitude,
                                         double turbidity,
                                         double albedo,
                                         int    wavelength,
                                         double gamma,
                                         double alpha,
                                         double zero,
                                         int    gamma_segment,
                                         int    alpha_segment,
                                         int    zero_segment) const {
    const int    albedo_low = (int)albedo;
    const double factor     = albedo - (double)albedo_low;

    double res_low = interpolateVisibility(elevation,
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

    if (factor < 1e-6 || albedo_low >= (albedos - 1)) {
        return res_low;
    }

    double res_high = interpolateVisibility(elevation,
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

double PragueSkyModel::interpolateWavelength(double elevation,
                                             double altitude,
                                             double turbidity,
                                             double albedo,
                                             double wavelength,
                                             double gamma,
                                             double alpha,
                                             double zero,
                                             int    gamma_segment,
                                             int    alpha_segment,
                                             int    zero_segment) const {
    // Don't interpolate, use the bin it belongs to

    return interpolateAlbedo(elevation,
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

double PragueSkyModel::skyRadiance(const double theta,
                                   const double gamma,
                                   const double shadow,
                                   const double zero,
                                   const double elevation,
                                   const double altitude,
                                   const double turbidity,
                                   const double albedo,
                                   const double wavelength) const {
    (void)theta;

    // Translate parameter values to indices

    const double turbidity_control = mapParameter(turbidity, turbidities, turbidity_vals);
    const double albedo_control    = mapParameter(albedo, albedos, albedo_vals);
    const double altitude_control  = mapParameter(altitude, altitudes, altitude_vals);
    const double elevation_control = mapParameter(radiansToDegrees(elevation), elevations, elevation_vals);

    const double channel_control = (wavelength - channel_start) / channel_width;

    if (channel_control >= channels || channel_control < 0.)
        return 0.;

    // Get params corresponding to the indices, reconstruct result and interpolate

    const double alpha = elevation < 0 ? shadow : zero;

    const int gamma_segment = findSegment(gamma, sun_nbreaks, sun_breaks);
    const int alpha_segment = findSegment(alpha, zenith_nbreaks, zenith_breaks);
    const int zero_segment  = findSegment(zero, emph_nbreaks, emph_breaks);

    const double res = interpolateWavelength(elevation_control,
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

    // ASSERT_NONNEGATIVE_DOUBLE(res);

    return res;
}

///////////////////////////////////////////////
// Sun radiance
///////////////////////////////////////////////

double PragueSkyModel::sunRadiance(const double theta,
                                   const double gamma,
                                   const double shadow,
                                   const double zero,
                                   const double elevation,
                                   const double altitude,
                                   const double turbidity,
                                   const double albedo,
                                   const double wavelength) const {
    (void)gamma;
    (void)shadow;
    (void)zero;
    (void)elevation;
    (void)albedo;
    double idx         = (wavelength - SUN_RAD_START) / SUN_RAD_STEP;
    double sunRadiance = 0.0;

    if (idx >= 0.0) {
        int    lowIdx   = floor(idx);
        double idxFloat = idx - floor(idx);

        sunRadiance = SUN_RAD_TABLE[lowIdx] * (1.0 - idxFloat) + SUN_RAD_TABLE[lowIdx + 1] * idxFloat;
        // ASSERT_POSITIVE_DOUBLE(sunRadiance);
    }

    double tau = PragueSkyModel::transmittance(theta,
                                               altitude,
                                               turbidity,
                                               wavelength,
                                               std::numeric_limits<double>::max());
    // ASSERT_UNIT_RANGE_DOUBLE(tau);

    return sunRadiance * tau;
}

///////////////////////////////////////////////
// Transmittance
///////////////////////////////////////////////

int circleBounds2D(double x_v, double y_v, double y_c, double radius, double* d) {
    double qa = (x_v * x_v) + (y_v * y_v);
    double qb = 2.0 * y_c * y_v;
    double qc = (y_c * y_c) - (radius * radius);
    double n  = (qb * qb) - (4.0 * qa * qc);
    if (n <= 0) {
        return 0;
    }
    float d1;
    float d2;
    n  = sqrt(n);
    d1 = (-qb + n) / (2.0 * qa);
    d2 = (-qb - n) / (2.0 * qa);
    *d = (d1 > 0 && d2 > 0) ? (d1 < d2 ? d1 : d2) : (d1 > d2 ? d1 : d2); // It fits in one line.
    if (*d <= 0) {
        return 0;
    }
    return 1;
}

void scaleAD(double x_p, double y_p, double* a, double* d) {
    double n;
    n  = sqrt((x_p * x_p) + (y_p * y_p));
    *a = n - PLANET_RADIUS;
    *a = *a > 0 ? *a : 0;
    *a = pow(*a / 100000.0, 1.0 / 3.0);
    *d = acos(y_p / n) * PLANET_RADIUS;
    *d = *d / 1571524.413613; // Maximum distance to the edge of the atmosphere in
                              // the transmittance model
    *d = pow(*d, 0.25);
    *d = *d > 1.0 ? 1.0 : *d;
}

void toAD(double theta, double distance, double altitude, double* a, double* d) {
    // Ray circle intersection
    double x_v       = sin(theta);
    double y_v       = cos(theta);
    double y_c       = PLANET_RADIUS + altitude;
    double atmo_edge = PLANET_RADIUS + 90000;
    double n;
    if (altitude < 0.001) // Handle altitudes close to 0 separately to avoid reporting
                          // intersections on the other side of the planet
    {
        if (theta <= 0.5 * PI) {
            if (circleBounds2D(x_v, y_v, y_c, atmo_edge, &n) == 0) {
                // Then we have a problem!
                // Return something, but this should never happen so long as the camera
                // is inside the atmosphere Which it should be in this work
                *a = 0;
                *d = 0;
                return;
            }
        } else {
            n = 0;
        }
    } else {
        if (circleBounds2D(x_v, y_v, y_c, PLANET_RADIUS, &n) == 1) // Check for planet intersection
        {
            if (n <= distance) // We do intersect the planet so return a and d at the
                               // surface
            {
                double x_p = x_v * n;
                double y_p = (y_v * n) + PLANET_RADIUS + altitude;
                scaleAD(x_p, y_p, a, d);
                return;
            }
        }
        if (circleBounds2D(x_v, y_v, y_c, atmo_edge, &n) == 0) {
            // Then we have a problem!
            // Return something, but this should never happen so long as the camera is
            // inside the atmosphere Which it should be in this work
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
    double y_p = (y_v * distance_corrected) + PLANET_RADIUS + altitude;
    scaleAD(x_p, y_p, a, d);
}

const float* PragueSkyModel::transmittanceCoefsIndex(const int turbidity,
                                                     const int altitude,
                                                     const int wavelength) const {
    int transmittance_values_per_turbidity = trans_rank * 11 * trans_altitudes;
    return &transmission_dataset_V[(turbidity * transmittance_values_per_turbidity) +
                                   (((altitude * 11) + wavelength) * trans_rank)];
}

void PragueSkyModel::transmittanceInterpolateWaveLength(const int    turbidity,
                                                        const int    altitude,
                                                        const int    wavelength_low,
                                                        const int    wavelength_inc,
                                                        const double wavelength_w,
                                                        double*      coefficients) const {
    const float* wll = transmittanceCoefsIndex(turbidity, altitude, wavelength_low);
    const float* wlu = transmittanceCoefsIndex(turbidity, altitude, wavelength_low + wavelength_inc);
    for (int i = 0; i < trans_rank; i++) {
        coefficients[i] = lerp(wll[i], wlu[i], wavelength_w);
    }
}

double PragueSkyModel::calcTransmittanceSVDAltitude(const int    turbidity,
                                                    const int    altitude,
                                                    const int    wavelength_low,
                                                    const int    wavelength_inc,
                                                    const double wavelength_factor,
                                                    const int    a_int,
                                                    const int    d_int,
                                                    const int    a_inc,
                                                    const int    d_inc,
                                                    const double wa,
                                                    const double wd) const {
    float  t[4] = { 0.0, 0.0, 0.0, 0.0 };
    double interpolated_coefficients[SVD_RANK];
    transmittanceInterpolateWaveLength(turbidity,
                                       altitude,
                                       wavelength_low,
                                       wavelength_inc,
                                       wavelength_factor,
                                       interpolated_coefficients);
    int index = 0;
    // Calculate pow space values
    for (int al = a_int; al <= a_int + a_inc; al++) {
        for (int dl = d_int; dl <= d_int + d_inc; dl++) {
            for (int i = 0; i < trans_rank; i++) {
                t[index] =
                    t[index] + (transmission_dataset_U[(altitude * trans_n_a * trans_n_d * trans_rank) +
                                                       (((dl * trans_n_a) + al) * trans_rank) + i] *
                                interpolated_coefficients[i]);
            }
            index++;
        }
    }
    if (d_inc == 1) {
        t[0] = lerp(t[0], t[1], wd);
        t[1] = lerp(t[2], t[3], wd);
    }
    if (a_inc == 1) {
        t[0] = lerp(t[0], t[1], wa);
    }
    return t[0];
}

double PragueSkyModel::calcTransmittanceSVD(const double a,
                                            const double d,
                                            const int    turbidity,
                                            const int    wavelength_low,
                                            const int    wavelength_inc,
                                            const double wavelength_factor,
                                            const int    altitude_low,
                                            const int    altitude_inc,
                                            const double altitude_factor) const {
    int    a_int = (int)floor(a * (double)trans_n_a);
    int    d_int = (int)floor(d * (double)trans_n_d);
    int    a_inc = 0;
    int    d_inc = 0;
    double wa    = (a * (double)trans_n_a) - (double)a_int;
    double wd    = (d * (double)trans_n_d) - (double)d_int;
    if (a_int < (trans_n_a - 1)) {
        a_inc = 1;
        wa    = nonlinlerp((double)a_int / (double)trans_n_a,
                        (double)(a_int + a_inc) / (double)trans_n_a,
                        a,
                        3.0);
    } else {
        a_int = trans_n_a - 1;
        wa    = 0;
    }
    if (d_int < (trans_n_d - 1)) {
        d_inc = 1;
        wd    = nonlinlerp((double)d_int / (double)trans_n_d,
                        (double)(d_int + d_inc) / (double)trans_n_d,
                        d,
                        4.0);
    } else {
        d_int = trans_n_d - 1;
        wd    = 0;
    }
    wa = wa < 0 ? 0 : wa;
    wa = wa > 1.0 ? 1.0 : wa;
    wd = wd < 0 ? 0 : wd;
    wd = wd > 1.0 ? 1.0 : wd;
    double trans[2];
    trans[0] = calcTransmittanceSVDAltitude(turbidity,
                                            altitude_low,
                                            wavelength_low,
                                            wavelength_inc,
                                            wavelength_factor,
                                            a_int,
                                            d_int,
                                            a_inc,
                                            d_inc,
                                            wa,
                                            wd);
    if (altitude_inc == 1) {
        trans[1] = calcTransmittanceSVDAltitude(turbidity,
                                                altitude_low + altitude_inc,
                                                wavelength_low,
                                                wavelength_inc,
                                                wavelength_factor,
                                                a_int,
                                                d_int,
                                                a_inc,
                                                d_inc,
                                                wa,
                                                wd);
        trans[0] = lerp(trans[0], trans[1], altitude_factor);
    }
    return trans[0];
}

double PragueSkyModel::transmittance(const double theta,
                                     const double altitude,
                                     const double turbidity,
                                     const double wavelength,
                                     const double distance) const {
    /*ASSERT_DOUBLE_WITHIN_RANGE(theta, 0.0, PI);
    ASSERT_DOUBLE_WITHIN_RANGE(altitude, 0.0, 15000.0);
    ASSERT_NONNEGATIVE_DOUBLE(turbidity);
    ASSERT_POSITIVE_DOUBLE(wavelength);
    ASSERT_POSITIVE_DOUBLE(distance);*/

    const double wavelength_norm = (wavelength - channel_start) / channel_width;
    if (wavelength_norm >= channels || wavelength_norm < 0.)
        return 0.;
    const int    wavelength_low    = (int)wavelength_norm;
    const double wavelength_factor = 0.0;
    const int    wavelength_inc    = wavelength_low < 10 ? 1 : 0;
    /*ASSERT_INTEGER_WITHIN_RANGE(wavelength_low, 0, 10);
    ASSERT_DOUBLE_WITHIN_RANGE(wavelength_factor, 0.0, 1.0);
    ASSERT_INTEGER_WITHIN_RANGE(wavelength_inc, 0, 1);*/

    int    altitude_low;
    double altitude_factor;
    int    altitude_inc;
    findInArray(transmission_altitudes,
                trans_altitudes,
                altitude,
                &altitude_low,
                &altitude_inc,
                &altitude_factor);
    /*ASSERT_INTEGER_WITHIN_RANGE(altitude_low, 0, 21);
  ASSERT_DOUBLE_WITHIN_RANGE(altitude_factor, 0.0, 1.0);
  ASSERT_INTEGER_WITHIN_RANGE(altitude_inc, 0, 1);*/

    int    turb_low;
    double turb_w;
    int    turb_inc;
    findInArray(transmission_turbities, trans_turbidities, turbidity, &turb_low, &turb_inc, &turb_w);
    /*ASSERT_INTEGER_WITHIN_RANGE(turb_low, 0, 2);
  ASSERT_DOUBLE_WITHIN_RANGE(turb_w, 0.0, 1.0);
  ASSERT_INTEGER_WITHIN_RANGE(turb_inc, 0, 1);*/

    // Calculate normalized and non-linearly scaled position in the atmosphere
    double a;
    double d;
    toAD(theta, distance, altitude, &a, &d);
    /*ASSERT_NONNEGATIVE_DOUBLE(a);
    ASSERT_NONNEGATIVE_DOUBLE(d);*/

    // Evaluate basis at low turbidity
    double trans_low = calcTransmittanceSVD(a,
                                            d,
                                            turb_low,
                                            wavelength_low,
                                            wavelength_inc,
                                            wavelength_factor,
                                            altitude_low,
                                            altitude_inc,
                                            altitude_factor);

    // Evaluate basis at high turbidity
    double trans_high = calcTransmittanceSVD(a,
                                             d,
                                             turb_low + turb_inc,
                                             wavelength_low,
                                             wavelength_inc,
                                             wavelength_factor,
                                             altitude_low,
                                             altitude_inc,
                                             altitude_factor);

    // Return interpolated transmittance values
    double trans = lerp(trans_low, trans_high, turb_w);
    // ASSERT_VALID_DOUBLE(trans);

    trans = clamp01(trans);
    trans = trans * trans;
    // ASSERT_UNIT_RANGE_DOUBLE(trans);

    return trans;
}

///////////////////////////////////////////////
// Polarisation mono version
///////////////////////////////////////////////

double PragueSkyModel::interpolateElevationPol(double elevation,
                                               int    altitude,
                                               int    turbidity,
                                               int    albedo,
                                               int    wavelength,
                                               double gamma,
                                               double alpha,
                                               int    gamma_segment,
                                               int    alpha_segment) const {
    const int    elevation_low = (int)elevation;
    const double factor        = elevation - (double)elevation_low;

    const double* control_params_low = controlParams(polarisation_dataset,
                                                     total_coefs_single_config_pol,
                                                     elevation_low,
                                                     altitude,
                                                     turbidity,
                                                     albedo,
                                                     wavelength);

    double res_low = reconstructPol(gamma, alpha, gamma_segment, alpha_segment, control_params_low);

    if (factor < 1e-6 || elevation_low >= (elevations - 1)) {
        return res_low;
    }

    const double* control_params_high = controlParams(polarisation_dataset,
                                                      total_coefs_single_config_pol,
                                                      elevation_low + 1,
                                                      altitude,
                                                      turbidity,
                                                      albedo,
                                                      wavelength);

    double res_high = reconstructPol(gamma, alpha, gamma_segment, alpha_segment, control_params_high);

    return lerp(res_low, res_high, factor);
}

double PragueSkyModel::interpolateAltitudePol(double elevation,
                                              double altitude,
                                              int    turbidity,
                                              int    albedo,
                                              int    wavelength,
                                              double gamma,
                                              double alpha,
                                              int    gamma_segment,
                                              int    alpha_segment) const {
    const int    altitude_low = (int)altitude;
    const double factor       = altitude - (double)altitude_low;

    double res_low = interpolateElevationPol(elevation,
                                             altitude_low,
                                             turbidity,
                                             albedo,
                                             wavelength,
                                             gamma,
                                             alpha,
                                             gamma_segment,
                                             alpha_segment);

    if (factor < 1e-6 || altitude_low >= (altitudes - 1)) {
        return res_low;
    }

    double res_high = interpolateElevationPol(elevation,
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

double PragueSkyModel::interpolateVisibilityPol(double elevation,
                                                double altitude,
                                                double turbidity,
                                                int    albedo,
                                                int    wavelength,
                                                double gamma,
                                                double alpha,
                                                int    gamma_segment,
                                                int    alpha_segment) const {
    // Ignore turbidity

    return interpolateAltitudePol(elevation,
                                  altitude,
                                  (int)turbidity,
                                  albedo,
                                  wavelength,
                                  gamma,
                                  alpha,
                                  gamma_segment,
                                  alpha_segment);
}

double PragueSkyModel::interpolateAlbedoPol(double elevation,
                                            double altitude,
                                            double turbidity,
                                            double albedo,
                                            int    wavelength,
                                            double gamma,
                                            double alpha,
                                            int    gamma_segment,
                                            int    alpha_segment) const {
    const int    albedo_low = (int)albedo;
    const double factor     = albedo - (double)albedo_low;

    double res_low = interpolateVisibilityPol(elevation,
                                              altitude,
                                              turbidity,
                                              albedo_low,
                                              wavelength,
                                              gamma,
                                              alpha,
                                              gamma_segment,
                                              alpha_segment);

    if (factor < 1e-6 || albedo_low >= (albedos - 1)) {
        return res_low;
    }

    double res_high = interpolateVisibilityPol(elevation,
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

double PragueSkyModel::interpolateWavelengthPol(double elevation,
                                                double altitude,
                                                double turbidity,
                                                double albedo,
                                                double wavelength,
                                                double gamma,
                                                double alpha,
                                                int    gamma_segment,
                                                int    alpha_segment) const {
    // Don't interpolate, use the bin it belongs to

    return interpolateAlbedoPol(elevation,
                                altitude,
                                turbidity,
                                albedo,
                                (int)wavelength,
                                gamma,
                                alpha,
                                gamma_segment,
                                alpha_segment);
}

double PragueSkyModel::polarisation(const double theta,
                                    const double gamma,
                                    const double elevation,
                                    const double altitude,
                                    const double turbidity,
                                    const double albedo,
                                    const double wavelength) const {
    // If no polarisation data available
    if (tensor_components_pol == 0) {
        return 0.0;
    }

    // Translate parameter values to indices

    const double turbidity_control = mapParameter(turbidity, turbidities, turbidity_vals);
    const double albedo_control    = mapParameter(albedo, albedos, albedo_vals);
    const double altitude_control  = mapParameter(altitude, altitudes, altitude_vals);
    const double elevation_control = mapParameter(radiansToDegrees(elevation), elevations, elevation_vals);

    const double channel_control = (wavelength - channel_start) / channel_width;
    if (channel_control >= channels || channel_control < 0.)
        return 0.;

    // Get params corresponding to the indices, reconstruct result and interpolate

    const int gamma_segment = findSegment(gamma, sun_nbreaks_pol, sun_breaks_pol);
    const int theta_segment = findSegment(theta, zenith_nbreaks_pol, zenith_breaks_pol);

    return -interpolateWavelengthPol(elevation_control,
                                     altitude_control,
                                     turbidity_control,
                                     albedo_control,
                                     channel_control,
                                     gamma,
                                     theta,
                                     gamma_segment,
                                     theta_segment);
}
