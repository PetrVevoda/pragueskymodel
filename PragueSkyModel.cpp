#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <limits>
#include "PragueSkyModel.h"

/////////////////////////////////////////////////////////////////////////////////////
// Constants
/////////////////////////////////////////////////////////////////////////////////////

constexpr double PI              = 3.141592653589793;
constexpr double PLANET_RADIUS   = 6378000.0;
constexpr double SAFETY_ALTITUDE = 50.0;
constexpr double SUN_RADIUS      = 0.004654793; // = 0.2667 degrees
constexpr double DIST_TO_EDGE = 1571524.413613; // Maximum distance to the edge of the atmosphere in the transmittance model
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



/////////////////////////////////////////////////////////////////////////////////////
// Conversion functions
/////////////////////////////////////////////////////////////////////////////////////

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
    const uint64 dbits = uint64(hi) << 32;
    double out;
    std::memcpy(&out, &dbits, sizeof(double));
    return out;
}



/////////////////////////////////////////////////////////////////////////////////////
// Vector3 operations
/////////////////////////////////////////////////////////////////////////////////////

float dot(const PragueSkyModel::Vector3& v1, const PragueSkyModel::Vector3& v2) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

double magnitude(const PragueSkyModel::Vector3& vector) {
    return std::sqrt(dot(vector, vector));
}

PragueSkyModel::Vector3 normalize(const PragueSkyModel::Vector3& vector) {
    return vector / magnitude(vector);
}



/////////////////////////////////////////////////////////////////////////////////////
// Helper functions
/////////////////////////////////////////////////////////////////////////////////////

double lerp(const double from, const double to, const double factor) {
    return (1.0 - factor) * from + factor * to;
}

double nonlinlerp(const double a, const double b, const double w, const double p) {
    const double c1 = pow(a, p);
    const double c2 = pow(b, p);
    return ((pow(w, p) - c1) / (c2 - c1));
}

double clamp01(const double x) {
    return (x < 0.0 ? 0.0 : (x > 1.0 ? 1.0 : x));
}



/////////////////////////////////////////////////////////////////////////////////////
// Data reading
/////////////////////////////////////////////////////////////////////////////////////

int unpackCoefsFromHalf(const std::vector<double>& breaks,
                        const std::vector<uint16>& values,
                        std::vector<double>&       coefs,
                        const int                  offset,
                        const double               scale) {
    for (int i = 0; i < breaks.size() - 1; ++i) {
        const double val1 = doubleFromHalf(values[i + 1]) / scale;
        const double val2 = doubleFromHalf(values[i]) / scale;
        const double diff = val1 - val2;

        coefs[offset + 2 * i]     = diff / (breaks[i + 1] - breaks[i]);
        coefs[offset + 2 * i + 1] = val2;
    }
    return 2 * breaks.size() - 2;
}

int unpackCoefsFromFloat(const std::vector<double>& breaks,
                         const std::vector<float>&  values,
                         std::vector<double>&       coefs,
                         const int                  offset) {
    for (int i = 0; i < breaks.size() - 1; ++i) {
        coefs[offset + 2 * i]     = (double(values[i + 1]) - double(values[i])) / (breaks[i + 1] - breaks[i]);
        coefs[offset + 2 * i + 1] = double(values[i]);
    }
    return 2 * breaks.size() - 2;
}

void PragueSkyModel::readRadiance(FILE* handle) {
    // Read metadata

    // Structure of the metadata part of the data file:
    // visibilityCount   (1 * int),  visibilities (visibilityCount * double),
    // albedoCount      (1 * int),  albedos    (albedoCount * double),
    // altitudeCount    (1 * int),  altitudes  (altitudeCount * double),
    // elevationCount   (1 * int),  elevations (elevationCount * double),
    // channels          (1 * int),  channelStart  (1 * double), channelWidth (1
    // * double), rankRad (1 * int), sunBreaksCountRad       (1 * int),
    // sunBreaksRad     (sunBreaksCountRad * double), zenithBreaksCountRad    (1 * int),
    // zenithBreaksRad  (zenithBreaksCountRad * double), emphBreaksCountRad      (1 * int),
    // emphBreaksRad    (emphBreaksCountRad * double)

    int valsRead;

    int visibilityCount = 0;
    valsRead            = fread(&visibilityCount, sizeof(int), 1, handle);
    if (valsRead != 1 || visibilityCount < 1)
        throw DatasetReadException("visibilityCountRad");

    visibilitiesRad.resize(visibilityCount);
    valsRead = fread(visibilitiesRad.data(), sizeof(double), visibilityCount, handle);
    if (valsRead != visibilityCount)
        throw DatasetReadException("visibilitesRad");

    int albedoCount = 0;
    valsRead = fread(&albedoCount, sizeof(int), 1, handle);
    if (valsRead != 1 || albedoCount < 1)
        throw DatasetReadException("albedoCountRad");

    albedosRad.resize(albedoCount);
    valsRead    = fread(albedosRad.data(), sizeof(double), albedoCount, handle);
    if (valsRead != albedoCount)
        throw DatasetReadException("albedosRad");

    int altitudeCount = 0;
    valsRead = fread(&altitudeCount, sizeof(int), 1, handle);
    if (valsRead != 1 || altitudeCount < 1)
        throw DatasetReadException("altitudeCountRad");

    altitudesRad.resize(altitudeCount);
    valsRead      = fread(altitudesRad.data(), sizeof(double), altitudeCount, handle);
    if (valsRead != altitudeCount)
        throw DatasetReadException("altitudesRad");

    int elevationCount = 0;
    valsRead = fread(&elevationCount, sizeof(int), 1, handle);
    if (valsRead != 1 || elevationCount < 1)
        throw DatasetReadException(" elevationCountRad");

    elevationsRad.resize(elevationCount);
    valsRead       = fread(elevationsRad.data(), sizeof(double), elevationCount, handle);
    if (valsRead != elevationCount)
        throw DatasetReadException("elevationsRad");

    valsRead = fread(&channels, sizeof(int), 1, handle);
    if (valsRead != 1 || channels < 1)
        throw DatasetReadException("channels");

    valsRead = fread(&channelStart, sizeof(double), 1, handle);
    if (valsRead != 1 || channelStart < 0)
        throw DatasetReadException("channelStart");

    valsRead = fread(&channelWidth, sizeof(double), 1, handle);
    if (valsRead != 1 || channelWidth <= 0)
        throw DatasetReadException("channelWidth");

    valsRead = fread(&rankRad, sizeof(int), 1, handle);
    if (valsRead != 1 || rankRad < 1)
        throw DatasetReadException("rankRad");

    int sunBreaksCount = 0;
    valsRead = fread(&sunBreaksCount, sizeof(int), 1, handle);
    if (valsRead != 1 || sunBreaksCount < 2)
        throw DatasetReadException("sunBreaksCountRad");

    sunBreaksRad.resize(sunBreaksCount);
    valsRead   = fread(sunBreaksRad.data(), sizeof(double), sunBreaksCount, handle);
    if (valsRead != sunBreaksCount)
        throw DatasetReadException("sunBreaksRad");

    int zenitBreaksCount = 0;
    valsRead = fread(&zenitBreaksCount, sizeof(int), 1, handle);
    if (valsRead != 1 || zenitBreaksCount < 2)
        throw DatasetReadException("zenitBreaksCountRad");

    zenithBreaksRad.resize(zenitBreaksCount);
    valsRead      = fread(zenithBreaksRad.data(), sizeof(double), zenitBreaksCount, handle);
    if (valsRead != zenitBreaksCount)
        throw DatasetReadException("zenithBreaksRad");

    int emphBreaksCount = 0;
    valsRead = fread(&emphBreaksCount, sizeof(int), 1, handle);
    if (valsRead != 1 || emphBreaksCount < 2)
        throw DatasetReadException("emphBreaksCountRad");

    emphBreaksRad.resize(emphBreaksCount);
    valsRead    = fread(emphBreaksRad.data(), sizeof(double), emphBreaksCount, handle);
    if (valsRead != emphBreaksCount)
        throw DatasetReadException("emphBreaksRad");

    // Calculate offsets and strides

    sunOffsetRad = 0;
    sunStrideRad = 2 * sunBreaksRad.size() - 2 + 2 * zenithBreaksRad.size() - 2;

    zenithOffsetRad = sunOffsetRad + 2 * sunBreaksRad.size() - 2;
    zenithStrideRad = sunStrideRad;

    emphOffsetRad = sunOffsetRad + rankRad * sunStrideRad;

    totalCoefsSingleConfigRad = emphOffsetRad + 2 * emphBreaksRad.size() - 2; // this is for one specific configuration
    totalConfigsRad = channels * elevationsRad.size() * altitudesRad.size() * albedosRad.size() * visibilitiesRad.size();
    totalCoefsAllConfigsRad = totalCoefsSingleConfigRad * totalConfigsRad;

    // Read data

    // Structure of the data part of the data file:
    // [[[[[[ sunCoefsRad (sunBreaksCountRad * half), zenithScale (1 * double),
    // zenithCoefsRad (zenithBreaksCountRad * half) ] * rankRad, emphCoefsRad
    // (emphBreaksCountRad * half) ]
    //   * channels ] * elevationCount ] * altitudeCount ] * albedoCount ] * visibilityCount

    int offset       = 0;
    datasetRad.resize(totalCoefsAllConfigsRad);

    std::vector<uint16> radianceTemp;
    radianceTemp.resize(std::max(sunBreaksRad.size(), std::max(zenithBreaksRad.size(), emphBreaksRad.size())));

    for (int con = 0; con < totalConfigsRad; ++con) {
        for (int tc = 0; tc < rankRad; ++tc) {
            valsRead = fread(radianceTemp.data(), sizeof(unsigned short), sunBreaksRad.size(), handle);
            if (valsRead != sunBreaksRad.size())
                throw DatasetReadException("sunCoefsRad");
            offset += unpackCoefsFromHalf(sunBreaksRad, radianceTemp, datasetRad, offset, 1.0);

            double zenithScale;
            valsRead = fread(&zenithScale, sizeof(double), 1, handle);
            if (valsRead != 1)
                throw DatasetReadException("zenithScaleRad");

            valsRead = fread(radianceTemp.data(), sizeof(unsigned short), zenithBreaksRad.size(), handle);
            if (valsRead != zenithBreaksRad.size())
                throw DatasetReadException("zenithCoefsRad");
            offset += unpackCoefsFromHalf(zenithBreaksRad, radianceTemp, datasetRad, offset, zenithScale);
        }

        valsRead = fread(radianceTemp.data(), sizeof(unsigned short), emphBreaksRad.size(), handle);
        if (valsRead != emphBreaksRad.size())
            throw DatasetReadException("emphCoefsRad");
        offset += unpackCoefsFromHalf(emphBreaksRad, radianceTemp, datasetRad, offset, 1.0);
    }
}

void PragueSkyModel::readTransmittance(FILE* handle) {
    // Read metadata

    int valsRead;

    valsRead = fread(&dDim, sizeof(int), 1, handle);
    if (valsRead != 1 || dDim < 1)
        throw DatasetReadException("dDim");

    valsRead = fread(&aDim, sizeof(int), 1, handle);
	if (valsRead != 1 || aDim < 1)
        throw DatasetReadException("aDim");

    int visibilitiesCount = 0;
    valsRead = fread(&visibilitiesCount, sizeof(int), 1, handle);
    if (valsRead != 1 || visibilitiesCount < 1)
        throw DatasetReadException("visibilitiesCountTrans");

    int altitudesCount = 0;
    valsRead = fread(&altitudesCount, sizeof(int), 1, handle);
    if (valsRead != 1 || altitudesCount < 1)
        throw DatasetReadException("altitudesCountTrans");

    valsRead = fread(&rankTrans, sizeof(int), 1, handle);
    if (valsRead != 1 || rankTrans < 1)
        throw DatasetReadException("rankTrans");

    altitudesTrans.resize(altitudesCount);
    valsRead               = fread(altitudesTrans.data(), sizeof(float), altitudesCount, handle);
    if (valsRead != altitudesCount)
        throw DatasetReadException("altitudesTrans");

    visibilitiesTrans.resize(visibilitiesCount);
    valsRead               = fread(visibilitiesTrans.data(), sizeof(float), visibilitiesCount, handle);
    if (valsRead != visibilitiesCount)
        throw DatasetReadException("visibilitiesTrans");

    const int totalCoefsU = dDim * aDim * rankTrans * altitudesTrans.size();
    const int totalCoefsV = visibilitiesTrans.size() * rankTrans * 11 * altitudesTrans.size();

    // Read data

    datasetTransU.resize(totalCoefsU);
    valsRead               = fread(datasetTransU.data(), sizeof(float), totalCoefsU, handle);
    if (valsRead != totalCoefsU)
        throw DatasetReadException("datasetTransU");

    datasetTransV.resize(totalCoefsV);
    valsRead               = fread(datasetTransV.data(), sizeof(float), totalCoefsV, handle);
    if (valsRead != totalCoefsV)
        throw DatasetReadException("datasetTransV");
}

void PragueSkyModel::readPolarisation(FILE* handle) {
    // Read metadata

    // Structure of the metadata part of the data file:
    // rankPol (1 * int), sunBreaksCountPol (1 * int), sunBreaksPol (sunBreaksCountPol * double), zenithBreaksCountPol (1
    // * int), zenithBreaksPol (zenithBreaksCountPol * double), empBreaksCountPol (1 * int), emphBreaksPol (empBreaksCountPol
    // * double)

    int valsRead;

    valsRead = fread(&rankPol, sizeof(int), 1, handle);
    if (valsRead != 1) {
        // Polarisation dataset not present
        rankPol = 0;
        return;
    }

    int sunBreaksCount = 0;
    valsRead = fread(&sunBreaksCount, sizeof(int), 1, handle);
    if (valsRead != 1 || sunBreaksCount < 1)
        throw DatasetReadException("sunBreaksCountPol");

    sunBreaksPol.resize(sunBreaksCount);
    valsRead       = fread(sunBreaksPol.data(), sizeof(double), sunBreaksCount, handle);
    if (valsRead != sunBreaksCount)
        throw DatasetReadException("sunBreaksPol");

    int zenithBreaksCount = 0;
    valsRead = fread(&zenithBreaksCount, sizeof(int), 1, handle);
    if (valsRead != 1 || zenithBreaksCount < 1)
        throw DatasetReadException("zenithBreaksCountPol");

    zenithBreaksPol.resize(zenithBreaksCount);
    valsRead          = fread(zenithBreaksPol.data(), sizeof(double), zenithBreaksCount, handle);
    if (valsRead != zenithBreaksCount)
        throw DatasetReadException("zenithBreaksPol");

    // Calculate offsets and strides

    sunOffsetPol = 0;
    sunStridePol = 2 * sunBreaksPol.size() - 2 + 2 * zenithBreaksPol.size() - 2;

    zenithOffsetPol = sunOffsetPol + 2 * sunBreaksPol.size() - 2;
    zenithStridePol = sunStridePol;

    totalCoefsSingleConfigPol =
        sunOffsetPol + rankPol * sunStridePol; // this is for one specific configuration
    totalCoefsAllConfigsPol = totalCoefsSingleConfigPol * totalConfigsRad;

    // Read data

    // Structure of the data part of the data file:
    // [[[[[[ sunCoefsPol (sunBreaksCountPol * float), zenithCoefsPol
    // (zenithBreaksCountPol * float) ] * rankPol]
    //   * channels ] * elevationCount ] * altitudeCount ] * albedoCount ] * visibilityCount

    int offset               = 0;
    datasetPol.resize(totalCoefsAllConfigsPol);
    
    std::vector<float> polarisationTemp;
    polarisationTemp.resize(std::max(sunBreaksPol.size(), zenithBreaksPol.size()));

    for (int con = 0; con < totalConfigsRad; ++con) {
        for (int tc = 0; tc < rankPol; ++tc) {
            valsRead = fread(polarisationTemp.data(), sizeof(float), sunBreaksPol.size(), handle);
            if (valsRead != sunBreaksPol.size())
                throw DatasetReadException("sunCoefsPol");
            offset += unpackCoefsFromFloat(sunBreaksPol,
                                           polarisationTemp,
                                           datasetPol,
                                           offset);

            valsRead = fread(polarisationTemp.data(), sizeof(float), zenithBreaksPol.size(), handle);
            if (valsRead != zenithBreaksPol.size())
                throw DatasetReadException("zenithCoefsPol");
            offset += unpackCoefsFromFloat(zenithBreaksPol,
                                           polarisationTemp,
                                           datasetPol,
                                           offset);
        }
    }
}



/////////////////////////////////////////////////////////////////////////////////////
// Constructor
/////////////////////////////////////////////////////////////////////////////////////

PragueSkyModel::PragueSkyModel(const std::string& filename) {
    if (FILE* handle = fopen(filename.c_str(), "rb")) {
        // Read data
        readRadiance(handle);
        readTransmittance(handle);
        readPolarisation(handle);
        fclose(handle);
    } else {
        throw DatasetNotFoundException(filename);
    }
}



/////////////////////////////////////////////////////////////////////////////////////
// Angles
/////////////////////////////////////////////////////////////////////////////////////

void computeAltitudeAndElevation(const PragueSkyModel::Vector3& viewpoint,
                                 const double                   groundLevelSolarElevationAtOrigin,
                                 const double                   groundLevelSolarAzimuthAtOrigin,
                                 double&                        solarElevationAtViewpoint,
                                 double&                        altitudeOfViewpoint,
                                 double&                        distanceToView,
                                 PragueSkyModel::Vector3&       directionToZenithN,
                                 PragueSkyModel::Vector3&       directionToSunN) {
    // Direction to zenith

    PragueSkyModel::Vector3 centerOfTheEarth  = PragueSkyModel::Vector3(0.0, 0.0, -PLANET_RADIUS);
    PragueSkyModel::Vector3 directionToZenith = viewpoint - centerOfTheEarth;
    directionToZenithN                       = normalize(directionToZenith);

    // Altitude of viewpoint

    distanceToView = magnitude(directionToZenith);
    assert(distanceToView > -0.0001);
    distanceToView = std::max(distanceToView, 0.0);

    altitudeOfViewpoint = distanceToView - PLANET_RADIUS;
    altitudeOfViewpoint = std::max(altitudeOfViewpoint, 0.0);

    // Direction to sun

    directionToSunN.x = cos(groundLevelSolarAzimuthAtOrigin) * cos(groundLevelSolarElevationAtOrigin);
    directionToSunN.y = sin(groundLevelSolarAzimuthAtOrigin) * cos(groundLevelSolarElevationAtOrigin);
    directionToSunN.z = sin(groundLevelSolarElevationAtOrigin);

    // Solar elevation at viewpoint (more precisely, solar elevation at the point
    // on the ground directly below viewpoint)

    const double dotZenithSun = dot(directionToZenithN, directionToSunN);

    solarElevationAtViewpoint = 0.5 * PI - acos(dotZenithSun);
}

PragueSkyModel::Parameters PragueSkyModel::computeParameters(const Vector3& viewpoint,
                                                             const Vector3& viewDirection,
                                                             const double   groundLevelSolarElevationAtOrigin,
                                                             const double   groundLevelSolarAzimuthAtOrigin,
                                                             const double   visibility,
                                                             const double   albedo) const {
    assert(viewpoint.z >= 0.0);
    assert(magnitude(viewDirection) > 0.0);
    assert(visibility >= 0.0);
    assert(albedo >= 0.0 && albedo <= 1.0);

    Parameters params;
    params.visibility = visibility;
    params.albedo = albedo;

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
                                params.elevation,
                                params.altitude,
                                distanceToView,
                                directionToZenithN,
                                directionToSunN);

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
    params.gamma = acos(dotProductSun);

    // Shadow angle - requires correction

    const double effectiveElevation = groundLevelSolarElevationAtOrigin;
    const double effectiveAzimuth   = groundLevelSolarAzimuthAtOrigin;
    const double shadowAngle       = effectiveElevation + PI * 0.5;

    Vector3 shadowDirectionN = Vector3(cos(shadowAngle) * cos(effectiveAzimuth),
                                       cos(shadowAngle) * sin(effectiveAzimuth),
                                       sin(shadowAngle));

    const double dotProductShadow = dot(correctViewN, shadowDirectionN);
    params.shadow = acos(dotProductShadow);

    // Zenith angle (theta) - corrected version stored in otherwise unused zero
    // angle

    double cosThetaCor = dot(correctViewN, directionToZenithN);
    params.zero = acos(cosThetaCor);

    // Zenith angle (theta) - uncorrected version goes outside

    double cosTheta = dot(viewDirectionN, directionToZenithN);
    params.theta = acos(cosTheta);

    return params;
}



/////////////////////////////////////////////////////////////////////////////////////
// Parameters
/////////////////////////////////////////////////////////////////////////////////////

void findInArray(const std::vector<float>& arr, const double value, int& index, int& inc, double& w) {
    inc = 0.0;
    if (value <= arr.front()) {
        index = 0.0;
        w     = 1.0;
    } else if (value >= arr.back()) {
        index = arr.size() - 1;
        w     = 0.0;
    } else {
        for (int i = 1; i < arr.size(); i++) {
            if (value < arr[i]) {
                index = i - 1;
                inc   = 1;
                w     = (value - arr[i - 1]) / (arr[i] - arr[i - 1]); // Assume linear
                assert(w >= 0.0 && w <= 1.0);
                return;
            }
        }
    }
}

double mapParameter(const double param,  const std::vector<double>& values) {
    double mapped = 0.0;

    if (param < values.front()) {
        mapped = 0.0;
    } else if (param > values.back()) {
        mapped = (double)values.size() - 1.0;
    } else {
        for (int v = 0; v < values.size(); ++v) {
            const double val = values[v];
            if (std::abs(val - param) < 1e-6) {
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

int findSegment(const double x, const std::vector<double>& breaks) {
    int segment = 0;
    for (segment = 0; segment < breaks.size(); ++segment) {
        if (breaks[segment + 1] >= x)
            break;
    }
    return segment;
}

double evalPP(const double x, const int segment, const std::vector<double>& breaks, const std::vector<double>::const_iterator coefs) {
    const double  x0 = x - breaks[segment];
    const std::vector<double>::const_iterator sc = coefs + 2 * segment; // segment coefs
    return *sc * x0 + *(sc + 1);
}

std::vector<double>::const_iterator PragueSkyModel::controlParams(const std::vector<double>& dataset,
                                                                  const int totalCoefsSingleConfig,
                                                                  const int elevation,
                                                                  const int altitude,
                                                                  const int visibility,
                                                                  const int albedo,
                                                                  const int wavelength) const {
    return dataset.cbegin() +
           (totalCoefsSingleConfig *
            (wavelength + channels * elevation + channels * elevationsRad.size() * altitude +
             channels * elevationsRad.size() * altitudesRad.size() * albedo +
             channels * elevationsRad.size() * altitudesRad.size() * albedosRad.size() * visibility));
}

double PragueSkyModel::reconstruct(const double                              gamma,
                                   const double                              alpha,
                                   const double                              zero,
                                   const int                                 gammaSegment,
                                   const int                                 alphaSegment,
                                   const int                                 zeroSegment,
                                   const std::vector<double>::const_iterator controlParams) const {
    double res = 0.0;
    for (int t = 0; t < rankRad; ++t) {
        const double sunVal =
            evalPP(gamma, gammaSegment, sunBreaksRad, controlParams + sunOffsetRad + t * sunStrideRad);
        const double zenithVal = evalPP(alpha,
                                        alphaSegment,
                                        zenithBreaksRad,
                                        controlParams + zenithOffsetRad + t * zenithStrideRad);
        res += sunVal * zenithVal;
    }
    const double emphVal = evalPP(zero, zeroSegment, emphBreaksRad, controlParams + emphOffsetRad);
    res *= emphVal;

    return std::max(res, 0.0);
}

double PragueSkyModel::reconstructPol(const double                              gamma,
                                      const double                              alpha,
                                      const int                                 gammaSegment,
                                      const int                                 alphaSegment,
                                      const std::vector<double>::const_iterator controlParams) const {
    double res = 0;
    for (int t = 0; t < rankPol; ++t) {
        const double sunVal =
            evalPP(gamma, gammaSegment, sunBreaksPol, controlParams + sunOffsetPol + t * sunStridePol);
        const double zenithVal = evalPP(alpha,
                                        alphaSegment,
                                        zenithBreaksPol,
                                        controlParams + zenithOffsetPol + t * zenithStridePol);
        res += sunVal * zenithVal;
    }
    return res;
}



/////////////////////////////////////////////////////////////////////////////////////
// Sky radiance
/////////////////////////////////////////////////////////////////////////////////////

double PragueSkyModel::interpolateElevation(double elevation,
                                            int    altitude,
                                            int    visibility,
                                            int    albedo,
                                            int    wavelength,
                                            double gamma,
                                            double alpha,
                                            double zero,
                                            int    gammaSegment,
                                            int    alphaSegment,
                                            int    zeroSegment) const {
    const int    elevationLow = (int)elevation;
    const double factor       = elevation - (double)elevationLow;

    const std::vector<double>::const_iterator controlParamsLow = controlParams(datasetRad,
                                                                               totalCoefsSingleConfigRad,
                                                                               elevationLow,
                                                                               altitude,
                                                                               visibility,
                                                                               albedo,
                                                                               wavelength);

    double resLow =
        reconstruct(gamma, alpha, zero, gammaSegment, alphaSegment, zeroSegment, controlParamsLow);

    if (factor < 1e-6 || elevationLow >= (elevationsRad.size() - 1)) {
        return resLow;
    }

    const std::vector<double>::const_iterator controlParamsHigh = controlParams(datasetRad,
                                                                                totalCoefsSingleConfigRad,
                                                                                elevationLow + 1,
                                                                                altitude,
                                                                                visibility,
                                                                                albedo,
                                                                                wavelength);

    double resHigh =
        reconstruct(gamma, alpha, zero, gammaSegment, alphaSegment, zeroSegment, controlParamsHigh);

    return lerp(resLow, resHigh, factor);
}

double PragueSkyModel::interpolateAltitude(double elevation,
                                           double altitude,
                                           int    visibility,
                                           int    albedo,
                                           int    wavelength,
                                           double gamma,
                                           double alpha,
                                           double zero,
                                           int    gammaSegment,
                                           int    alphaSegment,
                                           int    zeroSegment) const {
    const int    altitudeLow = (int)altitude;
    const double factor      = altitude - (double)altitudeLow;

    double resLow = interpolateElevation(elevation,
                                         altitudeLow,
                                         visibility,
                                         albedo,
                                         wavelength,
                                         gamma,
                                         alpha,
                                         zero,
                                         gammaSegment,
                                         alphaSegment,
                                         zeroSegment);

    if (factor < 1e-6 || altitudeLow >= (altitudesRad.size() - 1)) {
        return resLow;
    }

    double resHigh = interpolateElevation(elevation,
                                          altitudeLow + 1,
                                          visibility,
                                          albedo,
                                          wavelength,
                                          gamma,
                                          alpha,
                                          zero,
                                          gammaSegment,
                                          alphaSegment,
                                          zeroSegment);

    return lerp(resLow, resHigh, factor);
}

double PragueSkyModel::interpolateVisibility(double elevation,
                                             double altitude,
                                             double visibility,
                                             int    albedo,
                                             int    wavelength,
                                             double gamma,
                                             double alpha,
                                             double zero,
                                             int    gammaSegment,
                                             int    alphaSegment,
                                             int    zeroSegment) const {
    const int    visibilityLow = (int)visibility;
    const double factor        = visibility - (double)visibilityLow;

    double resLow = interpolateAltitude(elevation,
                                        altitude,
                                        visibilityLow,
                                        albedo,
                                        wavelength,
                                        gamma,
                                        alpha,
                                        zero,
                                        gammaSegment,
                                        alphaSegment,
                                        zeroSegment);

    if (factor < 1e-6 || visibilityLow >= (visibilitiesRad.size() - 1)) {
        return resLow;
    }

    double resHigh = interpolateAltitude(elevation,
                                         altitude,
                                         visibilityLow + 1,
                                         albedo,
                                         wavelength,
                                         gamma,
                                         alpha,
                                         zero,
                                         gammaSegment,
                                         alphaSegment,
                                         zeroSegment);

    return lerp(resLow, resHigh, factor);
}

double PragueSkyModel::interpolateAlbedo(double elevation,
                                         double altitude,
                                         double visibility,
                                         double albedo,
                                         int    wavelength,
                                         double gamma,
                                         double alpha,
                                         double zero,
                                         int    gammaSegment,
                                         int    alphaSegment,
                                         int    zeroSegment) const {
    const int    albedoLow = (int)albedo;
    const double factor    = albedo - (double)albedoLow;

    double resLow = interpolateVisibility(elevation,
                                          altitude,
                                          visibility,
                                          albedoLow,
                                          wavelength,
                                          gamma,
                                          alpha,
                                          zero,
                                          gammaSegment,
                                          alphaSegment,
                                          zeroSegment);

    if (factor < 1e-6 || albedoLow >= (albedosRad.size() - 1)) {
        return resLow;
    }

    double resHigh = interpolateVisibility(elevation,
                                           altitude,
                                           visibility,
                                           albedoLow + 1,
                                           wavelength,
                                           gamma,
                                           alpha,
                                           zero,
                                           gammaSegment,
                                           alphaSegment,
                                           zeroSegment);

    return lerp(resLow, resHigh, factor);
}

double PragueSkyModel::interpolateWavelength(double elevation,
                                             double altitude,
                                             double visibility,
                                             double albedo,
                                             double wavelength,
                                             double gamma,
                                             double alpha,
                                             double zero,
                                             int    gammaSegment,
                                             int    alphaSegment,
                                             int    zeroSegment) const {
    // Don't interpolate, use the bin it belongs to

    return interpolateAlbedo(elevation,
                             altitude,
                             visibility,
                             albedo,
                             (int)wavelength,
                             gamma,
                             alpha,
                             zero,
                             gammaSegment,
                             alphaSegment,
                             zeroSegment);
}

double PragueSkyModel::skyRadiance(const Parameters& params, const double wavelength) const {
    assert(wavelength > 0);

    // Translate parameter values to indices

    const double visibilityControl = mapParameter(params.visibility, visibilitiesRad);
    const double albedoControl     = mapParameter(params.albedo, albedosRad);
    const double altitudeControl   = mapParameter(params.altitude, altitudesRad);
    const double elevationControl  = mapParameter(radiansToDegrees(params.elevation), elevationsRad);

    const double channelControl = (wavelength - channelStart) / channelWidth;

    if (channelControl >= channels || channelControl < 0.0)
        return 0.0;

    // Get params corresponding to the indices, reconstruct result and interpolate

    const double alpha = params.elevation < 0.0 ? params.shadow : params.zero;

    const int gammaSegment = findSegment(params.gamma, sunBreaksRad);
    const int alphaSegment = findSegment(alpha, zenithBreaksRad);
    const int zeroSegment  = findSegment(params.zero, emphBreaksRad);

    const double res = interpolateWavelength(elevationControl,
                                             altitudeControl,
                                             visibilityControl,
                                             albedoControl,
                                             channelControl,
                                             params.gamma,
                                             alpha,
                                             params.zero,
                                             gammaSegment,
                                             alphaSegment,
                                             zeroSegment);

    assert(res >= 0.0);
    return res;
}



/////////////////////////////////////////////////////////////////////////////////////
// Sun radiance
/////////////////////////////////////////////////////////////////////////////////////

double PragueSkyModel::sunRadiance(const Parameters& params, const double wavelength) const {
    assert(wavelength > 0);

    if (params.gamma > SUN_RADIUS) {
        return 0.0;
    }

    double idx         = (wavelength - SUN_RAD_START) / SUN_RAD_STEP;
    double sunRadiance = 0.0;

    if (idx >= 0.0) {
        int    lowIdx   = floor(idx);
        double idxFloat = idx - floor(idx);

        sunRadiance = SUN_RAD_TABLE[lowIdx] * (1.0 - idxFloat) + SUN_RAD_TABLE[lowIdx + 1] * idxFloat;
        assert(sunRadiance > 0.0);
    }

    double tau = PragueSkyModel::transmittance(params, wavelength, std::numeric_limits<double>::max());
    assert(tau >= 0.0 && tau <= 1.0);

    return sunRadiance * tau;
}



/////////////////////////////////////////////////////////////////////////////////////
// Transmittance
/////////////////////////////////////////////////////////////////////////////////////

bool circleBounds2D(double xV, double yV, double yC, double radius, double& d) {
    const double qa = (xV * xV) + (yV * yV);
    const double qb = 2.0 * yC * yV;
    const double qc = (yC * yC) - (radius * radius);
    double n  = (qb * qb) - (4.0 * qa * qc);
    if (n <= 0) {
        return false;
    }
    n  = sqrt(n);
    const float d1 = (-qb + n) / (2.0 * qa);
    const float d2 = (-qb - n) / (2.0 * qa);
    d = (d1 > 0 && d2 > 0) ? (d1 < d2 ? d1 : d2) : (d1 > d2 ? d1 : d2);
    if (d <= 0) {
        return false;
    }
    return true;
}

void scaleAD(double xP, double yP, double& a, double& d) {
    const double n = sqrt((xP * xP) + (yP * yP));
    a = n - PLANET_RADIUS;
    a = a > 0.0 ? a : 0.0;
    a = pow(a / 100000.0, 1.0 / 3.0);
    d = acos(yP / n) * PLANET_RADIUS;
    d = d / DIST_TO_EDGE;
    d = pow(d, 0.25);
    d = d > 1.0 ? 1.0 : d;
	assert(a >= 0.0);
	assert(d >= 0.0);
}

void toAD(double theta, double distance, double altitude, double& a, double& d) {
    // Ray circle intersection
    double xV       = sin(theta);
    double yV       = cos(theta);
    double yC       = PLANET_RADIUS + altitude;
    double atmoEdge = PLANET_RADIUS + 90000;
    double n;
    if (altitude < 0.001) // Handle altitudes close to 0 separately to avoid reporting
                          // intersections on the other side of the planet
    {
        if (theta <= 0.5 * PI) {
            if (!circleBounds2D(xV, yV, yC, atmoEdge, n)) {
                // Then we have a problem!
                // Return something, but this should never happen so long as the camera
                // is inside the atmosphere Which it should be in this work
                a = 0;
                d = 0;
                return;
            }
        } else {
            n = 0;
        }
    } else {
        if (circleBounds2D(xV, yV, yC, PLANET_RADIUS, n)) // Check for planet intersection
        {
            if (n <= distance) // We do intersect the planet so return a and d at the
                               // surface
            {
                const double xP = xV * n;
                const double yP = (yV * n) + PLANET_RADIUS + altitude;
                scaleAD(xP, yP, a, d);
                return;
            }
        }
        if (!circleBounds2D(xV, yV, yC, atmoEdge, n)) {
            // Then we have a problem!
            // Return something, but this should never happen so long as the camera is
            // inside the atmosphere Which it should be in this work
            a = 0;
            d = 0;
            return;
        }
    }
    double distanceCorrected = n;
    // Use the smaller of the distances
    distanceCorrected = distance < distanceCorrected ? distance : distanceCorrected;
    // Points in world space
    const double xP = xV * distanceCorrected;
    const double yP = (yV * distanceCorrected) + PLANET_RADIUS + altitude;
    scaleAD(xP, yP, a, d);
}

std::vector<float>::const_iterator PragueSkyModel::transmittanceCoefsIndex(const int visibility,
                                                                           const int altitude,
                                                                           const int wavelength) const {
    int transmittanceValuesPerVisibility = rankTrans * 11 * altitudesTrans.size();
    return datasetTransV.cbegin() + (visibility * transmittanceValuesPerVisibility) +
           (((altitude * 11) + wavelength) * rankTrans);
}

void PragueSkyModel::transmittanceInterpolateWaveLength(const int    visibility,
                                                        const int    altitude,
                                                        const int    wavelengthLow,
                                                        const int    wavelengthInc,
                                                        const double wavelengthW,
                                                        double*      coefficients) const {
    const std::vector<float>::const_iterator wll =
        transmittanceCoefsIndex(visibility, altitude, wavelengthLow);
    const std::vector<float>::const_iterator wlu =
        transmittanceCoefsIndex(visibility, altitude, wavelengthLow + wavelengthInc);
    for (int i = 0; i < rankTrans; i++) {
        coefficients[i] = lerp(wll[i], wlu[i], wavelengthW);
    }
}

double PragueSkyModel::calcTransmittanceSVDAltitude(const int    visibility,
                                                    const int    altitude,
                                                    const int    wavelengthLow,
                                                    const int    wavelengthInc,
                                                    const double wavelengthFactor,
                                                    const int    aInt,
                                                    const int    dInt,
                                                    const int    aInc,
                                                    const int    dInc,
                                                    const double wa,
                                                    const double wd) const {
    std::array<float, 4>         t = { 0.0, 0.0, 0.0, 0.0 };
    std::array<double, SVD_RANK> interpolatedCoefficients;
    transmittanceInterpolateWaveLength(visibility,
                                       altitude,
                                       wavelengthLow,
                                       wavelengthInc,
                                       wavelengthFactor,
                                       interpolatedCoefficients.data());
    int index = 0;
    // Calculate pow space values
    for (int al = aInt; al <= aInt + aInc; al++) {
        for (int dl = dInt; dl <= dInt + dInc; dl++) {
            for (int i = 0; i < rankTrans; i++) {
                t[index] =
                    t[index] + (datasetTransU[(altitude * aDim * dDim * rankTrans) +
                                                       (((dl * aDim) + al) * rankTrans) + i] *
                        interpolatedCoefficients[i]);
            }
            index++;
        }
    }
    if (dInc == 1) {
        t[0] = lerp(t[0], t[1], wd);
        t[1] = lerp(t[2], t[3], wd);
    }
    if (aInc == 1) {
        t[0] = lerp(t[0], t[1], wa);
    }
    return t[0];
}

double PragueSkyModel::calcTransmittanceSVD(const double a,
                                            const double d,
                                            const int    visibility,
                                            const int    wavelengthLow,
                                            const int    wavelengthInc,
                                            const double wavelengthFactor,
                                            const int    altitudeLow,
                                            const int    altitudeInc,
                                            const double altitudeFactor) const {
    int    aInt = int(floor(a * double(aDim)));
    int    dInt = int(floor(d * double(dDim)));
    int    aInc = 0;
    int    dInc = 0;
    double wa   = (a * double(aDim)) - double(aInt);
    double wd   = (d * double(dDim)) - double(dInt);
    if (aInt < (aDim - 1)) {
        aInc = 1;
        wa   = nonlinlerp(double(aInt) / double(aDim), double(aInt + aInc) / double(aDim), a, 3.0);
    } else {
        aInt = aDim - 1;
        wa   = 0;
    }
    if (dInt < (dDim - 1)) {
        dInc = 1;
        wd   = nonlinlerp(double(dInt) / double(dDim), double(dInt + dInc) / double(dDim), d, 4.0);
    } else {
        dInt = dDim - 1;
        wd   = 0;
    }
    wa = clamp01(wa);
    wd = clamp01(wd);

    double trans = calcTransmittanceSVDAltitude(visibility,
                                                altitudeLow,
                                                wavelengthLow,
                                                wavelengthInc,
                                                wavelengthFactor,
                                                aInt,
                                                dInt,
                                                aInc,
                                                dInc,
                                                wa,
                                                wd);
    if (altitudeInc == 1) {
        const double transHigh = calcTransmittanceSVDAltitude(visibility,
                                                              altitudeLow + altitudeInc,
                                                              wavelengthLow,
                                                              wavelengthInc,
                                                              wavelengthFactor,
                                                              aInt,
                                                              dInt,
                                                              aInc,
                                                              dInc,
                                                              wa,
                                                              wd);
        trans                  = lerp(trans, transHigh, altitudeFactor);
    }
    return trans;
}

double PragueSkyModel::transmittance(const Parameters& params,
                                     const double      wavelength,
                                     const double      distance) const {
    assert(wavelength > 0);
    assert(distance > 0);

    const double wavelengthNorm = (wavelength - channelStart) / channelWidth;
    if (wavelengthNorm >= channels || wavelengthNorm < 0.)
        return 0.;
    const int    wavelengthLow    = (int)wavelengthNorm;
    const double wavelengthFactor = 0.0;
    const int    wavelengthInc    = wavelengthLow < 10 ? 1 : 0;

    int    altitudeLow;
    double altitudeFactor;
    int    altitudeInc;
    findInArray(altitudesTrans, params.altitude, altitudeLow, altitudeInc, altitudeFactor);

    int    visibilityLow;
    double visibilityW;
    int    visibilityInc;
    findInArray(visibilitiesTrans, params.visibility, visibilityLow, visibilityInc, visibilityW);

    // Calculate normalized and non-linearly scaled position in the atmosphere
    double a;
    double d;
    toAD(params.theta, distance, params.altitude, a, d);

    // Evaluate basis at low visibility
    double transLow = calcTransmittanceSVD(a,
                                           d,
                                           visibilityLow,
                                           wavelengthLow,
                                           wavelengthInc,
                                           wavelengthFactor,
                                           altitudeLow,
                                           altitudeInc,
                                           altitudeFactor);

    // Evaluate basis at high visibility
    double transHigh = calcTransmittanceSVD(a,
                                            d,
                                            visibilityLow + visibilityInc,
                                            wavelengthLow,
                                            wavelengthInc,
                                            wavelengthFactor,
                                            altitudeLow,
                                            altitudeInc,
                                            altitudeFactor);

    // Return interpolated transmittance values
    double trans = lerp(transLow, transHigh, visibilityW);
    trans        = clamp01(trans);
    trans        = trans * trans;

    assert(trans >= 0.0 && trans <= 1.0);
    return trans;
}


/////////////////////////////////////////////////////////////////////////////////////
// Polarisation
/////////////////////////////////////////////////////////////////////////////////////

double PragueSkyModel::interpolateElevationPol(double elevation,
                                               int    altitude,
                                               int    visibility,
                                               int    albedo,
                                               int    wavelength,
                                               double gamma,
                                               double alpha,
                                               int    gammaSegment,
                                               int    alphaSegment) const {
    const int    elevationLow = (int)elevation;
    const double factor        = elevation - (double)elevationLow;

    const std::vector<double>::const_iterator controlParamsLow = controlParams(datasetPol,
                                                     totalCoefsSingleConfigPol,
                                                     elevationLow,
                                                     altitude,
                                                     visibility,
                                                     albedo,
                                                     wavelength);

    double resLow = reconstructPol(gamma, alpha, gammaSegment, alphaSegment, controlParamsLow);

    if (factor < 1e-6 || elevationLow >= (elevationsRad.size() - 1)) {
        return resLow;
    }

    const std::vector<double>::const_iterator controlParamsHigh = controlParams(datasetPol,
                                                      totalCoefsSingleConfigPol,
                                                      elevationLow + 1,
                                                      altitude,
                                                      visibility,
                                                      albedo,
                                                      wavelength);

    double resHigh = reconstructPol(gamma, alpha, gammaSegment, alphaSegment, controlParamsHigh);

    return lerp(resLow, resHigh, factor);
}

double PragueSkyModel::interpolateAltitudePol(double elevation,
                                              double altitude,
                                              int    visibility,
                                              int    albedo,
                                              int    wavelength,
                                              double gamma,
                                              double alpha,
                                              int    gammaSegment,
                                              int    alphaSegment) const {
    const int    altitudeLow = (int)altitude;
    const double factor       = altitude - (double)altitudeLow;

    double resLow = interpolateElevationPol(elevation,
                                             altitudeLow,
                                             visibility,
                                             albedo,
                                             wavelength,
                                             gamma,
                                             alpha,
                                             gammaSegment,
                                             alphaSegment);

    if (factor < 1e-6 || altitudeLow >= (altitudesRad.size() - 1)) {
        return resLow;
    }

    double resHigh = interpolateElevationPol(elevation,
                                              altitudeLow + 1,
                                              visibility,
                                              albedo,
                                              wavelength,
                                              gamma,
                                              alpha,
                                              gammaSegment,
                                              alphaSegment);

    return lerp(resLow, resHigh, factor);
}

double PragueSkyModel::interpolateVisibilityPol(double elevation,
                                                double altitude,
                                                double visibility,
                                                int    albedo,
                                                int    wavelength,
                                                double gamma,
                                                double alpha,
                                                int    gammaSegment,
                                                int    alphaSegment) const {
    // Ignore visibility

    return interpolateAltitudePol(elevation,
                                  altitude,
                                  (int)visibility,
                                  albedo,
                                  wavelength,
                                  gamma,
                                  alpha,
                                  gammaSegment,
                                  alphaSegment);
}

double PragueSkyModel::interpolateAlbedoPol(double elevation,
                                            double altitude,
                                            double visibility,
                                            double albedo,
                                            int    wavelength,
                                            double gamma,
                                            double alpha,
                                            int    gammaSegment,
                                            int    alphaSegment) const {
    const int    albedoLow = (int)albedo;
    const double factor     = albedo - (double)albedoLow;

    double resLow = interpolateVisibilityPol(elevation,
                                              altitude,
                                              visibility,
                                              albedoLow,
                                              wavelength,
                                              gamma,
                                              alpha,
                                              gammaSegment,
                                              alphaSegment);

    if (factor < 1e-6 || albedoLow >= (albedosRad.size() - 1)) {
        return resLow;
    }

    double resHigh = interpolateVisibilityPol(elevation,
                                               altitude,
                                               visibility,
                                               albedoLow + 1,
                                               wavelength,
                                               gamma,
                                               alpha,
                                               gammaSegment,
                                               alphaSegment);

    return lerp(resLow, resHigh, factor);
}

double PragueSkyModel::interpolateWavelengthPol(double elevation,
                                                double altitude,
                                                double visibility,
                                                double albedo,
                                                double wavelength,
                                                double gamma,
                                                double alpha,
                                                int    gammaSegment,
                                                int    alphaSegment) const {
    // Don't interpolate, use the bin it belongs to

    return interpolateAlbedoPol(elevation,
                                altitude,
                                visibility,
                                albedo,
                                (int)wavelength,
                                gamma,
                                alpha,
                                gammaSegment,
                                alphaSegment);
}

double PragueSkyModel::polarisation(const Parameters& params, const double wavelength) const {
    // If no polarisation data available
    if (rankPol == 0) {
        throw NoPolarisationException();
    }

    // Translate parameter values to indices

    const double visibilityControl = mapParameter(params.visibility, visibilitiesRad);
    const double albedoControl     = mapParameter(params.albedo, albedosRad);
    const double altitudeControl   = mapParameter(params.altitude, altitudesRad);
    const double elevationControl  = mapParameter(radiansToDegrees(params.elevation), elevationsRad);

    const double channelControl = (wavelength - channelStart) / channelWidth;
    if (channelControl >= channels || channelControl < 0.)
        return 0.;

    // Get params corresponding to the indices, reconstruct result and interpolate

    const int gammaSegment = findSegment(params.gamma, sunBreaksPol);
    const int thetaSegment = findSegment(params.theta, zenithBreaksPol);

    return -interpolateWavelengthPol(elevationControl,
                                     altitudeControl,
                                     visibilityControl,
                                     albedoControl,
                                     channelControl,
                                     params.gamma,
                                     params.theta,
                                     gammaSegment,
                                     thetaSegment);
}
