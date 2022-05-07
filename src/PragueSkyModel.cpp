// Copyright 2022 Charles University
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <array>
#include <cassert>
#include <limits>
#include <cstring>

#include "PragueSkyModel.h"

/////////////////////////////////////////////////////////////////////////////////////
// Constants
/////////////////////////////////////////////////////////////////////////////////////

constexpr double PI               = 3.141592653589793;
constexpr double RAD_TO_DEG       = 180.0 / PI;
constexpr double PLANET_RADIUS    = 6378000.0;
constexpr double SAFETY_ALTITUDE  = 50.0;
constexpr double SUN_RADIUS       = 0.004654793; // = 0.2667 degrees
constexpr double DIST_TO_EDGE     = 1571524.413613; // Maximum distance to the edge of the atmosphere in the transmittance model
constexpr double ATMOSPHERE_WIDTH = 100000.0;
constexpr double SUN_RAD_START    = 310;
constexpr double SUN_RAD_STEP     = 1;
constexpr double SUN_RAD_TABLE[]  = {
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
constexpr double SUN_RAD_END = SUN_RAD_START + SUN_RAD_STEP * std::size(SUN_RAD_TABLE);


/////////////////////////////////////////////////////////////////////////////////////
// Conversion functions
/////////////////////////////////////////////////////////////////////////////////////

double radiansToDegrees(const double radians) {
    return radians * RAD_TO_DEG;
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

double dot(const PragueSkyModel::Vector3& v1, const PragueSkyModel::Vector3& v2) {
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

double nonlerp(const double a, const double b, const double w, const double p) {
    const double c1 = pow(a, p);
    const double c2 = pow(b, p);
    return ((pow(w, p) - c1) / (c2 - c1));
}


/////////////////////////////////////////////////////////////////////////////////////
// Data reading
/////////////////////////////////////////////////////////////////////////////////////

void PragueSkyModel::readRadiance(FILE* handle, const double singleVisibility) {
    // Read metadata.

    // Structure of the metadata part of the data file:
    // visibilityCount      (1 * int), visibilities         (visibilityCount * double),
    // albedoCount          (1 * int), albedos                  (albedoCount * double),
    // altitudeCount        (1 * int), altitudes              (altitudeCount * double),
    // elevationCount       (1 * int), elevations            (elevationCount * double),
    // channels             (1 * int), channelStart                       (1 * double), channelWidth (1 * double), 
    // rankRad              (1 * int), 
    // sunBreaksCountRad    (1 * int), sunBreaksRad       (sunBreaksCountRad * double), 
    // zenithBreaksCountRad (1 * int), zenithBreaksRad (zenithBreaksCountRad * double), 
    // emphBreaksCountRad   (1 * int), emphBreaksRad     (emphBreaksCountRad * double)

    size_t valsRead;

    int visibilityCount = 0;
    valsRead            = fread(&visibilityCount, sizeof(int), 1, handle);
    if (valsRead != 1 || visibilityCount < 1)
        throw DatasetReadException("visibilityCountRad");

    std::vector<double> visibilitiesRadInFile;
    visibilitiesRadInFile.resize(visibilityCount);
    valsRead = fread(visibilitiesRadInFile.data(), sizeof(double), visibilityCount, handle);
    if (valsRead != visibilityCount)
        throw DatasetReadException("visibilitesRad");

    int skippedVisibilities = 0;
    if (singleVisibility <= 0.0 || visibilityCount <= 1) {
        // If the given single visibility value is not valid or there is just one visibility present (so there
        // is nothing to save), all visibilities present are loaded.
        visibilitiesRad.resize(visibilityCount);
        visibilitiesRad = visibilitiesRadInFile;
    } else {
        // Otherwise, only the two visibilities closest to the single visibility value are loaded and the rest
        // is skipped.
        if (singleVisibility <= visibilitiesRadInFile.front()) {
            visibilitiesRad.resize(1);
            visibilitiesRad[0] = visibilitiesRadInFile.front();
        } else if (singleVisibility >= visibilitiesRadInFile.back()) {
            visibilitiesRad.resize(1);
            visibilitiesRad[0] = visibilitiesRadInFile.back();
            skippedVisibilities = visibilitiesRadInFile.size() - 1;
        }
        else {
            int visIdx = 0;
            while (singleVisibility >= visibilitiesRadInFile[visIdx]) {
                ++visIdx;
            }
            visibilitiesRad.resize(2);
            visibilitiesRad[0] = visibilitiesRadInFile[visIdx - 1];
            visibilitiesRad[1] = visibilitiesRadInFile[visIdx];
            skippedVisibilities = visIdx - 1;
        }        
    }

    int albedoCount = 0;
    valsRead        = fread(&albedoCount, sizeof(int), 1, handle);
    if (valsRead != 1 || albedoCount < 1)
        throw DatasetReadException("albedoCountRad");

    albedosRad.resize(albedoCount);
    valsRead = fread(albedosRad.data(), sizeof(double), albedoCount, handle);
    if (valsRead != albedoCount)
        throw DatasetReadException("albedosRad");

    int altitudeCount = 0;
    valsRead          = fread(&altitudeCount, sizeof(int), 1, handle);
    if (valsRead != 1 || altitudeCount < 1)
        throw DatasetReadException("altitudeCountRad");

    altitudesRad.resize(altitudeCount);
    valsRead = fread(altitudesRad.data(), sizeof(double), altitudeCount, handle);
    if (valsRead != altitudeCount)
        throw DatasetReadException("altitudesRad");

    int elevationCount = 0;
    valsRead           = fread(&elevationCount, sizeof(int), 1, handle);
    if (valsRead != 1 || elevationCount < 1)
        throw DatasetReadException(" elevationCountRad");

    elevationsRad.resize(elevationCount);
    valsRead = fread(elevationsRad.data(), sizeof(double), elevationCount, handle);
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

    totalConfigs =
        channels * elevationsRad.size() * altitudesRad.size() * albedosRad.size() * visibilitiesRad.size();

    skippedConfigsBegin =
        channels * elevationsRad.size() * altitudesRad.size() * albedosRad.size() * skippedVisibilities;

    skippedConfigsEnd = channels * elevationsRad.size() * altitudesRad.size() * albedosRad.size() *
                        (visibilitiesRadInFile.size() - skippedVisibilities - visibilitiesRad.size());

    valsRead = fread(&metadataRad.rank, sizeof(int), 1, handle);
    if (valsRead != 1 || metadataRad.rank < 1)
        throw DatasetReadException("rankRad");

    int sunBreaksCount = 0;
    valsRead           = fread(&sunBreaksCount, sizeof(int), 1, handle);
    if (valsRead != 1 || sunBreaksCount < 2)
        throw DatasetReadException("sunBreaksCountRad");

    metadataRad.sunBreaks.resize(sunBreaksCount);
    valsRead = fread(metadataRad.sunBreaks.data(), sizeof(double), sunBreaksCount, handle);
    if (valsRead != sunBreaksCount)
        throw DatasetReadException("sunBreaksRad");

    int zenitBreaksCount = 0;
    valsRead             = fread(&zenitBreaksCount, sizeof(int), 1, handle);
    if (valsRead != 1 || zenitBreaksCount < 2)
        throw DatasetReadException("zenitBreaksCountRad");

    metadataRad.zenithBreaks.resize(zenitBreaksCount);
    valsRead = fread(metadataRad.zenithBreaks.data(), sizeof(double), zenitBreaksCount, handle);
    if (valsRead != zenitBreaksCount)
        throw DatasetReadException("zenithBreaksRad");

    int emphBreaksCount = 0;
    valsRead            = fread(&emphBreaksCount, sizeof(int), 1, handle);
    if (valsRead != 1 || emphBreaksCount < 2)
        throw DatasetReadException("emphBreaksCountRad");

    metadataRad.emphBreaks.resize(emphBreaksCount);
    valsRead = fread(metadataRad.emphBreaks.data(), sizeof(double), emphBreaksCount, handle);
    if (valsRead != emphBreaksCount)
        throw DatasetReadException("emphBreaksRad");

    // Calculate offsets and strides.

    metadataRad.sunOffset = 0;
    metadataRad.sunStride = metadataRad.sunBreaks.size() + metadataRad.zenithBreaks.size();

    metadataRad.zenithOffset = metadataRad.sunOffset + metadataRad.sunBreaks.size();
    metadataRad.zenithStride = metadataRad.sunStride;

    metadataRad.emphOffset = metadataRad.sunOffset + metadataRad.rank * metadataRad.sunStride;

    metadataRad.totalCoefsSingleConfig =
        metadataRad.emphOffset + metadataRad.emphBreaks.size(); // this is for one specific configuration
    metadataRad.totalCoefsAllConfigs = metadataRad.totalCoefsSingleConfig * totalConfigs;

    // Read data.

    // Structure of the data part of the data file:
    // [[[[[[ sunCoefsRad       (sunBreaksCountRad * half), zenithScale (1 * double), 
    //        zenithCoefsRad (zenithBreaksCountRad * half) ] * rankRad, 
    //        emphCoefsRad     (emphBreaksCountRad * half) ]
    //  * channels ] * elevationCount ] * altitudeCount ] * albedoCount ] * visibilityCount

    int offset = 0;
    dataRad.resize(metadataRad.totalCoefsAllConfigs);

    std::vector<uint16> radianceTemp;
    radianceTemp.resize(std::max(metadataRad.sunBreaks.size(),
                                 std::max(metadataRad.zenithBreaks.size(), metadataRad.emphBreaks.size())));

    size_t oneConfigByteCount =
        ((metadataRad.sunBreaks.size() + metadataRad.zenithBreaks.size()) * sizeof(uint16) + sizeof(double)) *
            metadataRad.rank +
        metadataRad.emphBreaks.size() * sizeof(uint16);

    // If a single visibility was requested, skip all configurations from the beginning till those needed for
    // the requested visibility.
    fseek(handle, oneConfigByteCount * skippedConfigsBegin, SEEK_CUR);

    // Read configurations needed for the requested visibility (or all if none requested).
    for (int con = 0; con < totalConfigs; ++con) {
        for (int r = 0; r < metadataRad.rank; ++r) {
            // Read sun params.
            valsRead = fread(radianceTemp.data(), sizeof(uint16), metadataRad.sunBreaks.size(), handle);
            if (valsRead != metadataRad.sunBreaks.size())
                throw DatasetReadException("sunCoefsRad");

            // Unpack sun params from half.
            for (int i = 0; i < metadataRad.sunBreaks.size(); ++i) {
                dataRad[offset++] = float(doubleFromHalf(radianceTemp[i]));
            }

            // Read scaling factor for zenith params.
            double zenithScale;
            valsRead = fread(&zenithScale, sizeof(double), 1, handle);
            if (valsRead != 1)
                throw DatasetReadException("zenithScaleRad");

            // Read zenith params.
            valsRead = fread(radianceTemp.data(), sizeof(uint16), metadataRad.zenithBreaks.size(), handle);
            if (valsRead != metadataRad.zenithBreaks.size())
                throw DatasetReadException("zenithCoefsRad");

            // Unpack zenith params from half (these need additional rescaling).
            for (int i = 0; i < metadataRad.zenithBreaks.size(); ++i) {
                dataRad[offset++] = float(doubleFromHalf(radianceTemp[i]) / zenithScale);
            }
        }

        // Read emphasize params.
        valsRead = fread(radianceTemp.data(), sizeof(uint16), metadataRad.emphBreaks.size(), handle);
        if (valsRead != metadataRad.emphBreaks.size())
            throw DatasetReadException("emphCoefsRad");

        // Unpack emphasize params from half.
        for (int i = 0; i < metadataRad.emphBreaks.size(); ++i) {
            dataRad[offset++] = float(doubleFromHalf(radianceTemp[i]));
        }
    }

    // Skip remaining configurations till the end.
    fseek(handle, oneConfigByteCount * skippedConfigsEnd, SEEK_CUR);
}

void PragueSkyModel::readTransmittance(FILE* handle) {
    // Read metadata

    size_t valsRead;

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

    std::vector<float> temp;
    temp.resize(std::max(altitudesCount, visibilitiesCount));

    altitudesTrans.resize(altitudesCount);
    valsRead = fread(temp.data(), sizeof(float), altitudesCount, handle);
    if (valsRead != altitudesCount)
        throw DatasetReadException("altitudesTrans");
    for (int i = 0; i < altitudesCount; i++) {
        altitudesTrans[i] = double(temp[i]);
    }

    visibilitiesTrans.resize(visibilitiesCount);
    valsRead = fread(temp.data(), sizeof(float), visibilitiesCount, handle);
    if (valsRead != visibilitiesCount)
        throw DatasetReadException("visibilitiesTrans");
    for (int i = 0; i < visibilitiesCount; i++) {
        visibilitiesTrans[i] = double(temp[i]);
    }

    const size_t totalCoefsU = dDim * aDim * rankTrans * altitudesTrans.size();
    const size_t totalCoefsV = visibilitiesTrans.size() * rankTrans * 11 * altitudesTrans.size();

    // Read data

    dataTransU.resize(totalCoefsU);
    valsRead = fread(dataTransU.data(), sizeof(float), totalCoefsU, handle);
    if (valsRead != totalCoefsU)
        throw DatasetReadException("datasetTransU");

    dataTransV.resize(totalCoefsV);
    valsRead = fread(dataTransV.data(), sizeof(float), totalCoefsV, handle);
    if (valsRead != totalCoefsV)
        throw DatasetReadException("datasetTransV");
}

void PragueSkyModel::readPolarisation(FILE* handle) {
    // Read metadata.

    // Structure of the metadata part of the data file:
    // rankPol              (1 * int), 
    // sunBreaksCountPol    (1 * int), sunBreaksPol    (sunBreaksCountPol * double),
    // zenithBreaksCountPol (1 * int), zenithBreaksPol (zenithBreaksCountPol * double), 
    // empBreaksCountPol    (1 * int), emphBreaksPol   (empBreaksCountPol * double)

    size_t valsRead;

    valsRead = fread(&metadataPol.rank, sizeof(int), 1, handle);
    if (valsRead != 1) {
        // Polarisation dataset not present
        metadataPol.rank = 0;
        return;
    }

    int sunBreaksCount = 0;
    valsRead           = fread(&sunBreaksCount, sizeof(int), 1, handle);
    if (valsRead != 1 || sunBreaksCount < 1)
        throw DatasetReadException("sunBreaksCountPol");

    metadataPol.sunBreaks.resize(sunBreaksCount);
    valsRead = fread(metadataPol.sunBreaks.data(), sizeof(double), sunBreaksCount, handle);
    if (valsRead != sunBreaksCount)
        throw DatasetReadException("sunBreaksPol");

    int zenithBreaksCount = 0;
    valsRead              = fread(&zenithBreaksCount, sizeof(int), 1, handle);
    if (valsRead != 1 || zenithBreaksCount < 1)
        throw DatasetReadException("zenithBreaksCountPol");

    metadataPol.zenithBreaks.resize(zenithBreaksCount);
    valsRead = fread(metadataPol.zenithBreaks.data(), sizeof(double), zenithBreaksCount, handle);
    if (valsRead != zenithBreaksCount)
        throw DatasetReadException("zenithBreaksPol");

    // Calculate offsets and strides.

    metadataPol.sunOffset = 0;
    metadataPol.sunStride = metadataPol.sunBreaks.size() + metadataPol.zenithBreaks.size();

    metadataPol.zenithOffset = metadataPol.sunOffset + metadataPol.sunBreaks.size();
    metadataPol.zenithStride = metadataPol.sunStride;

    metadataPol.totalCoefsSingleConfig =
        metadataPol.sunOffset +
        metadataPol.rank * metadataPol.sunStride; // this is for one specific configuration
    metadataPol.totalCoefsAllConfigs = metadataPol.totalCoefsSingleConfig * totalConfigs;

    // Read data.

    // Structure of the data part of the data file:
    // [[[[[[ sunCoefsPol       (sunBreaksCountPol * float), 
    //        zenithCoefsPol (zenithBreaksCountPol * float) ] * rankPol] 
    // * channels ] * elevationCount ] * altitudeCount ] * albedoCount ] * visibilityCount

    size_t offset = 0;
    dataPol.resize(metadataPol.totalCoefsAllConfigs);

    size_t oneConfigByteCount =
        (metadataPol.sunBreaks.size() + metadataPol.zenithBreaks.size()) * sizeof(float) * metadataPol.rank;

    // If a single visibility was requested, skip all configurations from the beginning till those needed for
    // the requested visibility.
    fseek(handle, oneConfigByteCount * skippedConfigsBegin, SEEK_CUR);

    for (int con = 0; con < totalConfigs; ++con) {
        for (int r = 0; r < metadataPol.rank; ++r) {
            // Read sun params.
            valsRead = fread(dataPol.data() + offset, sizeof(float), metadataPol.sunBreaks.size(), handle);
            if (valsRead != metadataPol.sunBreaks.size())
                throw DatasetReadException("sunCoefsPol");
            offset += metadataPol.sunBreaks.size();

            // Read zenith params.
            valsRead = fread(dataPol.data() + offset, sizeof(float), metadataPol.zenithBreaks.size(), handle);
            if (valsRead != metadataPol.zenithBreaks.size())
                throw DatasetReadException("zenithCoefsPol");
            offset += metadataPol.zenithBreaks.size();
        }
    }
}


/////////////////////////////////////////////////////////////////////////////////////
// Initialization
/////////////////////////////////////////////////////////////////////////////////////

void PragueSkyModel::initialize(const std::string& filename, const double singleVisibility) {
    if (FILE* handle = fopen(filename.c_str(), "rb")) {
        initialized = false;
        // Read data
        readRadiance(handle, singleVisibility);
        readTransmittance(handle);
        readPolarisation(handle);
        fclose(handle);
        initialized = true;
    } else {
        throw DatasetNotFoundException(filename);
    }
}

PragueSkyModel::AvailableData PragueSkyModel::getAvailableData() const {
    if (!initialized) {
        throw NotInitializedException();
    }

    AvailableData available;
    available.albedoMin     = albedosRad.front();
    available.albedoMax     = albedosRad.back();
    available.altitudeMin   = altitudesRad.front();
    available.altitudeMax   = altitudesRad.back();
    available.elevationMin  = elevationsRad.front();
    available.elevationMax  = elevationsRad.back();
    available.visibilityMin = visibilitiesRad.front();
    available.visibilityMax = visibilitiesRad.back();
    available.polarisation  = (metadataPol.rank != 0);

    return available;
}


/////////////////////////////////////////////////////////////////////////////////////
// Filling the Parameters structure
/////////////////////////////////////////////////////////////////////////////////////

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
    params.albedo     = albedo;

    // Shift viewpoint about safety altitude up

    const Vector3 centerOfTheEarth   = Vector3(0.0, 0.0, -PLANET_RADIUS);
    Vector3       toViewpoint        = viewpoint - centerOfTheEarth;
    Vector3       toViewpointN       = normalize(toViewpoint);
    const double  distanceToView     = magnitude(toViewpoint) + SAFETY_ALTITUDE;
    Vector3       toShiftedViewpoint = toViewpointN * distanceToView;
    Vector3       shiftedViewpoint   = centerOfTheEarth + toShiftedViewpoint;

    Vector3 viewDirectionN = normalize(viewDirection);

    // Compute altitude of viewpoint

    params.altitude = distanceToView - PLANET_RADIUS;
    params.altitude = std::max(params.altitude, 0.0);

    // Direction to sun

    Vector3 directionToSunN;
    directionToSunN.x = cos(groundLevelSolarAzimuthAtOrigin) * cos(groundLevelSolarElevationAtOrigin);
    directionToSunN.y = sin(groundLevelSolarAzimuthAtOrigin) * cos(groundLevelSolarElevationAtOrigin);
    directionToSunN.z = sin(groundLevelSolarElevationAtOrigin);

    // Solar elevation at viewpoint (more precisely, solar elevation at the point
    // on the ground directly below viewpoint)

    const double dotZenithSun = dot(toViewpointN, directionToSunN);
    params.elevation          = 0.5 * PI - acos(dotZenithSun);

    // Altitude-corrected view direction

    Vector3 correctViewN;
    if (distanceToView > PLANET_RADIUS) {
        Vector3 lookAt = shiftedViewpoint + viewDirectionN;

        const double correction =
            sqrt(distanceToView * distanceToView - PLANET_RADIUS * PLANET_RADIUS) / distanceToView;

        Vector3 toNewOrigin = toViewpointN * (distanceToView - correction);
        Vector3 newOrigin   = centerOfTheEarth + toNewOrigin;
        Vector3 correctView = lookAt - newOrigin;

        correctViewN = normalize(correctView);
    } else {
        correctViewN = viewDirectionN;
    }

    // Sun angle (gamma) - no correction

    double dotProductSun = dot(viewDirectionN, directionToSunN);
    params.gamma         = acos(dotProductSun);

    // Shadow angle - requires correction

    const double effectiveElevation = groundLevelSolarElevationAtOrigin;
    const double effectiveAzimuth   = groundLevelSolarAzimuthAtOrigin;
    const double shadowAngle        = effectiveElevation + PI * 0.5;

    Vector3 shadowDirectionN = Vector3(cos(shadowAngle) * cos(effectiveAzimuth),
                                       cos(shadowAngle) * sin(effectiveAzimuth),
                                       sin(shadowAngle));

    const double dotProductShadow = dot(correctViewN, shadowDirectionN);
    params.shadow                 = acos(dotProductShadow);

    // Zenith angle (theta) - corrected version stored in otherwise unused zero
    // angle

    double cosThetaCor = dot(correctViewN, toViewpointN);
    params.zero        = acos(cosThetaCor);

    // Zenith angle (theta) - uncorrected version goes outside

    double cosTheta = dot(viewDirectionN, toViewpointN);
    params.theta    = acos(cosTheta);

    return params;
}


/////////////////////////////////////////////////////////////////////////////////////
// Model evaluation for sky radiance and polarisation
/////////////////////////////////////////////////////////////////////////////////////

/// Evaluates piecewise linear approximation
double evalPL(const std::vector<float>::const_iterator coefs, const double factor) {
    return (double(coefs[1]) - double(coefs[0])) * factor + double(coefs[0]);
}

std::vector<float>::const_iterator PragueSkyModel::getCoefficients(const std::vector<float>& dataset,
                                                                   const int totalCoefsSingleConfig,
                                                                   const int elevation,
                                                                   const int altitude,
                                                                   const int visibility,
                                                                   const int albedo,
                                                                   const int wavelength) const {
    return dataset.cbegin() +
           (totalCoefsSingleConfig *
            (wavelength + size_t(channels) * elevation + channels * elevationsRad.size() * altitude +
             channels * elevationsRad.size() * altitudesRad.size() * albedo +
             channels * elevationsRad.size() * altitudesRad.size() * albedosRad.size() * visibility));
}

PragueSkyModel::InterpolationParameter PragueSkyModel::getInterpolationParameter(
    const double               queryVal,
    const std::vector<double>& breaks) const {
    // Clamp the value to the valid range.
    const double clamped = std::clamp(queryVal, breaks.front(), breaks.back());

    // Get the nearest greater parameter value.
    const auto next = std::upper_bound(breaks.cbegin() + 1, breaks.cend(), clamped);

    // Compute the index and float factor.
    InterpolationParameter parameter;
    parameter.index = int(next - breaks.cbegin()) - 1;
    if (next == breaks.cend()) {
        parameter.factor = 0.0;
    } else {
        parameter.factor = (clamped - *(next - 1)) / (*next - *(next - 1));
    }

    assert(0 <= parameter.index && parameter.index < breaks.size() &&
           (parameter.index < breaks.size() - 1 || parameter.factor == 0.f));
    assert(0 <= parameter.factor && parameter.factor <= 1);
    return parameter;
}

double PragueSkyModel::reconstruct(const AngleParameters&                   radianceParameters,
                                   const std::vector<float>::const_iterator channelParameters,
                                   const Metadata&                          metadata) const {
    // The original image was emphasized (for radiance only), re-parametrized into gamma-alpha space,
    // decomposed into sum of rankRad outer products of 'sun' and 'zenith' vectors and these vectors were
    // stored as piece-wise polynomial approximation. Here the process is reversed (for one point).

    double result = 0.0;
    for (int r = 0; r < metadata.rank; ++r) {
        // Restore the right value in the 'sun' vector
        const double sunParam = evalPL(channelParameters + metadata.sunOffset + r * metadata.sunStride +
                                           radianceParameters.gamma.index,
                                       radianceParameters.gamma.factor);

        // Restore the right value in the 'zenith' vector
        const double zenithParam = evalPL(channelParameters + metadata.zenithOffset +
                                              r * metadata.zenithStride + radianceParameters.alpha.index,
                                          radianceParameters.alpha.factor);

        // Accumulate their "outer" product
        result += sunParam * zenithParam;
    }

    // De-emphasize (for radiance only)
    if (!metadata.emphBreaks.empty()) {
        const double emphParam =
            evalPL(channelParameters + metadata.emphOffset + radianceParameters.zero.index,
                   radianceParameters.zero.factor);
        result *= emphParam;
        result = std::max(result, 0.0);
    }

    return result;
}

double PragueSkyModel::evaluateModel(const Parameters&         params,
                                     const double              wavelength,
                                     const std::vector<float>& data,
                                     const Metadata&           metadata) const {
    // Ignore wavelengths outside the dataset range.
    if (wavelength < channelStart || wavelength >= (channelStart + channels * channelWidth)) {
        return 0.0;
    }
    // Don't interpolate wavelengths inside the dataset range.
    const int channelIndex = int(floor((wavelength - channelStart) / channelWidth));

    // Translate angle values to indices and interpolation factors.
    AngleParameters angleParameters;
    angleParameters.gamma = getInterpolationParameter(params.gamma, metadata.sunBreaks);
    if (!metadata.emphBreaks.empty()) { // for radiance
        angleParameters.alpha =
            getInterpolationParameter(params.elevation < 0.0 ? params.shadow : params.zero,
                                      metadata.zenithBreaks);
        angleParameters.zero = getInterpolationParameter(params.zero, metadata.emphBreaks);
    } else { // for polarisation
        angleParameters.alpha = getInterpolationParameter(params.theta, metadata.zenithBreaks);
    }

    // Translate configuration values to indices and interpolation factors.
    const InterpolationParameter visibilityParam =
        getInterpolationParameter(params.visibility, visibilitiesRad);
    const InterpolationParameter albedoParam   = getInterpolationParameter(params.albedo, albedosRad);
    const InterpolationParameter altitudeParam = getInterpolationParameter(params.altitude, altitudesRad);
    const InterpolationParameter elevationParam =
        getInterpolationParameter(radiansToDegrees(params.elevation), elevationsRad);

    // Prepare parameters controlling the interpolation.
    ControlParameters controlParameters;
    for (int i = 0; i < 16; ++i) {
        const int visibilityIndex = std::min(visibilityParam.index + i / 8, int(visibilitiesRad.size() - 1));
        const int albedoIndex     = std::min(albedoParam.index + (i % 8) / 4, int(albedosRad.size() - 1));
        const int altitudeIndex   = std::min(altitudeParam.index + (i % 4) / 2, int(altitudesRad.size() - 1));
        const int elevationIndex  = std::min(elevationParam.index + i % 2, int(elevationsRad.size() - 1));

        controlParameters.coefficients[i] = getCoefficients(data,
                                                            metadata.totalCoefsSingleConfig,
                                                            elevationIndex,
                                                            altitudeIndex,
                                                            visibilityIndex,
                                                            albedoIndex,
                                                            channelIndex);
    }
    controlParameters.interpolationFactor[0] = visibilityParam.factor;
    controlParameters.interpolationFactor[1] = albedoParam.factor;
    controlParameters.interpolationFactor[2] = altitudeParam.factor;
    controlParameters.interpolationFactor[3] = elevationParam.factor;

    // Interpolate.
    const double result = interpolate<0, 0>(angleParameters, controlParameters, metadata);
    assert(metadata.emphBreaks.empty() || result >= 0.0); // polarisation can be negative

    return result;
}


/////////////////////////////////////////////////////////////////////////////////////
// Sky radiance
/////////////////////////////////////////////////////////////////////////////////////

double PragueSkyModel::skyRadiance(const Parameters& params, const double wavelength) const {
    if (!initialized) {
        throw NotInitializedException();
    }

	return evaluateModel(params, wavelength, dataRad, metadataRad);
}


/////////////////////////////////////////////////////////////////////////////////////
// Sun radiance
/////////////////////////////////////////////////////////////////////////////////////

double PragueSkyModel::sunRadiance(const Parameters& params, const double wavelength) const {
    if (!initialized) {
        throw NotInitializedException();
    }

    // Ignore wavelengths outside the dataset range.
    if (wavelength < SUN_RAD_START || wavelength >= SUN_RAD_END) {
        return 0.0;
    }

    // Return zero for rays not hitting the sun.
    if (params.gamma > SUN_RADIUS) {
        return 0.0;
    }

    // Compute index into the sun radiance table.
    const double idx = (wavelength - SUN_RAD_START) / SUN_RAD_STEP;
    assert(idx >= 0 && idx < std::size(SUN_RAD_TABLE) - 1);
    const int    idxInt   = int(floor(idx));
    const double idxFloat = idx - floor(idx);

    // Interpolate between the two closest values in the sun radiance table.
    const double sunRadiance =
        SUN_RAD_TABLE[idxInt] * (1.0 - idxFloat) + SUN_RAD_TABLE[idxInt + 1] * idxFloat;
    assert(sunRadiance > 0.0);

    // Compute transmittance towards the sun.
    const double tau = PragueSkyModel::transmittance(params, wavelength, std::numeric_limits<double>::max());
    assert(tau >= 0.0 && tau <= 1.0);

    // Combine.
    return sunRadiance * tau;
}


/////////////////////////////////////////////////////////////////////////////////////
// Polarisation
/////////////////////////////////////////////////////////////////////////////////////

double PragueSkyModel::polarisation(const Parameters& params, const double wavelength) const {
    if (!initialized) {
        throw NotInitializedException();
    }

    // If no polarisation data available
    if (metadataPol.rank == 0) {
        throw NoPolarisationException();
    }

    return -evaluateModel(params, wavelength, dataPol, metadataPol);
}


/////////////////////////////////////////////////////////////////////////////////////
// Transmittance
/////////////////////////////////////////////////////////////////////////////////////

std::vector<float>::const_iterator PragueSkyModel::getCoefficientsTrans(const int visibility,
                                                                        const int altitude,
                                                                        const int wavelength) const {
    return dataTransV.cbegin() +
           ((visibility * altitudesTrans.size() + altitude) * channels + wavelength) * rankTrans;
}

std::vector<float>::const_iterator PragueSkyModel::getCoefficientsTransBase(const int altitude,
                                                                            const int a,
                                                                            const int d) const {
    return dataTransU.cbegin() + size_t(altitude) * aDim * dDim * rankTrans +
           (size_t(d) * aDim + a) * rankTrans;
}

/// Intersects the given ray (assuming rayPosX == 0) with a circle at origin with the given radius.
/// \return In the case the ray intersects the circle, distance to the circle is returned, otherwise the
///         function returns negative number.
double intersectRayWithCircle2D(const double rayDirX,
                                const double rayDirY,
                                const double rayPosY,
                                const double circleRadius) {
    assert(rayPosY > 0.0);
    assert(circleRadius > 0.0);

    // Compute discriminant
    const double qA      = rayDirX * rayDirX + rayDirY * rayDirY;
    const double qB      = 2.0 * rayPosY * rayDirY;
    const double qC      = rayPosY * rayPosY - circleRadius * circleRadius;
    double       discrim = qB * qB - 4.0 * qA * qC;

    // No intersection or touch only
    if (discrim <= 0.0) {
        return -1.0;
    }

    discrim = std::sqrt(discrim);

    // Compute distances to both intersections
    const double d1 = (-qB + discrim) / (2.0 * qA);
    const double d2 = (-qB - discrim) / (2.0 * qA);

    // Try to take the nearest positive one
    const double distToIsect = (d1 > 0.0 && d2 > 0.0) ? std::min(d1, d2) : std::max(d1, d2);
    return distToIsect;
}

/// Auxiliary function for \ref toTransmittanceParams. Computes [altitude, distance] parameters from
/// coordinates of intersection of view ray and ground or atmosphere edge. Note: using floats here instead of
/// doubles will cause banding artifacts.
std::tuple<double, double> isectToAltitudeDistance(const double isectX, const double isectY) {
    // Distance to the intersection from world origin (not along ray as distToIsect in the calling method).
    const double isectDist = std::sqrt(isectX * isectX + isectY * isectY);
    assert(isectDist > 0.0);

    // Compute normalized and non-linearly scaled position in the atmosphere
    double altitude = std::clamp(isectDist - PLANET_RADIUS, 0.0, ATMOSPHERE_WIDTH);
    altitude        = std::pow(altitude / ATMOSPHERE_WIDTH, 1.0 / 3.0);
    double distance = std::acos(isectY / isectDist) * PLANET_RADIUS;
    distance        = std::sqrt(distance / DIST_TO_EDGE);
    distance        = std::sqrt(distance); // Calling twice sqrt, since it is faster than pow(...,0.25)
    distance        = std::min(1.0, distance);
    assert(0.0 <= altitude && altitude <= 1.0);
    assert(0.0 <= distance && distance <= 1.0);

    return std::make_tuple(altitude, distance);
}

PragueSkyModel::InterpolationParameter PragueSkyModel::getInterpolationParameterTrans(const double value,
                                                                                      const int    paramCount,
                                                                                      const int power) const {
    InterpolationParameter param;
    param.index  = std::min(int(value * paramCount), paramCount - 1);
    param.factor = 0.0;
    if (param.index < paramCount - 1) {
        param.factor =
            nonlerp(double(param.index) / paramCount, double(param.index + 1) / paramCount, value, power);
        param.factor = std::clamp(param.factor, 0.0, 1.0);
    }
    return param;
}

PragueSkyModel::TransmittanceParameters PragueSkyModel::toTransmittanceParams(const double theta,
                                                                              const double distance,
                                                                              const double altitude) const {
    assert(theta >= 0.0 && theta <= PI);
    assert(distance >= 0.0);
    assert(altitude >= 0.0);

    const double rayDirX = sin(theta);
    const double rayDirY = cos(theta);
    const double rayPosY = PLANET_RADIUS + altitude; // rayPosX == 0

    constexpr double ATMOSPHERE_EDGE = PLANET_RADIUS + ATMOSPHERE_WIDTH;

    // Find intersection of the ground-to-sun ray with edge of the atmosphere (in 2D)
    double           distToIsect  = -1.0;
    constexpr double LOW_ALTITUDE = 0.3;
    if (altitude < LOW_ALTITUDE) {
        // Special handling of almost zero altitude case to avoid numerical issues.
        if (theta <= 0.5 * PI) {
            distToIsect = intersectRayWithCircle2D(rayDirX, rayDirY, rayPosY, ATMOSPHERE_EDGE);
        } else {
            distToIsect = 0.0;
        }
    } else {
        distToIsect = intersectRayWithCircle2D(rayDirX, rayDirY, rayPosY, PLANET_RADIUS);
        if (distToIsect < 0.0) {
            distToIsect = intersectRayWithCircle2D(rayDirX, rayDirY, rayPosY, ATMOSPHERE_EDGE);
        }
    }
    // The ray should always hit either the edge of the atmosphere or the planet (we are starting inside the
    // atmosphere).
    assert(distToIsect >= 0.0);

    distToIsect = std::min(distToIsect, distance);

    // Compute intersection coordinates
    const double isectX = rayDirX * distToIsect;
    const double isectY = rayDirY * distToIsect + rayPosY;

    // Get the internal [altitude, distance] parameters
    const auto [altitudeParam, distanceParam] = isectToAltitudeDistance(isectX, isectY);

    // Convert to interpolation parameters
    TransmittanceParameters params;
    params.altitude = getInterpolationParameterTrans(altitudeParam, aDim, 3);
    params.distance = getInterpolationParameterTrans(distanceParam, dDim, 4);

    return params;
}

double PragueSkyModel::reconstructTrans(const int                     visibilityIndex,
                                        const int                     altitudeIndex,
                                        const TransmittanceParameters transParams,
                                        const int                     channelIndex) const {
    const std::vector<float>::const_iterator coefs =
        getCoefficientsTrans(visibilityIndex, altitudeIndex, channelIndex);

    // Load transmittance values for bi-linear interpolation
    std::array<double, 4> transmittance = { 0.0, 0.0, 0.0, 0.0 };
    int                   index         = 0;
    for (int a = transParams.altitude.index; a <= transParams.altitude.index + 1; ++a) {
        if (a < aDim) {
            for (int d = transParams.distance.index; d <= transParams.distance.index + 1; ++d) {
                if (d < dDim) {
                    const std::vector<float>::const_iterator baseCoefs =
                        getCoefficientsTransBase(altitudeIndex, a, d);
                    for (int i = 0; i < rankTrans; ++i) {
                        // Reconstruct transmittance value
                        transmittance[index] += double(baseCoefs[i]) * double(coefs[i]);
                    }
                    index++;
                }
            }
        }
    }

    // Perform bi-linear interpolation
    if (transParams.distance.factor > 0.f) {
        transmittance[0] = lerp(transmittance[0], transmittance[1], transParams.distance.factor);
        transmittance[1] = lerp(transmittance[2], transmittance[3], transParams.distance.factor);
    }
    transmittance[0] = std::max(transmittance[0], 0.0);

    if (transParams.altitude.factor > 0.f) {
        transmittance[1] = std::max(transmittance[1], 0.0);
        transmittance[0] = lerp(transmittance[0], transmittance[1], transParams.altitude.factor);
    }

    assert(transmittance[0] >= 0.0);
    return transmittance[0];
}

double PragueSkyModel::interpolateTrans(const int                     visibilityIndex,
                                        const InterpolationParameter  altitudeParam,
                                        const TransmittanceParameters transParams,
                                        const int                     channelIndex) const {
    // Get transmittance for the nearest lower altitude.
    double trans = reconstructTrans(visibilityIndex, altitudeParam.index, transParams, channelIndex);

    // Interpolate with transmittance for the nearest higher altitude if needed.
    if (altitudeParam.factor > 0.0) {
        const double transHigh =
            reconstructTrans(visibilityIndex, altitudeParam.index + 1, transParams, channelIndex);
        trans = lerp(trans, transHigh, altitudeParam.factor);
    }

    return trans;
}

double PragueSkyModel::transmittance(const Parameters& params,
                                     const double      wavelength,
                                     const double      distance) const {
    assert(distance > 0.0);

    if (!initialized) {
        throw NotInitializedException();
    }

    // Ignore wavelengths outside the dataset range.
    if (wavelength < channelStart || wavelength >= (channelStart + channels * channelWidth)) {
        return 0.0;
    }
    // Don't interpolate wavelengths inside the dataset range.
    const int channelIndex = int(floor((wavelength - channelStart) / channelWidth));

    // Translate configuration values to indices and interpolation factors.
    const InterpolationParameter visibilityParam =
        getInterpolationParameter(params.visibility, visibilitiesTrans);
    const InterpolationParameter altitudeParam = getInterpolationParameter(params.altitude, altitudesTrans);

    // Calculate position in the atmosphere.
    const TransmittanceParameters transParams = toTransmittanceParams(params.theta, distance, params.altitude);

    // Get transmittance for the nearest lower visibility.
    double trans = interpolateTrans(visibilityParam.index, altitudeParam, transParams, channelIndex);

    // Interpolate with transmittance for the nearest higher visibility if needed.
    if (visibilityParam.factor > 0.0) {
        const double transHigh =
            interpolateTrans(visibilityParam.index + 1, altitudeParam, transParams, channelIndex);

        trans = lerp(trans, transHigh, visibilityParam.factor);
    }

    // Transmittance is stored as a square root. Needs to be restored here.
    trans = trans * trans;
    assert(trans >= 0.0 && trans <= 1.0);

    return trans;
}