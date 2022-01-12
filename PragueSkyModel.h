#include <exception>
#include <string>
#include <vector>

/// Physically-based sky model by Wilkie et al. [2020]. Improves on previous work especially in accuracy of
/// sunset scenarios. Based on reconstruction of radiance from a small dataset fitted to a large set of images
/// obtained by brute force atmosphere simulation.
///
/// Provides evaluation of spectral sky radiance, sun radiance, transmittance and polarisation for observer at
/// a specific altitude above ground. The range of configurations depends on supplied dataset. The full
/// version models atmosphere of visibility (meteorological range) from 20 km to 131.8 km for sun elevations
/// from -4.2 degrees to 90 degrees, observer altitudes from 0 km to 15 km and ground albedo from 0 to 1, and
/// provides results for wavelengths from 320 nm to 760 nm.
///
/// Usage:
/// 1. Create PragueSkyModel object by calling its constructor with a path to the dataset file.
/// 2. The model is parametrized by several values that are gathered in PragueSkyModel::Parameters structure.
/// You can either fill this structure manually (in that case see its description) or just call
/// computeParameters, which will compute it for you based on a few basic parameters.
/// 3. Use the Parameters structure when calling skyRadiance, sunRadiance, transmittance, or polarisation
/// methods to obtain the respective quantities.
///
/// Throws:
/// - DatasetNotFoundException: if the specified dataset file could not be found
/// - DatasetReadException: if an error occurred while reading the dataset file
/// - NoPolarisationException: if the polarisation method is called but the model does not contain
/// polarisation data
///
/// Note:
/// The entire model is written in a single class and does not depend on anything except of STL. It defines a
/// simple Vector3 class to simplify working with points and directions and expects using this class when
/// passing viewing point and direction to the computeParameters method.
class PragueSkyModel {



/////////////////////////////////////////////////////////////////////////////////////
// Public types
/////////////////////////////////////////////////////////////////////////////////////
public:
    /// Exception thrown by the constructor if the passed dataset file could not be found.
    class DatasetNotFoundException : public std::exception {
    private:
        const std::string message;

    public:
        DatasetNotFoundException(const std::string& filename)
            : message(std::string("Dataset file ") + filename + std::string(" not found")) {}

        virtual const char* what() const throw() { return message.c_str(); }
    };

    /// Exception thrown by the constructor if an error occurred while reading the passed dataset file.
    class DatasetReadException : public std::exception {
    private:
        const std::string message;

    public:
        DatasetReadException(const std::string& parameterName)
            : message(std::string("Dataset reading failed at ") + parameterName) {}

        virtual const char* what() const throw() { return message.c_str(); }
    };

    /// Exception thrown by the polarisation method if the dataset passed to the constructor does not contain
    /// polarisation data.
    class NoPolarisationException : public std::exception {
    private:
        const std::string message;

    public:
        NoPolarisationException()
            : message(std::string("The supplied dataset does not contain polarisation data")) {}

        virtual const char* what() const throw() { return message.c_str(); }
    };

    /// A simple 3D vector implementation. Provides some basic operations.
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

        Vector3 operator*(const double factor) const { return Vector3(x * factor, y * factor, z * factor); }

        Vector3 operator/(const double factor) const { return Vector3(x / factor, y / factor, z / factor); }

        bool isZero() const { return x == 0.0 && y == 0.0 && z == 0.0; }
    };

    /// Structure holding all parameters necessary for querying the model.
    struct Parameters {
        /// Angle between view direction and direction to zenith in radians, supported values in range [0,
        /// PI].
        double theta;

        /// Angle between view direction and direction to sun in radians, supported values in range [0, PI].
        double gamma;

        /// Altitude-corrected angle between view direction and direction perpendicular to a shadow plane (=
        /// direction to sun rotated PI / 2 towards direction to zenith) in radians, used for negative solar
        /// elevations only, supported values in range [0, PI]
        double shadow;

        /// Altitude-corrected version of the theta angle in radians, supported values in range [0, PI].
        double zero;

        /// Sun elevation at view point in radians, supported values in range [-0.073, PI/2] (for full
        /// dataset). For view points above ground differs from the ground level sun elevation expected by the
        /// computeParameters method.
        double elevation;

        /// Altitude of view point in meters, supported values in range [0, 15000] (for full dataset).
        double altitude;

        /// Horizontal visibility (meteorological range) at ground level in kilometers, supported values in
        /// range [20, 131.8] (for full dataset).
        double visibility;

        /// Ground albedo, supported values in range [0, 1] (for full dataset).
        double albedo;
    };



/////////////////////////////////////////////////////////////////////////////////////
// Private data
/////////////////////////////////////////////////////////////////////////////////////
private:
    int    channels;
    double channelStart;
    double channelWidth;

    // Radiance metadata

    std::vector<double> visibilitiesRad;
    std::vector<double> albedosRad;
    std::vector<double> altitudesRad;
    std::vector<double> elevationsRad;

    int rankRad;

    int                 sunOffsetRad;
    int                 sunStrideRad;
    std::vector<double> sunBreaksRad;

    int                 zenithOffsetRad;
    int                 zenithStrideRad;
    std::vector<double> zenithBreaksRad;

    int                 emphOffsetRad;
    std::vector<double> emphBreaksRad;

    int totalCoefsSingleConfigRad;
    int totalCoefsAllConfigsRad;
    int totalConfigsRad;

    // Radiance data

    std::vector<float> datasetRad;

    // Tranmittance metadata

    int                aDim;
    int                dDim;
    int                rankTrans;
    std::vector<float> altitudesTrans;
    std::vector<float> visibilitiesTrans;

    // Tranmittance data

    std::vector<float> datasetTransU;
    std::vector<float> datasetTransV;

    // Polarisation metadata

    int rankPol;

    int                 sunOffsetPol;
    int                 sunStridePol;
    std::vector<double> sunBreaksPol;

    int                 zenithOffsetPol;
    int                 zenithStridePol;
    std::vector<double> zenithBreaksPol;

    int totalCoefsSingleConfigPol;
    int totalCoefsAllConfigsPol;

    // Polarisation data

    std::vector<float> datasetPol;



/////////////////////////////////////////////////////////////////////////////////////
// Public methods
/////////////////////////////////////////////////////////////////////////////////////
public:
    /// Prepares the model and loads the given dataset file into memory.
    ///
	/// Throws:
    /// - DatasetNotFoundException: if the specified dataset file could not be found
    /// - DatasetReadException: if an error occurred while reading the dataset file
    PragueSkyModel(const std::string& filename);

    /// Computes all the parameters in the Parameters structure necessary for querying the model.
    ///
    /// Expects view point and direction, sun elevation and azimuth at origin, ground level visibility and
    /// ground albedo. Assumes origin at [0,0,0] with Z axis pointing up. Thus view point [0, 0, 100] defines
    /// observer altitude 100 m. Range of available values depends on the used dataset. The full version
    /// supports altitude from [0, 15000], elevation from [-0.073, PI/2], azimuth from [0, PI], visibility
    /// from [20, 131.8], and albedo from [0, 1]. Values outside range of the used dataset are clamped to the
    /// nearest supported value.
    Parameters computeParameters(const Vector3& viewPoint,
                                 const Vector3& viewDirection,
                                 const double   groundLevelSolarElevationAtOrigin,
                                 const double   groundLevelSolarAzimuthAtOrigin,
                                 const double   visibility,
                                 const double   albedo) const;

    /// Computes sky radiance only (without direct sun contribution) for given parameters and wavelength (full
    /// dataset supports wavelengths from 320 nm to 760 nm).
    double skyRadiance(const Parameters& params, const double wavelength) const;

    /// Computes sun radiance only (without radiance inscattered from the sky) for given parameters and
    /// wavelength (full dataset supports wavelengths from 320 nm to 760 nm).
    ///
    /// Checks whether the parameters correspond to view direction hitting the sun and returns 0 if not.
    double sunRadiance(const Parameters& params, const double wavelength) const;

    /// Computes degree of polarisation for given parameters and wavelength (full
    /// dataset supports wavelengths from 320 nm to 760 nm). Can be negative.
    ///
    /// Throws NoPolarisationException if the polarisation method is called but the model does not contain
    /// polarisation data.
    double polarisation(const Parameters& params, const double wavelength) const;

    /// Computes transmittance between view point and a point certain distance away from it along view
    /// direction.
    ///
    /// Expects the Parameters structure, wavelength (full dataset supports wavelengths from 320 nm
    /// to 760 nm) and the distance (any positive number, use std::numeric_limits<double>::max() for
    /// infinity).
    double transmittance(const Parameters& params, const double wavelength, const double distance) const;


/////////////////////////////////////////////////////////////////////////////////////
// Private methods
/////////////////////////////////////////////////////////////////////////////////////
private:
    void readRadiance(FILE* handle);
    void readTransmittance(FILE* handle);
    void readPolarisation(FILE* handle);

    std::vector<float>::const_iterator controlParams(const std::vector<float>& dataset,
                                                     const int                 totalCoefsSingleConfig,
                                                     const int                 elevation,
                                                     const int                 altitude,
                                                     const int                 visibility,
                                                     const int                 albedo,
                                                     const int                 wavelength) const;

    double reconstruct(const double                             gamma,
                       const double                             alpha,
                       const double                             zero,
                       const int                                gammaSegment,
                       const int                                alphaSegment,
                       const int                                zeroSegment,
                       const std::vector<float>::const_iterator controlParams) const;

    double reconstructPol(const double                             gamma,
                          const double                             alpha,
                          const int                                gammaSegment,
                          const int                                alphaSegment,
                          const std::vector<float>::const_iterator controlParams) const;

    double interpolateElevation(double elevation,
                                int    altitude,
                                int    visibility,
                                int    albedo,
                                int    wavelength,
                                double gamma,
                                double alpha,
                                double zero,
                                int    gammaSegment,
                                int    alphaSegment,
                                int    zeroSegment) const;

    double interpolateAltitude(double elevation,
                               double altitude,
                               int    visibility,
                               int    albedo,
                               int    wavelength,
                               double gamma,
                               double alpha,
                               double zero,
                               int    gammaSegment,
                               int    alphaSegment,
                               int    zeroSegment) const;

    double interpolateVisibility(double elevation,
                                 double altitude,
                                 double visibility,
                                 int    albedo,
                                 int    wavelength,
                                 double gamma,
                                 double alpha,
                                 double zero,
                                 int    gammaSegment,
                                 int    alphaSegment,
                                 int    zeroSegment) const;

    double interpolateAlbedo(double elevation,
                             double altitude,
                             double visibility,
                             double albedo,
                             int    wavelength,
                             double gamma,
                             double alpha,
                             double zero,
                             int    gammaSegment,
                             int    alphaSegment,
                             int    zeroSegment) const;

    double interpolateWavelength(double elevation,
                                 double altitude,
                                 double visibility,
                                 double albedo,
                                 double wavelength,
                                 double gamma,
                                 double alpha,
                                 double zero,
                                 int    gammaSegment,
                                 int    alphaSegment,
                                 int    zeroSegment) const;

    std::vector<float>::const_iterator transmittanceCoefsIndex(const int visibility, const int altitude, const int wavelength) const;

    void transmittanceInterpolateWaveLength(const int    visibility,
                                            const int    altitude,
                                            const int    wavelengthLow,
                                            const int    wavelengthInc,
                                            const double wavelengthW,
                                            double*      coefficients) const;

    double calcTransmittanceSVDAltitude(const int    visibility,
                                        const int    altitude,
                                        const int    wavelengthLow,
                                        const int    wavelengthInc,
                                        const double wavelengthFactor,
                                        const int    aInt,
                                        const int    dInt,
                                        const int    aInc,
                                        const int    dInc,
                                        const double wa,
                                        const double wd) const;

    double calcTransmittanceSVD(const double a,
                                const double d,
                                const int    visibility,
                                const int    wavelengthLow,
                                const int    wavelengthInc,
                                const double wavelengthFactor,
                                const int    altitudeLow,
                                const int    altitudeInc,
                                const double altitudeFactor) const;

    double interpolateElevationPol(double elevation,
                                   int    altitude,
                                   int    visibility,
                                   int    albedo,
                                   int    wavelength,
                                   double gamma,
                                   double alpha,
                                   int    gammaSegment,
                                   int    alphaSegment) const;

    double interpolateAltitudePol(double elevation,
                                  double altitude,
                                  int    visibility,
                                  int    albedo,
                                  int    wavelength,
                                  double gamma,
                                  double alpha,
                                  int    gammaSegment,
                                  int    alphaSegment) const;

    double interpolateVisibilityPol(double elevation,
                                    double altitude,
                                    double visibility,
                                    int    albedo,
                                    int    wavelength,
                                    double gamma,
                                    double alpha,
                                    int    gammaSegment,
                                    int    alphaSegment) const;

    double interpolateAlbedoPol(double elevation,
                                double altitude,
                                double visibility,
                                double albedo,
                                int    wavelength,
                                double gamma,
                                double alpha,
                                int    gammaSegment,
                                int    alphaSegment) const;

    double interpolateWavelengthPol(double elevation,
                                    double altitude,
                                    double visibility,
                                    double albedo,
                                    double wavelength,
                                    double gamma,
                                    double alpha,
                                    int    gammaSegment,
                                    int    alphaSegment) const;
};