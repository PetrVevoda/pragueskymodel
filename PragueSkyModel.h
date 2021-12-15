#include <exception>
#include <string>
#include <vector>

class PragueSkyModel {
public:
    class DatasetNotFoundException : public std::exception {
    private:
        const std::string message;

    public:
        DatasetNotFoundException(const std::string& filename)
            : message(std::string("Dataset file ") + filename + std::string(" not found")) {}

        virtual const char* what() const throw() { return message.c_str(); }
    };

    class DatasetReadException : public std::exception {
    private:
        const std::string message;

    public:
        DatasetReadException(const std::string& parameterName)
            : message(std::string("Dataset reading failed at ") + parameterName) {}

        virtual const char* what() const throw() { return message.c_str(); }
    };

	class NoPolarisationException : public std::exception {
	private:
		const std::string message;

	public:
        NoPolarisationException() : message(std::string("The supplied dataset does not contain polarisation data")) {}

		virtual const char* what() const throw() { return message.c_str(); }
	};

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

    struct Parameters {
        double theta;
        double gamma;
        double shadow;
        double zero;
        double elevation;
        double altitude;
        double visibility;
        double albedo;
    };

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

    int totalCoefsSingleConfigRad; // this is for one specific configuration
    int totalCoefsAllConfigsRad;
    int totalConfigsRad;

    // Radiance data

    std::vector<double> datasetRad;

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

    int totalCoefsSingleConfigPol; // this is for one specific configuration
    int totalCoefsAllConfigsPol;

    // Polarisation data

    std::vector<double> datasetPol;

public:
    PragueSkyModel(const std::string& filename);

    //   This computes the canonical angles of the model from
    //   a normalised view vector and solar elevation.
    Parameters computeParameters(const Vector3& viewpoint,
                                 const Vector3& viewDirection,
                                 const double   groundLevelSolarElevationAtOrigin,
                                 const double   groundLevelSolarAzimuthAtOrigin,
                                 const double   visibility,
                                 const double   albedo) const;

    //   theta  - zenith angle
    //   gamma  - sun angle
    //   shadow - angle from the shadow point, which is further 90 degrees above
    //   the sun zero   - angle from the zero point, which lies at the horizon
    //   below the sun altitude wavelength
    double skyRadiance(const Parameters& params, const double wavelength) const;

    /* ----------------------------------------------------------------------------

        arpragueskymodel_solar_radiance
        ---------------------------

        This computes transmittance between a point at 'altitude' and infinity in
        the direction 'theta' at a wavelength 'wavelength'.

    ----------------------------------------------------------------------------
  */
	double sunRadiance(const Parameters& params, const double wavelength) const;

	double polarisation(const Parameters& params, const double wavelength) const;

    /* ----------------------------------------------------------------------------

        arpragueskymodel_tau
        ------------------------------

        This computes transmittance between a point at 'altitude' and infinity in
        the direction 'theta' at a wavelength 'wavelength'.

    ----------------------------------------------------------------------------
  */
	double transmittance(const Parameters& params,
		const double wavelength,
		const double distance) const;

private:
    void readRadiance(FILE* handle);
    void readTransmittance(FILE* handle);
    void readPolarisation(FILE* handle);

    std::vector<double>::const_iterator controlParams(const std::vector<double>& dataset,
                                const int     totalCoefsSingleConfig,
                                const int     elevation,
                                const int     altitude,
                                const int     visibility,
                                const int     albedo,
                                const int     wavelength) const;

    double reconstruct(const double                              gamma,
                       const double                              alpha,
                       const double                              zero,
                       const int                                 gammaSegment,
                       const int                                 alphaSegment,
                       const int                                 zeroSegment,
                       const std::vector<double>::const_iterator controlParams) const;

    double reconstructPol(const double                              gamma,
                          const double                              alpha,
                          const int                                 gammaSegment,
                          const int                                 alphaSegment,
                          const std::vector<double>::const_iterator controlParams) const;

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