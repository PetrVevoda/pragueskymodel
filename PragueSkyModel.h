#include <vector>

class PragueSkyModel {
public:
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
    PragueSkyModel(const char* library_path);

    //   This computes the canonical angles of the model from
    //   a normalised view vector and solar elevation.
    void computeAngles(const Vector3& viewpoint,
                       const Vector3& viewDirection,
                       const double   groundLevelSolarElevationAtOrigin,
                       const double   groundLevelSolarAzimuthAtOrigin,
                       double*        solarElevationAtViewpoint,
                       double*        altitudeOfViewpoint,
                       double*        theta,
                       double*        gamma,
                       double*        shadow,
                       double*        zero) const;

    //   theta  - zenith angle
    //   gamma  - sun angle
    //   shadow - angle from the shadow point, which is further 90 degrees above
    //   the sun zero   - angle from the zero point, which lies at the horizon
    //   below the sun altitude wavelength
    double skyRadiance(const double theta,
                       const double gamma,
                       const double shadow,
                       const double zero,
                       const double elevation,
                       const double altitude,
                       const double turbidity,
                       const double albedo,
                       const double wavelength) const;

    /* ----------------------------------------------------------------------------

        arpragueskymodel_solar_radiance
        ---------------------------

        This computes transmittance between a point at 'altitude' and infinity in
        the direction 'theta' at a wavelength 'wavelength'.

    ----------------------------------------------------------------------------
  */
    double sunRadiance(const double theta,
                       const double gamma,
                       const double shadow,
                       const double zero,
                       const double elevation,
                       const double altitude,
                       const double turbidity,
                       const double albedo,
                       const double wavelength) const;

    double polarisation(const double theta,
                        const double gamma,
                        const double elevation,
                        const double altitude,
                        const double turbidity,
                        const double albedo,
                        const double wavelength) const;

    /* ----------------------------------------------------------------------------

        arpragueskymodel_tau
        ------------------------------

        This computes transmittance between a point at 'altitude' and infinity in
        the direction 'theta' at a wavelength 'wavelength'.

    ----------------------------------------------------------------------------
  */
    double transmittance(const double theta,
                         const double altitude,
                         const double turbidity,
                         const double wavelength,
                         const double distance) const;

private:
    void readRadiance(FILE* handle);
    void readTransmittance(FILE* handle);
    void readPolarisation(FILE* handle);

    std::vector<double>::const_iterator PragueSkyModel::controlParams(const std::vector<double>& dataset,
                                const int     totalCoefsSingleConfig,
                                const int     elevation,
                                const int     altitude,
                                const int     turbidity,
                                const int     albedo,
                                const int     wavelength) const;

    double reconstruct(const double                              gamma,
                       const double                              alpha,
                       const double                              zero,
                       const int                                 gamma_segment,
                       const int                                 alpha_segment,
                       const int                                 zero_segment,
                       const std::vector<double>::const_iterator controlParams) const;

    double reconstructPol(const double                              gamma,
                          const double                              alpha,
                          const int                                 gamma_segment,
                          const int                                 alpha_segment,
                          const std::vector<double>::const_iterator controlParams) const;

    double interpolateElevation(double elevation,
                                int    altitude,
                                int    turbidity,
                                int    albedo,
                                int    wavelength,
                                double gamma,
                                double alpha,
                                double zero,
                                int    gamma_segment,
                                int    alpha_segment,
                                int    zero_segment) const;

    double interpolateAltitude(double elevation,
                               double altitude,
                               int    turbidity,
                               int    albedo,
                               int    wavelength,
                               double gamma,
                               double alpha,
                               double zero,
                               int    gamma_segment,
                               int    alpha_segment,
                               int    zero_segment) const;

    double interpolateVisibility(double elevation,
                                 double altitude,
                                 double turbidity,
                                 int    albedo,
                                 int    wavelength,
                                 double gamma,
                                 double alpha,
                                 double zero,
                                 int    gamma_segment,
                                 int    alpha_segment,
                                 int    zero_segment) const;

    double interpolateAlbedo(double elevation,
                             double altitude,
                             double turbidity,
                             double albedo,
                             int    wavelength,
                             double gamma,
                             double alpha,
                             double zero,
                             int    gamma_segment,
                             int    alpha_segment,
                             int    zero_segment) const;

    double interpolateWavelength(double elevation,
                                 double altitude,
                                 double turbidity,
                                 double albedo,
                                 double wavelength,
                                 double gamma,
                                 double alpha,
                                 double zero,
                                 int    gamma_segment,
                                 int    alpha_segment,
                                 int    zero_segment) const;

    std::vector<float>::const_iterator transmittanceCoefsIndex(const int turbidity, const int altitude, const int wavelength) const;

    void transmittanceInterpolateWaveLength(const int    turbidity,
                                            const int    altitude,
                                            const int    wavelength_low,
                                            const int    wavelength_inc,
                                            const double wavelength_w,
                                            double*      coefficients) const;

    double calcTransmittanceSVDAltitude(const int    turbidity,
                                        const int    altitude,
                                        const int    wavelength_low,
                                        const int    wavelength_inc,
                                        const double wavelength_factor,
                                        const int    a_int,
                                        const int    d_int,
                                        const int    a_inc,
                                        const int    d_inc,
                                        const double wa,
                                        const double wd) const;

    double calcTransmittanceSVD(const double a,
                                const double d,
                                const int    turbidity,
                                const int    wavelength_low,
                                const int    wavelength_inc,
                                const double wavelength_factor,
                                const int    altitude_low,
                                const int    altitude_inc,
                                const double altitude_factor) const;

    double interpolateElevationPol(double elevation,
                                   int    altitude,
                                   int    turbidity,
                                   int    albedo,
                                   int    wavelength,
                                   double gamma,
                                   double alpha,
                                   int    gamma_segment,
                                   int    alpha_segment) const;

    double interpolateAltitudePol(double elevation,
                                  double altitude,
                                  int    turbidity,
                                  int    albedo,
                                  int    wavelength,
                                  double gamma,
                                  double alpha,
                                  int    gamma_segment,
                                  int    alpha_segment) const;

    double interpolateVisibilityPol(double elevation,
                                    double altitude,
                                    double turbidity,
                                    int    albedo,
                                    int    wavelength,
                                    double gamma,
                                    double alpha,
                                    int    gamma_segment,
                                    int    alpha_segment) const;

    double interpolateAlbedoPol(double elevation,
                                double altitude,
                                double turbidity,
                                double albedo,
                                int    wavelength,
                                double gamma,
                                double alpha,
                                int    gamma_segment,
                                int    alpha_segment) const;

    double interpolateWavelengthPol(double elevation,
                                    double altitude,
                                    double turbidity,
                                    double albedo,
                                    double wavelength,
                                    double gamma,
                                    double alpha,
                                    int    gamma_segment,
                                    int    alpha_segment) const;
};