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
    };

private:
    // Radiance metadata

    int     turbidities;
    double* turbidity_vals;

    int     albedos;
    double* albedo_vals;

    int     altitudes;
    double* altitude_vals;

    int     elevations;
    double* elevation_vals;

    int    channels;
    double channel_start;
    double channel_width;

    int tensor_components;

    int     sun_nbreaks;
    int     sun_offset;
    int     sun_stride;
    double* sun_breaks;

    int     zenith_nbreaks;
    int     zenith_offset;
    int     zenith_stride;
    double* zenith_breaks;

    int     emph_nbreaks;
    int     emph_offset;
    double* emph_breaks;

    int total_coefs_single_config; // this is for one specific configuration
    int total_coefs_all_configs;
    int total_configs;

    // Radiance data

    double* radiance_dataset;

    // Tranmittance metadata

    int    trans_n_a;
    int    trans_n_d;
    int    trans_turbidities;
    int    trans_altitudes;
    int    trans_rank;
    float* transmission_altitudes;
    float* transmission_turbities;

    // Tranmittance data

    float* transmission_dataset_U;
    float* transmission_dataset_V;

    // Polarisation metadata

    int tensor_components_pol;

    int     sun_nbreaks_pol;
    int     sun_offset_pol;
    int     sun_stride_pol;
    double* sun_breaks_pol;

    int     zenith_nbreaks_pol;
    int     zenith_offset_pol;
    int     zenith_stride_pol;
    double* zenith_breaks_pol;

    int total_coefs_single_config_pol; // this is for one specific configuration
    int total_coefs_all_configs_pol;

    // Polarisation data

    double* polarisation_dataset;

public:
    PragueSkyModel(const char* library_path);

    ~PragueSkyModel();

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

    const double* controlParams(const double* dataset,
                                const int     total_coefs_single_config,
                                const int     elevation,
                                const int     altitude,
                                const int     turbidity,
                                const int     albedo,
                                const int     wavelength) const;

    double reconstruct(const double  gamma,
                       const double  alpha,
                       const double  zero,
                       const int     gamma_segment,
                       const int     alpha_segment,
                       const int     zero_segment,
                       const double* control_params) const;

    double reconstructPol(const double  gamma,
                          const double  alpha,
                          const int     gamma_segment,
                          const int     alpha_segment,
                          const double* control_params) const;

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

    const float* transmittanceCoefsIndex(const int turbidity, const int altitude, const int wavelength) const;

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