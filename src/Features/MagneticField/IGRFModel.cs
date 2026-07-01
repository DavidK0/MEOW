using Brutal.Numerics;
using KSA;
using System.Globalization;

namespace MEOW;

public sealed class IGRFCoefficients {
    public const int MaxDegree = 13;

    // Coefficients indexed as [n, m].
    private readonly double[,] _g2025 = new double[MaxDegree + 1, MaxDegree + 1];
    private readonly double[,] _h2025 = new double[MaxDegree + 1, MaxDegree + 1];
    private readonly double[,] _gSv = new double[MaxDegree + 1, MaxDegree + 1];
    private readonly double[,] _hSv = new double[MaxDegree + 1, MaxDegree + 1];

    private IGRFCoefficients() {
    }

    public double GetG(int n, int m, double decimalYear) {
        double dt = decimalYear - 2025.0;
        return _g2025[n, m] + _gSv[n, m] * dt;
    }

    public double GetH(int n, int m, double decimalYear) {
        double dt = decimalYear - 2025.0;
        return _h2025[n, m] + _hSv[n, m] * dt;
    }

    public static IGRFCoefficients LoadFromFile(string path) {
        if(!File.Exists(path))
            throw new FileNotFoundException("Could not find IGRF coefficient file.", path);

        var coeffs = new IGRFCoefficients();

        foreach(string rawLine in File.ReadLines(path)) {
            string line = rawLine.Trim();

            if(line.Length == 0 || line.StartsWith("#"))
                continue;

            string[] parts = line.Split((char[]?)null, StringSplitOptions.RemoveEmptyEntries);

            if(parts.Length < 5)
                continue;

            // Data rows begin with "g" or "h".
            string kind = parts[0];

            if(kind != "g" && kind != "h")
                continue;

            int n = int.Parse(parts[1], CultureInfo.InvariantCulture);
            int m = int.Parse(parts[2], CultureInfo.InvariantCulture);

            if(n < 1 || n > MaxDegree || m < 0 || m > n)
                continue;

            // Your file has:
            // kind n m 1900 ... 2020 2025 SV
            //
            // Therefore:
            // parts[3]  = 1900.0
            // ...
            // parts[^2] = 2025.0
            // parts[^1] = 2025-2030 SV
            double value2025 = double.Parse(parts[^2], CultureInfo.InvariantCulture);
            double sv = double.Parse(parts[^1], CultureInfo.InvariantCulture);

            if(kind == "g") {
                coeffs._g2025[n, m] = value2025;
                coeffs._gSv[n, m] = sv;
            } else {
                coeffs._h2025[n, m] = value2025;
                coeffs._hSv[n, m] = sv;
            }
        }

        return coeffs;
    }

    public double3 GetDipoleNorthCcf(double decimalYear) {
        double g10 = GetG(1, 0, decimalYear);
        double g11 = GetG(1, 1, decimalYear);
        double h11 = GetH(1, 1, decimalYear);

        // KSA/CCF convention:
        // X = lon 0 equator
        // Y = lon 90 east equator
        // Z = geographic north
        return VectorMath.Normalize(new double3(
            -g11,
            -h11,
            -g10
        ));
    }
}

public sealed class IGRFModel : IBodyFrameFieldModel {
    private const int NMax = IGRFCoefficients.MaxDegree;

    // IGRF reference radius, in meters.
    // The IGRF spherical-harmonic expansion uses a = 6371.2 km.
    private const double ReferenceRadiusMeters = 6371200.0;

    private readonly IGRFCoefficients _coeffs;

    // Schmidt normalization and recursion constants.
    private readonly double[,] _snorm = new double[NMax + 1, NMax + 1];
    private readonly double[,] _k = new double[NMax + 1, NMax + 1];

    public IGRFModel(IGRFCoefficients coeffs) {
        _coeffs = coeffs;
        BuildSchmidtNormalization();
    }

    public double3 EvaluateField(double3 positionBodyFrame) {
        //SimTime time = Universe.GetElapsedSimTime();
        //double years = time.Years;
        double decimalYear = 2025.0;

        // KSA body frame -> IGRF frame
        double3 pIgrf = new double3(
            positionBodyFrame.X,
            positionBodyFrame.Z,
            positionBodyFrame.Y
        );

        double3 bIgrf = EvaluateFieldAtDecimalYear(pIgrf, decimalYear);

        // IGRF frame -> KSA body frame
        return new double3(
            bIgrf.X,
            bIgrf.Z,
            bIgrf.Y
        );
    }

    public double3 EvaluateFieldAtDecimalYear(double3 positionBodyFrame, double decimalYear) {
        double x = positionBodyFrame.X;
        double y = positionBodyFrame.Y;
        double z = positionBodyFrame.Z;

        double r = Math.Sqrt(x * x + y * y + z * z);

        if(r < 1.0)
            return new double3(0, 0, 0);

        // Body-frame convention:
        // y is geographic north.
        //
        // Spherical coordinates:
        // theta = geocentric colatitude, 0 at north pole.
        // phi   = longitude, positive east.
        double cosTheta = y / r;
        cosTheta = Math.Clamp(cosTheta, -1.0, 1.0);

        double sinTheta = Math.Sqrt(Math.Max(0.0, 1.0 - cosTheta * cosTheta));
        double phi = Math.Atan2(z, x);

        double sinPhi = Math.Sin(phi);
        double cosPhi = Math.Cos(phi);

        double[] sinMLon = new double[NMax + 1];
        double[] cosMLon = new double[NMax + 1];

        sinMLon[0] = 0.0;
        cosMLon[0] = 1.0;
        sinMLon[1] = sinPhi;
        cosMLon[1] = cosPhi;

        for(int m = 2; m <= NMax; m++) {
            sinMLon[m] = sinMLon[m - 1] * cosPhi + cosMLon[m - 1] * sinPhi;
            cosMLon[m] = cosMLon[m - 1] * cosPhi - sinMLon[m - 1] * sinPhi;
        }

        // Associated Legendre functions and their theta derivatives.
        //
        // This follows the common geomagnetic-model recursion style:
        // compute unnormalized P, multiply Gauss coefficients by Schmidt factors.
        double[,] p = new double[NMax + 1, NMax + 1];
        double[,] dp = new double[NMax + 1, NMax + 1];

        p[0, 0] = 1.0;
        dp[0, 0] = 0.0;

        for(int n = 1; n <= NMax; n++) {
            for(int m = 0; m <= n; m++) {
                if(n == m) {
                    p[n, m] = sinTheta * p[n - 1, m - 1];
                    dp[n, m] = sinTheta * dp[n - 1, m - 1] + cosTheta * p[n - 1, m - 1];
                } else if(n == 1 && m == 0) {
                    p[n, m] = cosTheta * p[n - 1, m];
                    dp[n, m] = cosTheta * dp[n - 1, m] - sinTheta * p[n - 1, m];
                } else {
                    p[n, m] = cosTheta * p[n - 1, m] - _k[n, m] * p[n - 2, m];
                    dp[n, m] = cosTheta * dp[n - 1, m]
                               - sinTheta * p[n - 1, m]
                               - _k[n, m] * dp[n - 2, m];
                }
            }
        }

        double br = 0.0;     // radial, positive outward
        double bTheta = 0.0; // positive toward increasing colatitude, i.e. southward
        double bPhi = 0.0;   // positive eastward

        double aOverR = ReferenceRadiusMeters / r;
        double arPower = aOverR * aOverR;

        for(int n = 1; n <= NMax; n++) {
            arPower *= aOverR;

            for(int m = 0; m <= n; m++) {
                double g = _coeffs.GetG(n, m, decimalYear) * _snorm[n, m];
                double h = _coeffs.GetH(n, m, decimalYear) * _snorm[n, m];

                double tmp = g * cosMLon[m] + h * sinMLon[m];

                br += arPower * (n + 1) * tmp * p[n, m];
                bTheta -= arPower * tmp * dp[n, m];

                if(m != 0) {
                    double tmpPhi = g * sinMLon[m] - h * cosMLon[m];

                    if(sinTheta > 1e-10) {
                        bPhi += arPower * m * tmpPhi * p[n, m] / sinTheta;
                    } else {
                        // Avoid singularity exactly at the geographic poles.
                        // Field-line tracing should rarely hit this exactly.
                        bPhi = 0.0;
                    }
                }
            }
        }

        // Convert spherical components back into your body frame.
        //
        // e_r     = outward
        // e_theta = southward
        // e_phi   = eastward
        double invHorizontal = sinTheta > 1e-10 ? 1.0 / sinTheta : 0.0;

        double3 eR = new double3(x / r, y / r, z / r);

        double3 eTheta;
        double3 ePhi;

        if(sinTheta > 1e-10) {
            eTheta = new double3(
                cosTheta * x / (r * sinTheta),
                -sinTheta,
                cosTheta * z / (r * sinTheta));

            ePhi = new double3(
                -z * invHorizontal / r,
                0.0,
                x * invHorizontal / r);
        } else {
            // Arbitrary stable basis at the pole.
            eTheta = new double3(1, 0, 0);
            ePhi = new double3(0, 0, 1);
        }

        // Result is in nanoTesla-scaled units.
        // For drawing field lines, only direction matters.
        return br * eR + bTheta * eTheta + bPhi * ePhi;
    }

    private void BuildSchmidtNormalization() {
        _snorm[0, 0] = 1.0;

        for(int n = 1; n <= NMax; n++) {
            _snorm[n, 0] = _snorm[n - 1, 0] * (2.0 * n - 1.0) / n;

            for(int m = 1; m <= n; m++) {
                double j = m == 1 ? 2.0 : 1.0;
                double factor = ((n - m + 1.0) * j) / (n + m);
                _snorm[n, m] = _snorm[n, m - 1] * Math.Sqrt(factor);
            }
        }

        for(int n = 2; n <= NMax; n++) {
            for(int m = 0; m <= n - 1; m++) {
                _k[n, m] =
                    (((n - 1.0) * (n - 1.0)) - (m * m))
                    / ((2.0 * n - 1.0) * (2.0 * n - 3.0));
            }
        }
    }
}
