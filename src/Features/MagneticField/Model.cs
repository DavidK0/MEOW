using Brutal.Numerics;
using KSA;

namespace MEOW;

public sealed class CompositeMagneticFieldModel : IFieldModel {
    private readonly IFieldModel[] _models;

    public CompositeMagneticFieldModel(params IFieldModel[] models) {
        _models = models;
    }

    public void BeginSimulationStep() {
        foreach(IFieldModel model in _models)
            model.BeginSimulationStep();
    }

    public double3 EvaluateField(double3 positionBodyFrame) {
        double3 total = new double3(0.0, 0.0, 0.0);

        foreach(IFieldModel model in _models) {
            total += model.EvaluateField(positionBodyFrame);
        }

        return total;
    }
}

public static class MagneticFieldModelFactory {
    public static BodyOverlayProfile CreateEarthProfile() {
        string userDocs = Environment.GetFolderPath(Environment.SpecialFolder.MyDocuments);
        string igrf14coeffsPath = Path.Combine(
            userDocs,
            "My Games",
            "Kitten Space Agency",
            "mods",
            "MEOW",
            "igrf14coeffs.txt");

        IGRFCoefficients igrf = IGRFCoefficients.LoadFromFile(igrf14coeffsPath);

        IGsmTransform gsm = new ApproxGsmTransform(igrf);

        IFieldModel field = new CompositeMagneticFieldModel(
            new IGRFModel(igrf),
            new T96Model(gsm));

        return new BodyOverlayProfile {
            BodyName = "Earth",
            MagneticFieldModel = field,
            GsmTransform = gsm
        };
    }
}

public interface IGsmTransform {
    void BeginSimulationStep() { }

    double3 BodyToGsmPosition(double3 bodyPosition);
    double3 GsmToBodyPosition(double3 gsmPosition);

    double3 BodyToGsmVector(double3 bodyVector);
    double3 GsmToBodyVector(double3 gsmVector);

    double GetDipoleTiltRadians();
}

public sealed class ApproxGsmTransform : IGsmTransform {
    private readonly IGRFCoefficients _igrf;

    private const double BasisUpdateAngleRadians = 1.0e-6;
    private static readonly double BasisUpdateDotThreshold =
        Math.Cos(BasisUpdateAngleRadians);

    private bool _hasBasis;
    private IParentBody? _cachedBody;
    private double3 _cachedSunDirection;
    private double3 _cachedDipoleDirection;
    private double3 _xGsmBody;
    private double3 _yGsmBody;
    private double3 _zGsmBody;

    public ApproxGsmTransform(IGRFCoefficients igrf) {
        _igrf = igrf;
    }

    public void BeginSimulationStep() {
        IParentBody earth = Program.ControlledVehicle.Orbit.Parent; // Assuming the controlled vehicle is orbiting Earth
        double3 earthPosEcl = earth.GetPositionEcl();
        Celestial earthCelestial = (Celestial)earth;
        IParentBody sun = earthCelestial.Parent;
        double3 sunPosEcl = sun.GetPositionEcl();
        double3 earthToSunEcl = sunPosEcl - earthPosEcl;

        // ECL/CCE direction -> CCI direction
        double3 earthToSunCci = earth.GetCce2Cci() * earthToSunEcl;
        double3 sunDirection = VectorMath.Normalize(earthToSunCci);

        double3 dipoleDirection = _igrf.GetDipoleNorthCcf(2026d); // Placeholder year

        if(_hasBasis &&
           ReferenceEquals(_cachedBody, earth) &&
           VectorMath.Dot(_cachedSunDirection, sunDirection) >= BasisUpdateDotThreshold &&
           VectorMath.Dot(_cachedDipoleDirection, dipoleDirection) >= BasisUpdateDotThreshold)
            return;

        _cachedBody = earth;
        _cachedSunDirection = sunDirection;
        _cachedDipoleDirection = dipoleDirection;

        _xGsmBody = sunDirection;
        _yGsmBody = VectorMath.Normalize(
            VectorMath.Cross(dipoleDirection, _xGsmBody));
        _zGsmBody = VectorMath.Normalize(
            VectorMath.Cross(_xGsmBody, _yGsmBody));
        _hasBasis = true;
    }

    public double3 BodyToGsmPosition(double3 pBody) {
        GetBasis(out var x, out var y, out var z);

        return new double3(
            VectorMath.Dot(pBody, x),
            VectorMath.Dot(pBody, y),
            VectorMath.Dot(pBody, z));
    }

    public double3 GsmToBodyPosition(double3 pGsm) {
        GetBasis(out var x, out var y, out var z);

        return pGsm.X * x +
               pGsm.Y * y +
               pGsm.Z * z;
    }

    public double3 BodyToGsmVector(double3 vBody) {
        GetBasis(out var x, out var y, out var z);

        return new double3(
            VectorMath.Dot(vBody, x),
            VectorMath.Dot(vBody, y),
            VectorMath.Dot(vBody, z));
    }

    public double3 GsmToBodyVector(double3 vGsm) {
        GetBasis(out var x, out var y, out var z);

        return vGsm.X * x +
               vGsm.Y * y +
               vGsm.Z * z;
    }

    public double GetDipoleTiltRadians() {
        GetBasis(out var x, out _, out _);

        return Math.Asin(Math.Clamp(
            VectorMath.Dot(_cachedDipoleDirection, x),
            -1.0,
            1.0));
    }

    private void GetBasis(
        out double3 xGsmBody,
        out double3 yGsmBody,
        out double3 zGsmBody) {

        if(!_hasBasis)
            BeginSimulationStep();

        xGsmBody = _xGsmBody;
        yGsmBody = _yGsmBody;
        zGsmBody = _zGsmBody;
    }
}
