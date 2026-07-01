using Brutal.Numerics;
using KSA;

namespace MEOW;

/// <summary>
/// Newtonian point-mass gravity evaluated in the absolute ECL frame.
/// </summary>
public sealed class GravitationalFieldModel : IFieldModel {
    private const double GravitationalConstant = 6.67430e-11;
    private readonly List<GravitySource> _sources = new();

    public void BeginSimulationStep() {
        _sources.Clear();

        CelestialSystem? system = Universe.CurrentSystem;
        if(system == null)
            return;

        foreach(Astronomical astronomical in system.All.AsSpan()) {
            if(astronomical is not Celestial body ||
               !double.IsFinite(body.Mass) ||
               body.Mass <= 0.0)
                continue;

            double minimumRadius = body is IParentBody parentBody
                ? Math.Max(parentBody.MeanRadius, 1.0)
                : 1.0;

            _sources.Add(new GravitySource(
                body.Mass,
                body.GetPositionEcl(),
                minimumRadius));
        }
    }

    public double3 EvaluateField(double3 positionEcl) {
        if(_sources.Count == 0)
            BeginSimulationStep();

        double3 acceleration = default;

        foreach(GravitySource source in _sources) {
            double3 displacement = source.PositionEcl - positionEcl;
            double distance = VectorMath.Length(displacement);
            if(!double.IsFinite(distance) || distance <= 0.0)
                continue;

            double clampedDistance = Math.Max(distance, source.MinimumRadius);
            double scale =
                GravitationalConstant * source.Mass /
                (clampedDistance * clampedDistance * distance);
            acceleration += displacement * scale;
        }

        return acceleration;
    }

    private readonly record struct GravitySource(
        double Mass,
        double3 PositionEcl,
        double MinimumRadius);
}
