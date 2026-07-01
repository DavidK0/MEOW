using Brutal.Numerics;
using KSA;

namespace MEOW;

/// <summary>
/// Defines the local sampling space used by a field renderer.
/// </summary>
public readonly struct FieldRenderDomain {
    public required double ReferenceLength { get; init; }
    public required double HalfExtent { get; init; }
    public required double GridSpacing { get; init; }
    public required Func<double3, double3> LocalToEclPoint { get; init; }
    public required Func<double3, double3> EclToLocalPoint { get; init; }
    public required Func<double3, double3> LocalToEclVector { get; init; }
    public required Func<double3, double3> EclToLocalVector { get; init; }

    public static FieldRenderDomain AroundBody(IParentBody body) {
        var localToEcl = body.GetCci2Cce();
        var eclToLocal = localToEcl.Inverse();

        return new FieldRenderDomain {
            ReferenceLength = body.MeanRadius,
            HalfExtent = body.MeanRadius * 4.0,
            GridSpacing = body.MeanRadius,
            LocalToEclPoint = point =>
                body.GetPositionEcl() + localToEcl * point,
            EclToLocalPoint = point =>
                eclToLocal * (point - body.GetPositionEcl()),
            LocalToEclVector = vector => localToEcl * vector,
            EclToLocalVector = vector => eclToLocal * vector
        };
    }
}
