using Brutal.Numerics;
using KSA;

namespace MEOW;

public readonly struct BodyOverlayContext {
    public required SimTime Time { get; init; }
    public required FieldRenderDomain Domain { get; init; }
    public required IFieldModel FieldModel { get; init; }
    public required BodyOverlayProfile Profile { get; init; }

    public double BodyRadius => Domain.ReferenceLength;
    public double3 BodyToWorld(double3 point) => Domain.LocalToEclPoint(point);
    public double3 WorldToBody(double3 point) => Domain.EclToLocalPoint(point);
}
