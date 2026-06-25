using Brutal.Numerics;
using KSA;

namespace MEOW;

public readonly struct BodyOverlayContext {
    public required SimTime Time { get; init; }
    public required double BodyRadius { get; init; }
    public required Func<double3, double3> BodyToWorld { get; init; }
    public required Func<double3, double3> WorldToBody { get; init; }
    public required BodyOverlayProfile Profile { get; init; }
}