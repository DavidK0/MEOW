using Brutal.Numerics;

namespace MEOW;

/// <summary>
/// Placeholder for a gravitational field model expressed in the body frame.
/// </summary>
public sealed class GravitationalFieldModel : IFieldModel {
    public double3 EvaluateField(double3 positionBodyFrame) {
        // TODO: Implement the selected gravity model.
        return default;
    }
}
