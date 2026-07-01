using Brutal.Numerics;

namespace MEOW;

/// <summary>
/// A vector field expressed in a body's local coordinate frame.
/// </summary>
public interface IFieldModel {
    void BeginSimulationStep() { }

    double3 EvaluateField(double3 positionBodyFrame);
}
