using Brutal.Numerics;

namespace MEOW;

/// <summary>
/// A vector field expressed in the solar system's absolute ECL frame.
/// </summary>
public interface IFieldModel {
    void BeginSimulationStep() { }

    double3 EvaluateField(double3 positionEcl);
}
