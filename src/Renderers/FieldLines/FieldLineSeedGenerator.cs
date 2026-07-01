using Brutal.Numerics;

namespace MEOW;

internal static class FieldLineSeedGenerator {
    public static IEnumerable<double3> LatticeAroundBody(
        double bodyRadius,
        double radiusMultiplier,
        int pointsPerAxis) {

        if(pointsPerAxis <= 0)
            yield break;

        double halfSideLength = bodyRadius * radiusMultiplier;

        if(pointsPerAxis == 1) {
            yield return default;
            yield break;
        }

        double spacing = 2.0 * halfSideLength / (pointsPerAxis - 1);

        for(int x = 0; x < pointsPerAxis; x++) {
            for(int y = 0; y < pointsPerAxis; y++) {
                for(int z = 0; z < pointsPerAxis; z++) {
                    yield return new double3(
                        -halfSideLength + x * spacing,
                        -halfSideLength + y * spacing,
                        -halfSideLength + z * spacing);
                }
            }
        }
    }
}
