using Brutal.Numerics;

namespace MEOW;

public static class VectorMath {
    private const float FloatEpsilon = 1e-6f;
    private const double DoubleEpsilon = 1e-12;

    public static float Length(float3 v) {
        return MathF.Sqrt(
            v.X * v.X +
            v.Y * v.Y +
            v.Z * v.Z
        );
    }

    public static double Length(double3 v) {
        return Math.Sqrt(
            v.X * v.X +
            v.Y * v.Y +
            v.Z * v.Z
        );
    }

    public static float3 Normalize(float3 v) {
        float len = Length(v);

        if(len <= FloatEpsilon) {
            return new float3(0.0f, 0.0f, 0.0f);
        }

        return new float3(
            v.X / len,
            v.Y / len,
            v.Z / len
        );
    }

    public static double3 Normalize(double3 v) {
        double len = Length(v);

        if(len <= DoubleEpsilon) {
            return new double3(0.0, 0.0, 0.0);
        }

        return new double3(
            v.X / len,
            v.Y / len,
            v.Z / len
        );
    }

    public static float3 NormalizeToFloat(double3 v) {
        double len = Length(v);

        if(len <= DoubleEpsilon) {
            return new float3(0.0f, 0.0f, 0.0f);
        }

        return new float3(
            (float)(v.X / len),
            (float)(v.Y / len),
            (float)(v.Z / len)
        );
    }

    public static bool IsZero(float3 v) {
        return Length(v) <= FloatEpsilon;
    }

    public static bool IsZero(double3 v) {
        return Length(v) <= DoubleEpsilon;
    }

    public static float3 Cross(float3 a, float3 b) {
        return new float3(
            a.Y * b.Z - a.Z * b.Y,
            a.Z * b.X - a.X * b.Z,
            a.X * b.Y - a.Y * b.X
        );
    }

    public static double3 Cross(double3 a, double3 b) {
        return new double3(
            a.Y * b.Z - a.Z * b.Y,
            a.Z * b.X - a.X * b.Z,
            a.X * b.Y - a.Y * b.X
        );
    }

    public static double Dot(double3 a, double3 b) {
        return
            a.X * b.X +
            a.Y * b.Y +
            a.Z * b.Z;
    }

    public static double3 RotateFromTo(double3 v, double3 from, double3 to) {
        from = Normalize(from);
        to = Normalize(to);

        double3 axis = Cross(from, to);
        double axisLength = Length(axis);

        double dot = Dot(from, to);

        if(axisLength < 1e-12) {
            return dot > 0.0 ? v : new double3(-v.X, -v.Y, -v.Z);
        }

        axis = axis / axisLength;

        double angle = Math.Atan2(axisLength, dot);

        return RotateAroundAxis(v, axis, angle);
    }

    public static double3 RotateAroundAxis(double3 v, double3 axis, double angle) {
        double cos = Math.Cos(angle);
        double sin = Math.Sin(angle);

        return v * cos
            + Cross(axis, v) * sin
            + axis * Dot(axis, v) * (1.0 - cos);
    }
}
