using Brutal.Numerics;
using KSA;

namespace MEOW;

public static class MEOWMath {
    public static double Lerp(double a, double b, double t) {
        return a + (b - a) * t;
    }
}