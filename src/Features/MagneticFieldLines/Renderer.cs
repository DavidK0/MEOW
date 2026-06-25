using Brutal.ImGuiApi;
using Brutal.Numerics;
using KSA;

namespace MEOW;

public sealed class FieldLine {
    public IReadOnlyList<double3> PointsBodyFrame { get; }

    public FieldLine(IReadOnlyList<double3> pointsBodyFrame) {
        PointsBodyFrame = pointsBodyFrame;
    }
}

public sealed class MagneticFieldOverlay : IOverlay {
    private readonly List<FieldLine> _fieldLines = new();

    private double _cachedRadius = -1;
    private int _cachedLineCount = -1;
    private double _cachedLineLength = -1;

    private const int DefaultLineCount = 96;
    private const double StartRadiusMultiplier = 1.05;
    private const double EndRadiusMultiplier = 3.0;

    private const int MaxTraceStepsPerDirection = 320;
    private const double TraceStepRadiusFraction = 0.035;
    private const double MinRadiusMultiplier = 1.00;

    public bool IsEnabled(MEOWSettings settings) {
        return settings.ShowFieldLines;
    }

    public void Update(BodyOverlayContext context, MEOWSettings settings) {
        if(!context.Profile.HasGlobalMagneticField)
            return;

        int lineCount = DefaultLineCount;
        double lineLength = context.BodyRadius * EndRadiusMultiplier;

        if(_cachedRadius == context.BodyRadius &&
           _cachedLineCount == lineCount &&
           Math.Abs(_cachedLineLength - lineLength) < 0.001)
            return;

        _cachedRadius = context.BodyRadius;
        _cachedLineCount = lineCount;
        _cachedLineLength = lineLength;

        RebuildFieldLines(context);
    }

    public void Draw(ImDrawListPtr draw_list, Camera camera, BodyOverlayContext context, MEOWSettings settings) {
        if(!context.Profile.HasGlobalMagneticField)
            return;

        foreach(var line in _fieldLines) {
            for(int i = 0; i < line.PointsBodyFrame.Count - 1; i++) {
                double3 a = context.BodyToWorld(line.PointsBodyFrame[i]);
                double3 b = context.BodyToWorld(line.PointsBodyFrame[i + 1]);

                ImDrawListExtensions.AddLine(
                    draw_list,
                    camera.EclToScreen(a),
                    camera.EclToScreen(b),
                    new ImColor8(255, 255, 255, 180),
                    2f);
            }
        }
    }


    private void RebuildFieldLines(BodyOverlayContext context) {
        _fieldLines.Clear();

        var model = context.Profile.MagneticFieldModel;

        if(model == null)
            return;

        double bodyRadius = context.BodyRadius;
        double startRadius = bodyRadius * StartRadiusMultiplier;
        double minRadius = bodyRadius * MinRadiusMultiplier;
        double maxRadius = bodyRadius * EndRadiusMultiplier;
        double stepLength = bodyRadius * TraceStepRadiusFraction;

        foreach(var seedDirection in GenerateFieldLineSeeds(DefaultLineCount)) {
            double3 seed = seedDirection * startRadius;

            List<double3> backward = TraceDirection(
                model,
                seed,
                -1.0,
                stepLength,
                minRadius,
                maxRadius);

            List<double3> forward = TraceDirection(
                model,
                seed,
                +1.0,
                stepLength,
                minRadius,
                maxRadius);

            var points = new List<double3>(backward.Count + 1 + forward.Count);

            backward.Reverse();
            points.AddRange(backward);
            points.Add(seed);
            points.AddRange(forward);

            if(points.Count >= 2)
                _fieldLines.Add(new FieldLine(points));
        }
    }

    private static List<double3> TraceDirection(
        IMagneticFieldModel model,
        double3 start,
        double sign,
        double stepLength,
        double minRadius,
        double maxRadius) {

        var points = new List<double3>();
        double3 p = start;

        for(int i = 0; i < MaxTraceStepsPerDirection; i++) {
            double3 next = Rk4Step(model, p, sign, stepLength);

            double r = VectorMath.Length(next);

            if(double.IsNaN(r) || double.IsInfinity(r))
                break;

            if(r < minRadius || r > maxRadius)
                break;

            // Prevent huge jumps near numerical singularities.
            if(VectorMath.Length(next - p) > stepLength * 2.0)
                break;

            points.Add(next);
            p = next;
        }

        return points;
    }

    private static double3 Rk4Step(
        IMagneticFieldModel model,
        double3 p,
        double sign,
        double h) {

        double3 k1 = FieldDirection(model, p, sign);
        double3 k2 = FieldDirection(model, p + k1 * (h * 0.5), sign);
        double3 k3 = FieldDirection(model, p + k2 * (h * 0.5), sign);
        double3 k4 = FieldDirection(model, p + k3 * h, sign);

        return p + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
    }

    private static double3 FieldDirection(
        IMagneticFieldModel model,
        double3 p,
        double sign) {

        double3 b = model.EvaluateField(p, new SimTime(0.0));
        double len = VectorMath.Length(b);

        if(len < 1e-12)
            return new double3(0, 0, 0);

        return sign * b / len;
    }

    private static IEnumerable<double3> GenerateFieldLineSeeds(int count) {
        if(count <= 0)
            yield break;

        double goldenAngle = Math.PI * (3.0 - Math.Sqrt(5.0));

        for(int i = 0; i < count; i++) {
            double y = 1.0 - (2.0 * i + 1.0) / count;
            double radius = Math.Sqrt(Math.Max(0.0, 1.0 - y * y));

            double theta = goldenAngle * i;

            double x = Math.Cos(theta) * radius;
            double z = Math.Sin(theta) * radius;

            yield return new double3(x, y, z);
        }
    }
}
