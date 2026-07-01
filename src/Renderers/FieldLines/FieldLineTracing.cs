using Brutal.Numerics;

namespace MEOW;

internal sealed class FieldLine {
    public FieldLine(IReadOnlyList<double3> pointsBodyFrame) {
        PointsBodyFrame = pointsBodyFrame;
    }

    public IReadOnlyList<double3> PointsBodyFrame { get; }
}

internal readonly record struct FieldTraceOptions(
    double StepLength,
    double MinRadius,
    double MaxRadius,
    int MaxStepsPerDirection,
    int MaxVisualPoints);

internal sealed class FieldLineBuilder {
    private const int MinimumPointCount = 8;

    private readonly FieldTracer _backward;
    private readonly FieldTracer _forward;
    private readonly double3 _seed;
    private readonly int _maxVisualPoints;

    public FieldLineBuilder(
        IFieldModel model,
        double3 seed,
        FieldTraceOptions options) {

        _seed = seed;
        _maxVisualPoints = options.MaxVisualPoints;
        _backward = new FieldTracer(model, seed, -1.0, options);
        _forward = new FieldTracer(model, seed, +1.0, options);
    }

    public bool IsAlive => _backward.IsAlive || _forward.IsAlive;

    public void Step() {
        if(_backward.IsAlive)
            _backward.Step();
        else
            _forward.Step();
    }

    public IReadOnlyList<double3> GetPoints() {
        var points = new List<double3>(
            _backward.Points.Count + 1 + _forward.Points.Count);

        for(int i = _backward.Points.Count - 1; i >= 0; i--)
            points.Add(_backward.Points[i]);

        points.Add(_seed);
        points.AddRange(_forward.Points);
        return points;
    }

    public FieldLine? Build() {
        IReadOnlyList<double3> points = GetPoints();
        if(points.Count < MinimumPointCount)
            return null;

        return new FieldLine(Decimate(points, _maxVisualPoints));
    }

    private static IReadOnlyList<double3> Decimate(
        IReadOnlyList<double3> points,
        int maxPoints) {

        if(points.Count <= maxPoints)
            return points;

        var result = new List<double3>(maxPoints);
        double stride = (points.Count - 1) / (double)(maxPoints - 1);

        for(int i = 0; i < maxPoints; i++) {
            int sourceIndex = Math.Clamp(
                (int)Math.Round(i * stride),
                0,
                points.Count - 1);
            result.Add(points[sourceIndex]);
        }

        return result;
    }
}

internal sealed class FieldTracer {
    private readonly IFieldModel _model;
    private readonly double _direction;
    private readonly FieldTraceOptions _options;

    private double3 _position;
    private int _steps;

    public FieldTracer(
        IFieldModel model,
        double3 start,
        double direction,
        FieldTraceOptions options) {

        _model = model;
        _position = start;
        _direction = direction;
        _options = options;
    }

    public List<double3> Points { get; } = new();
    public bool IsAlive { get; private set; } = true;

    public void Step() {
        if(!IsAlive)
            return;

        if(_steps >= _options.MaxStepsPerDirection) {
            IsAlive = false;
            return;
        }

        double3 next = Integrate(_position);
        double radius = VectorMath.Length(next);
        double movement = VectorMath.Length(next - _position);

        if(!double.IsFinite(radius) ||
           radius < _options.MinRadius ||
           radius > _options.MaxRadius ||
           movement < _options.StepLength * 1e-6 ||
           movement > _options.StepLength * 2.0) {

            IsAlive = false;
            return;
        }

        Points.Add(next);
        _position = next;
        _steps++;
    }

    private double3 Integrate(double3 position) {
        double h = _options.StepLength;
        double3 k1 = FieldDirection(position);
        double3 k2 = FieldDirection(position + k1 * (h * 0.5));
        double3 k3 = FieldDirection(position + k2 * (h * 0.5));
        double3 k4 = FieldDirection(position + k3 * h);

        return position + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
    }

    private double3 FieldDirection(double3 position) {
        double3 field = _model.EvaluateField(position);
        double length = VectorMath.Length(field);

        return length < 1e-12
            ? new double3(0.0, 0.0, 0.0)
            : _direction * field / length;
    }

}
