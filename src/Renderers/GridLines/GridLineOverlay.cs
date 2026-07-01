using Brutal.ImGuiApi;
using Brutal.Numerics;
using KSA;

namespace MEOW;

/// <summary>
/// Draws a body-frame lattice whose vertices are displaced by a vector field.
/// </summary>
public sealed class GridLineOverlay : IOverlay {
    private const double MinimumRadiusMultiplier = 1.01;
    private const float LineThickness = 1.0f;
    private const byte LineAlpha = 100;

    private readonly Func<BodyOverlayProfile, IFieldModel?> _getFieldModel;
    private readonly int _halfGridPointCount;
    private readonly double _gridSpacingRadiusMultiplier;
    private readonly double _dragStrength;
    private readonly List<GridPoint> _points = new();

    private int GridPointCount => _halfGridPointCount * 2 + 1;

    public GridLineOverlay(
        Func<BodyOverlayProfile, IFieldModel?> getFieldModel,
        int halfGridPointCount = 4,
        double gridSpacingRadiusMultiplier = 1.0,
        double dragStrength = 0.75) {

        ArgumentOutOfRangeException.ThrowIfLessThan(halfGridPointCount, 1);
        ArgumentOutOfRangeException.ThrowIfLessThanOrEqual(
            gridSpacingRadiusMultiplier,
            0.0);
        ArgumentOutOfRangeException.ThrowIfNegative(dragStrength);

        _getFieldModel = getFieldModel;
        _halfGridPointCount = halfGridPointCount;
        _gridSpacingRadiusMultiplier = gridSpacingRadiusMultiplier;
        _dragStrength = dragStrength;
    }

    public bool IsEnabled(MEOWSettings settings) => settings.ShowGridLines;

    public void Update(BodyOverlayContext context, MEOWSettings settings, double dt) {
        IFieldModel? model = _getFieldModel(context.Profile);
        if(model == null) {
            _points.Clear();
            return;
        }

        model.BeginSimulationStep();
        RebuildGrid(model, context.BodyRadius);
    }

    public void Draw(
        ImDrawListPtr drawList,
        Camera camera,
        BodyOverlayContext context,
        MEOWSettings settings) {

        if(_points.Count == 0)
            return;

        var color = new ImColor8(255, 255, 255, LineAlpha);
        int pointCount = GridPointCount;

        for(int x = 0; x < pointCount; x++) {
            for(int y = 0; y < pointCount; y++) {
                for(int z = 0; z < pointCount; z++) {
                    if(x + 1 < pointCount)
                        DrawSegment(drawList, camera, context, color, x, y, z, x + 1, y, z);
                    if(y + 1 < pointCount)
                        DrawSegment(drawList, camera, context, color, x, y, z, x, y + 1, z);
                    if(z + 1 < pointCount)
                        DrawSegment(drawList, camera, context, color, x, y, z, x, y, z + 1);
                }
            }
        }
    }

    private void RebuildGrid(IFieldModel model, double bodyRadius) {
        _points.Clear();

        int pointCount = GridPointCount;
        double spacing = bodyRadius * _gridSpacingRadiusMultiplier;
        double minimumRadius = bodyRadius * MinimumRadiusMultiplier;
        var samples = new List<FieldSample>(pointCount * pointCount * pointCount);
        double maximumFieldMagnitude = 0.0;

        for(int x = -_halfGridPointCount; x <= _halfGridPointCount; x++) {
            for(int y = -_halfGridPointCount; y <= _halfGridPointCount; y++) {
                for(int z = -_halfGridPointCount; z <= _halfGridPointCount; z++) {
                    var position = new double3(x * spacing, y * spacing, z * spacing);
                    bool isVisible = VectorMath.Length(position) >= minimumRadius;
                    double3 field = isVisible ? model.EvaluateField(position) : default;
                    double magnitude = VectorMath.Length(field);

                    if(double.IsFinite(magnitude))
                        maximumFieldMagnitude = Math.Max(maximumFieldMagnitude, magnitude);
                    else
                        field = default;

                    samples.Add(new FieldSample(position, field, isVisible));
                }
            }
        }

        double displacementScale = maximumFieldMagnitude > 0.0
            ? spacing * _dragStrength / maximumFieldMagnitude
            : 0.0;

        foreach(FieldSample sample in samples) {
            _points.Add(new GridPoint(
                sample.Position + sample.Field * displacementScale,
                sample.IsVisible));
        }
    }

    private void DrawSegment(
        ImDrawListPtr drawList,
        Camera camera,
        BodyOverlayContext context,
        ImColor8 color,
        int startX,
        int startY,
        int startZ,
        int endX,
        int endY,
        int endZ) {

        GridPoint start = _points[GetIndex(startX, startY, startZ)];
        GridPoint end = _points[GetIndex(endX, endY, endZ)];
        if(!start.IsVisible || !end.IsVisible)
            return;

        ImDrawListExtensions.AddLine(
            drawList,
            camera.EclToScreen(context.BodyToWorld(start.Position)),
            camera.EclToScreen(context.BodyToWorld(end.Position)),
            color,
            LineThickness);
    }

    private int GetIndex(int x, int y, int z) {
        int pointCount = GridPointCount;
        return (x * pointCount + y) * pointCount + z;
    }

    private readonly record struct FieldSample(
        double3 Position,
        double3 Field,
        bool IsVisible);

    private readonly record struct GridPoint(
        double3 Position,
        bool IsVisible);
}
