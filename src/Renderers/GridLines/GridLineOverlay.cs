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

    private readonly int? _halfGridPointCountOverride;
    private readonly double? _gridSpacingRadiusMultiplierOverride;
    private readonly double _dragStrength;
    private readonly List<GridPoint> _points = new();
    private int _gridPointCount;

    public GridLineOverlay(
        int? halfGridPointCount = null,
        double? gridSpacingRadiusMultiplier = null,
        double dragStrength = 0.75) {

        if(halfGridPointCount.HasValue)
            ArgumentOutOfRangeException.ThrowIfLessThan(
                halfGridPointCount.Value,
                1);
        if(gridSpacingRadiusMultiplier.HasValue)
            ArgumentOutOfRangeException.ThrowIfLessThanOrEqual(
                gridSpacingRadiusMultiplier.Value,
                0.0);
        ArgumentOutOfRangeException.ThrowIfNegative(dragStrength);

        _halfGridPointCountOverride = halfGridPointCount;
        _gridSpacingRadiusMultiplierOverride = gridSpacingRadiusMultiplier;
        _dragStrength = dragStrength;
    }

    public bool IsEnabled(MEOWSettings settings) => settings.ShowGridLines;

    public void Update(BodyOverlayContext context, MEOWSettings settings, double dt) {
        IFieldModel model = context.FieldModel;
        model.BeginSimulationStep();
        RebuildGrid(model, context.Domain);
    }

    public void Draw(
        ImDrawListPtr drawList,
        Camera camera,
        BodyOverlayContext context,
        MEOWSettings settings) {

        if(_points.Count == 0)
            return;

        var color = new ImColor8(255, 255, 255, LineAlpha);
        int pointCount = _gridPointCount;

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

    private void RebuildGrid(IFieldModel model, FieldRenderDomain domain) {
        _points.Clear();

        double spacing = _gridSpacingRadiusMultiplierOverride.HasValue
            ? domain.ReferenceLength * _gridSpacingRadiusMultiplierOverride.Value
            : domain.GridSpacing;
        int halfGridPointCount = _halfGridPointCountOverride ??
            Math.Max(1, (int)Math.Floor(domain.HalfExtent / spacing));
        int pointCount = halfGridPointCount * 2 + 1;
        _gridPointCount = pointCount;
        double minimumRadius = domain.ReferenceLength * MinimumRadiusMultiplier;
        var samples = new List<FieldSample>(pointCount * pointCount * pointCount);
        double maximumFieldMagnitude = 0.0;

        for(int x = -halfGridPointCount; x <= halfGridPointCount; x++) {
            for(int y = -halfGridPointCount; y <= halfGridPointCount; y++) {
                for(int z = -halfGridPointCount; z <= halfGridPointCount; z++) {
                    var position = new double3(x * spacing, y * spacing, z * spacing);
                    bool isVisible = VectorMath.Length(position) >= minimumRadius;
                    double3 field = isVisible
                        ? domain.EclToLocalVector(
                            model.EvaluateField(domain.LocalToEclPoint(position)))
                        : default;
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
        return (x * _gridPointCount + y) * _gridPointCount + z;
    }

    private readonly record struct FieldSample(
        double3 Position,
        double3 Field,
        bool IsVisible);

    private readonly record struct GridPoint(
        double3 Position,
        bool IsVisible);
}
