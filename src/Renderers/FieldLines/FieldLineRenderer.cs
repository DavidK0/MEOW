using Brutal.ImGuiApi;
using Brutal.Numerics;
using KSA;

namespace MEOW;

public sealed class FieldLineOverlay : IOverlay {
    private const int SeedLatticePointsPerAxis = 4;
    private const int FieldLineCount =
        SeedLatticePointsPerAxis *
        SeedLatticePointsPerAxis *
        SeedLatticePointsPerAxis;

    private const double SeedRadiusMultiplier = 4.0;
    private const double MinTraceRadiusMultiplier = 1.01;
    private const double MaxTraceRadiusMultiplier = 90.0;
    private const double TraceStepRadiusFraction = 0.04;

    private const int MaxTraceStepsPerDirection = 2500;
    private const int MaxTraceStepsPerFrame = 256;
    private const int MaxVisualPointsPerLine = 512;

    private const byte BaseAlpha = 180;
    private readonly List<FieldLineSlot> _lineSlots = new();

    private FieldLineCacheKey? _cachedKey;
    private double _seedBodyRadius = -1.0;
    private long _refreshGeneration;
    private int _nextSlotIndex;
    private bool _initialSetPublished;

    public bool IsEnabled(MEOWSettings settings) => settings.ShowFieldLines;

    public void Update(BodyOverlayContext context, MEOWSettings settings, double dt) {
        IFieldModel model = context.FieldModel;

        model.BeginSimulationStep();

        EnsureLineSlots(context.BodyRadius);

        FieldLineCacheKey cacheKey = FieldLineCacheKey.Create(context, model);
        if(_cachedKey != cacheKey) {
            _cachedKey = cacheKey;
            RequestRefresh();
        }

        FieldTraceOptions options = CreateTraceOptions(context);
        AdvanceTracing(model, options);
        PublishInitialSetWhenReady();
    }

    public void Draw(
        ImDrawListPtr drawList,
        Camera camera,
        BodyOverlayContext context,
        MEOWSettings settings) {

        foreach(FieldLineSlot slot in _lineSlots) {
            FieldLine? line = slot.CurrentLine;
            if(line == null)
                continue;

            DrawPolyline(
                drawList,
                camera,
                context,
                line.PointsBodyFrame,
                BaseAlpha);
        }

    }

    private static FieldTraceOptions CreateTraceOptions(BodyOverlayContext context) {
        double bodyRadius = context.BodyRadius;
        return new FieldTraceOptions(
            StepLength: bodyRadius * TraceStepRadiusFraction,
            MinRadius: bodyRadius * MinTraceRadiusMultiplier,
            MaxRadius: bodyRadius * MaxTraceRadiusMultiplier,
            MaxStepsPerDirection: MaxTraceStepsPerDirection,
            MaxVisualPoints: MaxVisualPointsPerLine,
            Domain: context.Domain);
    }

    private void EnsureLineSlots(double bodyRadius) {
        if(_lineSlots.Count == FieldLineCount &&
           _seedBodyRadius == bodyRadius)
            return;

        _lineSlots.Clear();
        foreach(double3 seed in FieldLineSeedGenerator.LatticeAroundBody(
            bodyRadius,
            SeedRadiusMultiplier,
            SeedLatticePointsPerAxis)) {

            _lineSlots.Add(new FieldLineSlot(seed));
        }

        _seedBodyRadius = bodyRadius;
        _initialSetPublished = false;
        _nextSlotIndex = 0;
        _cachedKey = null;
    }

    private void RequestRefresh() {
        _refreshGeneration++;

        foreach(FieldLineSlot slot in _lineSlots)
            slot.DesiredGeneration = _refreshGeneration;
    }

    private void AdvanceTracing(IFieldModel model, FieldTraceOptions options) {
        for(int step = 0; step < MaxTraceStepsPerFrame; step++) {
            FieldLineSlot? slot = FindNextSlotNeedingWork();
            if(slot == null)
                return;

            if(slot.Builder == null) {
                slot.Builder = new FieldLineBuilder(model, slot.Seed, options);
                slot.BuilderGeneration = slot.DesiredGeneration;
            }

            slot.Builder.Step();
            if(slot.Builder.IsAlive)
                continue;

            FieldLine? replacement = slot.Builder.Build();
            slot.Builder = null;

            if(replacement == null)
                continue;

            if(_initialSetPublished)
                slot.CurrentLine = replacement;
            else
                slot.PendingInitialLine = replacement;

            slot.CompletedGeneration = slot.BuilderGeneration;
        }
    }

    private FieldLineSlot? FindNextSlotNeedingWork() {
        for(int offset = 0; offset < _lineSlots.Count; offset++) {
            int index = (_nextSlotIndex + offset) % _lineSlots.Count;
            FieldLineSlot slot = _lineSlots[index];

            if(slot.Builder == null &&
               slot.CompletedGeneration >= slot.DesiredGeneration)
                continue;

            _nextSlotIndex = (index + 1) % _lineSlots.Count;
            return slot;
        }

        return null;
    }

    private void PublishInitialSetWhenReady() {
        if(_initialSetPublished ||
           _lineSlots.Any(slot => slot.PendingInitialLine == null))
            return;

        foreach(FieldLineSlot slot in _lineSlots) {
            slot.CurrentLine = slot.PendingInitialLine;
            slot.PendingInitialLine = null;
        }

        _initialSetPublished = true;
    }

    private static void DrawPolyline(
        ImDrawListPtr drawList,
        Camera camera,
        BodyOverlayContext context,
        IReadOnlyList<double3> points,
        byte alpha) {

        var color = new ImColor8(255, 255, 255, alpha);

        for(int i = 0; i < points.Count - 1; i++) {
            double3 start = context.BodyToWorld(points[i]);
            double3 end = context.BodyToWorld(points[i + 1]);

            ImDrawListExtensions.AddLine(
                drawList,
                camera.EclToScreen(start),
                camera.EclToScreen(end),
                color,
                2f);
        }
    }

    private readonly record struct FieldLineCacheKey(
        double BodyRadius,
        IFieldModel Model) {

        public static FieldLineCacheKey Create(
            BodyOverlayContext context,
            IFieldModel model) =>
            new(context.BodyRadius, model);
    }

    private sealed class FieldLineSlot {
        public FieldLineSlot(double3 seed) {
            Seed = seed;
        }

        public double3 Seed { get; }
        public FieldLine? CurrentLine { get; set; }
        public FieldLine? PendingInitialLine { get; set; }
        public FieldLineBuilder? Builder { get; set; }
        public long BuilderGeneration { get; set; }
        public long CompletedGeneration { get; set; }
        public long DesiredGeneration { get; set; }
    }
}
