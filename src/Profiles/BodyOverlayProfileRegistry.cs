namespace MEOW;

public sealed class BodyOverlayProfile {
    public required string BodyName { get; init; }

    public IFieldModel? MagneticFieldModel { get; init; }
    public bool HasGlobalMagneticField => MagneticFieldModel != null;
    public IFieldModel? GravitationalFieldModel { get; init; }
    public bool HasGlobalGravitationalField => GravitationalFieldModel != null;
    public IGsmTransform? GsmTransform { get; init; }
}

public static class BodyOverlayProfileRegistry {
    private static readonly Dictionary<string, BodyOverlayProfile> Profiles = new();

    public static void Register(BodyOverlayProfile profile) {
        Profiles[profile.BodyName] = profile;
    }

    public static BodyOverlayProfile? Get(string bodyName) {
        return Profiles.TryGetValue(bodyName, out var profile) ? profile : null;
    }
}
