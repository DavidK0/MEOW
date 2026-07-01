using Brutal.Numerics;

namespace MEOW;

public sealed class MEOWSettings {
    public bool ShowFieldLines = true;
    public bool ShowGridLines;
    public bool ShowRadiationBelts = true;
    public bool ShowMagnetopause = true;
    public bool ShowMagnetotail = true;
    public bool ShowMagneticCusps = true;
    public bool ShowBowShock = true;
    public bool ShowMagnetosheath = true;

    public MEOWSettings Clone() {
        return new MEOWSettings {
            ShowFieldLines = ShowFieldLines,
            ShowGridLines = ShowGridLines,
            ShowRadiationBelts = ShowRadiationBelts,
            ShowMagnetopause = ShowMagnetopause,
            ShowMagnetotail = ShowMagnetotail,
            ShowMagneticCusps = ShowMagneticCusps,
            ShowBowShock = ShowBowShock,
            ShowMagnetosheath = ShowMagnetosheath
        };
    }
}
internal static class MEOWSettingsStore {
    private static SaveScopedSettingsStore<MEOWSettings>? _store;
    private static MEOWSettings _current = new();

    public static MEOWSettings Current {
        get {
            EnsureInitialized();
            return _current;
        }
    }

    public static void LoadForSave(string saveId) {
        EnsureInitialized();
        
        if(string.IsNullOrEmpty(saveId)) {
            _current = new MEOWSettings();
            return;
        }
        
        _store.Load();
        _current = _store.GetCurrent(saveId).Clone();
    }

    public static void SaveForSave(string saveId) {
        EnsureInitialized();

        if(string.IsNullOrEmpty(saveId))
            return;

        _store.Set(saveId, _current.Clone());
        _store.Save(saveId);
    }

    public static void SetCurrentFromDefaults() {
        _current = new MEOWSettings();
    }

    private static void EnsureInitialized() {
        if(_store == null)
            throw new InvalidOperationException("StellariumCatalogSettingsStore.Init() must be called before use.");
    }

    public static void Init() {
        string userDocs = Environment.GetFolderPath(Environment.SpecialFolder.MyDocuments);

        string savesDir = Path.Combine(
            userDocs,
            "My Games",
            "Kitten Space Agency",
            "saves");

        _store = new SaveScopedSettingsStore<MEOWSettings>(
            savesDir,
            "StellariumCatalog_settings.toml",
            () => new MEOWSettings(),
            StellariumCatalogSettingsToml.Read,
            StellariumCatalogSettingsToml.Write);
    }

    public static void Load() {
        EnsureInitialized();
        _store.Load();
    }

    public static void Save() {
        EnsureInitialized();
        _store.Save();
    }
}

internal static class StellariumCatalogSettingsToml {
    public static MEOWSettings Read(SettingsBlock block) {
        var s = new MEOWSettings();

        return s;
    }

    public static void Write(
        SettingsBlockWriter writer,
        string saveId,
        MEOWSettings s) {

        writer.EndBlock();
    }
}
