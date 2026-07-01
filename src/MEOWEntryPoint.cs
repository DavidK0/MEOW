using Brutal.ImGuiApi;
using Brutal.Numerics;
using HarmonyLib;
using KSA;
using ModMenu;
using StarMap.API;

namespace MEOW;

[StarMapMod]
public class MEOWEntryPoint {
    private static Harmony? _harmony;

    [StarMapAllModsLoaded]
    public static void OnFullyLoaded() {
        _harmony = new Harmony("dejvid.meow");

        MEOWSettingsStore.Init();
        MEOWSettingsStore.Load();
        SaveLoadObserver.ApplyPatches(_harmony);

        MEOWRenderer.Init();

        BodyOverlayProfile earthProfile = MagneticFieldModelFactory.CreateEarthProfile();

        BodyOverlayProfileRegistry.Register(earthProfile);
    }

    [ModMenuEntry("MEOW")]
    public static void DrawMenu() {
        ImGui.Text("Planetary Magnetic Fields");
        ImGui.Checkbox("Field Lines", ref MEOWSettingsStore.Current.ShowFieldLines);

        // Not implemented yet, but leaving this here for future reference
        //ImGui.Checkbox("Radiation Belts", ref MEOWSettingsStore.Current.ShowRadiationBelts);
        //ImGui.Separator();
        //
        //ImGui.Text("Magnetospheres");
        //ImGui.Checkbox("Magnetopause", ref MEOWSettingsStore.Current.ShowMagnetopause);
        //ImGui.Checkbox("Magnetotail", ref MEOWSettingsStore.Current.ShowMagnetotail);
        //ImGui.Checkbox("Magnetic Cusps", ref MEOWSettingsStore.Current.ShowMagneticCusps);
        //
        //ImGui.Separator();
        //
        //ImGui.Text("Solar Wind Interaction");
        //ImGui.Checkbox("Bow Shock", ref MEOWSettingsStore.Current.ShowBowShock);
        //ImGui.Checkbox("Magnetosheath", ref MEOWSettingsStore.Current.ShowMagnetosheath);

        // Later I will also want to add, hadley/farrel/polar cells, ocean currents, lagrange point, and more
    }

    [StarMapAfterGui]
    public static void OnAfterUi(double dt) {
        // get the position of Earth like this
        double3 earthPosition = Program.ControlledVehicle.Orbit.Parent.GetPositionEcl();

        MEOWRenderer.Draw();
    }
}
