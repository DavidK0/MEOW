using Brutal.ImGuiApi;
using Brutal.Numerics;
using KSA;

namespace MEOW;

public unsafe static class MEOWRenderer {
    private static readonly OverlayManager _overlayManager = new();

    public static void Init() {
        _overlayManager.Register(new FieldLineOverlay(
            profile => profile.MagneticFieldModel));
    }

    public static void Draw() {
        var vehicle = Program.ControlledVehicle;
        var camera = Program.GetMainCamera();

        if(vehicle == null || camera == null)
            return;

        var parentBody = vehicle.Orbit?.Parent;
        if(parentBody == null)
            return;

        ImGuiViewport* viewport = ImGui.GetMainViewport();
        ImDrawListPtr? draw_list = CreateWindow(viewport);
        if(draw_list == null)
            return;

        BodyOverlayProfile? profile = BodyOverlayProfileRegistry.Get("Earth");

        if(profile == null)
            return;

        var cciToCce = parentBody.GetCci2Cce();
        var cceToCci = cciToCce.Inverse();

        var context = new BodyOverlayContext {
            Time = Universe.GetElapsedSimTime(),
            BodyRadius = parentBody.MeanRadius,

            BodyToWorld = pCci => parentBody.GetPositionEcl() + cciToCce * pCci,
            WorldToBody = pEcl => cceToCci * (pEcl - parentBody.GetPositionEcl()),

            Profile = profile
        };

        _overlayManager.Update(context, MEOWSettingsStore.Current, 1d/60d);
        _overlayManager.Draw(draw_list.Value, camera, context, MEOWSettingsStore.Current);

        ImGui.End();
    }

    public static ImDrawListPtr? CreateWindow(ImGuiViewport* viewport) {
        float2 window_size = viewport->Size;
        ImGui.SetNextWindowPos(viewport->Pos);
        ImGui.SetNextWindowSize(window_size);
        ImGui.SetNextWindowViewport(viewport->ID);
        ImGuiWindowFlags flags =
            ImGuiWindowFlags.NoTitleBar |
            ImGuiWindowFlags.NoResize |
            ImGuiWindowFlags.NoMove |
            ImGuiWindowFlags.NoCollapse |
            ImGuiWindowFlags.NoBringToFrontOnFocus |
            ImGuiWindowFlags.NoBackground |
            ImGuiWindowFlags.NoFocusOnAppearing |
            ImGuiWindowFlags.NoInputs |
            ImGuiWindowFlags.NoNavFocus;
        if(!ImGui.Begin("MEOWWindow", flags)) {
            ImGui.End();
            return null;
        }
        return ImGui.GetWindowDrawList();
    }
}
