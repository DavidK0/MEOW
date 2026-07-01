using Brutal.ImGuiApi;
using KSA;

namespace MEOW;

public sealed class OverlayManager {
    private readonly List<IOverlay> _overlays = new();

    public void Register(IOverlay overlay) {
        _overlays.Add(overlay);
    }

    public void Update(BodyOverlayContext context, MEOWSettings settings, double dt) {
        foreach(var overlay in _overlays) {
            if(overlay.IsEnabled(settings))
                overlay.Update(context, settings, dt);
        }
    }

    public void Draw(ImDrawListPtr draw_list, Camera camera, BodyOverlayContext context, MEOWSettings settings) {
        foreach(var overlay in _overlays) {
            if(overlay.IsEnabled(settings))
                overlay.Draw(draw_list, camera, context, settings);
        }
    }
}