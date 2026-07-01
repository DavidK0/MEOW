using Brutal.ImGuiApi;
using KSA;

namespace MEOW;

public interface IOverlay {
    bool IsEnabled(MEOWSettings settings);
    void Update(BodyOverlayContext context, MEOWSettings settings, double dt);
    void Draw(ImDrawListPtr draw_list, Camera camera, BodyOverlayContext context, MEOWSettings settings);
}