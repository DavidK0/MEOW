using KSA;

namespace MEOW;

public static class CelestialBodyResolver {
    public static IParentBody? FindParentBody(string name) {
        CelestialSystem? system = Universe.CurrentSystem;
        if(system == null)
            return null;

        foreach(Astronomical body in system.All.AsSpan()) {
            if (body is IParentBody parentBody && parentBody.Id == name)
                return parentBody;
        }

        return null;
    }
}
