using Brutal.Numerics;
using KSA;

namespace MEOW;

public sealed class T96Model : IBodyFrameFieldModel {
    private const double EarthRadiusMeters = 6371200.0;

    private readonly IGsmTransform _gsm;

    public double PdynNpa { get; set; } = 2.0;
    public double DstNt { get; set; } = -5.0;
    public double ByImfNt { get; set; } = 0.0;
    public double BzImfNt { get; set; } = 0.0;

    public T96Model(IGsmTransform gsm) {
        _gsm = gsm;
    }

    public void BeginSimulationStep() {
        _gsm.BeginSimulationStep();
    }

    public double3 EvaluateField(double3 positionBodyFrame) {
        double3 posGsmMeters = _gsm.BodyToGsmPosition(positionBodyFrame);

        double xRe = posGsmMeters.X / EarthRadiusMeters;
        double yRe = posGsmMeters.Y / EarthRadiusMeters;
        double zRe = posGsmMeters.Z / EarthRadiusMeters;

        double ps = _gsm.GetDipoleTiltRadians();

        double3 bGsmNt = T96.Evaluate(
            xRe,
            yRe,
            zRe,
            PdynNpa,
            DstNt,
            ByImfNt,
            BzImfNt,
            ps);

        return _gsm.GsmToBodyVector(bGsmNt);
    }
}

internal static class T96 {
    private const double Pi = 3.141592654;

    public static double3 Evaluate(double x, double y, double z, double pdyn, double dst, double byImf, double bzImf, double ps) {
        double pdyn0 = 2.0, eps10 = 3630.7;
        double[] a = OneBased(1.162, 22.344, 18.50, 2.602, 6.903, 5.287, 0.5790, 0.4462, 0.7850);
        double am0 = 70.0, s0 = 1.08, x00 = 5.48, dsig = 0.005;
        double delimfx = 20.0, delimfy = 10.0;

        double depr = 0.8 * dst - 13.0 * Math.Sqrt(pdyn);
        double bt = Math.Sqrt(byImf * byImf + bzImf * bzImf);
        double theta = 0.0;
        if (byImf != 0.0 || bzImf != 0.0) {
            theta = Math.Atan2(byImf, bzImf);
            if (theta <= 0.0) theta += 6.2831853;
        }

        double ct = Math.Cos(theta);
        double st = Math.Sin(theta);
        double eps = 718.5 * Math.Sqrt(pdyn) * bt * Math.Sin(theta / 2.0);
        double facteps = eps / eps10 - 1.0;
        double factpd = Math.Sqrt(pdyn / pdyn0) - 1.0;
        double rcampl = -a[1] * depr;
        double tampl2 = a[2] + a[3] * factpd + a[4] * facteps;
        double tampl3 = a[5] + a[6] * factpd;
        double b1ampl = a[7] + a[8] * facteps;
        double b2ampl = 20.0 * b1ampl;
        double reconn = a[9];
        double xappa = Math.Pow(pdyn / pdyn0, 0.14);
        double xappa3 = xappa * xappa * xappa;
        double ys = y * ct - z * st;
        double zs = z * ct + y * st;
        double factimf = Math.Exp(x / delimfx - Math.Pow(ys / delimfy, 2.0));
        double oimfx = 0.0;
        double oimfy = reconn * byImf * factimf;
        double oimfz = reconn * bzImf * factimf;
        double rimfampl = reconn * bt;
        double xx = x * xappa;
        double yy = y * xappa;
        double zz = z * xappa;

        double x0 = x00 / xappa;
        double am = am0 / xappa;
        double rho2 = y * y + z * z;
        double asq = am * am;
        double xmxm = am + x - x0;
        if (xmxm < 0.0) xmxm = 0.0;
        double axx0 = xmxm * xmxm;
        double aro = asq + rho2;
        double sigma = Math.Sqrt((aro + axx0 + Math.Sqrt(Math.Pow(aro + axx0, 2.0) - 4.0 * asq * axx0)) / (2.0 * asq));

        if (sigma < s0 + dsig) {
            var cf = DipShld(ps, xx, yy, zz);
            var tr = TailRc96(Math.Sin(ps), xx, yy, zz);
            var r1 = Birk1Tot02(ps, xx, yy, zz);
            var r2 = Birk2Tot02(ps, xx, yy, zz);
            var rimfs = Intercon(xx, ys * xappa, zs * xappa);
            double rimfy = rimfs.by * ct + rimfs.bz * st;
            double rimfz = rimfs.bz * ct - rimfs.by * st;
            double fx = cf.bx * xappa3 + rcampl * tr.bxrc + tampl2 * tr.bxt2 + tampl3 * tr.bxt3 + b1ampl * r1.bx + b2ampl * r2.bx + rimfampl * rimfs.bx;
            double fy = cf.by * xappa3 + rcampl * tr.byrc + tampl2 * tr.byt2 + tampl3 * tr.byt3 + b1ampl * r1.by + b2ampl * r2.by + rimfampl * rimfy;
            double fz = cf.bz * xappa3 + rcampl * tr.bzrc + tampl2 * tr.bzt2 + tampl3 * tr.bzt3 + b1ampl * r1.bz + b2ampl * r2.bz + rimfampl * rimfz;
            if (sigma < s0 - dsig) return new double3(fx, fy, fz);

            double fint = 0.5 * (1.0 - (sigma - s0) / dsig);
            double fext = 0.5 * (1.0 + (sigma - s0) / dsig);
            var q = Dipole(ps, x, y, z);
            return new double3((fx + q.bx) * fint + oimfx * fext - q.bx,
                               (fy + q.by) * fint + oimfy * fext - q.by,
                               (fz + q.bz) * fint + oimfz * fext - q.bz);
        }

        var d = Dipole(ps, x, y, z);
        return new double3(oimfx - d.bx, oimfy - d.by, oimfz - d.bz);
    }

    private static (double bx, double by, double bz) DipShld(double ps, double x, double y, double z) {
        double[] a1 = OneBased(.24777, -27.003, -.46815, 7.0637, -1.5918, -.90317E-01, 57.522, 13.757, 2.0100, 10.458, 4.5798, 2.1695);
        double[] a2 = OneBased(-.65385, -18.061, -.40457, -5.0995, 1.2846, .78231E-01, 39.592, 13.291, 1.9970, 10.062, 4.5140, 2.1558);
        double cps = Math.Cos(ps), sps = Math.Sin(ps);
        var h = CylHarm(a1, x, y, z);
        var f = CylHar1(a2, x, y, z);
        return (h.bx * cps + f.bx * sps, h.by * cps + f.by * sps, h.bz * cps + f.bz * sps);
    }

    private static (double bx, double by, double bz) CylHarm(double[] a, double x, double y, double z) {
        double rho = Math.Sqrt(y * y + z * z);
        double sinfi, cosfi;
        if (rho < 1e-8) { sinfi = 1.0; cosfi = 0.0; rho = 1e-8; }
        else { sinfi = z / rho; cosfi = y / rho; }
        double sinfi2 = sinfi * sinfi;
        double si2co2 = sinfi2 - cosfi * cosfi;
        double bx = 0.0, by = 0.0, bz = 0.0;
        for (int i = 1; i <= 3; i++) {
            double dzeta = rho / a[i + 6];
            double xj0 = Bes(dzeta, 0);
            double xj1 = Bes(dzeta, 1);
            double xexp = Math.Exp(x / a[i + 6]);
            bx -= a[i] * xj1 * xexp * sinfi;
            by += a[i] * (2.0 * xj1 / dzeta - xj0) * xexp * sinfi * cosfi;
            bz += a[i] * (xj1 / dzeta * si2co2 - xj0 * sinfi2) * xexp;
        }
        for (int i = 4; i <= 6; i++) {
            double dzeta = rho / a[i + 6];
            double xksi = x / a[i + 6];
            double xj0 = Bes(dzeta, 0);
            double xj1 = Bes(dzeta, 1);
            double xexp = Math.Exp(xksi);
            double brho = (xksi * xj0 - (dzeta * dzeta + xksi - 1.0) * xj1 / dzeta) * xexp * sinfi;
            double bphi = (xj0 + xj1 / dzeta * (xksi - 1.0)) * xexp * cosfi;
            bx += a[i] * (dzeta * xj0 + xksi * xj1) * xexp * sinfi;
            by += a[i] * (brho * cosfi - bphi * sinfi);
            bz += a[i] * (brho * sinfi + bphi * cosfi);
        }
        return (bx, by, bz);
    }

    private static (double bx, double by, double bz) CylHar1(double[] a, double x, double y, double z) {
        double rho = Math.Sqrt(y * y + z * z);
        double sinfi, cosfi;
        if (rho < 1e-10) { sinfi = 1.0; cosfi = 0.0; }
        else { sinfi = z / rho; cosfi = y / rho; }
        double bx = 0.0, by = 0.0, bz = 0.0;
        for (int i = 1; i <= 3; i++) {
            double dzeta = rho / a[i + 6], xksi = x / a[i + 6];
            double xj0 = Bes(dzeta, 0), xj1 = Bes(dzeta, 1), xexp = Math.Exp(xksi);
            double brho = xj1 * xexp;
            bx -= a[i] * xj0 * xexp;
            by += a[i] * brho * cosfi;
            bz += a[i] * brho * sinfi;
        }
        for (int i = 4; i <= 6; i++) {
            double dzeta = rho / a[i + 6], xksi = x / a[i + 6];
            double xj0 = Bes(dzeta, 0), xj1 = Bes(dzeta, 1), xexp = Math.Exp(xksi);
            double brho = (dzeta * xj0 + xksi * xj1) * xexp;
            bx += a[i] * (dzeta * xj1 - xj0 * (xksi + 1.0)) * xexp;
            by += a[i] * brho * cosfi;
            bz += a[i] * brho * sinfi;
        }
        return (bx, by, bz);
    }

    private static double Bes(double x, int k) {
        if (k == 0) return Bes0(x);
        if (k == 1) return Bes1(x);
        if (x == 0.0) return 0.0;
        double g = 2.0 / x;
        if (x > k) {
            int n = 1;
            double xjn = Bes1(x), xjnm1 = Bes0(x), xjnp1;
            while (true) {
                xjnp1 = g * n * xjn - xjnm1;
                n++;
                if (n >= k) return xjnp1;
                xjnm1 = xjn;
                xjn = xjnp1;
            }
        } else {
            int n = 24;
            double xjn = 1.0, xjnp1 = 0.0, sum = 0.0, bes = 0.0;
            while (true) {
                if (n % 2 == 0) sum += xjn;
                double xjnm1 = g * n * xjn - xjnp1;
                n--;
                xjnp1 = xjn;
                xjn = xjnm1;
                if (n == k) bes = xjn;
                if (Math.Abs(xjn) > 1e5) {
                    xjnp1 *= 1e-5;
                    xjn *= 1e-5;
                    sum *= 1e-5;
                    if (n <= k) bes *= 1e-5;
                }
                if (n == 0) {
                    sum = xjn + 2.0 * sum;
                    return bes / sum;
                }
            }
        }
    }

    private static double Bes0(double x) {
        if (Math.Abs(x) < 3.0) {
            double x32 = Math.Pow(x / 3.0, 2.0);
            return 1.0 - x32 * (2.2499997 - x32 * (1.2656208 - x32 * (0.3163866 - x32 * (0.0444479 - x32 * (0.0039444 - x32 * 0.00021)))));
        }
        double xd3 = 3.0 / x;
        double f0 = 0.79788456 - xd3 * (0.00000077 + xd3 * (0.00552740 + xd3 * (0.00009512 - xd3 * (0.00137237 - xd3 * (0.00072805 - xd3 * 0.00014476)))));
        double t0 = x - 0.78539816 - xd3 * (0.04166397 + xd3 * (0.00003954 - xd3 * (0.00262573 - xd3 * (0.00054125 + xd3 * (0.00029333 - xd3 * 0.00013558)))));
        return f0 / Math.Sqrt(x) * Math.Cos(t0);
    }

    private static double Bes1(double x) {
        if (Math.Abs(x) < 3.0) {
            double x32 = Math.Pow(x / 3.0, 2.0);
            double bes1xm1 = 0.5 - x32 * (0.56249985 - x32 * (0.21093573 - x32 * (0.03954289 - x32 * (0.00443319 - x32 * (0.00031761 - x32 * 0.00001109)))));
            return bes1xm1 * x;
        }
        double xd3 = 3.0 / x;
        double f1 = 0.79788456 + xd3 * (0.00000156 + xd3 * (0.01659667 + xd3 * (0.00017105 - xd3 * (0.00249511 - xd3 * (0.00113653 - xd3 * 0.00020033)))));
        double t1 = x - 2.35619449 + xd3 * (0.12499612 + xd3 * (0.0000565 - xd3 * (0.00637879 - xd3 * (0.00074348 + xd3 * (0.00079824 - xd3 * 0.00029166)))));
        return f1 / Math.Sqrt(x) * Math.Cos(t1);
    }

    private static (double bx, double by, double bz) Intercon(double x, double y, double z) {
        double[] a = OneBased(-8.411078731, 5932254.951, -9073284.93, -11.68794634, 6027598.824, -9218378.368, -6.508798398, -11824.42793, 18015.66212, 7.99754043, 13.9669886, 90.24475036, 16.75728834, 1015.645781, 1553.493216);
        double[] p = OneBased(a[10], a[11], a[12]);
        double[] r = OneBased(a[13], a[14], a[15]);
        double[] rp = OneBased(1.0 / p[1], 1.0 / p[2], 1.0 / p[3]);
        double[] rr = OneBased(1.0 / r[1], 1.0 / r[2], 1.0 / r[3]);
        int l = 0;
        double bx = 0.0, by = 0.0, bz = 0.0;
        for (int i = 1; i <= 3; i++) {
            double cypi = Math.Cos(y * rp[i]);
            double sypi = Math.Sin(y * rp[i]);
            for (int k = 1; k <= 3; k++) {
                double szrk = Math.Sin(z * rr[k]);
                double czrk = Math.Cos(z * rr[k]);
                double sqpr = Math.Sqrt(rp[i] * rp[i] + rr[k] * rr[k]);
                double epr = Math.Exp(x * sqpr);
                double hx = -sqpr * epr * cypi * szrk;
                double hy = rp[i] * epr * sypi * szrk;
                double hz = -rr[k] * epr * cypi * czrk;
                l++;
                bx += a[l] * hx;
                by += a[l] * hy;
                bz += a[l] * hz;
            }
        }
        return (bx, by, bz);
    }

    private sealed class WarpState {
        public double Cpss, Spss, Dpsrr, Rps, Warp, D, Xs, Zs, Dxsx, Dxsy, Dxsz, Dzsx, Dzsy, Dzsz, Dzetas, Ddzetadx, Ddzetady, Ddzetadz, Zsww;
    }

    private static (double bxrc, double byrc, double bzrc, double bxt2, double byt2, double bzt2, double bxt3, double byt3, double bzt3) TailRc96(double sps, double x, double y, double z) {
        double[] arc = OneBased(-3.087699646, 3.516259114, 18.81380577, -13.95772338, -5.497076303, 0.1712890838, 2.392629189, -2.728020808, -14.79349936, 11.08738083, 4.388174084, 0.2492163197E-01, 0.7030375685, -.7966023165, -3.835041334, 2.642228681, -0.2405352424, -0.7297705678, -0.3680255045, 0.1333685557, 2.795140897, -1.078379954, 0.8014028630, 0.1245825565, 0.6149982835, -0.2207267314, -4.424578723, 1.730471572, -1.716313926, -0.2306302941, -0.2450342688, 0.8617173961E-01, 1.54697858, -0.6569391113, -0.6537525353, 0.2079417515, 12.75434981, 11.37659788, 636.4346279, 1.752483754, 3.604231143, 12.83078674, 7.412066636, 9.434625736, 676.7557193, 1.701162737, 3.580307144, 14.64298662);
        double[] atail2 = OneBased(.8747515218, -.9116821411, 2.209365387, -2.159059518, -7.059828867, 5.924671028, -1.916935691, 1.996707344, -3.877101873, 3.947666061, 11.38715899, -8.343210833, 1.194109867, -1.244316975, 3.73895491, -4.406522465, -20.66884863, 3.020952989, .2189908481, -.09942543549, -.927225562, .1555224669, .6994137909, -.08111721003, -.7565493881, .4686588792, 4.266058082, -.3717470262, -3.920787807, .02298569870, .7039506341, -.5498352719, -6.675140817, .8279283559, -2.234773608, -1.622656137, 5.187666221, 6.802472048, 39.13543412, 2.784722096, 6.979576616, 25.71716760, 4.495005873, 8.068408272, 93.47887103, 4.158030104, 9.313492566, 57.18240483);
        double[] atail3 = OneBased(-19091.95061, -3011.613928, 20582.16203, 4242.918430, -2377.091102, -1504.820043, 19884.04650, 2725.150544, -21389.04845, -3990.475093, 2401.610097, 1548.171792, -946.5493963, 490.1528941, 986.9156625, -489.3265930, -67.99278499, 8.711175710, -45.15734260, -10.76106500, 210.7927312, 11.41764141, -178.0262808, .7558830028, 339.3806753, 9.904695974, 69.50583193, -118.0271581, 22.85935896, 45.91014857, -425.6607164, 15.47250738, 118.2988915, 65.58594397, -201.4478068, -14.57062940, 19.69877970, 20.30095680, 86.45407420, 22.50403727, 23.41617329, 48.48140573, 24.61031329, 123.5395974, 223.5367692, 39.50824342, 65.83385762, 266.2948657);
        double rh = 9.0, dr = 4.0, g = 10.0, d0 = 2.0, deltady = 10.0;
        var w = new WarpState();
        double dr2 = dr * dr;
        double c11 = Math.Sqrt((1.0 + rh) * (1.0 + rh) + dr2);
        double c12 = Math.Sqrt((1.0 - rh) * (1.0 - rh) + dr2);
        double c1 = c11 - c12;
        double spsc1 = sps / c1;
        w.Rps = 0.5 * (c11 + c12) * sps;
        double r = Math.Sqrt(x * x + y * y + z * z);
        double sq1 = Math.Sqrt((r + rh) * (r + rh) + dr2);
        double sq2 = Math.Sqrt((r - rh) * (r - rh) + dr2);
        double c = sq1 - sq2;
        double cs = (r + rh) / sq1 - (r - rh) / sq2;
        w.Spss = spsc1 / r * c;
        w.Cpss = Math.Sqrt(1.0 - w.Spss * w.Spss);
        w.Dpsrr = sps / (r * r) * (cs * r - c) / Math.Sqrt(Math.Pow(r * c1, 2.0) - Math.Pow(c * sps, 2.0));
        double wfac = y / (Math.Pow(y, 4.0) + 1.0e4);
        double ww = wfac * y * y * y;
        double ws = 4.0e4 * y * wfac * wfac;
        w.Warp = g * sps * ww;
        w.Xs = x * w.Cpss - z * w.Spss;
        w.Zsww = z * w.Cpss + x * w.Spss;
        w.Zs = w.Zsww + w.Warp;
        w.Dxsx = w.Cpss - x * w.Zsww * w.Dpsrr;
        w.Dxsy = -y * w.Zsww * w.Dpsrr;
        w.Dxsz = -w.Spss - z * w.Zsww * w.Dpsrr;
        w.Dzsx = w.Spss + x * w.Xs * w.Dpsrr;
        w.Dzsy = w.Xs * y * w.Dpsrr + g * sps * ws;
        w.Dzsz = w.Cpss + w.Xs * z * w.Dpsrr;
        w.D = d0 + deltady * Math.Pow(y / 20.0, 2.0);
        double dddy = deltady * y * 0.005;
        w.Dzetas = Math.Sqrt(w.Zs * w.Zs + w.D * w.D);
        w.Ddzetadx = w.Zs * w.Dzsx / w.Dzetas;
        w.Ddzetady = (w.Zs * w.Dzsy + w.D * dddy) / w.Dzetas;
        w.Ddzetadz = w.Zs * w.Dzsz / w.Dzetas;

        var sh = ShlCar3x3(arc, x, y, z, sps);
        var rc = RingCurr96(x, y, z, w);
        var sh2 = ShlCar3x3(atail2, x, y, z, sps);
        var td = TailDisk(x, y, z, w);
        var sh3 = ShlCar3x3(atail3, x, y, z, sps);
        var t87 = Tail87(x, z, w);
        return (sh.bx + rc.bx, sh.by + rc.by, sh.bz + rc.bz,
                sh2.bx + td.bx, sh2.by + td.by, sh2.bz + td.bz,
                sh3.bx + t87.bx, sh3.by, sh3.bz + t87.bz);
    }

    private static (double bx, double by, double bz) ShlCar3x3(double[] a, double x, double y, double z, double sps) {
        double cps = Math.Sqrt(1.0 - sps * sps);
        double s3ps = 4.0 * cps * cps - 1.0;
        double hx = 0.0, hy = 0.0, hz = 0.0;
        int l = 0;
        for (int m = 1; m <= 2; m++)
        for (int i = 1; i <= 3; i++) {
            double p = a[36 + i], q = a[42 + i];
            double cypi = Math.Cos(y / p), cyqi = Math.Cos(y / q);
            double sypi = Math.Sin(y / p), syqi = Math.Sin(y / q);
            for (int k = 1; k <= 3; k++) {
                double r = a[39 + k], s = a[45 + k];
                double szrk = Math.Sin(z / r), czsk = Math.Cos(z / s), czrk = Math.Cos(z / r), szsk = Math.Sin(z / s);
                double sqpr = Math.Sqrt(1.0 / (p * p) + 1.0 / (r * r));
                double sqqs = Math.Sqrt(1.0 / (q * q) + 1.0 / (s * s));
                double epr = Math.Exp(x * sqpr), eqs = Math.Exp(x * sqqs);
                double dx = 0.0, dy = 0.0, dz = 0.0;
                for (int n = 1; n <= 2; n++) {
                    l++;
                    if (m == 1) {
                        if (n == 1) {
                            dx = -sqpr * epr * cypi * szrk;
                            dy = epr / p * sypi * szrk;
                            dz = -epr / r * cypi * czrk;
                        } else { dx *= cps; dy *= cps; dz *= cps; }
                    } else {
                        if (n == 1) {
                            dx = -sps * sqqs * eqs * cyqi * czsk;
                            dy = sps * eqs / q * syqi * czsk;
                            dz = sps * eqs / s * cyqi * szsk;
                        } else { dx *= s3ps; dy *= s3ps; dz *= s3ps; }
                    }
                    hx += a[l] * dx;
                    hy += a[l] * dy;
                    hz += a[l] * dz;
                }
            }
        }
        return (hx, hy, hz);
    }

    private static (double bx, double by, double bz) RingCurr96(double x, double y, double z, WarpState w) {
        double[] f = OneBased(569.895366, -1603.386993);
        double[] beta = OneBased(2.722188, 3.766875);
        double d0 = 2.0, deltadx = 0.0, xd = 0.0, xldx = 4.0;
        double dzsy = w.Xs * y * w.Dpsrr;
        double xxd = x - xd;
        double fdx = 0.5 * (1.0 + xxd / Math.Sqrt(xxd * xxd + xldx * xldx));
        double dddx = deltadx * 0.5 * xldx * xldx / Math.Pow(Math.Sqrt(xxd * xxd + xldx * xldx), 3.0);
        double d = d0 + deltadx * fdx;
        double dzetas = Math.Sqrt(w.Zsww * w.Zsww + d * d);
        double rhos = Math.Sqrt(w.Xs * w.Xs + y * y);
        double ddzetadx = (w.Zsww * w.Dzsx + d * dddx) / dzetas;
        double ddzetady = w.Zsww * dzsy / dzetas;
        double ddzetadz = w.Zsww * w.Dzsz / dzetas;
        double drhosdx, drhosdy, drhosdz;
        if (rhos < 1e-5) { drhosdx = 0.0; drhosdy = Sign1(y); drhosdz = 0.0; }
        else { drhosdx = w.Xs * w.Dxsx / rhos; drhosdy = (w.Xs * w.Dxsy + y) / rhos; drhosdz = w.Xs * w.Dxsz / rhos; }
        return DiskLoopField(2, f, beta, x, w.Xs, y, z, dzetas, rhos, ddzetadx, ddzetady, ddzetadz, drhosdx, drhosdy, drhosdz, w, 0.0, w.Zsww);
    }

    private static (double bx, double by, double bz) TailDisk(double x, double y, double z, WarpState w) {
        double[] f = OneBased(-745796.7338, 1176470.141, -444610.529, -57508.01028);
        double[] beta = OneBased(7.9250000, 8.0850000, 8.4712500, 27.89500);
        double xshift = 4.5;
        double rhos = Math.Sqrt(Math.Pow(w.Xs - xshift, 2.0) + y * y);
        double drhosdx, drhosdy, drhosdz;
        if (rhos < 1e-5) { drhosdx = 0.0; drhosdy = Sign1(y); drhosdz = 0.0; }
        else { drhosdx = (w.Xs - xshift) * w.Dxsx / rhos; drhosdy = ((w.Xs - xshift) * w.Dxsy + y) / rhos; drhosdz = (w.Xs - xshift) * w.Dxsz / rhos; }
        return DiskLoopField(4, f, beta, x, w.Xs, y, z, w.Dzetas, rhos, w.Ddzetadx, w.Ddzetady, w.Ddzetadz, drhosdx, drhosdy, drhosdz, w, xshift, w.Zsww);
    }

    private static (double bx, double by, double bz) DiskLoopField(int nmax, double[] f, double[] beta, double x, double xs, double y, double z, double dzetas, double rhos, double ddzetadx, double ddzetady, double ddzetadz, double drhosdx, double drhosdy, double drhosdz, WarpState w, double xshift, double zsww) {
        double bx = 0.0, by = 0.0, bz = 0.0;
        for (int i = 1; i <= nmax; i++) {
            double bi = beta[i];
            double s1 = Math.Sqrt(Math.Pow(dzetas + bi, 2.0) + Math.Pow(rhos + bi, 2.0));
            double s2 = Math.Sqrt(Math.Pow(dzetas + bi, 2.0) + Math.Pow(rhos - bi, 2.0));
            double ds1ddz = (dzetas + bi) / s1, ds2ddz = (dzetas + bi) / s2;
            double ds1drhos = (rhos + bi) / s1, ds2drhos = (rhos - bi) / s2;
            double ds1dx = ds1ddz * ddzetadx + ds1drhos * drhosdx;
            double ds1dy = ds1ddz * ddzetady + ds1drhos * drhosdy;
            double ds1dz = ds1ddz * ddzetadz + ds1drhos * drhosdz;
            double ds2dx = ds2ddz * ddzetadx + ds2drhos * drhosdx;
            double ds2dy = ds2ddz * ddzetady + ds2drhos * drhosdy;
            double ds2dz = ds2ddz * ddzetadz + ds2drhos * drhosdz;
            double s1ts2 = s1 * s2, s1ps2 = s1 + s2, s1ps2sq = s1ps2 * s1ps2;
            double fac1 = Math.Sqrt(s1ps2sq - Math.Pow(2.0 * bi, 2.0));
            double ass = fac1 / (s1ts2 * s1ps2sq);
            double term1 = 1.0 / (s1ts2 * s1ps2 * fac1);
            double fac2 = ass / s1ps2sq;
            double dasds1 = term1 - fac2 / s1 * (s2 * s2 + s1 * (3.0 * s1 + 4.0 * s2));
            double dasds2 = term1 - fac2 / s2 * (s1 * s1 + s2 * (3.0 * s2 + 4.0 * s1));
            double dasdx = dasds1 * ds1dx + dasds2 * ds2dx;
            double dasdy = dasds1 * ds1dy + dasds2 * ds2dy;
            double dasdz = dasds1 * ds1dz + dasds2 * ds2dz;
            bx += f[i] * ((2.0 * ass + y * dasdy) * w.Spss - (xs - xshift) * dasdz + ass * w.Dpsrr * (y * y * w.Cpss + z * zsww));
            by -= f[i] * y * (ass * w.Dpsrr * xs + dasdz * w.Cpss + dasdx * w.Spss);
            bz += f[i] * ((2.0 * ass + y * dasdy) * w.Cpss + (xs - xshift) * dasdx - ass * w.Dpsrr * (x * zsww + y * y * w.Spss));
        }
        return (bx, by, bz);
    }

    private static (double bx, double bz) Tail87(double x, double z, WarpState w) {
        double dd = 3.0, hpi = 1.5707963, rt = 40.0, xn = -10.0, x1 = -1.261, x2 = -0.663, b0 = 0.391734, b1 = 5.89715, b2 = 24.6833, xn21 = 76.37, xnr = -0.1071, adln = 0.13238005;
        double zs = z - w.Rps + w.Warp, zp = z - rt, zm = z + rt;
        double xnx = xn - x, xnx2 = xnx * xnx, xc1 = x - x1, xc2 = x - x2, xc22 = xc2 * xc2, xr2 = xc2 * xnr, xc12 = xc1 * xc1, d2 = dd * dd;
        double b20 = zs * zs + d2, b2p = zp * zp + d2, b2m = zm * zm + d2;
        double b = Math.Sqrt(b20), bp = Math.Sqrt(b2p), bm = Math.Sqrt(b2m);
        double xa1 = xc12 + b20, xap1 = xc12 + b2p, xam1 = xc12 + b2m;
        double xa2 = 1.0 / (xc22 + b20), xap2 = 1.0 / (xc22 + b2p), xam2 = 1.0 / (xc22 + b2m);
        double xna = xnx2 + b20, xnap = xnx2 + b2p, xnam = xnx2 + b2m;
        double f = b20 - xc22, fp = b2p - xc22, fm = b2m - xc22;
        double xln1 = Math.Log(xn21 / xna), xlnp1 = Math.Log(xn21 / xnap), xlnm1 = Math.Log(xn21 / xnam);
        double xln2 = xln1 + adln, xlnp2 = xlnp1 + adln, xlnm2 = xlnm1 + adln;
        double aln = 0.25 * (xlnp1 + xlnm1 - 2.0 * xln1);
        double s0 = (Math.Atan(xnx / b) + hpi) / b, s0p = (Math.Atan(xnx / bp) + hpi) / bp, s0m = (Math.Atan(xnx / bm) + hpi) / bm;
        double s1 = (xln1 * .5 + xc1 * s0) / xa1, s1p = (xlnp1 * .5 + xc1 * s0p) / xap1, s1m = (xlnm1 * .5 + xc1 * s0m) / xam1;
        double s2 = (xc2 * xa2 * xln2 - xnr - f * xa2 * s0) * xa2, s2p = (xc2 * xap2 * xlnp2 - xnr - fp * xap2 * s0p) * xap2, s2m = (xc2 * xam2 * xlnm2 - xnr - fm * xam2 * s0m) * xam2;
        double g1 = (b20 * s0 - 0.5 * xc1 * xln1) / xa1, g1p = (b2p * s0p - 0.5 * xc1 * xlnp1) / xap1, g1m = (b2m * s0m - 0.5 * xc1 * xlnm1) / xam1;
        double g2 = ((0.5 * f * xln2 + 2.0 * s0 * b20 * xc2) * xa2 + xr2) * xa2, g2p = ((0.5 * fp * xlnp2 + 2.0 * s0p * b2p * xc2) * xap2 + xr2) * xap2, g2m = ((0.5 * fm * xlnm2 + 2.0 * s0m * b2m * xc2) * xam2 + xr2) * xam2;
        double bx = b0 * (zs * s0 - 0.5 * (zp * s0p + zm * s0m)) + b1 * (zs * s1 - 0.5 * (zp * s1p + zm * s1m)) + b2 * (zs * s2 - 0.5 * (zp * s2p + zm * s2m));
        double bz = b0 * aln + b1 * (g1 - 0.5 * (g1p + g1m)) + b2 * (g2 - 0.5 * (g2p + g2m));
        return (bx, bz);
    }

    private static (double bx, double by, double bz) Birk1Tot02(double ps, double x, double y, double z) {
        double[] c1 = OneBased(-0.911582E-03, -0.376654E-02, -0.727423E-02, -0.270084E-02, -0.123899E-02, -0.154387E-02, -0.340040E-02, -0.191858E-01, -0.518979E-01, 0.635061E-01, 0.440680, -0.396570, 0.561238E-02, 0.160938E-02, -0.451229E-02, -0.251810E-02, -0.151599E-02, -0.133665E-02, -0.962089E-03, -0.272085E-01, -0.524319E-01, 0.717024E-01, 0.523439, -0.405015, -89.5587, 23.2806);
        double[] c2 = OneBased(6.04133, .305415, .606066E-02, .128379E-03, -.179406E-04, 1.41714, -27.2586, -4.28833, -1.30675, 35.5607, 8.95792, .961617E-03, -.801477E-03, -.782795E-03, -1.65242, -16.5242, -5.33798, .424878E-03, .331787E-03, -.704305E-03, .844342E-03, .953682E-04, .886271E-03, 25.1120, 20.9299, 5.14569, -44.1670, -51.0672, -1.87725, 20.2998, 48.7505, -2.97415, 3.35184, -54.2921, -.838712, -10.5123, 70.7594, -4.94104, .106166E-03, .465791E-03, -.193719E-03, 10.8439, -29.7968, 8.08068, .463507E-03, -.224475E-04, .177035E-03, -.317581E-03, -.264487E-03, .102075E-03, 7.71390, 10.1915, -4.99797, -23.1114, -29.2043, 12.2928, 10.9542, 33.6671, -9.3851, .174615E-03, -.789777E-06, .686047E-03, .460104E-04, -.345216E-02, .221871E-02, .110078E-01, -.661373E-02, .249201E-02, .343978E-01, -.193145E-05, .493963E-05, -.535748E-04, .191833E-04, -.100496E-03, -.210103E-03, -.232195E-02, .315335E-02, -.134320E-01, -.263222E-01);
        double xltd = 78.0, xltnght = 70.0, dtet0 = 0.034906, rh = 9.0, dr = 4.0;
        double tnoonn = (90.0 - xltd) * 0.01745329;
        double tnoons = Pi - tnoonn;
        double dtetdn = (xltd - xltnght) * 0.01745329;
        double sps = Math.Sin(ps);
        double r2 = x * x + y * y + z * z, r = Math.Sqrt(r2), r3 = r * r2;
        double c = Math.Sqrt(Math.Pow(r + rh, 2.0) + dr * dr) - Math.Sqrt(Math.Pow(r - rh, 2.0) + dr * dr);
        double q = Math.Sqrt(Math.Pow(rh + 1.0, 2.0) + dr * dr) - Math.Sqrt(Math.Pow(rh - 1.0, 2.0) + dr * dr);
        double spsas = sps / r * c / q, cpsas = Math.Sqrt(1.0 - spsas * spsas);
        double xas = x * cpsas - z * spsas, zas = x * spsas + z * cpsas;
        double pas = xas != 0.0 || y != 0.0 ? Math.Atan2(y, xas) : 0.0;
        double tas = Math.Atan2(Math.Sqrt(xas * xas + y * y), zas);
        double stas = Math.Sin(tas);
        double ff = stas / Math.Pow(Math.Pow(stas, 6.0) * (1.0 - r3) + r3, 0.1666666667);
        double tet0 = Math.Asin(ff);
        if (tas > 1.5707963) tet0 = Pi - tet0;
        double dtet = dtetdn * Math.Pow(Math.Sin(pas * 0.5), 2.0);
        double tetr1n = tnoonn + dtet, tetr1s = tnoons - dtet;
        int loc = 0;
        if (tet0 < tetr1n - dtet0 || tet0 > tetr1s + dtet0) loc = 1;
        if (tet0 > tetr1n + dtet0 && tet0 < tetr1s - dtet0) loc = 2;
        if (tet0 >= tetr1n - dtet0 && tet0 <= tetr1n + dtet0) loc = 3;
        if (tet0 >= tetr1s - dtet0 && tet0 <= tetr1s + dtet0) loc = 4;

        (double bx, double by, double bz) b;
        if (loc == 1) b = EvalDipLoop1(c1, x, y, z, ps);
        else if (loc == 2) b = EvalConDip1(c2, x, y, z, ps);
        else {
            double t01 = loc == 3 ? tetr1n - dtet0 : tetr1s - dtet0;
            double t02 = loc == 3 ? tetr1n + dtet0 : tetr1s + dtet0;
            double sqr = Math.Sqrt(r);
            double st01as = sqr / Math.Pow(r3 + 1.0 / Math.Pow(Math.Sin(t01), 6.0) - 1.0, 0.1666666667);
            double st02as = sqr / Math.Pow(r3 + 1.0 / Math.Pow(Math.Sin(t02), 6.0) - 1.0, 0.1666666667);
            double ct01as = Math.Sqrt(1.0 - st01as * st01as), ct02as = Math.Sqrt(1.0 - st02as * st02as);
            if (loc == 4) { ct01as = -ct01as; ct02as = -ct02as; }
            var p1 = BoundaryPoint(r, st01as, ct01as, pas, cpsas, spsas);
            var p2 = BoundaryPoint(r, st02as, ct02as, pas, cpsas, spsas);
            var b1 = loc == 3 ? EvalDipLoop1(c1, p1.x, p1.y, p1.z, ps) : EvalConDip1(c2, p1.x, p1.y, p1.z, ps);
            var b2 = loc == 3 ? EvalConDip1(c2, p2.x, p2.y, p2.z, ps) : EvalDipLoop1(c1, p2.x, p2.y, p2.z, ps);
            double ss = Math.Sqrt(Math.Pow(p2.x - p1.x, 2.0) + Math.Pow(p2.y - p1.y, 2.0) + Math.Pow(p2.z - p1.z, 2.0));
            double ds = Math.Sqrt(Math.Pow(x - p1.x, 2.0) + Math.Pow(y - p1.y, 2.0) + Math.Pow(z - p1.z, 2.0));
            double frac = ds / ss;
            b = (b1.bx * (1.0 - frac) + b2.bx * frac, b1.by * (1.0 - frac) + b2.by * frac, b1.bz * (1.0 - frac) + b2.bz * frac);
        }
        var sh = Birk1Shld(ps, x, y, z);
        return (b.bx + sh.bx, b.by + sh.by, b.bz + sh.bz);
    }

    private static (double x, double y, double z) BoundaryPoint(double r, double st, double ct, double pas, double cpsas, double spsas) {
        double xas = r * st * Math.Cos(pas), y = r * st * Math.Sin(pas), zas = r * ct;
        return (xas * cpsas + zas * spsas, y, -xas * spsas + zas * cpsas);
    }

    private static (double bx, double by, double bz) EvalDipLoop1(double[] c, double x, double y, double z, double ps) {
        var d = DipLoop1(x, y, z, ps);
        double bx = 0, by = 0, bz = 0;
        for (int i = 1; i <= 26; i++) { bx += c[i] * d[1, i]; by += c[i] * d[2, i]; bz += c[i] * d[3, i]; }
        return (bx, by, bz);
    }

    private static (double bx, double by, double bz) EvalConDip1(double[] c, double x, double y, double z, double ps) {
        var d = ConDip1(x, y, z, ps);
        double bx = 0, by = 0, bz = 0;
        for (int i = 1; i <= 79; i++) { bx += c[i] * d[1, i]; by += c[i] * d[2, i]; bz += c[i] * d[3, i]; }
        return (bx, by, bz);
    }

    private static double[,] DipLoop1(double x, double y, double z, double ps) {
        double[] xx = OneBased(-11.0, -7.0, -7.0, -3.0, -3.0, 1.0, 1.0, 1.0, 5.0, 5.0, 9.0, 9.0);
        double[] yy = OneBased(2.0, 0.0, 4.0, 2.0, 6.0, 0.0, 4.0, 8.0, 2.0, 6.0, 0.0, 4.0);
        double tilt = 1.00891, dipx = 1.12541, dipy = 0.945719, rh = 9.0, dr = 4.0;
        double[] xcentre = OneBased(2.28397, -5.60831), radius = OneBased(1.86106, 7.83281);
        double[,] d = new double[4, 27];
        double sps = Math.Sin(ps), dr2 = dr * dr, q = Math.Sqrt(Math.Pow(rh + 1.0, 2.0) + dr2) - Math.Sqrt(Math.Pow(rh - 1.0, 2.0) + dr2);
        for (int i = 1; i <= 12; i++) {
            double r = Math.Sqrt(Math.Pow(xx[i] * dipx, 2.0) + Math.Pow(yy[i] * dipy, 2.0));
            double c = Math.Sqrt(Math.Pow(r + rh, 2.0) + dr2) - Math.Sqrt(Math.Pow(r - rh, 2.0) + dr2);
            double spsas = sps / r * c / q, cpsas = Math.Sqrt(1.0 - spsas * spsas);
            double xd = xx[i] * dipx * cpsas, yd = yy[i] * dipy, zd = -xx[i] * dipx * spsas;
            var b1 = DipXyz(x - xd, y - yd, z - zd);
            var b2 = Math.Abs(yd) > 1e-10 ? DipXyz(x - xd, y + yd, z - zd) : ZeroDip();
            d[1, i] = b1.bxz + b2.bxz; d[2, i] = b1.byz + b2.byz; d[3, i] = b1.bzz + b2.bzz;
            d[1, i + 12] = (b1.bxx + b2.bxx) * sps; d[2, i + 12] = (b1.byx + b2.byx) * sps; d[3, i + 12] = (b1.bzx + b2.bzx) * sps;
        }
        for (int j = 1; j <= 2; j++) {
            double rr = j == 1 ? Math.Abs(xcentre[1] + radius[1]) : Math.Abs(radius[2] - xcentre[2]);
            double c = Math.Sqrt(Math.Pow(rr + rh, 2.0) + dr2) - Math.Sqrt(Math.Pow(rr - rh, 2.0) + dr2);
            double spsas = sps / rr * c / q, cpsas = Math.Sqrt(1.0 - spsas * spsas);
            if (j == 1) {
                double xo = x * cpsas - z * spsas, zo = x * spsas + z * cpsas;
                var b = CrossLp(xo, y, zo, xcentre[1], radius[1], tilt);
                d[1, 25] = b.bx * cpsas + b.bz * spsas; d[2, 25] = b.by; d[3, 25] = -b.bx * spsas + b.bz * cpsas;
            } else {
                double xo = x * cpsas - z * spsas - xcentre[2], zo = x * spsas + z * cpsas;
                var b = Circle(xo, y, zo, radius[2]);
                d[1, 26] = b.bx * cpsas + b.bz * spsas; d[2, 26] = b.by; d[3, 26] = -b.bx * spsas + b.bz * cpsas;
            }
        }
        return d;
    }

    private static double[,] ConDip1(double x, double y, double z, double ps) {
        double dx = -0.16, scalein = 0.08, scaleout = 0.4;
        double[] xx = OneBased(-10.0, -7.0, -4.0, -4.0, 0.0, 4.0, 4.0, 7.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        double[] yy = OneBased(3.0, 6.0, 3.0, 9.0, 6.0, 3.0, 9.0, 6.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        double[] zz = OneBased(20.0, 20.0, 4.0, 20.0, 4.0, 4.0, 20.0, 20.0, 20.0, 2.0, 3.0, 4.5, 7.0, 10.0);
        double[,] d = new double[4, 80];
        double sps = Math.Sin(ps), cps = Math.Cos(ps);
        double xsm = x * cps - z * sps - dx, zsm = z * cps + x * sps;
        double ro2 = xsm * xsm + y * y, ro = Math.Sqrt(ro2);
        double[] cf = new double[6], sf = new double[6];
        cf[1] = xsm / ro; sf[1] = y / ro;
        for (int i = 2; i <= 5; i++) { cf[i] = cf[i - 1] * cf[1] - sf[i - 1] * sf[1]; sf[i] = sf[i - 1] * cf[1] + cf[i - 1] * sf[1]; }
        double r = Math.Sqrt(ro2 + zsm * zsm), c = zsm / r, s = ro / r, ch = Math.Sqrt(0.5 * (1.0 + c)), sh = Math.Sqrt(0.5 * (1.0 - c)), tnh = sh / ch, cnh = 1.0 / tnh;
        for (int m = 1; m <= 5; m++) {
            double bt = m * cf[m] / (r * s) * (Math.Pow(tnh, m) + Math.Pow(cnh, m));
            double bf = -0.5 * m * sf[m] / r * (Math.Pow(tnh, m - 1) / (ch * ch) - Math.Pow(cnh, m - 1) / (sh * sh));
            double bxsm = bt * c * cf[1] - bf * sf[1], by = bt * c * sf[1] + bf * cf[1], bzsm = -bt * s;
            d[1, m] = bxsm * cps + bzsm * sps; d[2, m] = by; d[3, m] = -bxsm * sps + bzsm * cps;
        }
        xsm = x * cps - z * sps; zsm = z * cps + x * sps;
        for (int i = 1; i <= 9; i++) {
            double sc = (i == 3 || i == 5 || i == 6) ? scalein : scaleout;
            double xd = xx[i] * sc, yd = yy[i] * sc, zd = zz[i];
            var b1 = DipXyz(xsm - xd, y - yd, zsm - zd);
            var b2 = DipXyz(xsm - xd, y + yd, zsm - zd);
            var b3 = DipXyz(xsm - xd, y - yd, zsm + zd);
            var b4 = DipXyz(xsm - xd, y + yd, zsm + zd);
            int ix = i * 3 + 3, iy = ix + 1, iz = iy + 1;
            PutDipCombo(d, ix, cps, sps, b1.bxx + b2.bxx - b3.bxx - b4.bxx, b1.byx + b2.byx - b3.byx - b4.byx, b1.bzx + b2.bzx - b3.bzx - b4.bzx);
            PutDipCombo(d, iy, cps, sps, b1.bxy - b2.bxy - b3.bxy + b4.bxy, b1.byy - b2.byy - b3.byy + b4.byy, b1.bzy - b2.bzy - b3.bzy + b4.bzy);
            PutDipCombo(d, iz, cps, sps, b1.bxz + b2.bxz + b3.bxz + b4.bxz, b1.byz + b2.byz + b3.byz + b4.byz, b1.bzz + b2.bzz + b3.bzz + b4.bzz);
            ix += 27; iy += 27; iz += 27;
            PutDipCombo(d, ix, cps, sps, sps * (b1.bxx + b2.bxx + b3.bxx + b4.bxx), sps * (b1.byx + b2.byx + b3.byx + b4.byx), sps * (b1.bzx + b2.bzx + b3.bzx + b4.bzx));
            PutDipCombo(d, iy, cps, sps, sps * (b1.bxy - b2.bxy + b3.bxy - b4.bxy), sps * (b1.byy - b2.byy + b3.byy - b4.byy), sps * (b1.bzy - b2.bzy + b3.bzy - b4.bzy));
            PutDipCombo(d, iz, cps, sps, sps * (b1.bxz + b2.bxz - b3.bxz - b4.bxz), sps * (b1.byz + b2.byz - b3.byz - b4.byz), sps * (b1.bzz + b2.bzz - b3.bzz - b4.bzz));
        }
        for (int i = 1; i <= 5; i++) {
            double zd = zz[i + 9];
            var b1 = DipXyz(xsm, y, zsm - zd);
            var b2 = DipXyz(xsm, y, zsm + zd);
            int ix = 58 + i * 2, iz = ix + 1;
            PutDipCombo(d, ix, cps, sps, b1.bxx - b2.bxx, b1.byx - b2.byx, b1.bzx - b2.bzx);
            PutDipCombo(d, iz, cps, sps, b1.bxz + b2.bxz, b1.byz + b2.byz, b1.bzz + b2.bzz);
            ix += 10; iz += 10;
            PutDipCombo(d, ix, cps, sps, sps * (b1.bxx + b2.bxx), sps * (b1.byx + b2.byx), sps * (b1.bzx + b2.bzx));
            PutDipCombo(d, iz, cps, sps, sps * (b1.bxz - b2.bxz), sps * (b1.byz - b2.byz), sps * (b1.bzz - b2.bzz));
        }
        return d;
    }

    private static void PutDipCombo(double[,] d, int i, double cps, double sps, double bxsm, double by, double bzsm) {
        d[1, i] = bxsm * cps + bzsm * sps;
        d[2, i] = by;
        d[3, i] = bzsm * cps - bxsm * sps;
    }

    private static (double bxx, double byx, double bzx, double bxy, double byy, double bzy, double bxz, double byz, double bzz) DipXyz(double x, double y, double z) {
        double x2 = x * x, y2 = y * y, z2 = z * z, r2 = x2 + y2 + z2;
        double xmr5 = 30574.0 / (r2 * r2 * Math.Sqrt(r2)), xmr53 = 3.0 * xmr5;
        double bxx = xmr5 * (3.0 * x2 - r2), byx = xmr53 * x * y, bzx = xmr53 * x * z;
        double byy = xmr5 * (3.0 * y2 - r2), bzy = xmr53 * y * z, bzz = xmr5 * (3.0 * z2 - r2);
        return (bxx, byx, bzx, byx, byy, bzy, bzx, bzy, bzz);
    }

    private static (double bxx, double byx, double bzx, double bxy, double byy, double bzy, double bxz, double byz, double bzz) ZeroDip() => (0, 0, 0, 0, 0, 0, 0, 0, 0);

    private static (double bx, double by, double bz) Birk1Shld(double ps, double x, double y, double z) {
        double[] a = OneBased(1.174198045, -1.463820502, 4.840161537, -3.674506864, 82.18368896, -94.94071588, -4122.331796, 4670.278676, -21.54975037, 26.72661293, -72.81365728, 44.09887902, 40.08073706, -51.23563510, 1955.348537, -1940.971550, 794.0496433, -982.2441344, 1889.837171, -558.9779727, -1260.543238, 1260.063802, -293.5942373, 344.7250789, -773.7002492, 957.0094135, -1824.143669, 520.7994379, 1192.484774, -1192.184565, 89.15537624, -98.52042999, -0.8168777675E-01, 0.4255969908E-01, 0.3155237661, -0.3841755213, 2.494553332, -0.6571440817E-01, -2.765661310, 0.4331001908, 0.1099181537, -0.6154126980E-01, -0.3258649260, 0.6698439193, -5.542735524, 0.1604203535, 5.854456934, -0.8323632049, 3.732608869, -3.130002153, 107.0972607, -32.28483411, -115.2389298, 54.45064360, -0.5826853320, -3.582482231, -4.046544561, 3.311978102, -104.0839563, 30.26401293, 97.29109008, -50.62370872, -296.3734955, 127.7872523, 5.303648988, 10.40368955, 69.65230348, 466.5099509, 1.645049286, 3.825838190, 11.66675599, 558.9781177, 1.826531343, 2.066018073, 25.40971369, 990.2795225, 2.319489258, 4.555148484, 9.691185703, 591.8280358);
        double cps = Math.Cos(ps), sps = Math.Sin(ps), s3ps = 4.0 * cps * cps - 1.0;
        double bx = 0, by = 0, bz = 0; int l = 0;
        for (int m = 1; m <= 2; m++)
        for (int i = 1; i <= 4; i++) {
            double rp = 1.0 / a[64 + i], rr = 1.0 / a[68 + i], rq = 1.0 / a[72 + i], rs = 1.0 / a[76 + i];
            double cypi = Math.Cos(y * rp), cyqi = Math.Cos(y * rq), sypi = Math.Sin(y * rp), syqi = Math.Sin(y * rq);
            for (int k = 1; k <= 4; k++) {
                double rrk = 1.0 / a[68 + k], rsk = 1.0 / a[76 + k];
                double szrk = Math.Sin(z * rrk), czsk = Math.Cos(z * rsk), czrk = Math.Cos(z * rrk), szsk = Math.Sin(z * rsk);
                double sqpr = Math.Sqrt(rp * rp + rrk * rrk), sqqs = Math.Sqrt(rq * rq + rsk * rsk);
                double epr = Math.Exp(x * sqpr), eqs = Math.Exp(x * sqqs);
                double hx = 0, hy = 0, hz = 0;
                for (int n = 1; n <= 2; n++) {
                    if (m == 1) {
                        if (n == 1) { hx = -sqpr * epr * cypi * szrk; hy = rp * epr * sypi * szrk; hz = -rrk * epr * cypi * czrk; }
                        else { hx *= cps; hy *= cps; hz *= cps; }
                    } else {
                        if (n == 1) { hx = -sps * sqqs * eqs * cyqi * czsk; hy = sps * rq * eqs * syqi * czsk; hz = sps * rsk * eqs * cyqi * szsk; }
                        else { hx *= s3ps; hy *= s3ps; hz *= s3ps; }
                    }
                    l++; bx += a[l] * hx; by += a[l] * hy; bz += a[l] * hz;
                }
            }
        }
        return (bx, by, bz);
    }

    private static (double bx, double by, double bz) Birk2Tot02(double ps, double x, double y, double z) {
        var w = Birk2Shl(x, y, z, ps);
        var h = R2Birk(x, y, z, ps);
        return (w.bx + h.bx, w.by + h.by, w.bz + h.bz);
    }

    private static (double bx, double by, double bz) Birk2Shl(double x, double y, double z, double ps) {
        double[] a = OneBased(-111.6371348, 124.5402702, 110.3735178, -122.0095905, 111.9448247, -129.1957743, -110.7586562, 126.5649012, -0.7865034384, -0.2483462721, 0.8026023894, 0.2531397188, 10.72890902, 0.8483902118, -10.96884315, -0.8583297219, 13.85650567, 14.90554500, 10.21914434, 10.09021632, 6.340382460, 14.40432686, 12.71023437, 12.83966657);
        double cps = Math.Cos(ps), sps = Math.Sin(ps), s3ps = 4.0 * cps * cps - 1.0;
        double hx = 0.0, hy = 0.0, hz = 0.0;
        int l = 0;
        for (int m = 1; m <= 2; m++)
        for (int i = 1; i <= 2; i++) {
            double p = a[16 + i];
            double q = a[20 + i];
            double cypi = Math.Cos(y / p), cyqi = Math.Cos(y / q);
            double sypi = Math.Sin(y / p), syqi = Math.Sin(y / q);
            for (int k = 1; k <= 2; k++) {
                double r = a[18 + k], s = a[22 + k];
                double szrk = Math.Sin(z / r), czsk = Math.Cos(z / s), czrk = Math.Cos(z / r), szsk = Math.Sin(z / s);
                double sqpr = Math.Sqrt(1.0 / (p * p) + 1.0 / (r * r));
                double sqqs = Math.Sqrt(1.0 / (q * q) + 1.0 / (s * s));
                double epr = Math.Exp(x * sqpr), eqs = Math.Exp(x * sqqs);
                double dx = 0.0, dy = 0.0, dz = 0.0;
                for (int n = 1; n <= 2; n++) {
                    l++;
                    if (m == 1) {
                        if (n == 1) { dx = -sqpr * epr * cypi * szrk; dy = epr / p * sypi * szrk; dz = -epr / r * cypi * czrk; }
                        else { dx *= cps; dy *= cps; dz *= cps; }
                    } else {
                        if (n == 1) { dx = -sps * sqqs * eqs * cyqi * czsk; dy = sps * eqs / q * syqi * czsk; dz = sps * eqs / s * cyqi * szsk; }
                        else { dx *= s3ps; dy *= s3ps; dz *= s3ps; }
                    }
                    hx += a[l] * dx; hy += a[l] * dy; hz += a[l] * dz;
                }
            }
        }
        return (hx, hy, hz);
    }

    private static (double bx, double by, double bz) R2Birk(double x, double y, double z, double ps) {
        const double delarg = 0.030, delarg1 = 0.015;
        double cps = Math.Cos(ps), sps = Math.Sin(ps);
        double xsm = x * cps - z * sps;
        double zsm = z * cps + x * sps;
        double xks = Xksi(xsm, y, zsm);
        double bxsm = 0.0, bysm = 0.0, bzsm = 0.0;
        if (xks < -(delarg + delarg1)) {
            var outer = R2Outer(xsm, y, zsm);
            bxsm = -outer.bx * 0.02; bysm = -outer.by * 0.02; bzsm = -outer.bz * 0.02;
        } else if (xks >= -(delarg + delarg1) && xks < -delarg + delarg1) {
            var outer = R2Outer(xsm, y, zsm);
            var sheet = R2Sheet(xsm, y, zsm);
            double f2 = -0.02 * Tksi(xks, -delarg, delarg1);
            double f1 = -0.02 - f2;
            bxsm = outer.bx * f1 + sheet.bx * f2; bysm = outer.by * f1 + sheet.by * f2; bzsm = outer.bz * f1 + sheet.bz * f2;
        } else if (xks >= -delarg + delarg1 && xks < delarg - delarg1) {
            var sheet = R2Sheet(xsm, y, zsm);
            bxsm = -sheet.bx * 0.02; bysm = -sheet.by * 0.02; bzsm = -sheet.bz * 0.02;
        } else if (xks >= delarg - delarg1 && xks < delarg + delarg1) {
            var inner = R2Inner(xsm, y, zsm);
            var sheet = R2Sheet(xsm, y, zsm);
            double f1 = -0.02 * Tksi(xks, delarg, delarg1);
            double f2 = -0.02 - f1;
            bxsm = inner.bx * f1 + sheet.bx * f2; bysm = inner.by * f1 + sheet.by * f2; bzsm = inner.bz * f1 + sheet.bz * f2;
        } else {
            var inner = R2Inner(xsm, y, zsm);
            bxsm = -inner.bx * 0.02; bysm = -inner.by * 0.02; bzsm = -inner.bz * 0.02;
        }
        return (bxsm * cps + bzsm * sps, bysm, bzsm * cps - bxsm * sps);
    }

    private static (double bx, double by, double bz) R2Inner(double x, double y, double z) {
        double[] pl = OneBased(154.185, -2.12446, .601735E-01, -.153954E-02, .355077E-04, 29.9996, 262.886, 99.9132);
        double[] pn = OneBased(-8.1902, 6.5239, 5.504, 7.7815, .8573, 3.0986, .0774, -.038);
        var c = BConic(x, y, z, 5);
        var db8 = Loops4(x, y, z, pn[1], pn[2], pn[3], pn[4], pn[5], pn[6]);
        var db6 = DipDistr(x - pn[7], y, z, 0);
        var db7 = DipDistr(x - pn[8], y, z, 1);
        double bx = 0, by = 0, bz = 0;
        for (int i = 1; i <= 5; i++) { bx += pl[i] * c.bx[i]; by += pl[i] * c.by[i]; bz += pl[i] * c.bz[i]; }
        bx += pl[6] * db6.bx + pl[7] * db7.bx + pl[8] * db8.bx;
        by += pl[6] * db6.by + pl[7] * db7.by + pl[8] * db8.by;
        bz += pl[6] * db6.bz + pl[7] * db7.bz + pl[8] * db8.bz;
        return (bx, by, bz);
    }

    private static (double bx, double by, double bz) R2Outer(double x, double y, double z) {
        double[] pl = OneBased(-34.105, -2.00019, 628.639, 73.4847, 12.5162);
        double[] pn = OneBased(.55, .694, .0031, 1.55, 2.8, .1375, -.7, .2, .9625, -2.994, 2.925, -1.775, 4.3, -.275, 2.7, .4312, 1.55);
        var db1 = CrossLp(x, y, z, pn[1], pn[2], pn[3]);
        var db2 = CrossLp(x, y, z, pn[4], pn[5], pn[6]);
        var db3 = CrossLp(x, y, z, pn[7], pn[8], pn[9]);
        var db4 = Circle(x - pn[10], y, z, pn[11]);
        var db5 = Loops4(x, y, z, pn[12], pn[13], pn[14], pn[15], pn[16], pn[17]);
        return (pl[1] * db1.bx + pl[2] * db2.bx + pl[3] * db3.bx + pl[4] * db4.bx + pl[5] * db5.bx,
                pl[1] * db1.by + pl[2] * db2.by + pl[3] * db3.by + pl[4] * db4.by + pl[5] * db5.by,
                pl[1] * db1.bz + pl[2] * db2.bz + pl[3] * db3.bz + pl[4] * db4.bz + pl[5] * db5.bz);
    }

    private static (double bx, double by, double bz) R2Sheet(double x, double y, double z) {
        double[] pnx = OneBased(-19.0969, -9.28828, -0.129687, 5.58594, 22.5055, 0.483750E-01, 0.396953E-01, 0.579023E-01);
        double[] pny = OneBased(-13.6750, -6.70625, 2.31875, 11.4062, 20.4562, 0.478750E-01, 0.363750E-01, 0.567500E-01);
        double[] pnz = OneBased(-16.7125, -16.4625, -0.1625, 5.1, 23.7125, 0.355625E-01, 0.318750E-01, 0.538750E-01);
        double[] aa = OneBased(8.07190, -7.39582, -7.62341, 0.684671, -13.5672, 11.6681, 13.1154, -0.890217, 7.78726, -5.38346, -8.08738, 0.609385, -2.70410, 3.53741, 3.15549, -1.11069, -8.47555, 0.278122, 2.73514, 4.55625, 13.1134, 1.15848, -3.52648, -8.24698, -6.85710, -2.81369, 2.03795, 4.64383, 2.49309, -1.22041, -1.67432, -0.422526, -5.39796, 7.10326, 5.53730, -13.1918, 4.67853, -7.60329, -2.53066, 7.76338, 5.60165, 5.34816, -4.56441, 7.05976, -2.62723, -0.529078, 1.42019, -2.93919, 55.6338, -1.55181, 39.8311, -80.6561, -46.9655, 32.8925, -6.32296, 19.7841, 124.731, 10.4347, -30.7581, 102.680, -47.4037, -3.31278, 9.37141, -50.0268, -533.319, 110.426, 1000.20, -1051.40, 1619.48, 589.855, -1462.73, 1087.10, -1994.73, -1654.12, 1263.33, -260.210, 1424.84, 1255.71, -956.733, 219.946);
        double[] bb = OneBased(-9.08427, 10.6777, 10.3288, -0.969987, 6.45257, -8.42508, -7.97464, 1.41996, -1.92490, 3.93575, 2.83283, -1.48621, 0.244033, -0.757941, -0.386557, 0.344566, 9.56674, -2.5365, -3.32916, -5.86712, -6.19625, 1.83879, 2.52772, 4.34417, 1.87268, -2.13213, -1.69134, -.176379, -.261359, .566419, 0.3138, -0.134699, -3.83086, -8.4154, 4.77005, -9.31479, 37.5715, 19.3992, -17.9582, 36.4604, -14.9993, -3.1442, 6.17409, -15.5519, 2.28621, -0.891549E-2, -.462912, 2.47314, 41.7555, 208.614, -45.7861, -77.8687, 239.357, -67.9226, 66.8743, 238.534, -112.136, 16.2069, -40.4706, -134.328, 21.56, -0.201725, 2.21, 32.5855, -108.217, -1005.98, 585.753, 323.668, -817.056, 235.750, -560.965, -576.892, 684.193, 85.0275, 168.394, 477.776, -289.253, -123.216, 75.6501, -178.605);
        double[] cc = OneBased(1167.61, -917.782, -1253.2, -274.128, -1538.75, 1257.62, 1745.07, 113.479, 393.326, -426.858, -641.1, 190.833, -29.9435, -1.04881, 117.125, -25.7663, -1168.16, 910.247, 1239.31, 289.515, 1540.56, -1248.29, -1727.61, -131.785, -394.577, 426.163, 637.422, -187.965, 30.0348, 0.221898, -116.68, 26.0291, 12.6804, 4.84091, 1.18166, -2.75946, -17.9822, -6.80357, -1.47134, 3.02266, 4.79648, 0.665255, -0.256229, -0.857282E-1, -0.588997, 0.634812E-1, 0.164303, -0.15285, 22.2524, -22.4376, -3.85595, 6.07625, -105.959, -41.6698, 0.378615, 1.55958, 44.3981, 18.8521, 3.19466, 5.89142, -8.63227, -2.36418, -1.027, -2.31515, 1035.38, 2040.66, -131.881, -744.533, -3274.93, -4845.61, 482.438, 1567.43, 1354.02, 2040.47, -151.653, -845.012, -111.723, -265.343, -26.1171, 216.632);
        double xks = Xksi(x, y, z);
        double[] tx = R2SheetT(xks, pnx), ty = R2SheetT(xks, pny), tz = R2SheetT(xks, pnz);
        double rho2 = x * x + y * y, r = Math.Sqrt(rho2 + z * z), rho = Math.Sqrt(rho2);
        double c1p = x / rho, s1p = y / rho, s2p = 2.0 * s1p * c1p, c2p = c1p * c1p - s1p * s1p, s3p = s2p * c1p + c2p * s1p, c3p = c2p * c1p - s2p * s1p, s4p = s3p * c1p + c3p * s1p;
        double ct = z / r;
        double[] hx = OneBased(1.0, c1p, c2p, c3p);
        double[] hy = OneBased(s1p, s2p, s3p, s4p);
        double[] sx = OneBased(Fexp(ct, pnx[1]), Fexp(ct, pnx[2]), Fexp(ct, pnx[3]), Fexp(ct, pnx[4]), Fexp(ct, pnx[5]));
        double[] sy = OneBased(Fexp(ct, pny[1]), Fexp(ct, pny[2]), Fexp(ct, pny[3]), Fexp(ct, pny[4]), Fexp(ct, pny[5]));
        double[] sz = OneBased(Fexp1(ct, pnz[1]), Fexp1(ct, pnz[2]), Fexp1(ct, pnz[3]), Fexp1(ct, pnz[4]), Fexp1(ct, pnz[5]));
        return (EvalSheet(aa, sx, hx, tx), EvalSheet(bb, sy, hy, ty), EvalSheet(cc, sz, hx, tz));
    }

    private static double[] R2SheetT(double xks, double[] p) => OneBased(
        1.0,
        xks / Math.Sqrt(xks * xks + p[6] * p[6]),
        Math.Pow(p[7], 3.0) / Math.Pow(Math.Sqrt(xks * xks + p[7] * p[7]), 3.0),
        xks / Math.Pow(Math.Sqrt(xks * xks + p[8] * p[8]), 5.0) * 3.493856 * Math.Pow(p[8], 4.0));

    private static double EvalSheet(double[] coeff, double[] s, double[] h, double[] t) {
        double result = 0.0;
        for (int isv = 1; isv <= 5; isv++)
        for (int ih = 1; ih <= 4; ih++)
        for (int it = 1; it <= 4; it++) {
            int index = (isv - 1) * 16 + (ih - 1) * 4 + it;
            result += s[isv] * h[ih] * coeff[index] * t[it];
        }
        return result;
    }

    private static (double[] bx, double[] by, double[] bz) BConic(double x, double y, double z, int nmax) {
        double[] cbx = new double[nmax + 1], cby = new double[nmax + 1], cbz = new double[nmax + 1];
        double ro2 = x * x + y * y, ro = Math.Sqrt(ro2);
        double cf = x / ro, sf = y / ro, cfm1 = 1.0, sfm1 = 0.0;
        double r2 = ro2 + z * z, r = Math.Sqrt(r2), c = z / r, s = ro / r;
        double ch = Math.Sqrt(0.5 * (1.0 + c)), sh = Math.Sqrt(0.5 * (1.0 - c));
        double tnhm1 = 1.0, cnhm1 = 1.0, tnh = sh / ch, cnh = 1.0 / tnh;
        for (int m = 1; m <= nmax; m++) {
            double cfm = cfm1 * cf - sfm1 * sf, sfm = cfm1 * sf + sfm1 * cf;
            cfm1 = cfm; sfm1 = sfm;
            double tnhm = tnhm1 * tnh, cnhm = cnhm1 * cnh;
            double bt = m * cfm / (r * s) * (tnhm + cnhm);
            double bf = -0.5 * m * sfm / r * (tnhm1 / (ch * ch) - cnhm1 / (sh * sh));
            tnhm1 = tnhm; cnhm1 = cnhm;
            cbx[m] = bt * c * cf - bf * sf;
            cby[m] = bt * c * sf + bf * cf;
            cbz[m] = -bt * s;
        }
        return (cbx, cby, cbz);
    }

    private static (double bx, double by, double bz) DipDistr(double x, double y, double z, int mode) {
        double x2 = x * x, rho2 = x2 + y * y, r2 = rho2 + z * z, r3 = r2 * Math.Sqrt(r2);
        if (mode == 0) return (z / Math.Pow(rho2, 2.0) * (r2 * (y * y - x2) - rho2 * x2) / r3,
                               -x * y * z / Math.Pow(rho2, 2.0) * (2.0 * r2 + rho2) / r3,
                               x / r3);
        return (z / Math.Pow(rho2, 2.0) * (y * y - x2),
                -2.0 * x * y * z / Math.Pow(rho2, 2.0),
                x / rho2);
    }

    private static (double bx, double by, double bz) Circle(double x, double y, double z, double rl) {
        double pi = 3.141592654;
        double rho2 = x * x + y * y, rho = Math.Sqrt(rho2);
        double r22 = z * z + Math.Pow(rho + rl, 2.0), r2 = Math.Sqrt(r22), r12 = r22 - 4.0 * rho * rl, r32 = 0.5 * (r12 + r22);
        double xk2 = 1.0 - r12 / r22, xk2s = 1.0 - xk2, dl = Math.Log(1.0 / xk2s);
        double k = 1.38629436112 + xk2s * (0.09666344259 + xk2s * (0.03590092383 + xk2s * (0.03742563713 + xk2s * 0.01451196212))) + dl * (0.5 + xk2s * (0.12498593597 + xk2s * (0.06880248576 + xk2s * (0.03328355346 + xk2s * 0.00441787012))));
        double e = 1.0 + xk2s * (0.44325141463 + xk2s * (0.0626060122 + xk2s * (0.04757383546 + xk2s * 0.01736506451))) + dl * xk2s * (0.2499836831 + xk2s * (0.09200180037 + xk2s * (0.04069697526 + xk2s * 0.00526449639)));
        double brho = rho > 1e-6 ? z / (rho2 * r2) * (r32 / r12 * e - k) : pi * rl / r2 * (rl - rho) / r12 * z / (r32 - rho2);
        return (brho * x, brho * y, (k - e * (r32 - 2.0 * rl * rl) / r12) / r2);
    }

    private static (double bx, double by, double bz) CrossLp(double x, double y, double z, double xc, double rl, double al) {
        double cal = Math.Cos(al), sal = Math.Sin(al);
        double y1 = y * cal - z * sal, z1 = y * sal + z * cal;
        double y2 = y * cal + z * sal, z2 = -y * sal + z * cal;
        var b1 = Circle(x - xc, y1, z1, rl);
        var b2 = Circle(x - xc, y2, z2, rl);
        return (b1.bx + b2.bx, (b1.by + b2.by) * cal + (b1.bz - b2.bz) * sal, -(b1.by - b2.by) * sal + (b1.bz + b2.bz) * cal);
    }

    private static (double bx, double by, double bz) Loops4(double x, double y, double z, double xc, double yc, double zc, double r, double theta, double phi) {
        double ct = Math.Cos(theta), st = Math.Sin(theta), cp = Math.Cos(phi), sp = Math.Sin(phi);
        (double bx, double by, double bz) One(double xs, double yss, double zs, int q) {
            double xss = xs * ct - zs * st, zss = zs * ct + xs * st;
            var b = Circle(xss, yss, zss, r);
            double bxs = b.bx * ct + b.bz * st;
            double bzq = b.bz * ct - b.bx * st;
            if (q == 1) return (bxs * cp - b.by * sp, bxs * sp + b.by * cp, bzq);
            if (q == 2) return (bxs * cp + b.by * sp, -bxs * sp + b.by * cp, bzq);
            if (q == 3) return (-bxs * cp - b.by * sp, bxs * sp - b.by * cp, bzq);
            return (-bxs * cp + b.by * sp, -bxs * sp - b.by * cp, bzq);
        }
        var b1 = One((x - xc) * cp + (y - yc) * sp, (y - yc) * cp - (x - xc) * sp, z - zc, 1);
        var b2 = One((x - xc) * cp - (y + yc) * sp, (y + yc) * cp + (x - xc) * sp, z - zc, 2);
        var b3 = One(-(x - xc) * cp + (y + yc) * sp, -(y + yc) * cp - (x - xc) * sp, z + zc, 3);
        var b4 = One(-(x - xc) * cp - (y - yc) * sp, -(y - yc) * cp + (x - xc) * sp, z + zc, 4);
        return (b1.bx + b2.bx + b3.bx + b4.bx, b1.by + b2.by + b3.by + b4.by, b1.bz + b2.bz + b3.bz + b4.bz);
    }

    private static double Xksi(double x, double y, double z) {
        double a11a12 = 0.305662, a21a22 = -0.383593, a41a42 = 0.2677733, a51a52 = -0.097656, a61a62 = -0.636034, b11b12 = -0.359862, b21b22 = 0.424706, c61c62 = -0.126366, c71c72 = 0.292578, r0 = 1.21563, dr = 7.50937;
        double tnoon = 0.3665191, dteta = 0.09599309;
        double x2 = x * x, y2 = y * y, z2 = z * z, r2 = x2 + y2 + z2, r = Math.Sqrt(r2);
        double xr = x / r, yr = y / r, zr = z / r;
        double pr = r < r0 ? 0.0 : Math.Sqrt(Math.Pow(r - r0, 2.0) + dr * dr) - dr;
        double f = x + pr * (a11a12 + a21a22 * xr + a41a42 * xr * xr + a51a52 * yr * yr + a61a62 * zr * zr);
        double g = y + pr * (b11b12 * yr + b21b22 * xr * yr);
        double h = z + pr * (c61c62 * zr + c71c72 * xr * zr);
        double fchsg2 = f * f + g * g;
        if (fchsg2 < 1e-5) return -1.0;
        double fgh = fchsg2 + h * h;
        double alpha = fchsg2 / Math.Pow(Math.Sqrt(fgh), 3.0);
        double theta = tnoon + 0.5 * dteta * (1.0 - f / Math.Sqrt(fchsg2));
        return alpha - Math.Pow(Math.Sin(theta), 2.0);
    }

    private static double Fexp(double s, double a) {
        const double e = 2.718281828459;
        return a < 0.0 ? Math.Sqrt(-2.0 * a * e) * s * Math.Exp(a * s * s) : s * Math.Exp(a * (s * s - 1.0));
    }

    private static double Fexp1(double s, double a) => a <= 0.0 ? Math.Exp(a * s * s) : Math.Exp(a * (s * s - 1.0));

    private static double Tksi(double xksi, double xks0, double dxksi) {
        double tdz3 = 2.0 * Math.Pow(dxksi, 3.0);
        if (xksi - xks0 < -dxksi) return 0.0;
        if (xksi - xks0 >= dxksi) return 1.0;
        if (xksi >= xks0 - dxksi && xksi < xks0) {
            double br3 = Math.Pow(xksi - xks0 + dxksi, 3.0);
            return 1.5 * br3 / (tdz3 + br3);
        }
        double br32 = Math.Pow(xksi - xks0 - dxksi, 3.0);
        return 1.0 + 1.5 * br32 / (tdz3 - br32);
    }

    private static (double bx, double by, double bz) Dipole(double ps, double x, double y, double z) {
        double sps = Math.Sin(ps), cps = Math.Cos(ps);
        double p = x * x, u = z * z, v = 3.0 * z * x, t = y * y;
        double q = 30574.0 / Math.Pow(Math.Sqrt(p + t + u), 5.0);
        double bx = q * ((t + u - 2.0 * p) * sps - v * cps);
        double by = -3.0 * y * q * (x * sps + z * cps);
        double bz = q * ((p + t - 2.0 * u) * cps - v * sps);
        return (bx, by, bz);
    }

    private static double[] OneBased(params double[] values) {
        double[] result = new double[values.Length + 1];
        Array.Copy(values, 0, result, 1, values.Length);
        return result;
    }

    private static double Sign1(double value) => value >= 0.0 ? 1.0 : -1.0;
}
