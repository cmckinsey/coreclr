// Licensed to the .NET Foundation under one or more agreements.
// The .NET Foundation licenses this file to you under the MIT license.
// See the LICENSE file in the project root for more information.
//
// C# adaptation of C implementation of Livermore Loops Fortran benchmark.

/* Livermore Loops coded in C        Latest File Modification  20 Oct 92,
 *  by Tim Peters, Kendall Square Res. Corp. tim@ksr.com, ksr!tim@uunet.uu.net
 *     SUBROUTINE KERNEL( TK)  replaces the Fortran routine in LFK Test program.
 ************************************************************************
 *                                                                      *
 *            KERNEL     executes 24 samples of "C" computation         *
 *                                                                      *
 *                TK(1) - total cpu time to execute only the 24 kernels.*
 *                TK(2) - total Flops executed by the 24 Kernels        *
 *                                                                      *
 ************************************************************************
 *                                                                      *
 *     L. L. N. L.   " C "   K E R N E L S:   M F L O P S               *
 *                                                                      *
 *     These kernels measure   " C "   numerical computation            *
 *     rates for  a  spectrum  of  cpu-limited computational            *
 *     structures or benchmarks.   Mathematical  through-put            *
 *     is measured  in  units  of millions of floating-point            *
 *     operations executed per second, called Megaflops/sec.            *
 *                                                                      *
 *     Fonzi's Law: There is not now and there never will be a language *
 *                  in which it is the least bit difficult to write     *
 *                  bad programs.                                       *
 *                                                    F.H.MCMAHON  1972 *
 ************************************************************************
 *Originally from  Greg Astfalk, AT&T, P.O.Box 900, Princeton, NJ. 08540*
 *               by way of Frank McMahon (LLNL).                        *
 *                                                                      *
 *                               REFERENCE                              *
 *                                                                      *
 *              F.H.McMahon,   The Livermore Fortran Kernels:           *
 *              A Computer Test Of The Numerical Performance Range,     *
 *              Lawrence Livermore National Laboratory,                 *
 *              Livermore, California, UCRL-53745, December 1986.       *
 *                                                                      *
 *       from:  National Technical Information Service                  *
 *              U.S. Department of Commerce                             *
 *              5285 Port Royal Road                                    *
 *              Springfield, VA.  22161                                 *
 *                                                                      *
 *    Changes made to correct many array subscripting problems,         *
 *      make more readable (added #define's), include the original      *
 *      FORTRAN versions of the runs as comments, and make more         *
 *      portable by Kelly O'Hair (LLNL) and Chuck Rasbold (LLNL).       *
 *                                                                      *
 ************************************************************************
 */

using Microsoft.Xunit.Performance;
using System;
using System.Runtime.CompilerServices;
using Xunit;

[assembly: OptimizeForBenchmarks]
[assembly: MeasureInstructionsRetired]

public class LLoops
{

#if DEBUG
    public const int Iterations = 1;
#else
    public const int Iterations = 4000;
#endif

    const double MaxErr = 1.0e-6;

    double[] m_x = new double[1002];
    double[] m_y = new double[1002];
    double[] m_z = new double[1002];
    double[] m_u = new double[501];
    double[][] m_px;
    double[][] m_cx;
    double[][][] m_u1;
    double[][][] m_u2;
    double[][][] m_u3;
    double[][] m_b;
    double[] m_bnk1 = new double[6];
    double[][] m_c;
    double[] m_bnk2 = new double[6];
    double[][] m_p;
    double[] m_bnk3 = new double[6];
    double[][] m_h;
    double[] m_bnk4 = new double[6];
    double[] m_bnk5 = new double[6];
    double[] m_ex = new double[68];
    double[] m_rh = new double[68];
    double[] m_dex = new double[68];
    double[] m_vx = new double[151];
    double[] m_xx = new double[151];
    double[] m_grd = new double[151];
    int[] m_e = new int[193];
    int[] m_f = new int[193];
    int[] m_nrops = { 0, 5, 10, 2, 2, 2, 2, 16, 36, 17, 9, 1, 1, 7, 11 };
    int[] m_loops = { 0, 400, 200, 1000, 510, 1000, 1000, 120, 40, 100, 100, 1000, 1000, 128, 150 };
    double[] m_checks = {
        0, 0.811986948148e+07, 0.356310000000e+03, 0.356310000000e+03, -0.402412007078e+05,
        0.136579037764e+06, 0.419716278716e+06,
        0.429449847526e+07, 0.314064400000e+06,
        0.182709000000e+07, -0.140415250000e+09,
        0.374895020500e+09, 0.000000000000e+00,
        0.171449024000e+06, -0.510829560800e+07
    };

    public static volatile object VolatileObject;

    [MethodImpl(MethodImplOptions.NoInlining)]
    static void Escape(object obj) {
        VolatileObject = obj;
    }

    static T[][] AllocArray<T>(int n1, int n2) {
        T[][] a = new T[n1][];
        for (int i = 0; i < n1; ++i) {
            a[i] = new T[n2];
        }
        return a;
    }

    static T[][][] AllocArray<T>(int n1, int n2, int n3)
    {
        T[][][] a = new T[n1][][];
        for (int i = 0; i < n1; ++i) {
            a[i] = new T[n2][];
            for (int j = 0; j < n2; j++) {
                a[i][j] = new T[n3];
            }
        }

        return a;
    }

    [MethodImpl(MethodImplOptions.NoInlining)]
    bool Bench()
    {
        m_px = AllocArray<double>(16, 101);
        m_cx = AllocArray<double>(16, 101);

        m_u1 = AllocArray<double>(6, 23, 3);
        m_u2 = AllocArray<double>(6, 23, 3);
        m_u3 = AllocArray<double>(6, 23, 3);

        m_b = AllocArray<double>(65, 9);
        m_c = AllocArray<double>(65, 9);
        m_h = AllocArray<double>(65, 9);

        m_p = AllocArray<double>(5, 513);

        for (int i = 0; i < Iterations; i++) {
           Main1(i < Iterations - 1 ? 0 : 1);
        }

        return true;
    }

    static int Clock()
    {
        return 0;
    }

    void Main1(int output)
    {
        int nt, lw, nl1, nl2;
        int i, i1, i2, ip, ir, ix, j, j1, j2, k, kx, ky, l, m;
        double[] ts = new double[21];
        double[] rt = new double[21];
        double[] rpm = new double[21];
        double[] cksum = new double[21];
        double r, t, a11, a12, a13, sig, a21, a22, a23, a31, a32, a33;
        double b28, b27, b26, b25, b24, b23, b22, c0, flx, rx1;
        double q, s, scale, uu, du1, du2, du3, ar, br, cr, xi, ri;
        int[] mops = new int[20];

        for (i = 1; i <= 20; i++) {
            cksum[i] = 0.0;
        }

        r = 4.86;
        t = 276.0;
        a11 = 0.5;
        a12 = 0.33;
        a13 = 0.25;
        sig = 0.8;
        a21 = 0.20;
        a22 = 0.167;
        a23 = 0.141;
        a31 = 0.125;
        a32 = 0.111;
        a33 = 0.10;
        b28 = 0.1;
        b27 = 0.2;
        b26 = 0.3;
        b25 = 0.4;
        b24 = 0.5;
        b23 = 0.6;
        b22 = 0.7;
        c0 = 0.8;
        flx = 4.689;
        rx1 = 64.0;

        /*
         *     end of initialization -- begin timing
         */

        /* loop 1      hydro excerpt */

        Init();
        ts[1] = (double)Clock();
        q = 0.0;
        for (k = 1; k <= 400; k++) {
            m_x[k] = q + m_y[k] * (r * m_z[k + 10] + t * m_z[k + 11]);
        }
        ts[1] = (double)Clock() - ts[1];
        for (k = 1; k <= 400; k++) {
            cksum[1] += (double)k * m_x[k];
        }

        /* loop 2      mlr, inner product */

        Init();
        ts[2] = (double)Clock();
        q = 0.0;
        for (k = 1; k <= 996; k += 5) {
            q += m_z[k] * m_x[k] + m_z[k + 1] * m_x[k + 1] + m_z[k + 2] * m_x[k + 2] + m_z[k + 3] * m_x[k + 3] + m_z[k + 4] * m_x[k + 4];
        }
        ts[2] = (double)Clock() - ts[2];
        cksum[2] = q;

        /* loop 3      inner prod */

        Init();
        ts[3] = (double)Clock();
        q = 0.0;
        for (k = 1; k <= 1000; k++) {
            q += m_z[k] * m_x[k];
        }
        ts[3] = (double)Clock() - ts[3];
        cksum[3] = q;

        /* loop 4      banded linear equarions */

        Init();
        ts[4] = (double)Clock();
        for (l = 7; l <= 107; l += 50) {
            lw = l;
            for (j = 30; j <= 870; j += 5) {
                m_x[l - 1] -= m_x[lw++] * m_y[j];
            }
            m_x[l - 1] = m_y[5] * m_x[l - 1];
        }
        ts[4] = (double)Clock() - ts[4];
        for (l = 7; l <= 107; l += 50) {
            cksum[4] += (double)l * m_x[l - 1];
        }

        /* loop 5      tri-diagonal elimination, below diagonal */

        Init();
        ts[5] = (double)Clock();
        for (i = 2; i <= 998; i += 3) {
            m_x[i] = m_z[i] * (m_y[i] - m_x[i - 1]);
            m_x[i + 1] = m_z[i + 1] * (m_y[i + 1] - m_x[i]);
            m_x[i + 2] = m_z[i + 2] * (m_y[i + 2] - m_x[i + 1]);
        }
        ts[5] = (double)Clock() - ts[5];
        for (i = 2; i <= 1000; i++) {
            cksum[5] += (double)i * m_x[i];
        }

        /* loop 6      tri-diagonal elimination, above diagonal */

        Init();
        ts[6] = (double)Clock();
        for (j = 3; j <= 999; j += 3) {
            i = 1003 - j;
            m_x[i] = m_x[i] - m_z[i] * m_x[i + 1];
            m_x[i - 1] = m_x[i - 1] - m_z[i - 1] * m_x[i];
            m_x[i - 2] = m_x[i - 2] - m_z[i - 2] * m_x[i - 1];
        }
        ts[6] = (double)Clock() - ts[6];
        for (j = 1; j <= 999; j++) {
            l = 1001 - j;
            cksum[6] += (double)j * m_x[l];
        }

        /* loop 7      equation of state excerpt */

        Init();
        ts[7] = (double)Clock();
        for (m = 1; m <= 120; m++) {
            m_x[m] = m_u[m] + r * (m_z[m] + r * m_y[m]) + t * (m_u[m + 3] + r * (m_u[m + 2] + r * m_u[m + 1]) + t * (m_u[m + 6] + r * (m_u[m + 5] + r * m_u[m + 4])));
        }
        ts[7] = (double)Clock() - ts[7];
        for (m = 1; m <= 120; m++) {
            cksum[7] += (double)m * m_x[m];
        }

        /* loop 8      p.d.e. integration */

        Init();
        ts[8] = (double)Clock();
        nl1 = 1; 
        nl2 = 2;
        for (kx = 2; kx <= 3; kx++) {
            for (ky = 2; ky <= 21; ky++) {
                du1 = m_u1[kx][ky + 1][nl1] - m_u1[kx][ky - 1][nl1];
                du2 = m_u2[kx][ky + 1][nl1] - m_u2[kx][ky - 1][nl1];
                du3 = m_u3[kx][ky + 1][nl1] - m_u3[kx][ky - 1][nl1];
                m_u1[kx][ky][nl2] = m_u1[kx][ky][nl1] + a11 * du1 + a12 * du2 + a13 * du3 + sig * (m_u1[kx + 1][ky][nl1]
                   - 2.0 * m_u1[kx][ky][nl1] + m_u1[kx - 1][ky][nl1]);
                m_u2[kx][ky][nl2] = m_u2[kx][ky][nl1] + a21 * du1 + a22 * du2 + a23 * du3 + sig * (m_u2[kx + 1][ky][nl1]
                   - 2.0 * m_u2[kx][ky][nl1] + m_u2[kx - 1][ky][nl1]);
                m_u3[kx][ky][nl2] = m_u3[kx][ky][nl1] + a31 * du1 + a32 * du2 + a33 * du3 + sig * (m_u3[kx + 1][ky][nl1]
                   - 2.0 * m_u3[kx][ky][nl1] + m_u3[kx - 1][ky][nl1]);
            }
        }
        ts[8] = (double)Clock() - ts[8];
        for (i = 1; i <= 2; i++) {
            for (kx = 2; kx <= 3; kx++) {
                for (ky = 2; ky <= 21; ky++) {
                    cksum[8] += (double)kx * (double)ky * (double)i * (m_u1[kx][ky][i] + m_u2[kx][ky][i] + m_u3[kx][ky][i]);
              }
           }
        }

        /* loop 9      integrate predictors */

        Init();
        ts[9] = (double)Clock();
        for (i = 1; i <= 100; i++) {
            m_px[1][i] = b28 * m_px[13][i] + b27 * m_px[12][i] + b26 * m_px[11][i] + b25 * m_px[10][i] + b24 * m_px[9][i] +
               b23 * m_px[8][i] + b22 * m_px[7][i] + c0 * (m_px[5][i] + m_px[6][i]) + m_px[3][i];
        }
        ts[9] = (double)Clock() - ts[9];
        for (i = 1; i <= 100; i++) {
            cksum[9] += (double)i * m_px[1][i];
        }

        /* loop 10     difference predictors */

        Init();
        ts[10] = (double)Clock();
        for (i = 1; i <= 100; i++) {
            ar = m_cx[5][i];
            br = ar - m_px[5][i];
            m_px[5][i] = ar;
            cr = br - m_px[6][i];
            m_px[6][i] = br;
            ar = cr - m_px[7][i];
            m_px[7][i] = cr;
            br = ar - m_px[8][i];
            m_px[8][i] = ar;
            cr = br - m_px[9][i];
            m_px[9][i] = br;
            ar = cr - m_px[10][i];
            m_px[10][i] = cr;
            br = ar - m_px[11][i];
            m_px[11][i] = ar;
            cr = br - m_px[12][i];
            m_px[12][i] = br;
            m_px[14][i] = cr - m_px[13][i];
            m_px[13][i] = cr;
        }
        ts[10] = (double)Clock() - ts[10];
        for (i = 1; i <= 100; i++) {
            for (k = 5; k <= 14; k++) {
                cksum[10] += (double)k * (double)i * m_px[k][i];
           }
        }

        /* loop 11     first sum. */

        Init();
        ts[11] = (double)Clock();
        m_x[1] = m_y[1];
        for (k = 2; k <= 1000; k++) {
            m_x[k] = m_x[k - 1] + m_y[k];
        }
        ts[11] = (double)Clock() - ts[11];
        for (k = 1; k <= 1000; k++) {
            cksum[11] += (double)k * m_x[k];
        }

        /* loop 12     first diff. */

        Init();
        ts[12] = (double)Clock();
        for (k = 1; k <= 999; k++) {
            m_x[k] = m_y[k + 1] - m_y[k];
        }
        ts[12] = (double)Clock() - ts[12];
        for (k = 1; k <= 999; k++) {
            cksum[12] += (double)k * m_x[k];
        }

        /* loop 13      2-d particle pusher */

        Init();
        ts[13] = (double)Clock();
        for (ip = 1; ip <= 128; ip++) {
            i1 = (int)m_p[1][ip];
            j1 = (int)m_p[2][ip];
            m_p[3][ip] += m_b[i1][j1];
            m_p[4][ip] += m_c[i1][j1];
            m_p[1][ip] += m_p[3][ip];
            m_p[2][ip] += m_p[4][ip];
            // Each element of m_p, m_b and m_c is initialized to 1.00025 in Init().
            // From the assignments above,
            // i2 = m_p[1][ip] = m_p[1][ip] + m_p[3][ip] = m_p[1][ip] + m_p[3][ip] + m_b[i1][j1] = 1 + 1 + 1 = 3
            // j2 = m_p[2][ip] = m_p[2][ip] + m_p[4][ip] = m_p[2][ip] + m_p[4][ip] + m_c[i1][j1] = 1 + 1 + 1 = 3
            i2 = (int)m_p[1][ip];
            j2 = (int)m_p[2][ip];
            // Accessing m_y, m_z upto 35
            m_p[1][ip] += m_y[i2 + 32];
            m_p[2][ip] += m_z[j2 + 32];

            i2 += m_e[i2 + 32];
            j2 += m_f[j2 + 32];
            m_h[i2][j2] += 1.0;
        }
        ts[13] = (double)Clock() - ts[13];
        for (ip = 1; ip <= 128; ip++) {
            cksum[13] += (double)ip * (m_p[3][ip] + m_p[4][ip] + m_p[1][ip] + m_p[2][ip]);
        }
        for (k = 1; k <= 64; k++) {
            for (ix = 1; ix <= 8; ix++) {
                cksum[13] += (double)k * (double)ix * m_h[k][ix];
            }
        }

        /* loop 14      1-d particle pusher */

        Init();
        ts[14] = (double)Clock();
        for (k = 1; k <= 150; k++) {
            // m_grd[150] = 13.636
            // Therefore ix <= 13
            ix = (int)m_grd[k];
            xi = (double)ix;
            m_vx[k] += m_ex[ix] + (m_xx[k] - xi) * m_dex[ix];
            m_xx[k] += m_vx[k] + flx;
            ir = (int)m_xx[k];
            ri = (double)ir;
            rx1 = m_xx[k] - ri;
            ir = System.Math.Abs(ir % 64);
            m_xx[k] = ri + rx1;
            // ir < 64 since ir = ir % 64
            // So m_rh is accessed upto 64
            m_rh[ir] += 1.0 - rx1;
            m_rh[ir + 1] += rx1;
        }
        ts[14] = (double)Clock() - ts[14];
        for (k = 1; k <= 150; k++) {
            cksum[14] += (double)k * (m_vx[k] + m_xx[k]);
        }
        for (k = 1; k <= 67; k++) {
            cksum[14] += (double)k * m_rh[k];
        }

        /* time the clock call */

        ts[15] = (double)Clock();
        ts[15] = (double)Clock() - ts[15];

        /* scale= set to convert time to micro-seconds */

        scale = 1.0;
        rt[15] = ts[15] * scale;

        nt = 14;
        t = s = uu = 0.0;
        for (k = 1; k <= nt; k++) {
            rt[k] = (ts[k] - ts[15]) * scale;
            t += rt[k];
            mops[k] = m_nrops[k] * m_loops[k];
            s += (double)mops[k];
            rpm[k] = 0.0;
            if (rt[k] != 0.0) {
                rpm[k] = (double)mops[k] / rt[k];
            }
            uu += rpm[k];
        }
        uu /= (double)nt;
        s /= t;

        // Ensure that the array elements are live-out
        Escape(ts);
        Escape(rt);
        Escape(rpm);
        Escape(cksum);
        Escape(mops);
    }

    void Init()
    {
        int j, k, l;

        for (k = 1; k <= 1000; k++) {
            m_x[k] = 1.11;
            m_y[k] = 1.123;
            m_z[k] = 0.321;
        }

        for (k = 1; k <= 500; k++) {
            m_u[k] = 0.00025;
        }

        for (k = 1; k <= 15; k++) {
            for (l = 1; l <= 100; l++) {
               m_px[k][l] = l;
               m_cx[k][l] = l;
            }
        }

        for (j = 1; j < 6; j++) {
            for (k = 1; k < 23; k++) {
                for (l = 1; l < 3; l++) {
                    m_u1[j][k][l] = k;
                    m_u2[j][k][l] = k + k;
                    m_u3[j][k][l] = k + k + k;
                }
            }
        }

        for (j = 1; j < 65; j++) {
            for (k = 1; k < 9; k++) {
                m_b[j][k] = 1.00025;
                m_c[j][k] = 1.00025;
                m_h[j][k] = 1.00025;
            }
        }

        for (j = 1; j < 6; j++) {
            m_bnk1[j] = j * 100;
            m_bnk2[j] = j * 110;
            m_bnk3[j] = j * 120;
            m_bnk4[j] = j * 130;
            m_bnk5[j] = j * 140;
        }

        for (j = 1; j < 5; j++) {
            for (k = 1; k < 513; k++) {
                m_p[j][k] = 1.00025;
            }
        }

        for (j = 1; j < 193; j++) {
            m_e[j] = m_f[j] = 1;
        }

        for (j = 1; j < 68; j++) {
            m_ex[j] = m_rh[j] = m_dex[j] = (double)j;
        }

        for (j = 1; j < 151; j++) {
            m_vx[j] = 0.001;
            m_xx[j] = 0.001;
            m_grd[j] = (double)(j / 8 + 3);
        }
    }

    [Benchmark]
    public static void Test() {
        var lloops = new LLoops();
        foreach (var iteration in Benchmark.Iterations) {
            using (iteration.StartMeasurement()) {
                lloops.Bench();
            }
        }
    }

    bool TestBase() {
        bool result = Bench();
        return result;
    }

    public static int Main() {
        var lloops = new LLoops();
        bool result = lloops.TestBase();
        return (result ? 100 : -1);
    }

}
