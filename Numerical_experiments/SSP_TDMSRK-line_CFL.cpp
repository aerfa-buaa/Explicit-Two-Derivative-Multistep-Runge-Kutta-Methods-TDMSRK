#include <math.h>
#include <iostream>
#include <ctime>
#include <cstring>
#include <string>
#include <fstream>
#include <iterator>
#include <valarray>
using namespace std;
class cfda
{

private:
    int w;

public:
    double PI = 3.141592653589793;
    // double PI = 3.1415926535897932385;

    int n;                                         // 网格结点数（有限差分）
    int nb;                                        // 边界网格结点数（有限差分）
    double t, tt, CFL;                             // 计算总时 间t,实际计算时间tt;
    int nt, kt, ij, ijj, onet_out = 2;             // 计算时间步长;
    double nx_begin, nx_end;                       // 网格长度;
    int n1;                                        // 总的网格结点数（有限差分）
    double *Ltt_n0, *x;                            // 声明变长数组
    double *L_n0, *L_n1, *L_n2, *L_n3, *L_n4;      // f=L
    double *Lt_n0, *Lt_n1, *Lt_n2, *Lt_n3, *Lt_n4; // L_t

    double *TL_n0, *TL_n1, *TL_n2, *TL_n3, *TL_n4;      // Tf=L
    double *TLt_n0, *TLt_n1, *TLt_n2, *TLt_n3, *TLt_n4; // TL_t

    double *u_n0, *u_n1, *u_n2, *u_n3, *u_n4, *u_n5, *u_n6, *u_n7, *u_n8, *u_n9; // 推进步中间时刻
    double *f_n0, *f_n1, *f_n2, *f_n3, *f_n4, *f_n5, *f_n6, *f_n7, *f_n8, *f_n9; // 推进步中间时刻
    double *u_nn, *f_nn, *Lttx_n0;                                               // 下一时刻（n+1)
    double *Tu1, *Tu2;
    // double  u_excat[513 + 6 * 2];
    double *df, *u_excat; // du[513 + 6 * 2],df[513 + 6 * 2],
    ////-----------------------------------------
    double dx;
    double dt;
    double *dd; // ss
    double A1, A2, A3, A4, A5, B0, B1, B2, B3, C0, C1, C2, C3, D0, D1, D2, D3, maxuu, TVmax, k, b;
    double a21, a32, a31, aa12, aa21, aa31, aa32, w1, w2, v1, v2, v3, vv1, vv2, vv3, ww1, ww2;
    double put_obj[2][1300];   //
    double put_one_n[2][1300]; // onet_out = 10;
    ////-----------------------------------------
    void intc()
    {
        dx = (nx_end - nx_begin) / (1.0 * (n - 1));

        // for (int i = 0; i < n1; i++)
        // {
        //     x[i] = nx_begin + (i - nb) * dx;
        //     u_n0[i] = intc_fun(x[i]);
        //     df[i] = 1;
        // }

        for (int i = 0; i < n1; i++) // 间断
        {
            x[i] = nx_begin + (i - nb) * dx;
            u_n0[i] = 0.0;
            if (x[i] >= -0.5 && x[i] <= 0.0)
            {
                u_n0[i] = 1.0;
            }
            // df[i] = 1;
        }
        // Write_obj(0);
    }

    double intc_fun(double xxf)
    {
        double yyf;
        // yyf = (sin(2 * PI * xxf) / 2.0 + 0.5) * 1.0; //  / PIbergure
        // yyf = sin(PI * xxf) / 2.0 + 0.25;//bergure
        yyf = sin(PI * xxf); // line
        return yyf;
    }

    // //-----------------------------------------
    void border()
    {
        borderfun(u_n0), borderfun(u_n1), borderfun(u_n2), borderfun(u_n3), borderfun(u_n4), borderfun(u_n5);
        borderfun(u_n6), borderfun(u_n7), borderfun(u_n8), borderfun(u_n9), borderfun(u_nn);
        borderfun(f_n0), borderfun(f_n1), borderfun(f_n2), borderfun(f_n3), borderfun(f_n4), borderfun(f_n5);
        borderfun(f_n6), borderfun(f_n7), borderfun(f_n8), borderfun(f_n9), borderfun(f_nn);
        borderfun(L_n0), borderfun(L_n1), borderfun(L_n2), borderfun(L_n3), borderfun(L_n4);
        borderfun(Lt_n0), borderfun(Lt_n1), borderfun(Lt_n2), borderfun(Lt_n3), borderfun(Lt_n4);
        borderfun(Ltt_n0), borderfun(Lttx_n0);
    }

    void borderfun(double *ffbc)
    {

        for (int i = 0; i < nb; i++)
        {

            *(ffbc + i) = *(ffbc + n1 - 2 * nb + i - 1);  // 周期边条
            *(ffbc + n1 - nb + i) = *(ffbc + nb + i + 1); // 周期边条

            // *(ffbc + i) = 0.;
            // *(ffbc + n1 - nb + i) = 0.;

            // f_n0[i] = f_n0[nb];                   //恒定边条
            // f_n0[n1 - 1 - i] = f_n0[n1 - 1 - nb]; //恒定边条

            // *(ffbc + i) = *(ffbc + nb);                   //恒定边条
            // *(ffbc + n1 - 1 - i) = *(ffbc + n1 - 1 - nb); //恒定边条

            // *(ffbc + i) = 0.0;                   //恒定边条
            // *(ffbc + n1 - 1 - i) = 0.0; //恒定边条
            // ffbc++;
        }
    }

    // //-----------------------------------------
    void carr(double *p1, double *p2, int i)
    {
        while (i-- > 0)
        {
            *p1++ = *p2++;
        }
    }
    // //-----------------------------------------
    void compt_2_stage_int()
    {
        int i, j;
        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput L(n0)
        {
            j = i + nb;
            L_n0[j] = -1. * dfdx(f_n0, j);
            Lt_n0[j] = dfdxx(f_n0, j);

            Tu2[j] = Tu1[j], TL_n1[j] = TL_n0[j], TLt_n1[j] = TLt_n0[j];
            Tu1[j] = u_n0[j], TL_n0[j] = L_n0[j], TLt_n0[j] = Lt_n0[j];
        }
    }

    void compt_SGLM_int()
    {
        int i, j;
        a21 = 0.85; // 0.532  0.468; //0.765;   // 1. / 210. * (93. + sqrt(2139.));
        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput L(n0)
        {
            j = i + nb;
            L_n0[j] = -1. * dfdx(f_n0, j);
            Lt_n0[j] = dfdxx(f_n0, j);
            // Tu1[j] = u_n0[j] + 0.270774947096390 * dt * L_n0[j] + 0.104435936229161 * dt * dt * Lt_n0[j];
            // Tu1[j] = u_n0[j] + 0.527015480634761 * dt * L_n0[j] + 0.027140536548443 * dt * dt * Lt_n0[j];
            // Tu1[j] = u_n0[j] + 0.032910530894256 * dt * L_n0[j] + 0.000541551521871 * dt * dt * Lt_n0[j];
            Tu1[j] = u_n0[j] + 0.093858080521929 * dt * L_n0[j]; //- 0.008889300639912 * dt * dt * Lt_n0[j];
        }
    }

    void compt23_LLttnn_line_SSP()
    {

        int i, j;

        a21 = 0.594223212099088, v2 = 0.306027487008159;
        ww1 = 0.0, ww2 = 0.0, w1 = 0.0, w2 = 0.0;
        v1 = 1. - v2, vv1 = 1. / 6. * (3. - 1. / a21 - 3. * a21 * v2), vv2 = 1. / (6. * a21) - (a21 * v2) / 2.;

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput L(n0)
        {
            j = i + nb;
            L_n0[j] = -dfdx(f_n0, j);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput Lt(n0)
        {
            j = i + nb;
            Lt_n0[j] = dfdxx(f_n0, j);
            // Lt_n0[j] = -1. * dfdxcent(L_n0, j);
            u_n1[j] = u_n0[j] + a21 * dt * L_n0[j] + a21 * a21 / 2.0 * dt * dt * Lt_n0[j];
            f_n1[j] = f_u(u_n1[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n1[j] = -dfdx(f_n1, j);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            Lt_n1[j] = dfdxx(f_n1, j);
            // Lt_n1[j] = -1. * dfdxcent(L_n1, j);

            u_nn[j] = u_n0[j] + dt * (v1 * L_n0[j] + v2 * L_n1[j] + w1 * TL_n0[j] + w2 * TL_n1[j]) +
                      dt * dt * (vv1 * Lt_n0[j] + vv2 * Lt_n1[j] + ww1 * TLt_n0[j] + ww2 * TLt_n1[j]);

            // TL_n1[j] = L_n1[j], TL_n0[j] = L_n0[j];
            // TLt_n1[j] = Lt_n1[j], TLt_n0[j] = Lt_n0[j];
        }
    }

    void compt24_LLttnn_line_SSP()
    {
        int i, j;

        A1 = 0.5, B0 = 1., B1 = 0., C0 = 1. / 3., C1 = 2. / 3.;
        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput L(n0)
        {
            j = i + nb;
            L_n0[j] = -dfdx(f_n0, j);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput Lt(n0)
        {
            j = i + nb;
            Lt_n0[j] = dfdxx(f_n0, j);
            u_n1[j] = u_n0[j] + A1 * dt * L_n0[j] + A1 * A1 / 2.0 * dt * dt * Lt_n0[j];
            f_n1[j] = f_u(u_n1[j]);
        }
        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n1[j] = -dfdx(f_n1, j);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            Lt_n1[j] = dfdxx(f_n1, j);
            u_nn[j] = u_n0[j] + dt * (B0 * L_n0[j] + B1 * L_n1[j]) + dt * dt / 2.0 * (C0 * Lt_n0[j] + C1 * Lt_n1[j]);

            // TL_n1[j] = L_n1[j], TL_n0[j] = L_n0[j];
            // TLt_n1[j] = Lt_n1[j], TLt_n0[j] = Lt_n0[j];
        }
    }

    void compt34_LLttnn_line_SSP()
    {
        int i, j;
        // A1 = 1. / 4., A2 = 1. / 2., B0 = 1. / 3.; //稳定性足够  B0 = 1. / 5.;  B0 = 1. / 3.;
        // A1 = 1. / 4., A2 = 2. / 3., B0 = 1. / 3.;
        A1 = 0.5, A2 = 0.8, B0 = 0.6;
        B1 = -((pow(A2, 3) * (-1. + B0)) / (-pow(A1, 3) + pow(A2, 3)));
        B2 = -((pow(A1, 3) * (-1. + B0)) / (pow(A1, 3) - pow(A2, 3)));
        C0 = -((-1. + 2. * A1 + 2. * A2 - 6. * A1 * A2) / (6. * A1 * A2)) + (A2 * (-pow(A1, 3) + A1 * A2 * A2) * (-1. + B0)) / (-pow(A1, 3) + pow(A2, 3));
        C1 = -((-1. + 2. * A2) / (6. * A1 * (A1 - A2))) + (A1 * pow(A2, 3) * (-1. + B0)) / (-pow(A1, 3) + pow(A2, 3));
        C2 = -((1. - 2. * A1) / (6. * (A1 - A2) * A2)) - (pow(A1, 3) * A2 * (-1. + B0)) / (-pow(A1, 3) + pow(A2, 3));

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput L(n0)
        {
            j = i + nb;
            L_n0[j] = -dfdx(f_n0, j);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput Lt(n0)
        {
            j = i + nb;
            Lt_n0[j] = dfdxx(f_n0, j);
            u_n1[j] = u_n0[j] + A1 * dt * L_n0[j] + A1 * A1 / 2.0 * dt * dt * Lt_n0[j];
            u_n2[j] = u_n0[j] + A2 * dt * L_n0[j] + A2 * A2 / 2.0 * dt * dt * Lt_n0[j];
            f_n1[j] = f_u(u_n1[j]), f_n2[j] = f_u(u_n2[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n1[j] = -dfdx(f_n1, j);
            L_n2[j] = -dfdx(f_n2, j);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            Lt_n1[j] = dfdxx(f_n1, j);
            Lt_n2[j] = dfdxx(f_n2, j);

            u_nn[j] = u_n0[j] + dt * (B0 * L_n0[j] + B1 * L_n1[j] + B2 * L_n2[j]) + dt * dt / 2.0 * (C0 * Lt_n0[j] + C1 * Lt_n1[j] + C2 * Lt_n2[j]);

            // TL_n1[j] = L_n1[j], TL_n0[j] = L_n0[j];
            // TLt_n1[j] = Lt_n1[j], TLt_n0[j] = Lt_n0[j];
            // Tu2[j] = Tu1[j];
            // Tu1[j] = u_nn[j];
        }
    }

    void compt35_LLttnn_line_SSP()
    {
        int i, j;

        a21 = 0.7504;
        aa12 = -1. + 2. * a21;
        vv1 = (1. + 2. * a21 * (-4. + 5. * a21)) / (12. * a21 * (-3. + 5. * a21));
        vv2 = 1. / (12. * a21 * (3. + 10. * (-1. + a21) * a21));
        vv3 = (25. * pow(aa12, 3)) / (12. * (-3. + 5. * a21) * (3. + 10. * (-1 + a21) * a21));
        v1 = 1., v2 = 0., v3 = 0.;
        a31 = (3. - 5. * a21) / (5. - 10. * a21), a32 = 0.;
        aa31 = ((-3. + 5. * a21) * (-3. + 10. * a21) * (1. + 5. * (-1. + a21) * a21)) / (250. * a21 * pow(aa12, 3));
        aa32 = ((-3. + 5. * a21) * (3. + 10. * (-1. + a21) * a21)) / (250. * a21 * pow(aa12, 3));

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput L(n0)
        {
            j = i + nb;
            L_n0[j] = -1. * dfdx(f_n0, j);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput Lt(n0)
        {
            j = i + nb;
            Lt_n0[j] = dfdxx(f_n0, j);
            u_n1[j] = u_n0[j] + a21 * dt * L_n0[j] + a21 * a21 / 2.0 * dt * dt * Lt_n0[j];
            f_n1[j] = f_u(u_n1[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n1[j] = -1. * dfdx(f_n1, j);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            Lt_n1[j] = dfdxx(f_n1, j);
            u_n2[j] = u_n0[j] + dt * (a31 * L_n0[j] + a32 * L_n1[j]) + dt * dt * (aa31 * Lt_n0[j] + aa32 * Lt_n1[j]);
            f_n2[j] = f_u(u_n2[j]);
        }
        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n2[j] = -1. * dfdx(f_n2, j);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            Lt_n2[j] = dfdxx(f_n2, j);

            u_nn[j] = u_n0[j] + dt * (v1 * L_n0[j] + v2 * L_n1[j] + v3 * L_n2[j]) +
                      dt * dt * (vv1 * Lt_n0[j] + vv2 * Lt_n1[j] + vv3 * Lt_n2[j]);

            // TL_n1a[j] = L_n1[j], TL_n1[j] = L_n0[j];
            // TLt_n1a[j] = Lt_n1[j], TLt_n1[j] = Lt_n0[j];
        }
    }

    void compt46_LLttnn_line_SSP()
    {
        int i, j;

        double a41, aa42, aa43, vv4, a42, a43, aa41;
        // a21 = 0.188882174845445, aa41 = 0.106470264959241;
        // a31 = 0.478952613047175, aa42 = 0.069604114078399;
        // a41 = 0.793700870016536, aa43 = 0.138906156494863;
        // v1 = 1., vv1 = 0.059536043418564;
        // aa21 = 0.017838237987173, vv2 = 0.220386636642689;
        // aa31 = 0.001782796208054, vv3 = 0.157700637816082;
        // aa32 = 0.112915006564305, vv4 = 0.062376682122665;

        a21 = 0.235740260201908, aa41 = 0.0337154180289742;
        a31 = 0.766058541677957, aa42 = 0.0725363672383309;
        a41 = 0.504428678818826, aa43 = 0.0209723607400122;
        v1 = 1., vv1 = 0.0740479031154873;
        aa21 = 0.0277867351400587, vv2 = 0.253930692818331;
        aa31 = 0.00713461664453855, vv3 = 0.0765679011859131;
        aa32 = 0.286288227994331, vv4 = 0.0954535028801300;

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput Lt(n0)
        {
            j = i + nb;
            L_n0[j] = -1. * dfdx(f_n0, j);
            Lt_n0[j] = dfdxx(f_n0, j);
            u_n1[j] = u_n0[j] + a21 * dt * L_n0[j] + aa21 * dt * dt * Lt_n0[j];
            f_n1[j] = f_u(u_n1[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n1[j] = -1. * dfdx(f_n1, j);
            Lt_n1[j] = dfdxx(f_n1, j);
            u_n2[j] = u_n0[j] + dt * (a31 * L_n0[j] + a32 * L_n1[j]) + dt * dt * (aa31 * Lt_n0[j] + aa32 * Lt_n1[j]);
            f_n2[j] = f_u(u_n2[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n2[j] = -1. * dfdx(f_n2, j);
            Lt_n2[j] = dfdxx(f_n2, j);
            u_n3[j] = u_n0[j] + dt * (a41 * L_n0[j]) +
                      dt * dt * (aa41 * Lt_n0[j] + aa42 * Lt_n1[j] + aa43 * Lt_n2[j]);
            f_n3[j] = f_u(u_n3[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n3[j] = -1. * dfdx(f_n3, j);
            Lt_n3[j] = dfdxx(f_n3, j);

            u_nn[j] = u_n0[j] + dt * (v1 * L_n0[j]) +
                      dt * dt * (vv1 * Lt_n0[j] + vv2 * Lt_n1[j] + vv3 * Lt_n2[j] + vv4 * Lt_n3[j]);
        }
    }

    // -------------------------------------
    void compt_TDMSRK_223()
    {
        int i, j;
        double b1, b2, d31, d32;

        double arr[14] = {1.3073645988157900E-32, 1.0000000000000000E+00, 9.3046838672164100E-02, 9.0695316132783600E-01,
                          0.0000000000000000E+00, 5.4155814936085300E-01, 8.1533175638285000E-02, 4.9116787560567800E-01,
                          5.2034578742820200E-01, 1.2325951644078300E-32, 1.4664261456957600E-01, 1.2325951644078300E-32,
                          1.3299798286925600E-01, 1.2021427176561400E-01}; // -1.1412144558793900E+00

        d31 = arr[0], d32 = arr[1], b1 = arr[2], b2 = arr[3];
        a31 = arr[4], a32 = arr[5];
        v1 = arr[6], v2 = arr[7], v3 = arr[8];
        aa31 = arr[9], aa32 = arr[10];
        vv1 = arr[11], vv2 = arr[12], vv3 = arr[13];

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput L(n0)
        {
            j = i + nb;
            L_n0[j] = -dfdx(f_n0, j);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput Lt(n0)
        {
            j = i + nb;
            Lt_n0[j] = dfdxx(f_n0, j);
            // u_n1[j] = u_n0[j] + a32 * dt * L_n0[j] + aa32 * dt * dt * Lt_n0[j];
            u_n1[j] = d31 * Tu1[j] + d32 * u_n0[j] + dt * (a31 * TL_n0[j] + a32 * L_n0[j]) + dt * dt * (aa31 * TLt_n0[j] + aa32 * Lt_n0[j]);
            f_n1[j] = f_u(u_n1[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n1[j] = -dfdx(f_n1, j);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            Lt_n1[j] = dfdxx(f_n1, j);
            // Lt_n1[j] = -1. * dfdxcent(L_n1, j);

            u_nn[j] = b1 * Tu1[j] + b2 * u_n0[j] + dt * (v1 * TL_n0[j] + v2 * L_n0[j] + v3 * L_n1[j]) +
                      dt * dt * (vv1 * TLt_n0[j] + vv2 * Lt_n0[j] + vv3 * Lt_n1[j]);

            Tu1[j] = u_n0[j], TL_n0[j] = L_n0[j], TLt_n0[j] = Lt_n0[j];
        }
    }

    void compt_TDMSRK_224()
    {

        int i, j;
        double b1, b2, d31, d32;

        // double arr[14] = {1.1671925590488700E-01, 8.8328074409511300E-01, 9.6290355553452100E-02,
        //                   9.0370964444654800E-01, 1.1163876790256200E-01, 5.6522815642513400E-01, 7.1928985767073800E-02,
        //                   7.2213557460869500E-01, 3.0222579517768400E-01, 7.0736807855637400E-03, 2.0308816438431000E-01,
        //                   1.6624339991858600E-02, 1.3084993645048000E-01, 2.0701845704200000E-01}; // 9.3477573326746200E-01
        // double arr[14] = {0.091664742742202, 0.908335257257798, 0.081130199025293, 0.918869800974707, 0.094617906714575,
        //                   0.566238650266489, 0.060810431478413, 0.718757066125079, 0.301562701421802, 0.001297099946660, 0.209478096095077,
        //                   0.013290683328805, 0.134630471663731, 0.200677155844183}; //,(0.944338173269534)

        double arr[14] = {0.0859182382169719, 0.914081761783028, 0.0776580383519973, 0.922341961648003,
                          0.0907759259244949, 0.566387706162208,
                          0.0583024430067141, 0.717984222294227, 0.301371373051056,
                          -7.18464639633728E-17, 0.210977456824701,
                          0.0125444208848960, 0.135504000944643, 0.199267993301564};//-0.946487048652487
        d31 = arr[0], d32 = arr[1], b1 = arr[2], b2 = arr[3];
        a31 = arr[4], a32 = arr[5];
        v1 = arr[6], v2 = arr[7], v3 = arr[8];
        aa31 = arr[9], aa32 = arr[10];
        vv1 = arr[11], vv2 = arr[12], vv3 = arr[13];

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput L(n0)
        {
            j = i + nb;
            L_n0[j] = -dfdx(f_n0, j);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput Lt(n0)
        {
            j = i + nb;
            Lt_n0[j] = dfdxx(f_n0, j);
            u_n1[j] = d31 * Tu1[j] + d32 * u_n0[j] + dt * (a31 * TL_n0[j] + a32 * L_n0[j]) + dt * dt * (aa31 * TLt_n0[j] + aa32 * Lt_n0[j]);

            f_n1[j] = f_u(u_n1[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n1[j] = -dfdx(f_n1, j);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            Lt_n1[j] = dfdxx(f_n1, j);

            u_nn[j] = b1 * Tu1[j] + b2 * u_n0[j] + dt * (v1 * TL_n0[j] + v2 * L_n0[j] + v3 * L_n1[j]) +
                      dt * dt * (vv1 * TLt_n0[j] + vv2 * Lt_n0[j] + vv3 * Lt_n1[j]);

            // u_nn[j] = u_n0[j] + dt * (L_n0[j]) + 0.5 * dt * dt * (Lt_n0[j]);

            Tu1[j] = u_n0[j], TL_n0[j] = L_n0[j], TLt_n0[j] = Lt_n0[j];
        }
    }

    void compt_TDMSRK_235()
    {
        int i, j;
        double d31, d32, d33, b1, b2, b3, a33, v4, aa33, vv4;

        // double arr[20] = {3.0423035627486100E-02, 1.6806131761218300E-01, 8.0151564676033100E-01, 3.8790763256770900E-02,
        //                   1.3277234796796400E-01, 8.2843688877526500E-01, 1.2982328267980100E-02, 2.0203577550206300E-01,
        //                   5.8101105840906700E-01, 3.4236188650271200E-02, 1.0616396753415600E-01, 7.9489851109331000E-01,
        //                   2.7505520720376800E-01, 1.4003994043627200E-02, 0.0000000000000000E+00, 2.2993326081558200E-01,
        //                   7.3586945940719000E-03, 3.2126990022231500E-02, 1.2082329070476400E-01, 2.1436987214802000E-01}; // -8.3183939673330900E-01

        double arr[20] = {0.014767605406112, 0.167600981747149, 0.817631412846739, 0.030575274114486, 0.133045703046640,
                          0.836379022838873, 0.000013886700706, 0.199211349061870, 0.569735915417622, 0.027032987291601, 0.105732896401668,
                          0.780704117613404, 0.280726249968939, 0.010423447853488, 0.000000000000000, 0.238971864574849, 0.005532322014941,
                          0.031144942360066, 0.126836083983140, 0.208085846602162}; //,(0.841322457462487)

        d31 = arr[0],
        d32 = arr[1], d33 = arr[2], b1 = arr[3], b2 = arr[4], b3 = arr[5];
        a31 = arr[6], a32 = arr[7], a33 = arr[8];
        v1 = arr[9], v2 = arr[10], v3 = arr[11], v4 = arr[12];
        aa31 = arr[13], aa32 = arr[14], aa33 = arr[15];
        vv1 = arr[16], vv2 = arr[17], vv3 = arr[18], vv4 = arr[19];

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput L(n0)
        {
            j = i + nb;
            L_n0[j] = -dfdx(f_n0, j);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput Lt(n0)
        {
            j = i + nb;
            Lt_n0[j] = dfdxx(f_n0, j);
            u_n1[j] = d31 * Tu2[j] + d32 * Tu1[j] + d33 * u_n0[j] + dt * (a31 * TL_n1[j] + a32 * TL_n0[j] + a33 * L_n0[j]) +
                      dt * dt * (aa31 * TLt_n1[j] + aa32 * TLt_n0[j] + aa33 * Lt_n0[j]);
            f_n1[j] = f_u(u_n1[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n1[j] = -dfdx(f_n1, j);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            Lt_n1[j] = dfdxx(f_n1, j);

            u_nn[j] = b1 * Tu2[j] + b2 * Tu1[j] + b3 * u_n0[j] + dt * (v1 * TL_n1[j] + v2 * TL_n0[j] + v3 * L_n0[j] + v4 * L_n1[j]) +
                      dt * dt * (vv1 * TLt_n1[j] + vv2 * TLt_n0[j] + vv3 * Lt_n0[j] + vv4 * Lt_n1[j]);

            // u_nn[j] = u_n0[j] + dt * (L_n0[j]) + 0.5 * dt * dt * (Lt_n0[j]);
            Tu2[j] = Tu1[j], TL_n1[j] = TL_n0[j], TLt_n1[j] = TLt_n0[j];
            Tu1[j] = u_n0[j], TL_n0[j] = L_n0[j], TLt_n0[j] = Lt_n0[j];
        }
    }

    void compt_TDMSRK_326()
    {
        int i, j;

        double b1, b2, d31, d32, d41, d42, a41, a42, a43, aa41, aa42, aa43, v4, vv4;

        double arr[24] = {3.4192786601275000E-01, 9.2331503301445700E-02, 6.5807213398725000E-01, 9.0766849669855400E-01,
                          1.0426387973076500E-01, 8.9573612026923500E-01, 2.5790165593231500E-01, 4.9998434064656400E-01,
                          6.6076932775015300E-02, 7.2129498199869600E-01, 0.0000000000000000E+00, 7.5009584555929800E-02,
                          6.9357361470979600E-01, 4.4055647477858000E-33, 3.3568068046503900E-01, 5.9967923477784500E-02,
                          1.1348038264019000E-01, 1.2601479587349600E-02, 1.2645394691931600E-01, 1.2239634140739500E-01,
                          1.3819314026467700E-02, 1.3801717989153900E-01, 5.8941699256348000E-02, 7.8787813242481900E-02}; // 9.2651413378796700E-01
        d31 = arr[0], d41 = arr[1], d32 = arr[2], d42 = arr[3], b1 = arr[4], b2 = arr[5];
        a31 = arr[6], a32 = arr[7], a41 = arr[8], a42 = arr[9], a43 = arr[10];
        v1 = arr[11], v2 = arr[12], v3 = arr[13], v4 = arr[14];
        aa31 = arr[15], aa32 = arr[16], aa41 = arr[17], aa42 = arr[18], aa43 = arr[19];
        vv1 = arr[20], vv2 = arr[21], vv3 = arr[22], vv4 = arr[23];

        i = n1 - nb * 2;
        border();
        while (i-- > 0) // comput Lt(n0)
        {
            j = i + nb;
            L_n0[j] = -dfdx(f_n0, j);
            Lt_n0[j] = dfdxx(f_n0, j);
            u_n1[j] = d31 * Tu1[j] + d32 * u_n0[j] + dt * (a31 * TL_n0[j] + a32 * L_n0[j]) +
                      dt * dt * (aa31 * TLt_n0[j] + aa32 * Lt_n0[j]);
            f_n1[j] = f_u(u_n1[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n1[j] = -dfdx(f_n1, j);
            Lt_n1[j] = dfdxx(f_n1, j);
            u_n2[j] = d41 * Tu1[j] + d42 * u_n0[j] + dt * (a41 * TL_n0[j] + a42 * L_n0[j] + a43 * L_n1[j]) +
                      dt * dt * (aa41 * TLt_n0[j] + aa42 * Lt_n0[j] + aa43 * Lt_n1[j]);

            f_n2[j] = f_u(u_n2[j]);
        }

        i = n1 - nb * 2;
        border();
        while (i-- > 0)
        {
            j = i + nb;
            L_n2[j] = -dfdx(f_n2, j);
            Lt_n2[j] = dfdxx(f_n2, j);
            u_nn[j] = b1 * Tu1[j] + b2 * u_n0[j] + dt * (v1 * TL_n0[j] + v2 * L_n0[j] + v3 * L_n1[j] + v4 * L_n2[j]) +
                      dt * dt * (vv1 * TLt_n0[j] + vv2 * Lt_n0[j] + vv3 * Lt_n1[j] + vv4 * Lt_n2[j]);

            Tu1[j] = u_n0[j], TL_n0[j] = L_n0[j], TLt_n0[j] = Lt_n0[j];
        }
    }



    //------------------------------------------------
    void f_eq_u(double *ff, double *uu)
    {
        for (int i = 0; i < n1; i++)
        {

            *ff++ = f_u(*uu++);
        }
    }
    double f_u(double xx)
    {
        double y;
        // y = xx * xx / 2.0;//bergure
        y = xx; // line
        return y;
    }
    void Store_obj(double *ss, double *ff)
    {
        for (int ik = 0; ik < n1; ik++) //
        {
            *ss++ = *ff++;
        }
    }

    void Write_obj(int j)
    {
        if (j < -0.9)
        {
            // string title = "result_t_burgers_excat"; //+ to_string(t);
            string title = to_string(n);
            string Title = title + ".plt";
            ofstream ofs;
            ofs.open(Title, ios::out);
            ofs << " VARIABLES=x,u " << endl;

            for (int i = nb; i < n1 - nb; i++)
            {

                ofs.precision(10);
                ofs << x[i];
                for (int j = 0; j < kt; j++)
                {
                    ofs << "  ";
                    ofs.precision(18);
                    ofs << put_obj[j][i];
                }

                ofs << endl;
            }
            ofs.close();
        }

        if (j < 0.5 && j > -0.8)
        {
            // string title = "result_t_burgers_excat"; //+ to_string(t);
            string title = to_string(n);
            string Title = "excat.txt";
            ofstream ofs;
            ofs.open(Title, ios::out);
            ofs << "Title=" << title << endl;

            for (int i = nb; i < n1 - nb; i++)
            {

                ofs.precision(10);
                ofs << x[i];
                ofs << "  ";
                ofs.precision(18);
                ofs << *(u_excat + i);
                ofs << endl;
            }
            ofs.close();
        }

        if (j > 1.5)
        {
            string title = "result_one_t_burgers_every-" + to_string(j);
            string Title = title + ".plt";
            ofstream ofs;
            ofs.open(Title, ios::out);
            ofs << "Title=" << title << endl;
            for (int i = nb; i < n1 - nb; i++)
            {

                ofs.precision(10);
                ofs << x[i];
                for (j = 0; j < onet_out; j++)
                {
                    ofs << "  ";
                    ofs.precision(18);
                    ofs << put_one_n[j][i];
                }

                ofs << endl;
            }
            ofs.close();
        }
    }
    ////-----------------------------------------

    void RK33_compt()
    {
        // compt_t1(cff1.uu_t1, cff1.uu_t0,n_all);
        // uu_t1 = uu_t0 + dt * df;  df_dx(double *g, int i)
        int j, i;
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            df[nb + j] = -dfdx(f_n0, nb + j); //-uu_t0[nb + j];
            u_n1[nb + j] = u_n0[nb + j] + dt * df[nb + j];
            f_n1[nb + j] = f_u(u_n1[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            df[nb + j] = -dfdx(f_n1, nb + j); //-uu_t1[nb + j];
            u_n2[nb + j] = u_n0[nb + j] * 0.75 + 0.25 * u_n1[nb + j] + 0.25 * dt * df[nb + j];
            f_n2[nb + j] = f_u(u_n2[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            df[nb + j] = -dfdx(f_n2, nb + j);
            u_nn[nb + j] = u_n0[nb + j] * 1.0 / 3.0 + u_n2[nb + j] * 2.0 / 3.0 + 2.0 / 3.0 * dt * df[nb + j];
            j++;
        }
    }

    void SSPRK54_compt()
    {
        int j, i;
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            df[nb + j] = -dfdx(f_n0, nb + j); //-uu_t0[nb + j];
            u_n1[nb + j] = u_n0[nb + j] + 0.39175222700392 * dt * df[nb + j];
            f_n1[nb + j] = f_u(u_n1[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            df[nb + j] = -dfdx(f_n1, nb + j); //-uu_t1[nb + j];
            u_n2[nb + j] = 0.44437049406734 * u_n0[nb + j] + 0.55562950593266 * u_n1[nb + j] + 0.36841059262959 * dt * df[nb + j];
            f_n2[nb + j] = f_u(u_n2[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            df[nb + j] = -dfdx(f_n2, nb + j); //-uu_t1[nb + j];
            u_n3[nb + j] = 0.62010185138540 * u_n0[nb + j] + 0.37989814861460 * u_n2[nb + j] + 0.25189177424738 * dt * df[nb + j];
            f_n3[nb + j] = f_u(u_n3[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            df[nb + j] = -dfdx(f_n3, nb + j); //-uu_t1[nb + j];
            Ltt_n0[nb + j] = df[nb + j];
            u_n4[nb + j] = 0.17807995410773 * u_n0[nb + j] + 0.82192004589227 * u_n3[nb + j] + 0.54497475021237 * dt * df[nb + j];
            f_n4[nb + j] = f_u(u_n4[nb + j]);
            j++;
        }

        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            df[nb + j] = -dfdx(f_n4, nb + j);
            u_nn[nb + j] = 0.00683325884039 * u_n0[nb + j] + 0.51723167208978 * u_n2[nb + j] + 0.12759831133288 * u_n3[nb + j] +
                           0.34833675773694 * u_n4[nb + j] + 0.08460416338212 * dt * Ltt_n0[nb + j] + 0.22600748319395 * dt * df[nb + j];
            j++;
        }
    }

    void SSPRK54_compt_Shu()
    {
        int j, i;
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n0[nb + j] = -dfdx(f_n0, nb + j); //-uu_t0[nb + j];
            u_n1[nb + j] = u_n0[nb + j] + 0.391752226571890 * dt * L_n0[nb + j];
            f_n1[nb + j] = f_u(u_n1[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n1[nb + j] = -dfdx(f_n1, nb + j); //-uu_t1[nb + j];
            u_n2[nb + j] = 0.444370493651235 * u_n0[nb + j] + 0.555629506348765 * u_n1[nb + j] + 0.368410593050371 * dt * L_n1[nb + j];
            f_n2[nb + j] = f_u(u_n2[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n2[nb + j] = -dfdx(f_n2, nb + j); //-uu_t1[nb + j];
            u_n3[nb + j] = 0.620101851488403 * u_n0[nb + j] + 0.379898148511597 * u_n2[nb + j] + 0.251891774271694 * dt * L_n2[nb + j];
            f_n3[nb + j] = f_u(u_n3[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n3[nb + j] = -dfdx(f_n3, nb + j); //-uu_t1[nb + j];
            u_n4[nb + j] = 0.178079954393132 * u_n0[nb + j] + 0.821920045606868 * u_n3[nb + j] + 0.544974750228521 * dt * L_n3[nb + j];
            f_n4[nb + j] = f_u(u_n4[nb + j]);
            j++;
        }

        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            df[nb + j] = -dfdx(f_n4, nb + j);
            u_nn[nb + j] = 0.517231671970585 * u_n2[nb + j] + 0.096059710526147 * u_n3[nb + j] + 0.063692468666290 * dt * L_n3[nb + j] +
                           0.386708617503269 * u_n4[nb + j] + 0.226007483236906 * dt * df[nb + j];

            j++;
        }
    }

    void RK65LS_compt2() // Numerical Methods for Ordinary Differential Equations
    {
        int j, i;
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n0[nb + j] = -dfdx(f_n0, nb + j); //-uu_t0[nb + j];
            u_n1[nb + j] = u_n0[nb + j] + 1. / 4. * dt * L_n0[nb + j];
            f_n1[nb + j] = f_u(u_n1[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n1[nb + j] = -dfdx(f_n1, nb + j); //-uu_t1[nb + j];
            u_n2[nb + j] = u_n0[nb + j] + 1. / 8. * dt * L_n0[nb + j] + 1. / 8. * dt * L_n1[nb + j];
            f_n2[nb + j] = f_u(u_n2[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n2[nb + j] = -dfdx(f_n2, nb + j); //-uu_t1[nb + j];
            u_n3[nb + j] = u_n0[nb + j] + 0.0 * dt * L_n0[nb + j] + 0.0 * dt * L_n1[nb + j] +
                           1. / 2. * dt * L_n2[nb + j];
            f_n3[nb + j] = f_u(u_n3[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n3[nb + j] = -dfdx(f_n3, nb + j); //-uu_t1[nb + j];
            u_n4[nb + j] = u_n0[nb + j] + 3. / 16. * dt * L_n0[nb + j] - 3. / 8. * dt * L_n1[nb + j] +
                           3. / 8. * dt * L_n2[nb + j] + 9. / 16. * dt * L_n3[nb + j];
            f_n4[nb + j] = f_u(u_n4[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n4[nb + j] = -dfdx(f_n4, nb + j); //-uu_t1[nb + j];
            u_n5[nb + j] = u_n0[nb + j] - 3. / 7. * dt * L_n0[nb + j] + 8. / 7. * dt * L_n1[nb + j] + 6. / 7. * dt * L_n2[nb + j] -
                           12. / 7. * dt * L_n3[nb + j] + 8. / 7. * dt * L_n4[nb + j];
            f_n5[nb + j] = f_u(u_n5[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            df[nb + j] = -dfdx(f_n5, nb + j);
            u_nn[nb + j] = u_n0[nb + j] + 7. / 90. * dt * L_n0[nb + j] + 0. * dt * L_n1[nb + j] +
                           16. / 45. * dt * L_n2[nb + j] + 2. / 15. * dt * L_n3[nb + j] +
                           16. / 45. * dt * L_n4[nb + j] + 7. / 90. * dt * df[nb + j];
            j++;
        }
    }

    void RK76LS_compt()
    {
        int j, i;
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n0[nb + j] = -dfdx(f_n0, nb + j); //-uu_t0[nb + j];
            u_n1[nb + j] = u_n0[nb + j] + 1. / 6. * dt * L_n0[nb + j];
            f_n1[nb + j] = f_u(u_n1[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n1[nb + j] = -dfdx(f_n1, nb + j); //-uu_t1[nb + j];
            u_n2[nb + j] = u_n0[nb + j] + 12. / 169. * dt * L_n0[nb + j] + 27. / 169. * dt * L_n1[nb + j];
            f_n2[nb + j] = f_u(u_n2[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n2[nb + j] = -dfdx(f_n2, nb + j); //-uu_t1[nb + j];
            u_n3[nb + j] = u_n0[nb + j] + 107952. / 571787. * dt * L_n0[nb + j] - 406107. / 571787. * dt * L_n1[nb + j] + 566826. / 571787. * dt * L_n2[nb + j];
            f_n3[nb + j] = f_u(u_n3[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n3[nb + j] = -dfdx(f_n3, nb + j); //-uu_t1[nb + j];
            u_n4[nb + j] = u_n0[nb + j] + 1371923. / 6669000. * dt * L_n0[nb + j] - 819. / 3800. * dt * L_n1[nb + j] +
                           19411847. / 88236000. * dt * L_n2[nb + j] + 561846173. / 1147068000. * dt * L_n3[nb + j];
            f_n4[nb + j] = f_u(u_n4[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            L_n4[nb + j] = -dfdx(f_n4, nb + j); //-uu_t1[nb + j];
            u_n5[nb + j] = u_n0[nb + j] - 1563412. / 5835375. * dt * L_n0[nb + j] + 468. / 475. * dt * L_n1[nb + j] - 14488201. / 168199875. * dt * L_n2[nb + j] -
                           1711096709. / 6846562125. * dt * L_n3[nb + j] + 648832. / 1549583. * dt * L_n4[nb + j];
            f_n5[nb + j] = f_u(u_n5[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            Lt_n1[nb + j] = -dfdx(f_n5, nb + j); //-uu_t1[nb + j];
            u_n6[nb + j] = u_n0[nb + j] + 120115. / 277641. * dt * L_n0[nb + j] - 117. / 113. * dt * L_n1[nb + j] + 219237109. / 296102601. * dt * L_n2[nb + j] +
                           29855183083. / 44628054003. * dt * L_n3[nb + j] - 3009600. / 9215941. * dt * L_n4[nb + j] + 297825. / 572797. * dt * Lt_n1[nb + j];
            f_n6[nb + j] = f_u(u_n6[nb + j]);
            j++;
        }
        i = n1 - nb * 2, j = 0;
        border();
        while (i-- > 0)
        {
            df[nb + j] = -dfdx(f_n6, nb + j);
            u_nn[nb + j] = u_n0[nb + j] + 137. / 1872. * dt * L_n0[nb + j] + 371293. / 1164612. * dt * L_n2[nb + j] +
                           3939040643. / 23169727152. * dt * L_n3[nb + j] + 19000. / 104859. * dt * L_n4[nb + j] +
                           45125. / 243312. * dt * Lt_n1[nb + j] + 113. / 1584. * dt * df[nb + j];
            j++;
        }
    }

    //---------------------------------------------------------

    int Perimeter()
    {
        cout << " yes ";
        return 2 * t;
    }

    double dfdx(double *g, int i)
    {
        double y;
        // 1up
        y = (*(g + i) - *(g + i - 1)) / dx;

        return y;
    }

    double dfdxx(double *g, int i)
    {
        double y;
        // 2ord
        y = (*(g + i + 1) - *(g + i) * 2.0 + *(g + i - 1)) / dx / dx;
        // y = (*(g + i ) - *(g + i - 1) * 2.0 + *(g + i-2)) / dx / dx;
        return y;
    }

    void compt_CFL()
    {
        double *umax;
        // umax = max_element(u_n0 + 1, u_n0 + n1 - 1);
        // cout << *umax;
        // dt = CFL * dx / (*umax);
        // dt=0.02/1.0;
        dt = CFL * dx;
    }

    void maxTV()
    {
        double TV = 0.0;
        for (int i = 0; i < n1 - 1; i++) // 间断
        {
            TV = abs(u_nn[i + 1] - u_nn[i]) + TV;
        }
        TV = abs(TV - 2.0);
        // cout << " \n TV: " << TV;
        if (tt > 0.0)
        {
            if (TV > TVmax)
            {
                TVmax = TV;
            }
        }
    }

    void TVabs(int m)
    {
        double TV = 0.0;
        TV = TVmax;

        if (TV < 2.0E-16)
        {
            TV = 2.0E-16;
        }
        cout << "\n  CFL: " << CFL << ",  TV: " << TV << "   \n";
        string Title = "TV_CFL.plt";
        ofstream ofs;
        if (m < 2)
        {
            ofs.open(Title, ios::out);
            ofs << " VARIABLES=CFL,TV " << endl;
            ofs.precision(4), ofs << CFL, ofs << "  ";
            ofs.precision(10), ofs << TV << endl;
            ofs.close();
        }
        else
        {
            ofs.open(Title, ios::app);
            ofs.precision(4), ofs << CFL, ofs << "  ";
            ofs.precision(10), ofs << TV << endl;
            ofs.close();
        }
    }
};
////-----------------------------------------
void comput_tt(cfda &);
int comput_main(int, double, double, int);
int main()
{
    // int n, oder_n = 2, begin_cell = 1600;
    int n, oder_n = 300, begin_cell = 3200;
    double t = 100.0, CFL = 0.95; // 计算时间总长，CFL数  CFL = 0.66 1.45
    n = begin_cell + 1;
    for (int m = 1; m < oder_n; m++)
    {
        comput_main(n, t, CFL, m);
        CFL = CFL + 0.001;
    }
    // system("pause");
    return 2;
}
int comput_main(int n, double t, double CFL, int m)
{
    int nb = 6; // 边界网格结点数（有限差分）
    int n1 = n + 6 * 2;
    double main_df[n + 6 * 2], u_excat[n + 6 * 2];
    double Ltt_n0[n + 6 * 2], x[n + 6 * 2];                                                          // 声明变长数组
    double L_n0[n + 6 * 2], L_n1[n + 6 * 2], L_n2[n + 6 * 2], L_n3[n + 6 * 2], L_n4[n + 6 * 2];      // f=L
    double Lt_n0[n + 6 * 2], Lt_n1[n + 6 * 2], Lt_n2[n + 6 * 2], Lt_n3[n + 6 * 2], Lt_n4[n + 6 * 2]; // L_t
    double u_n0[n + 6 * 2], u_n1[n + 6 * 2], u_n2[n + 6 * 2], u_n3[n + 6 * 2], u_n4[n + 6 * 2];      // 推进步中间时刻
    double u_n5[n + 6 * 2], u_n6[n + 6 * 2], u_n7[n + 6 * 2], u_n8[n + 6 * 2], u_n9[n + 6 * 2];      // 推进步中间时刻
    double f_n0[n + 6 * 2], f_n1[n + 6 * 2], f_n2[n + 6 * 2], f_n3[n + 6 * 2], f_n4[n + 6 * 2];      // 推进步中间时刻
    double f_n5[n + 6 * 2], f_n6[n + 6 * 2], f_n7[n + 6 * 2], f_n8[n + 6 * 2], f_n9[n + 6 * 2];      // 推进步中间时刻
    double u_nn[n + 6 * 2], f_nn[n + 6 * 2], Lttx_n0[n + 6 * 2];                                     // 下一时刻（n+1)

    double TL_n0[n + 6 * 2], TL_n1[n + 6 * 2], TL_n2[n + 6 * 2], TL_n3[n + 6 * 2], TL_n4[n + 6 * 2];
    double TLt_n0[n + 6 * 2], TLt_n1[n + 6 * 2], TLt_n2[n + 6 * 2], TLt_n3[n + 6 * 2], TLt_n4[n + 6 * 2];
    double Tu1[n + 6 * 2], Tu2[n + 6 * 2];
    class cfda cff;
    cff.kt = 1;
    cff.A1 = 1 / 3.0;
    cff.A2 = 2 / 3.0;
    cff.t = t;           // 1.0 / cff.PI; //计算总时间
    cff.nx_begin = -0.7; // 网格右端点
    cff.nx_end = 1.3;    // 网格左端点
    cff.CFL = CFL, cff.TVmax = 0.0;

    cff.n = n, cff.nb = nb, cff.n1 = n1;
    cff.Ltt_n0 = Ltt_n0, cff.x = x;
    cff.u_n0 = u_n0, cff.u_n1 = u_n1, cff.u_n2 = u_n2, cff.u_n3 = u_n3, cff.u_n4 = u_n4;                     // 推进步中间时刻
    cff.u_n5 = u_n5, cff.u_n6 = u_n6, cff.u_n7 = u_n7, cff.u_n8 = u_n8, cff.u_n9 = u_n9;                     // 推进步中间时刻
    cff.f_n0 = f_n0, cff.f_n1 = f_n1, cff.f_n2 = f_n2, cff.f_n3 = f_n3, cff.f_n4 = f_n4;                     // 推进步中间时刻
    cff.f_n5 = f_n5, cff.f_n6 = f_n6, cff.f_n7 = f_n7, cff.f_n8 = f_n8, cff.f_n9 = f_n9;                     // 推进步中间时刻
                                                                                                             // 声明变长数组
    cff.u_nn = u_nn, cff.f_nn = f_nn, cff.Lttx_n0 = Lttx_n0;                                                 // 下一时刻（n+1)
    cff.L_n0 = L_n0, cff.L_n1 = L_n1, cff.L_n2 = L_n2, cff.L_n3 = L_n3, cff.L_n4 = L_n4;                     // f=L
    cff.Lt_n0 = Lt_n0, cff.Lt_n1 = Lt_n1, cff.Lt_n2 = Lt_n2, cff.Lt_n3 = Lt_n3, cff.Lt_n4 = Lt_n4;           // L_t
                                                                                                             // 存储一时刻（n-1)
    cff.TL_n0 = TL_n0, cff.TL_n1 = TL_n1, cff.TL_n2 = TL_n2, cff.TL_n3 = TL_n3, cff.TL_n4 = TL_n4;           // f=TL
    cff.TLt_n0 = TLt_n0, cff.TLt_n1 = TLt_n1, cff.TLt_n2 = TLt_n2, cff.TLt_n3 = TLt_n3, cff.TLt_n4 = TLt_n4; // TL_t

    cff.Tu1 = Tu1, cff.Tu2 = Tu2;
    cff.df = main_df, cff.u_excat = u_excat, cff.intc();
    for (cff.ij = 0; cff.ij < cff.kt; cff.ij++)
    {
        // cff.nt = 1000 * pow(2, cff.ij); //计算时间步数

        comput_tt(cff);

        cff.Store_obj(cff.put_obj[cff.ij], cff.u_nn); // 存储数据
        // cff.Write_obj(cff.ij + 1);                    //存储数据
    }
    cff.Write_obj(-1);
    cff.TVabs(m);
    return 2;
}
////-----------------------------------------
void comput_tt(cfda &cff1)
{
    int n_all, i3 = 0, tt_flag = 0;
    cff1.tt = 0, n_all = cff1.n1, cff1.border();
    cff1.f_eq_u(cff1.f_n0, cff1.u_n0), cff1.TVmax = 0.0;
    for (int i1 = 0; i1 < 100; i1++) // 500 200 50
    {
        cff1.compt_CFL();
        cff1.tt = cff1.tt + cff1.dt;
        if (i1 < 20000)
        {
            // if (i1 % 2 == 0)
            //     cff1.compt_2_stage_int(); //   cff1.compt_SGLM_int(); //  cff1.compt_2_stage_int();

            // cff1.dt = cff1.dt * 0.5;

            // cff1.compt23_LLttnn_line_SSP();
            // cff1.compt24_LLttnn_line_SSP(); // cff1.compt34_LLttnn_line_SSP();
            // cff1.compt35_LLttnn_line_SSP();
            // cff1.compt46_LLttnn_line_SSP();

            // cff1.compt_SGLM5_line_int();
            // cff1.compt_SGLM3_line_int();

            // cff1.RK33_compt(); // 3-RK
            // cff1.SSPRK54_compt_Shu();// cff1.SSPRK54_compt(); //
            // cff1.RK65LS_compt2();// cff1.RK65LS_compt();
            cff1.RK76LS_compt();// cff1.RK65LS_compt();
        }
        else
        {
            // cff1.compt_TDMSRK_223();
            // cff1.compt_TDMSRK_224();
            // cff1.compt_TDMSRK_235();
            cff1.compt_TDMSRK_326();
        }

        //___________________________________________________________________
        if (cff1.tt > cff1.t)
        {
            cff1.f_eq_u(cff1.f_nn, cff1.u_nn);
            cout.precision(18); // 精度为18，正常为6
            cout << " \n comput time: " << cff1.tt - cff1.dt << " comput time n:   " << i1;
            break;
        }
        cff1.carr(cff1.u_n0, cff1.u_nn, n_all);
        cff1.f_eq_u(cff1.f_n0, cff1.u_nn);

        // cout << " \n comput time: " << cff1.tt << " comput time n:   " << i1;

        cff1.border();
        cff1.maxTV();

        if (i3 * cff1.t / cff1.onet_out < cff1.tt + 0.00000001)
        {
            // cff1.Store_obj(cff1.put_one_n[i3], cff1.u_nn); //存储数据
            i3++;
        }
    }
}
