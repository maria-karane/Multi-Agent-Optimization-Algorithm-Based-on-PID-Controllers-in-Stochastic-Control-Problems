using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Управление
{
    public interface ITasks
    {
        int n();
        double t0();
        double t1();

        //для пучков
        double[] x0_a();
        double[] x1_a();
        int numb_tr();


        int[] LL(); 

        void setPosition(double[] x);

        double x_delta();

        //
        double c_down();
        double c_up();
        //
        //double u_cheb_min();
        //double u_cheb_max();

        double x0(int i);
        double v(int i, double u, double[] x);

        double[] Solve(U upr);

        double[] SolveT(U upr);

        void setUprType(double u0, bool f_basis, bool view_cos, bool chebyshev, int L0, int[] x_L); //, int fst = 1);

        void setIcount(bool I_count, bool I_count2);
        bool getIcount();
        bool getIcount2();

        bool getUprType();

        bool getViewCos();
        bool getChebyshev();

        //int getFST();

        int[] get_x_L();

        double[] x1();
        double[] x0();

        double I(double[] rx);



        double Dt();

        double U0();

        double U1();
        double U2();
    }
}
