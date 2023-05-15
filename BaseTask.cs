using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Управление
{
    public abstract class BaseTask : ITasks
    {

        public double NextGaussian(Random r, double mu = 0, double sigma = 1)
        {
            var u1 = r.NextDouble();
            var u2 = r.NextDouble();

            var rand_std_normal = Math.Sqrt(-2.0 * Math.Log(u1)) *
                                Math.Sin(2.0 * Math.PI * u2);

            var rand_normal = mu + sigma * rand_std_normal;

            return rand_normal;
        }

        public abstract int n();
        public abstract  double t0();
        public abstract  double t1();

        protected double[] _x0_a, _x1_a, position;
        protected int[] _LL;

        public double[] x0_a() { return _x0_a; }
        public double[] x1_a() { return _x1_a; }

        public double[] _x0, _x1;

        public double[] x0() { return _x0; }
        public double[] x1() { return _x1; }

        public int numb_tr()
        { return _x0_a.Length; }



        public int[] LL() { return _LL; }

        public void SetLL(int L0, int L1)
        {
            _LL = new int[2];
            _LL[0] = L0;
            _LL[1] = L1;
        }

        public void setPosition(double[] x) 
        {
            position = new double[x.Length];
            for (int i = 0; i < x.Length; i++)
            {
                position[i] = x[i];
            }

        }

        public abstract double x_delta();

        public abstract double c_down();
        public abstract double c_up();

        //public abstract double u_cheb_min();
        //public abstract double u_cheb_max();

        public virtual double x0(int i)
        {
            return position[i];
        }
        public abstract  double v(int i, double u, double[] x);

        public double U0()
        {
            return _u0;
        }

        public double U1()
        {
            return _u0 > 0 ? _u_pos : _u_neg;
        }

        public double U2()
        {
            return _u0 > 0 ? _u_neg : _u_pos;
        }

        public abstract double I(double[] rx);

        public abstract double Dt();

        public bool f_basis;
        public bool view_cos;
        public bool chebyshev;
        public int L0, L1, L2;
        public int[] x_L;


        protected double _u0 = 1;
        protected double _u_neg, _u_pos;

       // public int fst;

        public void setUprType(double u0, bool f_basis_, bool view_cos_, bool chebyshev_, int L0_, int[] x_L_) //, int fst_ = 1)
        {
            _u0 = u0;
            f_basis = f_basis_;
            view_cos = view_cos_;
            chebyshev = chebyshev_;
            //fst = fst_;
            L0 = L0_;
            x_L = x_L_;
        }


        public int[] get_x_L()
        {
            return x_L;
        }

        public double[] get_x0()
        {
            return _x0;
        }

        public double[] get_x1()
        {
            return _x1;
        }



        public bool I_count;
        public bool I_count2;
        public void setIcount(bool I_count_, bool I_count2_)
        {         
            I_count = I_count_;
            I_count2 = I_count2_;
        }
        public bool getIcount()
        {
            return I_count;
        }
        public bool getIcount2()
        {
            return I_count2;
        }


        public bool getUprType()
        {
            return f_basis;
        }

        public bool getViewCos()
        {
            return view_cos;
        }
        public bool getChebyshev()
        {
            return chebyshev;
        }
        //public bool getI_count()
        //{
        //    return I_count;
        //}

        //public int getFST()
        //{
        //    return fst;
        //}

        public double[] Solve(U upr)
        {
            double[] traj;
 
            return Solve_proc(upr, out traj, false); 
           

            
        }

        public double[] SolveT(U upr)
        {
            double[] traj;

            Solve_proc(upr, out traj, true); 
            
            return traj;
        }

        //public double[] Solve_proc2(U upr, out double[] traj, bool f0)
        //{
        //    return null;

        //}
        Random R = new Random();

        public double[] Solve_proc(U upr, out double[] traj, bool f0) 
        {
            //-------------------------------------------------------------------------------------
            //МЕТОД РУНГЕ-КУТТЫ
            double[] x = new double[2];
            double[] y = new double[2];
            double[] z = new double[2];
            int d = _x0_a.Length;

            double[] rx = new double[2*d];

            double[] tr = null;
            double dt = Dt();

            int m = (int)Math.Round((t1() - t0()) / dt);
            dt = (t1() - t0()) / m;

            if (f0)
            {
                tr = new double[d * 3 * (m + 1)];

            }

            //for (int s = 0; s < d; s++)
            //{                
            //    x[0] = _x0_a[s];
            //    x[1] = _x1_a[s];

            //    if (f0)
            //    {
            //        tr[s*3*(m+1)] = x[0];
            //        tr[s*3*(m+1) + m + 1] = x[1];
            //        tr[s*3*(m+1) + 2 * (m + 1)] = upr.value(t0(), x);
            //    }

            //    for (int k = 0; k < m; k++)
            //    {
            //        double t = t0() + k * dt;

            //        //double u = t < tp ? 1 : -1;

            //        var k10 = v(0, upr.value(t, x), x);
            //        var k11 = v(1, upr.value(t, x), x);

            //        z[0] = x[0] + k10 * dt / 2;
            //        z[1] = x[1] + k11 * dt / 2;

            //        var k20 = v(0, upr.value(t + dt / 2, z), z);
            //        var k21 = v(1, upr.value(t + dt / 2, z), z);

            //        z[0] = x[0] + k20 * dt / 2;
            //        z[1] = x[1] + k21 * dt / 2;

            //        var k30 = v(0, upr.value(t + dt / 2, z), z);
            //        var k31 = v(1, upr.value(t + dt / 2, z), z);

            //        z[0] = x[0] + k30 * dt;
            //        z[1] = x[1] + k31 * dt;

            //        var k40 = v(0, upr.value(t + dt, z), z);
            //        var k41 = v(1, upr.value(t + dt, z), z);

            //        y[0] = x[0] + (k10 + 2 * k20 + 2 * k30 + k40) * dt / 6;
            //        y[1] = x[1] + (k11 + 2 * k21 + 2 * k31 + k41) * dt / 6;

            //        //y[0] = x[0] + (6*k10 ) * dt / 6;
            //        //y[1] = x[1] + (6*k11 ) * dt / 6;

            //        x[0] = y[0];
            //        x[1] = y[1];
            //        if (f0)
            //        {
            //            tr[s*3*(m+1) + k + 1] = x[0];
            //            tr[s*3*(m+1) + m + 1 + k + 1] = x[1];
            //            tr[s*3*(m+1) + 2 * (m + 1) + k + 1] = upr.value(t, x);
            //        }

            //    }
            //    for (int i = 0; i < 2; i++)
            //    {
            //        rx[s * 2 + i] = x[i];
            //    }
            //}

            //traj = tr;
            //return rx;
            //-------------------------------------------------------------------------------------
            //МЕТОД ЭЙЛЕРА-МАРУЯМЫ


            for (int s = 0; s < d; s++)
            {                
                x[0] = _x0_a[s];
                x[1] = _x1_a[s];

                if (f0)
                {
                    tr[s*3*(m+1)] = x[0];
                    tr[s*3*(m+1) + m + 1] = x[1];
                    tr[s*3*(m+1) + 2 * (m + 1)] = upr.value(t0(), x);
                }

                for (int k = 0; k < m; k++)
                {                   
                    double sigma = 0.01;
                    
                    double t = t0() + k * dt;
                    
                    var k10 = v(0, upr.value(t, x), x);
                    var k11 = v(1, upr.value(t, x), x);

                    y[0] = x[0] + k10 * dt + sigma * NextGaussian(R, 0, 1);// вместо 1 можно Math.Sqrt(dt)
                    y[1] = x[1] + k11 * dt + sigma * NextGaussian(R, 0, 1);


                    x[0] = y[0];
                    x[1] = y[1];
                    if (f0)
                    {
                        tr[s*3*(m+1) + k + 1] = x[0];
                        tr[s*3*(m+1) + m + 1 + k + 1] = x[1];
                        tr[s*3*(m+1) + 2 * (m + 1) + k + 1] = upr.value(t, x);
                    }
                }
                for (int i = 0; i < 2; i++)
                {
                    rx[s * 2 + i] = x[i];
                }
            }

            traj = tr;
            return rx;
   
     //-------------------------------------------------------------------------------------
        }
    }
}
