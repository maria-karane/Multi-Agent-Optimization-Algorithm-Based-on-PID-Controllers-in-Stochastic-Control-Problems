using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Управление
{
    public class U
    {
        double t0, t1, p0, p1x, p1y, p2, w0, w1, w2;
        double[] x0;
        double[] x1;

        
        double[] t_k;


        int m, L, L0, x_dim;
        int[] x_L; 
        double u_min, u_max;

        public bool view_cos;
        public bool chebyshev;


        //public int fst;
        public double[] T, Y;   //T[i] - коэффициенты в разложении

        //
        public int factorial(int m)
        {
            int k = 1;
            int s;
            for (int i = 1; i <= m; i++)
            {
                s = k;
                k = s * i;
            }

            return k;
        }


        public void parse_U_index(int j, out int i, out int[] J)
        {
            if (x_L[1] == 0 || x_L[0] == 0) x_dim = 1;           
            int sigma = 0;
            
            J = new int[x_dim];
            i = j % L0;
            j = j / L0;
           

            for (int k = x_dim - 1; k >= 0; k--)
            {
                if (x_L[0] == 0) sigma = 1;
  
                J[k] = j % x_L[k + sigma];
                j = j / x_L[k + sigma];
            }
            
        }

        public int get_U_index(int i, int[] J)
        {
            int j = 0;

            //int coord = 0;
            //if (x_L[1] == 0) x_dim = x_dim - 1;
            //if (x_L[0] == 0) coord = 1;

            for (int k = 0; k < x_dim; k++)
            {
                j *= x_L[k];
                j += J[k];
            }
            j *= L0;

            return j + i;
        }

        public double value(double t, double[] X)
        {
            p0 = Math.Sqrt((1.0) / t1);
            w0 = Math.PI / t1;
           
            //        w2 = Math.PI / (x1[1] - x0[1]);

            //if (X.Length != x_dim)
            //    throw new InvalidOperationException("wrong dimension X");

            double _2 = Math.Sqrt(2.0/t1);
            if (view_cos)
            {
                var s = 0.0;

                if (x_L[0] == 0 && x_L[1] == 0)
                {
                    //throw new NotImplementedException();

                    for (int i = 0; i < m; i++)
                    {
                        //было      

                        //косинусоиды 1
                        //s += T[i] * Math.Cos(i * Math.PI * (2 * (t - t1) / (t1 - t0) - 1));
                        //

                        //стало

                        //Косинусоиды 2
                        if (i == 0)
                        {
                            s += T[i] * Math.Sqrt((1.0) / t1);
                        }
                        if (i > 0)
                        {
                            s += T[i] * Math.Sqrt((2.0) / t1) * Math.Cos(i * Math.PI * t / t1);
                        }
                        //

                        //Полиномы Лежандра
                        //if ((m - i) % 2 == 0)
                        //{
                        //    s += T[i] * Math.Sqrt((2.0 * m + 1) / t1)
                        //        * (factorial(m + i) / (factorial(i) * factorial(i) * factorial(m - i)))
                        //        * Math.Pow(t, i) / Math.Pow(t1, i);
                        //}
                        //if ((m - i) % 2 != 0)
                        //{
                        //    s += T[i] * Math.Sqrt((2.0 * m + 1) / t1)
                        //        * (-1) * (factorial(m + i) / (factorial(i) * factorial(i) * factorial(m - i)))
                        //        * Math.Pow(t, i) / Math.Pow(t1, i);
                        //}
                        //
                    }
                }
                else if (x_L[0] >= 1 && x_L[1] == 0)
                {
                    for (int j = 0; j < m; j++)
                    {
                        double w2 = 0, polynom_Cos = 0;
                        int i;
                        int[] J;
                        parse_U_index(j, out i, out J);

                        double p = T[j] * p0;
                        p *= _2 * Math.Cos(i * t * w0);
                        //for (int k = 0; k < x_dim; k++)
                        //{
                            int k = 0;
                            w1 = Math.PI / (x1[k] - x0[k]);

                            //когда L1 >= 1, L2 == 0 => k==0                           
                            if (k == 0)
                            {
                                polynom_Cos = Math.Cos(J[k] * (X[k] - x0[k]) * w1);
                                if (polynom_Cos == 1) w2 = Math.Sqrt(1 / (x1[k] - x0[k]));
                                else w2 = Math.Sqrt(2 / (x1[k] - x0[k]));
                            }
                            //else 
                            //{
                            //    polynom_Cos = 1;
                            //    w2 = 1;
                            //}

                            p *= w2 * polynom_Cos;
                        //}
                        s += p;
                    }
                }
                else if (x_L[0] == 0 && x_L[1] >= 1)
                {
                    for (int j = 0; j < m; j++)
                    {
                        double w2 = 0, polynom_Cos = 0;
                        int i;
                        int[] J;
                        parse_U_index(j, out i, out J);

                        double p = T[j] * p0;
                        p *= _2 * Math.Cos(i * t * w0);

                        //for (int k = 0; k < x_dim; k++)
                        //{
                            int k = 1;
                            w1 = Math.PI / (x1[k] - x0[k]);

                            

                            //когда L1 == 0, L2 >= 1 => k==1
                            if (k == 1)
                            {
                                polynom_Cos = Math.Cos(J[k-1] * (X[k] - x0[k]) * w1);
                                if (polynom_Cos == 1) w2 = Math.Sqrt(1 / (x1[k] - x0[k]));
                                else w2 = Math.Sqrt(2 / (x1[k] - x0[k]));
                            }
                            //else { w2 = 1; }
                            //w2 = 1;


                            p *= w2 * polynom_Cos;
                        //}
                        s += p;
                    }
                }

                else
                {
                    for (int j = 0; j < m; j++)
                    {
                        double w2 = 0, polynom_Cos = 0;
                        int i;
                        int[] J;
                        parse_U_index(j, out i, out J);



                        double p = T[j] * p0;
                        p *= _2 * Math.Cos(i * t * w0);

                        for (int k = 0; k < x_dim; k++)
                        {
                            w1 = Math.PI / (x1[k] - x0[k]);

                            polynom_Cos = Math.Cos(J[k] * (X[k] - x0[k]) * w1);

                            //когда L1 и L2 >= 1
                            if (polynom_Cos == 1) w2 = Math.Sqrt(1 / (x1[k] - x0[k]));
                            else w2 = Math.Sqrt(2 / (x1[k] - x0[k]));

                            //когда L1 >= 1, L2 == 0 => k==0                           
                            //if (k == 0)
                            //{
                            //    if (polynom_Cos == 1) w2 = Math.Sqrt(1 / (x1[k] - x0[k]));
                            //    else w2 = Math.Sqrt(2 / (x1[k] - x0[k]));
                            //}
                            //else { w2 = 1; }
                            //w2 = 1;

                            //когда L1 == 0, L2 >= 1 => k==1
                            //if (k == 1)
                            //{
                            //    if (polynom_Cos == 1) w2 = Math.Sqrt(1 / (x1[k] - x0[k]));
                            //    else w2 = Math.Sqrt(2 / (x1[k] - x0[k]));
                            //}
                            //else { w2 = 1; }
                            //w2 = 1;


                            p *= w2 * polynom_Cos;
                        }
                        s += p;
                    }
                }
                if (s < u_min)
                {
                    s = u_min;
                }
                if (s > u_max)
                {
                    s = u_max;
                }

                return s;
            }

            else if (chebyshev)
            {
                double polynom = 0;
                
                //разбиение отрезка t и подсчет узлов-------------------------------------------------
                int Nj = 8;
                Nj = m - 1;

                t_k = new double[Nj+1];


                for (int j = 0; j <= Nj; j++)
                {                    
                    t_k[j] = Math.Cos(Math.PI * j / Nj);                   
                }
                if (t_k[0] != t1 || t_k[Nj] != t0)
                    throw new InvalidOperationException("wrong t");

                //формула 
                for (int i = 0; i <= Nj; i++)
                {
                    //параметр с_i
                    int par_c_i;                  
                    if (i == 0 || i == Nj) par_c_i = 2;
                    else par_c_i = 1;
                    //
                    //коэффициент перед функцией
                    double param, product = 1, Func_Cheb = 0;
                    param = Math.Pow(2, Nj - 1) * Math.Pow((-1), i + 2) / (par_c_i * Nj); 
                    //
                    for (int j = 0; j <= Nj; j++)
                    {
                        if (i == j)
                        {
                            product *= 1;
                        }
                        else product *= t - t_k[j];
                    }
                    Func_Cheb = param * product;
       
                    polynom += T[i] * Func_Cheb;

                }
                //                        ОБРЕЗКА
                if (polynom < u_min)
                {
                    polynom = u_min;
                }
                if (polynom > u_max)
                {
                    polynom = u_max;
                }
                return polynom;
            }

            else
            {
                if (t <= t0) return Y[0];
                for (int i = 0; i < m; i++)
                {
                    if (t < T[i]) return Y[i];
                }
                return Y[m];
            }
        }

        //public double this[double t] 
        //{
        //    get
        //    {
        //        return value(t, null);
        //    }

        //}
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
        //public U(int m_, bool view_cos_)
        //{
        //    m = m_;
        //    view_cos = view_cos_;
        //    T = new double[m];
        //    Y = new double[m + 1];
        //}

        //public U(double[] tp, double u0, bool f_basis, bool view_cos_)
        //{
        //    m = tp.Length;
        //    view_cos = view_cos_;
        //    T = new double[m];
        //    Y = new double[m + 1];
        //    for (int i = 0; i < m; i++)
        //    {
        //        T[i] = tp[i];
        //        Y[i] = u0;
        //        u0 = -u0; 
        //    }
        //    Y[m] = u0;
        //}

        public U(double[] tp, int[] x_L, double u1, double u2, bool f_basis, bool view_cos_, bool chebyshev_, double t0_, double t1_, double[] x0_, double[] x1_)
        {
            this.x_L = x_L;
            x_dim = x_L.Length;
            if (x0_.Length != x_dim)
                throw new InvalidOperationException("wrong dimension x0");

            if (x1_.Length != x_dim)
                throw new InvalidOperationException("wrong dimension x1");

           
            t0 = t0_;
            t1 = t1_;
            m = tp.Length;

            L0 = m;

            x0 = x0_;
            x1 = x1_;

            if (x_L[0] == 0 && x_L[1] != 0) L0 /= x_L[1];
            else if (x_L[1] == 0 && x_L[0] != 0) L0 /= x_L[0];
            else if ((x_L[0] == 0) && (x_L[1] == 0)) { }
            else
            {
                for (int k = 0; k < x_dim; k++) L0 /= x_L[k];
            }
            

            view_cos = view_cos_;
            chebyshev = chebyshev_;
            T = new double[m];
            Y = new double[m + 1];
            for (int i = 0; i < m; i++)
            {
                T[i] = tp[i];
                Y[i] = ((i % 2) == 0) ? u1 : u2;

            }
            Y[m] = ((m % 2) == 0) ? u1 : u2;

            u_min = Math.Min(u1, u2);
            u_max = Math.Max(u1, u2);

        }

        //для пучков
        //public U(double[] tp, double u1, double u2, bool f_basis, bool view_cos_, int fst_,
        //    double t0_, double t1_, double[] x0_arr, double[] x1_arr, int[] LL)
        //{
        //    t0 = t0_;
        //    t1 = t1_;
        //    m = tp.Length;
        //    L = LL[1];
        //    L0 = LL[0];

        //    view_cos = view_cos_;
        //    fst = fst_;
        //    T = new double[m];
        //    Y = new double[m + 1];
        //    for (int i = 0; i < m; i++)
        //    {
        //        T[i] = tp[i];
        //        Y[i] = ((i % 2) == 0) ? u1 : u2;

        //    }
        //    Y[m] = ((m % 2) == 0) ? u1 : u2;

        //    u_min = Math.Min(u1, u2);
        //    u_max = Math.Max(u1, u2);

        //    int d = x0_arr.Length;
        //    x0 = new double[d];
        //    x1 = new double[d];

            

        //    for (int j = 0; j < d; j++) 
        //    {
        //        x0[j] = x0_arr[j];
        //        x1[j] = x1_arr[j];
        //    }

        //    p0 = Math.Sqrt((1.0) / t1);
        //    w0 = Math.PI / t1;
        //    if (d == 2)
        //    {
        //        p0 *= Math.Sqrt((1.0) / (x1[0] - x0[0]));
        //        p0 *= Math.Sqrt((1.0) / (x1[1] - x0[1]));
        //        w1 = Math.PI / (x1[0] - x0[0]);
        //        w2 = Math.PI / (x1[1] - x0[1]);
        //    }
        //    else { throw new Exception("wrong dimension " + d); }
        //}
   
    }

  


}
