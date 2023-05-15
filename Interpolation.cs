using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace Управление
{
    public class Interpolation_Search
    {
        public int Ip;
        public int K;
        public int N;
        public double[,] oblast;

        public ITasks task;

        public int z;
        public double thr;
        public double WScale;

        public double Mfend = 0;

        public List<double> AverF = new List<double>();
        public List<double> BestF = new List<double>();

        Random R = new Random();

        public int Nm;
        int m;
        int M1;
        int M2;
        double PRT;
        int nstep;
        int b2;
        bool f_use_spline, f_use_b_spline, f_use_besie3;


        public int Nf = 0;

        public Vector optimalStart;
        public Vector optimalD;

        public void Init_krill(int ip, int k, double[,] obl, int n, int z_, /*ITasks task_,*/
            /* параметры настройки алгоритма: */ double mu_, double Nmax_, double Vf_, double Dmax_, double ci_)
        {
            if (z_ >= 100) { f_interp = true; }
            K = k;
            N = n;
            oblast = obl;
            //----
            Ip = ip;
            z = z_;

            //-----------------------------------

            //task = task_;

            //oblast = new double[n, 2];
            //for (int i = 0; i < n; i++)
            //{
            //    if (task.getViewCos())
            //    {
            //        oblast[i, 0] = task.c_down() / Math.Sqrt(n);
            //        oblast[i, 1] = task.c_up() / Math.Sqrt(n);
            //    }
            //    else
            //    {
            //        oblast[i, 0] = task.t0();
            //        oblast[i, 1] = task.t1();
            //    }
            //}

            //-------------------------------------





            if (z_ >= 300)
            {
                z = z_ - 300;
                f_use_besie3 = true;
            }
            else if (z_ >= 200)
            {
                z = z_ - 200;
                f_use_b_spline = true;
            }

            else if (z_ >= 100)
            {
                z = z_ - 100;
                f_use_spline = true;
            }
            else
            {
                z = z_;
                f_use_spline = false;
                f_use_b_spline = false;
                f_use_besie3 = false;
            }





            //
            mu = mu_;
            Nmax = Nmax_;
            Vf = Vf_;
            Dmax = Dmax_;
            ci = ci_;
            //
            stepVolInitial = 0.001;
            stepVolFinal = 0.0001;
            //
            optimalStart = new Vector(n);
            optimalD = new Vector(n);
            for (int j = 0; j < n; j++) optimalStart.v[j] = -1;
        }

        public void Init(int ip, int k, int n, ITasks task_,
            /* параметры настройки алгоритма: */ int M1_, int M2_, double PRT_, int nstep_, int b2_
            )
        {
            z = 1;
            task = task_;
            K = k;
            N = n;
            oblast = new double[n, 2];
            for (int i = 0; i < n; i++)
            {
                if (task.getViewCos())
                {
                    oblast[i, 0] = task.c_down() / Math.Sqrt(n);
                    oblast[i, 1] = task.c_up() / Math.Sqrt(n);
                }
                else if (task.getChebyshev())
                {
                    oblast[i, 0] = task.c_down();
                    oblast[i, 1] = task.c_up();
                }
                else
                {
                    oblast[i, 0] = task.t0();
                    oblast[i, 1] = task.t1();
                }
            }
            //----
            Ip = ip;

            //~~~~~~~~~


            M1 = M1_;
            M2 = M2_;
            PRT = PRT_;
            nstep = nstep_;
            b2 = b2_;

            //~~~~~~~~~~~~~
            stepVolInitial = 0.001;
            stepVolFinal = 0.0001;
            //
            optimalStart = new Vector(n);
            optimalD = new Vector(n);
            for (int j = 0; j < n; j++) optimalStart.v[j] = -1;
        }

        class FishComparer : IComparer<Fish>
        {
            public int Compare(Fish x, Fish y)
            {
                if (x.f < y.f) return 1;
                if (x.f > y.f) return -1;
                return 0;
            }
        }


        private double Func(double[] ta, Fish fish)
        {
            if (f_interp)
            {
                int z = this.z;
                if (f_use_spline)
                {
                    z = 14;
                }
                if (f_use_b_spline)
                {
                    z = 15;
                }
                if (f_use_besie3)
                {
                    z = 16;
                }

                double func = 0;

                if (z == 14)
                {
                    double x1 = ta[0];
                    func = Spline(x1);
                }
                else if (z == 15)
                {
                    double x1 = ta[0];
                    func = B_Spline(x1);
                }
                else if (z == 16)
                {
                    double x1 = ta[0];
                    func = Besie3(x1);
                }
                else
                {
                    int stop = 1;
                }
                Nf++;
                return func;

            }


            else
            {
                Array.Sort(ta);
                // var tp = ta[0];

                var _x0 = task.x0();
                var _x1 = task.x1();

                var upr = new U(ta, task.get_x_L(), task.U1(), task.U2(), task.getUprType(), task.getViewCos(), task.getChebyshev(), task.t0(), task.t1(), _x0, _x1);

                var rx = task.Solve(upr);
                if (fish != null) fish.rx = rx;
                return -task.I(rx);
            }
        }

        public bool f_interp = false;

        private double Func1(double[] x, bool original = false)
        {
            int z = this.z;
            if (f_use_spline && !original)
            {
                z = 14;
            }
            if (f_use_b_spline && !original)
            {
                z = 15;
            }
            if (f_use_besie3 && !original)
            {
                z = 16;
            }

            double func = 0;

            if (z == 14)
            {
                double x1 = x[0];
                func = Spline(x1);
            }
            else if (z == 15)
            {
                double x1 = x[0];
                func = B_Spline(x1);
            }
            else if (z == 16)
            {
                double x1 = x[0];
                func = Besie3(x1);
            }
            else
            {
                int stop = 1;
            }
            Nf++;
            return func;
        }

        private static double[] Cpow(double x, double y, int p)
        {
            double[] Cp = new double[2];
            Cp[0] = x; Cp[1] = y;
            double x0 = 0;
            double y0 = 0;
            for (int i = 1; i < p; i++)
            {
                x0 = Cp[0] * x - Cp[1] * y;
                y0 = Cp[1] * x + Cp[0] * y;
                Cp[0] = x0; Cp[1] = y0;
            }
            return Cp;
        }
        //--------------------------------------------------------------------------------------------------

        public struct Vector
        {
            public double[] v;

            public Vector(int n)
            {
                v = new double[n];
            }

            public double distance(Vector y)
            {
                var S = 0.0;
                for (int i = 0; i < v.Length; i++)
                {
                    var d = y.v[i] - v[i];
                    S += d * d;
                }
                return Math.Sqrt(S);
            }

            public void normalize()
            {
                var s = 0.0;
                for (int i = 0; i < v.Length; i++)
                {
                    s += v[i] * v[i];
                }
                if (s == 0) return;
                s = Math.Sqrt(s);
                for (int i = 0; i < v.Length; i++)
                {
                    v[i] /= s;
                }
            }

            public Vector add(Vector y)
            {
                Vector z = new Vector(v.Length);
                for (int i = 0; i < v.Length; i++) z.v[i] = v[i] + y.v[i];
                return z;
            }

            public static Vector operator +(Vector x, Vector y)
            {
                return x.add(y);
            }

            public static Vector operator -(Vector x, Vector y)
            {
                Vector z = new Vector(x.v.Length);
                for (int i = 0; i < x.v.Length; i++) z.v[i] = x.v[i] - y.v[i];
                return z;
            }

            public static Vector operator /(Vector x, double k)
            {
                Vector z = new Vector(x.v.Length);
                for (int i = 0; i < x.v.Length; i++) z.v[i] = x.v[i] / k;
                return z;
            }

            public static Vector operator *(Vector x, double k)
            {
                Vector z = new Vector(x.v.Length);
                for (int i = 0; i < x.v.Length; i++) z.v[i] = x.v[i] * k;
                return z;
            }

            public static Vector operator *(double k, Vector x)
            {
                Vector z = new Vector(x.v.Length);
                for (int i = 0; i < x.v.Length; i++) z.v[i] = x.v[i] * k;
                return z;
            }

            public double this[int i]
            {
                get
                {
                    return v[i];
                }
            }

            public static Vector Zero(int n)
            {
                Vector z = new Vector(n);
                for (int i = 0; i < n; i++) z.v[i] = 0;
                return z;
            }
        }

        public class Fish
        {
            public Interpolation_Search A;

            public double[] rx;

            private Vector _x;
            public double[] dx;
            public double W;
            public double f;
            public double f2;
            public double normf;
            //public List<Fish> Kloni;

            public double df { get { return f2 - f; } }

            public Fish(Interpolation_Search alg)
            {
                A = alg;
                dx = new double[A.N];
                x2 = new double[A.N];
                W = 0.5 * A.WScale;
            }

            public Vector x
            {
                get { return _x; }
                set
                {
                    _x = value;
                    for (int i = 0; i < A.N; i++) x2[i] = _x.v[i];
                    f = A.Func(_x.v, this);
                    f2 = f;
                }
            }

            public double[] X { get { return x.v; } }

            public double[] x2;

            //public double[] X2 
            //{ 
            //    get 
            //    {
            //        var x2 = new double[A.N];
            //        for (int i = 0; i < A.N; i++) x2[i] = x.v[i] + dx[i];
            //        return x2; 
            //    } 
            //}

            public bool X2inOblast
            {
                get
                {
                    //var x2 = X2;
                    for (int i = 0; i < A.N; i++)
                        if (x2[i] < A.oblast[i, 0] || x2[i] > A.oblast[i, 1]) return false;
                    return true;
                }
            }

            public void CalcNextValueOfF()
            {
                f2 = A.Func(x2, this);
            }

            public void DropNexValueOfF()
            {
                for (int i = 0; i < A.N; i++)
                {
                    dx[i] = 0;
                    x2[i] = x.v[i];
                }
                f2 = f;
            }

            public void Move()
            {
                for (int i = 0; i < A.N; i++) x.v[i] = x.v[i] + dx[i];
                f = f2;
            }

            public int pairIndex;
        }

        public List<Fish> Population = new List<Fish>();

        private double maxOblSize
        {
            get
            {
                var m = oblast[0, 1] - oblast[0, 0];
                for (int i = 0; i < N; i++)
                {
                    if (oblast[i, 1] - oblast[i, 0] > m) m = oblast[i, 1] - oblast[i, 0];
                }
                return m;
            }
        }

        private double stepIndInitial;
        private double stepIndFinal;
        private double stepIndDelta { get { return (stepIndFinal - stepIndInitial) / K; } }
        private double stepInd { get { return stepIndInitial + stepIndDelta * k; } }

        private double stepVolInitial;
        private double stepVolFinal;
        private double stepVolDelta { get { return (stepVolFinal - stepVolInitial) / K; } }
        private double stepVol { get { return stepVolInitial + stepVolDelta * k; } }

        // Динамические параметры

        private int k;          // номер итерации
        private double dW;
        private double Mf;

        public int kpop
        {
            get
            {
                return k;
            }
            set
            {
                k = value;
            }
        }

        //private void CalcMf()
        //{
        //    Mf = 0;
        //    for (int i = 0; i < Population.Count; i++)
        //    {
        //        Fish member = Population[i];
        //        Mf += member.f;
        //    }
        //    Mf = Mf / (double)Population.Count;
        //}

        //создание начальной популяции

        public void Create_Pop()
        {
            fMin = 0;
            fMax = 0;
            for (int i = 0; i < Ip; i++)
            {
                Fish member = new Fish(this);
                //member.Kloni = new List<Fish>();

                Vector x0 = new Vector(N);
                for (int j = 0; j < N; j++)
                {
                    double u0 = oblast[j, 0];
                    double u1 = oblast[j, 1];
                    if (optimalStart.v[j] > 0)
                    {
                        u0 = optimalStart.v[j] - optimalD.v[j];
                        u1 = optimalStart.v[j] + optimalD.v[j];
                    }
                    x0.v[j] = R.NextDouble() * (u1 - u0) + u0;
                }
                member.x = x0;

                if (i == 0) fMin = fMax = member.f;
                else
                {
                    if (member.f < fMin) { fMin = member.f; jMin = i; }
                    if (member.f > fMax) { fMax = member.f; jMax = i; }
                }

                Population.Add(member);
                // Mf += member.f;     // ???
            }
            initAlgorithm();
            Analyze();
        }

        public double[] x_opt = new double[2];
        public double I_opt;

        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        public Vector P1, P2, P3, P4;

        double Spline(double t)//Besie 4         100
        {
            var P = Math.Pow(1 - t, 3) * P1 + 3 * t * Math.Pow(1 - t, 2) * P2 + 3 * (1 - t) * t * t * P3 + Math.Pow(t, 3) * P4;
            var y = Func1(P.v, true);

            return y;
        }

        double B_Spline(double t) //CR, B        200
        {
            var P = 0.5 * (-t * Math.Pow(1 - t, 2) * P1 + (2 - 5 * t * t + 3 * Math.Pow(t, 3)) * P2 + t * (1 + 4 * t - 3 * t * t) * P3 - t * t * (1 - t) * P4);
            var y = Func1(P.v, true);
            return y;
        }

        double Besie3(double t)                 //300
        {
            var P = Math.Pow(1 - t, 2) * P1 + 2 * (1 - t) * t * P2 + t * t * P3;
            var y = Func1(P.v, true);
            return y;
        }

        Vector getSplinePoint(double t)
        {
            if (f_use_b_spline)
            {
                return 0.5 * (-t * Math.Pow(1 - t, 2) * P1 + (2 - 5 * t * t + 3 * Math.Pow(t, 3)) * P2 + t * (1 + 4 * t - 3 * t * t) * P3 - t * t * (1 - t) * P4);
            }
            else if (f_use_besie3)
            {
                return Math.Pow(1 - t, 2) * P1 + 2 * (1 - t) * t * P2 + t * t * P3;
            }
            else if (f_use_spline)
            {
                return Math.Pow(1 - t, 3) * P1 + 3 * t * Math.Pow(1 - t, 2) * P2 + 3 * (1 - t) * t * t * P3 + Math.Pow(t, 3) * P4;
            }
            return new Vector(N);
        }





        //шаг 4. Интерполяционный поиск
        public void interp_search()
        {
            int NP = Population.Count;
            Population.Sort(new FishComparer());

            //4.1 Сплайн
            for (int k = 0; k < M1; k++)
            {
                int j1 = R.Next(3, NP - 1);
                int j2 = R.Next(2, j1 - 1);
                int j3 = R.Next(1, j2 - 1);


                const int Ip = 20;
                const int K = 100;
                const double mu = 0.5;
                const double Nmax = 0.01;
                const double Vf = 0.02;
                const double Dmax = 0.005;
                const double ci = 0.2;

                double[,] obl = new double[1, 2];
                obl[0, 0] = 0;
                obl[0, 1] = 1;
                var alg = new Interpolation_Search();
                alg.Init_krill(Ip, K, obl, 1, 100 + z, mu, Nmax, Vf, Dmax, ci);

                alg.P1 = Population[0].x;
                alg.P2 = Population[j1].x;
                alg.P3 = Population[j2].x;
                alg.P4 = Population[j3].x;

                alg.Create_Pop();
                alg.Work_krill();
                alg.Population.Sort(new FishComparer());
                Fish member = new Fish(this);
                member.x = alg.getSplinePoint(alg.Population[0].x.v[0]);
                Population.Add(member);
            }

            //4.2 B-сплайн
            for (int k = 0; k < M2; k++)
            {
                int j1 = R.Next(2, NP - 3);
                int j2 = R.Next(0, j1 - 1);
                int j3 = R.Next(j1 + 1, NP - 2);
                int j4 = R.Next(j3 + 1, NP - 1);


                const int Ip = 20;
                const int K = 100;
                const double mu = 0.5;
                const double Nmax = 0.01;
                const double Vf = 0.02;
                const double Dmax = 0.005;
                const double ci = 0.2;

                double[,] obl = new double[1, 2];
                obl[0, 0] = 0;
                obl[0, 1] = 1;
                var alg = new Interpolation_Search();
                alg.Init_krill(Ip, K, obl, 1, 200 + z, mu, Nmax, Vf, Dmax, ci);

                alg.P1 = Population[j1].x;
                alg.P2 = Population[j2].x;
                alg.P3 = Population[j3].x;
                alg.P4 = Population[j4].x;

                alg.Create_Pop();//?????????????????????????????????????????????????
                alg.Work_krill();
                alg.Population.Sort(new FishComparer());
                Fish member = new Fish(this);
                member.x = alg.getSplinePoint(alg.Population[0].x.v[0]);
                Population.Add(member);
            }

            //4.3
            Population.Sort(new FishComparer());
            Population.RemoveRange(NP, M1 + M2);
            //
        }

        //шаг 5. Миграция популяции
        public void migration()
        {

            for (int j = 1; j < NP; j++)
            {
                var x = Population[j].x;
                Vector PRTVector = new Vector(N);
                for (int i = 0; i < N; i++)
                {
                    PRTVector.v[i] = (R.NextDouble() < PRT) ? 1 : 0;
                }
                Vector u = new Vector(N);
                for (int i = 0; i < N; i++)
                {
                    u.v[i] = (Population[0].x.v[i] - x.v[i]) / (nstep / 2.0) * PRTVector.v[i];
                }
                double maxf = 0;
                Vector new_x = x;
                for (int m = 0; m <= nstep; m++)
                {
                    var x2 = x + m * u;
                    var f1 = Func(x2.v, null);
                    if (m == 0) { maxf = f1; }
                    else if (f1 > maxf)
                    {
                        maxf = f1;
                        new_x = x2;

                    }
                }


                //проверка на принадлежность области !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                for (int i = 0; i < N; i++)
                {
                    if (oblast[i, 0] > new_x.v[i])
                    {
                        new_x.v[i] = oblast[i, 0];
                    }
                    else if (new_x.v[i] > oblast[i, 1])
                    {
                        new_x.v[i] = oblast[i, 1];
                    }
                }
                //

                Population[j].x = new_x;


            }

            Population.Sort(new FishComparer());
        }


        //шаг 6. Фронтальный поиск
        public void frontal_search()
        {
            const int Ip = 20;
            const int K = 100;
            const double mu = 0.5;
            const double Nmax = 0.01;
            const double Vf = 0.02;
            const double Dmax = 0.005;
            const double ci = 0.2;

            double[,] obl = new double[1, 2];
            obl[0, 0] = 0;
            obl[0, 1] = 1;


            //6.1
            var alg = new Interpolation_Search();
            alg.Init_krill(Ip, K, obl, 1, 300 + z, mu, Nmax, Vf, Dmax, ci);

            alg.P1 = Population[0].x;
            alg.P2 = Population[2].x;
            alg.P3 = Population[1].x;

            alg.Create_Pop();
            alg.Work_krill();
            alg.Population.Sort(new FishComparer());
            Fish member1 = new Fish(this);
            member1.x = alg.getSplinePoint(alg.Population[0].x.v[0]);
            Population.Add(member1);



            //6.2
            alg = new Interpolation_Search();
            alg.Init_krill(Ip, K, obl, 1, 200 + z, mu, Nmax, Vf, Dmax, ci);

            alg.P1 = Population[2].x;
            alg.P2 = Population[0].x;
            alg.P3 = Population[1].x;
            alg.P4 = Population[3].x;

            alg.Create_Pop();
            alg.Work_krill();
            alg.Population.Sort(new FishComparer());
            Fish member2 = new Fish(this);
            member2.x = alg.getSplinePoint(alg.Population[0].x.v[0]);
            Population.Add(member2);



            //6.3
            alg = new Interpolation_Search();
            alg.Init_krill(Ip, K, obl, 1, 100 + z, mu, Nmax, Vf, Dmax, ci);

            alg.P1 = Population[0].x;
            alg.P2 = Population[2].x;
            alg.P3 = Population[3].x;
            alg.P4 = Population[1].x;

            alg.Create_Pop();
            alg.Work_krill();
            alg.Population.Sort(new FishComparer());
            Fish member3 = new Fish(this);
            member3.x = alg.getSplinePoint(alg.Population[0].x.v[0]);
            Population.Add(member3);



            //6.4
            alg = new Interpolation_Search();
            alg.Init_krill(Ip, K, obl, 1, 200 + z, mu, Nmax, Vf, Dmax, ci);

            alg.P1 = Population[0].x;
            alg.P2 = Population[2].x;
            alg.P3 = Population[3].x;
            alg.P4 = Population[1].x;

            alg.Create_Pop();
            alg.Work_krill();
            alg.Population.Sort(new FishComparer());
            Fish member4 = new Fish(this);
            member4.x = alg.getSplinePoint(alg.Population[0].x.v[0]);
            Population.Add(member4);


            //6.5

            Population.Sort(new FishComparer());
            Population.RemoveRange(NP - 1, 4);
        }


        //шаг 7. Сокращение популяции
        public void reduction_pop()
        {
            Population.RemoveRange(NP - b2 - 1, b2);

        }



        //шаг 8. Пополнение популяции
        public void refill_pop()
        {
            //Create_Pop
            fMin = 0;
            fMax = 0;
            for (int i = 0; i < b2; i++)
            {
                Fish member = new Fish(this);
                //member.Kloni = new List<Fish>();

                Vector x0 = new Vector(N);
                for (int j = 0; j < N; j++)
                {
                    double u0 = oblast[j, 0];
                    double u1 = oblast[j, 1];
                    if (optimalStart.v[j] > 0)
                    {
                        u0 = optimalStart.v[j] - optimalD.v[j];
                        u1 = optimalStart.v[j] + optimalD.v[j];
                    }
                    x0.v[j] = R.NextDouble() * (u1 - u0) + u0;
                }
                member.x = x0;

                if (i == 0) fMin = fMax = member.f;
                else
                {
                    if (member.f < fMin) { fMin = member.f; jMin = i; }
                    if (member.f > fMax) { fMax = member.f; jMax = i; }
                }

                Population.Add(member);
                // Mf += member.f;     // ???
            }



        }


        //шаг 10. Выбор решения из последней популяции.
        public void decision_choice()
        {
            //выбрать решение
            //Population[0].x;
        }



















        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // основной цикл

        /* F SRDF SRDF SRDF ... SRD   (F)
        */

        double fMin, fMax;
        int jMin, jMax;

        Fish fishMin { get { return Population[jMin]; } }
        Fish fishMax { get { return Population[jMax]; } }

        double[,] dist;
        double[] ds;
        Vector xfood;
        Vector[] alphaLoc;
        Vector[] alphaTarget;
        Vector[] alpha;
        double[] fbest;
        Vector[] Nold;
        Vector[] Nnew;
        Vector[] betafood;
        Vector[] betabest;
        Vector[] beta;
        Vector[] Fold;
        Vector[] Fnew;
        Vector[] xjbest;
        double[] fjbest;
        Vector[] xNew;
        Vector[] xCr;
        int[] neib;
        double mu = 3;
        double Nmax = 0.01;
        int NP { get { return Ip; } }
        double Vf = 0.02;
        double Dmax = 0.005;
        double ci = 0.2;
        double dt;

        public void initAlgorithm()
        {
            int NP = Population.Count;
            int NP0 = NP;
            //
            neib = new int[NP0];
            dist = new double[NP0, NP0];
            ds = new double[NP0];
            fbest = new double[NP0];
            alphaLoc = new Vector[NP0];
            alphaTarget = new Vector[NP0];
            alpha = new Vector[NP0];
            Nold = new Vector[NP0];
            Nnew = new Vector[NP0];
            betafood = new Vector[NP0];
            betabest = new Vector[NP0];
            beta = new Vector[NP0];
            Fold = new Vector[NP0];
            Fnew = new Vector[NP0];
            fjbest = new double[NP0];
            xjbest = new Vector[NP0];
            xNew = new Vector[NP0];
            xCr = new Vector[NP0];
            //
            dt = 0;
            for (int i = 0; i < N; i++)
            {
                dt += oblast[i, 1] - oblast[i, 0];
            }
            dt *= ci;
            //
            updateDistances();
            //
            for (int j = 0; j < NP; j++)
            {
                Nold[j] = Vector.Zero(N);
                Fold[j] = Vector.Zero(N);
                xjbest[j] = Vector.Zero(N);
                fjbest[j] = 0;
            }
        }

        public void updateDistances()
        {
            int NP = Population.Count;
            for (int i = 0; i < NP; i++)
            {
                dist[i, i] = 0;

                for (int j = 0; j < NP; j++)
                {
                    if (i != j)
                        dist[i, j] = Population[i].x.distance(Population[j].x);
                }
            }

            for (int i = 0; i < NP; i++)
            {

                double s = 0;
                for (int j = 0; j < NP; j++)
                {
                    s += dist[j, i];
                }
                ds[i] = s / (5.0 * NP);

                if (i == 0) fMin = fMax = Population[i].f;
                else
                {
                    if (Population[i].f < fMin) { fMin = Population[i].f; jMin = i; }
                    if (Population[i].f > fMax) { fMax = Population[i].f; jMax = i; }
                }

                if (Population[i].f > fjbest[i])
                {
                    xjbest[i] = Population[i].x;
                    fjbest[i] = Population[i].f;
                }
            }

            Analyze(false);
        }

        public void updateNeibours()
        {
            int NP = Population.Count;
            for (int j = 0; j < NP; j++)
            {
                int n1 = 0;
                Vector dxBar = Vector.Zero(N);
                double dfBar = 0;
                alphaLoc[j] = Vector.Zero(N);

                for (int k = 0; k < NP; k++)
                {
                    if (k != j)
                    {
                        if (dist[j, k] < ds[j])
                        {
                            n1++;
                            dxBar = (Population[k].x - Population[j].x) / (mu + dist[j, k]);
                            dfBar += (Population[j].f - Population[k].f) / (fMin - fMax);
                            alphaLoc[j] = alphaLoc[j] + dfBar * dxBar;
                        }
                    }
                }
                neib[j] = n1;
            }
        }

        public void updateAlphaTarget(int k)
        {
            double cbest = 2 * (R.NextDouble() + 1.0 * k / K);
            //
            int NP = Population.Count;
            for (int i = 0; i < NP; i++)
            {
                double dfBar = (Population[i].f - fMax) / (fMin - fMax);
                fbest[i] = dfBar;
                Vector dxBar = (fishMax.x - Population[i].x) / (mu + dist[i, jMax]);
                alphaTarget[i] = cbest * dfBar * dxBar;
                alpha[i] = alphaLoc[i] + alphaTarget[i];
            }
        }

        public void calcNewN(int k)
        {
            // ? // double omega = 0.8 * double(k)/K + 0.2 * R.NextDouble();
            //
            int NP = Population.Count;
            for (int i = 0; i < NP; i++)
            {
                double omega = 0.8 - (0.8 * k / K) + 0.2 * R.NextDouble();
                Nnew[i] = Nmax * alpha[i] + omega * Nold[i];
            }
        }

        private string xfoodLock = "";

        public void calcXFood()
        {
            xfood = Vector.Zero(N);
            double D = 0;
            //
            lock (xfoodLock)
            {
                int NP = Population.Count;
                for (int i = 0; i < NP; i++)
                {
                    double fi = 1.0 / Population[i].f;
                    D += fi;
                    xfood = xfood + fi * Population[i].x;
                }
                xfood = xfood / D;
            }
        }

        public void calcBetaFood(int k)
        {
            double ff = Func(xfood.v, null);
            //
            int NP = Population.Count;
            for (int i = 0; i < NP; i++)
            {
                Vector xj = Population[i].x;
                double dfBar = (Population[i].f - ff) / (fMin - ff);
                Vector dxBar = (xfood - xj) / (mu + xj.distance(xfood));
                double cfood = 2 * (1 - 1.0 * k / K);
                betafood[i] = cfood * dfBar * dxBar;
            }
        }

        public void calcBetaBest()
        {
            int NP = Population.Count;
            for (int i = 0; i < NP; i++)
            {
                Vector xj = Population[i].x;
                double dfBar = (Population[i].f - fjbest[i]) / (fMin - fMax); /// dF^{j,jbest}         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! // fbest[i];
                Vector dxBar = (xjbest[i] - Population[i].x) / (mu + xjbest[i].distance(xj));
                betabest[i] = dfBar * dxBar;
                beta[i] = betafood[i] + betabest[i];
            }
        }

        public void callsNewF(int k)
        {
            int NP = Population.Count;
            for (int j = 0; j < NP; j++)
            {
                double omegaf = 0.8 - (0.8 * k / K) + 0.2 * R.NextDouble();
                Fnew[j] = Vf * beta[j] + omegaf * Fold[j];
            }
        }

        public void steps567(int k)
        {
            int NP = Population.Count;
            for (int j = 0; j < NP; j++)
            {
                Vector D = new Vector(N);
                for (int i = 0; i < N; i++) D.v[i] = Dmax * (1 - 1.0 * k / K) * (-1 + 2 * R.NextDouble());
                //
                Vector V = D + Fnew[j] + Nnew[j];
                xNew[j] = Population[j].x + dt * V;
                //
                for (int i = 0; i < N; i++)
                {
                    if (xNew[j].v[i] < oblast[i, 0] || xNew[j].v[i] >= oblast[i, 1])
                        xNew[j].v[i] = oblast[i, 0] + (oblast[i, 1] - oblast[i, 0]) * R.NextDouble();
                }
            }
        }




        public void crossover()
        {
            int NP = Population.Count;





            for (int j = 0; j < NP; j++)
            {
                double Cr = 0.2 * (Population[j].f - fMax) / (fMin - fMax);
                //
                xCr[j] = Vector.Zero(N);
                if ((fMin - fMax) == 0)
                    continue;
                for (int i = 0; i < N; i++)
                {
                    if (R.NextDouble() < Cr)
                    {
                        int r = j;
                        while (r == j) r = R.Next() % NP;
                        xCr[j].v[i] = xNew[r].v[i];
                    }
                    else
                        xCr[j].v[i] = xNew[j].v[i];
                }
                int stop = 0;
                if (N > 1)
                    if (double.IsNaN(xCr[j].v[0]) || double.IsNaN(xCr[j].v[1]))
                        stop = 1;
            }
        }

        public void mutation()
        {
            int NP = Population.Count;
            for (int j = 0; j < NP; j++)
            {
                double mu = 0.05 * fbest[j];
                //
                Vector xMu = Vector.Zero(N);
                for (int i = 0; i < N; i++)
                {
                    double randi = R.NextDouble();
                    if (randi < mu)
                    {
                        int p = j;
                        while (p == j) p = R.Next() % NP;
                        int q = j;
                        while (q == j) q = R.Next() % NP;
                        //
                        double nu = R.NextDouble();
                        xMu.v[i] = Population[jMax].x.v[i] +
                            nu * (Population[p].x.v[i] - Population[q].x.v[i]);
                    }
                    else
                        xMu.v[i] = xCr[j].v[i];
                    //
                    if (xMu.v[i] < oblast[i, 0] || xMu.v[i] >= oblast[i, 1])
                        xMu.v[i] = oblast[i, 0] + (oblast[i, 1] - oblast[i, 0]) * R.NextDouble();
                }
                Population[j].x = xMu;
            }
        }
        public string toprotocol = "";
        public void Work()//работа алгоритма
        {
            initAlgorithm();

            for (k = 0; k < K; k++)
            {
                interp_search();
                migration();
                frontal_search();
                reduction_pop();
                refill_pop();
                Population.Sort(new FishComparer());
                //if (k > 0)
                //{
                //    Analyze();
                //    Console.WriteLine("k=" + k + " I=" + I_opt);
                //    if (I_opt < -0.13) break;
                //}
            }
            Analyze(true);
            if ((File.Exists("protocol.dt")))
            {
                FileStream fs = new FileStream("protocol.dt", FileMode.Append, FileAccess.Write);
                StreamWriter r1 = new StreamWriter(fs);
                //Population[0].x.v[0]

                string pr = "";
                pr += "\r\n" + "\r\n";
                pr += "Гибридный мультиагентный алгоритм интерполяционного поиска" + "\r\n" + "\r\n";
                pr += "КОНЕЧНАЯ ПОПУЛЯЦИЯ" + "\r\n" + "\r\n";
                pr += "|-----------------|--------------------------------------------------|----------------------------|\r\n";
                pr += "|   Номер агента  |                 Коэффициенты                     |         Критерий           |\r\n";
                pr += "|-----------------|--------------------------------------------------|----------------------------|\r\n";


                for (int j = 0; j < Ip; j++)
                {
                    pr += "|        " + (j + 1) + "        |";
                    for (int i = 0; i < N; i++)
                    {

                        pr += Math.Round(Population[j].x.v[i], 2) + " ";
                    }
                    pr += "|    " + (-Math.Round(Population[j].f, 4)) + "     |" + "\r\n";

                }

                toprotocol = pr;
                r1.Write(toprotocol);

                r1.Close();
                fs.Close();
            }
        }

        public void Work_krill()//работа алгоритма
        {
            initAlgorithm();

            for (k = 0; k < K; k++)
            {
                updateNeibours();
                updateAlphaTarget(k);
                calcNewN(k);
                calcXFood();
                calcBetaFood(k);
                calcBetaBest();
                callsNewF(k);
                steps567(k);
                crossover();
                mutation();
                updateDistances();



                // var fish = Population[0];
                /*
                LogIteration();
                if (k > 0) Food();//питание 
                Swiming();//плавание
                Reproduction();//размножение
                Delete();//удаление
                */
            }
            Analyze(true);
        }

        public double[] traj;
        public double[] ta;
        public int _numb_tr;

        void Analyze(bool final = false)
        {
            if (!f_interp)
            {
                double i0 = 1000;
                Fish f_opt = null;

                foreach (Fish member in Population)
                {
                    double v = task.I(member.rx);
                    if (v < i0)
                    {
                        I_opt = i0 = v;
                        x_opt = member.rx;
                        f_opt = member;
                    }

                }

                if (f_opt != null)
                {
                    ta = f_opt.x.v;
                    Array.Sort(ta);

                    var upr = new U(ta, task.get_x_L(), task.U1(), task.U2(), task.getUprType(), task.getViewCos(), task.getChebyshev(), task.t0(), task.t1(), task.x0(), task.x1());
                    traj = task.SolveT(upr);
                    _numb_tr = task.numb_tr();
                    // Console.WriteLine(ta);
                }
            }

            RecalcMaxDf();
            if (!final) return;

            Mfend = Math.Round(Mf, 6);


        }

        //операция плавания

        public void Swiming()
        {
            calcBetaBest();
            callsNewF(k);
            steps567(k);
        }
        //операция питания
        public void Food()
        {
            updateNeibours();
            updateAlphaTarget(k);
            calcNewN(k);
            calcXFood();
            calcBetaFood(k);
        }
        //операция разножения
        public void Reproduction()
        {
            crossover();
        }
        //удаление рыб с наим весом
        public void Delete()
        {
            mutation();
            updateDistances();
            RecalcMaxDf();
        }

        double maxDf = 0;

        private void RecalcMaxDf()
        {
            maxDf = 0;
            var S = 0.0;
            for (int j = 0; j < Population.Count; j++)
            {
                Fish member = Population[j];
                if (member.df > maxDf) maxDf = member.df;
                S += member.f;
            }
            BestF.Add(maxDf);
            Mf = S / Population.Count;
            AverF.Add(Mf);
        }

        private void LogIteration()
        {
            foreach (var fish in Population)
            {
                Console.Write(fish.W + "  ");
                Console.WriteLine();
            }
        }
    }
}
