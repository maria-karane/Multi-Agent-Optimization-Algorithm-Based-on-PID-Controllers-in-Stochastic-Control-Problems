using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using MathNet;
using MathNet.Numerics.LinearAlgebra.Double;

using System.IO;


namespace Управление
{
    public class PID_regulators
    {
        public int Ip;
        public int K;
        public int N;
        public double[,] oblast;

        public int z;
        public double thr;
        public double WScale;

        public double Mfend = 0;

        public List<double> AverF = new List<double>();
        public List<double> BestF = new List<double>();

        Random R = new Random();

        public ITasks task;

        public int Nm;
        int m;
        int M1;
        int M2;
        double PRT;
        int nstep;
        int b2;
        //bool f_use_spline, f_use_b_spline, f_use_besie3;


        public int Nf = 0;

        public Vector optimalStart;
        public Vector optimalD;


        public void Init(int ip, int k, double[,] obl, int n, int z_,
            /* параметры настройки алгоритма: */ int M1_, int M2_, double PRT_, int nstep_, int b2_
            )
        {
            K = k;
            N = n;
            oblast = obl;
            //----
            Ip = ip;
            z = z_;

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


        //
        public int Pmax;
        double kp;
        double kd1;
        double kd2;
        double ki;
        double h;
        //


        public void Init_NEW(int ip, int k, int n, ITasks task_,
            /* параметры настройки алгоритма: */ int Pmax_, double kp_, double kd1_, double kd2_, double ki_, double h_
    )
        {
            task = task_;

            K = k;//итерации
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
            Pmax = Pmax_;//число проходов
            kp = kp_;                               // kp
            kd1 = kd1_;                             // kd1
            kd2 = kd2_;                             // kd2
            ki = ki_;                               // ki
            h = h_;
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
            Array.Sort(ta);            
            var upr = new U(ta, task.get_x_L(), task.U1(), task.U2(), task.getUprType(), task.getViewCos(), task.getChebyshev(), task.t0(), task.t1(), task.x0(), task.x1());

            var rx = task.Solve(upr);
            if (fish != null) fish.rx = rx;
            return -task.I(rx);            
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
            public PID_regulators PIDR;
            public double[] rx;

            private Vector _x;
            public double[] dx;
            public double W;
            public double f;
            public double f2;
            public double normf;
            //public List<Fish> Kloni;

            public double df { get { return f2 - f; } }

            public Fish(PID_regulators pidr)
            {
                PIDR = pidr;
                dx = new double[PIDR.N];
                x2 = new double[PIDR.N];
                W = 0.5 * PIDR.WScale;
            }

            public Vector x
            {
                get { return _x; }
                set
                {
                    _x = value;
                    for (int i = 0; i < PIDR.N; i++) x2[i] = _x.v[i];
                    f = PIDR.Func(_x.v, this);
                    f2 = f;
                }
            }

            public double[] X { get { return x.v; } }

            public double[] x2;

            public bool X2inOblast
            {
                get
                {
                    for (int i = 0; i < PIDR.N; i++)
                        if (x2[i] < PIDR.oblast[i, 0] || x2[i] > PIDR.oblast[i, 1]) return false;
                    return true;
                }
            }

            public void CalcNextValueOfF()
            {
                f2 = PIDR.Func(x2, this);
            }

            public void UpdateF()
            {
                f = PIDR.Func(x.v, this);
            }

            public void DropNexValueOfF()
            {
                for (int i = 0; i < PIDR.N; i++)
                {
                    dx[i] = 0;
                    x2[i] = x.v[i];
                }
                f2 = f;
            }

            public void Move()
            {
                for (int i = 0; i < PIDR.N; i++) x.v[i] = x.v[i] + dx[i];
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


        //шаг 2. создание начальной популяции

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

            if ((File.Exists("protocol.dt")))
            {
                FileStream fs = new FileStream("protocol.dt", FileMode.Append, FileAccess.Write);
                StreamWriter r1 = new StreamWriter(fs);


                string pr = "";
                pr += "\r\n" + "\r\n";
                pr += "Мультиагентный алгоритм на основе ПИД-регуляторов" + "\r\n" + "\r\n";
                pr += "НАЧАЛЬНАЯ ПОПУЛЯЦИЯ" + "\r\n" + "\r\n";
                pr += "|-----------------|--------------------------------------------------|----------------------------|\r\n";
                pr += "|   Номер агента  |                 Коэффициенты                     |         Критерий           |\r\n";
                pr += "|-----------------|--------------------------------------------------|----------------------------|\r\n";

                toprotocol = pr;
                r1.Write(toprotocol);

                r1.Close();
                fs.Close();
            }

            Analyze();
        }


        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        //



        //шаг 5. Первая группа агентов
        public void move_1group()
        {

            int NP = Population.Count;
            Population.Sort(new FishComparer());
            var fish1 = Population[0];
            //var X_best = (MathNet.Numerics.LinearAlgebra.Matrix<double>)DenseMatrix.OfArray(new[,] { 
            //                                        { fish1.x.v[0]}, 
            //                                        { fish1.x.v[1]}, 
            //                                        { fish1.x.v[2]}, 
            //                                        { fish1.x.v[3]}
            //                                        //{ fish1.x.v[4]},
            //                                        //{ fish1.x.v[5]},
            //                                        //{ fish1.x.v[6]},
            //                                        //{ fish1.x.v[7]}
            //                                        });

            var X_best = (MathNet.Numerics.LinearAlgebra.Matrix<double>)DenseMatrix.Create(N, 1, 0);
            for (int i = 0; i < N; i++)
            {
                X_best[i, 0] = fish1.x.v[i];
            }

            for (int fish_i = 1; fish_i < NP; fish_i++)
            {
                var fish = Population[fish_i];
                var x = fish.x;
                int group = (int)(fish_i / (NP / 4.0)) + 1;
                if (group > 4) { group = 4; }


                //---------------------------------------------------------
                var A = DenseMatrix.Create(N, N, 0.0);
                for (int i = 0; i < N / 2; i++)
                {
                    for (int j = N / 2; j < N; j++)
                    {
                        if (i == (j - (N / 2))) { A[i, j] = 1; }
                        else { A[i, j] = 0; }
                    }
                }

                //-----------------------------------------------------------
                var B = DenseMatrix.Create(N, N / 2, 0.0);

                for (int i = N / 2; i < N; i++)
                {
                    for (int j = 0; j < N / 2; j++)
                    {
                        if (j == (i - N / 2)) { B[i, j] = 1; }
                        else { B[i, j] = 0; }
                    }
                }

                //KP N*2N
                //------------------------------------------------------------
                var KP = DenseMatrix.Create(N / 2, N, 0.0);
                for (int i = 0; i < N / 2; i++)
                {
                    for (int j = 0; j < N; j++)
                    {
                        if (j == i) { KP[i, j] = kp; }
                        else { KP[i, j] = 0; }
                    }
                }

                //KD1 N*2N
                //------------------------------------------------------------
                var KD1 = DenseMatrix.Create(N / 2, N, 0.0);
                for (int i = 0; i < N / 2; i++)
                {
                    for (int j = 0; j < N; j++)
                    {
                        if (j == i) { KD1[i, j] = kd1; }
                        else { KD1[i, j] = 0; }
                    }
                }

                //KD2 N*2N
                //------------------------------------------------------------
                var KD2 = DenseMatrix.Create(N / 2, N, 0.0);
                for (int i = 0; i < N / 2; i++)
                {
                    for (int j = 0; j < N; j++)
                    {
                        if (j == i) { KD2[i, j] = kd2; }
                        else { KD2[i, j] = 0; }
                    }
                }

                //KI N*2N
                //------------------------------------------------------------
                var KI = DenseMatrix.Create(N / 2, N, 0.0);

                for (int i = 0; i < N / 2; i++)
                {
                    for (int j = 0; j < N; j++)
                    {
                        if (j == i) { KI[i, j] = ki; }
                        else { KI[i, j] = 0; }
                    }
                }
                //------------------------------------------------------------
                //Решение ДУ dX/dt = (A - B * Q^-1 * B^T * P) * X => X

                //!!!!!! для n>2
                //var X = (MathNet.Numerics.LinearAlgebra.Matrix<double>)DenseMatrix.OfArray(new[,] { { x.v[0]}, 
                //                                    {x.v[1]}, 
                //                                    {x.v[2]}, 
                //                                    {x.v[3]}
                //                                    //{x.v[4]},
                //                                    //{x.v[5]}, 
                //                                    //{x.v[6]},
                //                                    //{x.v[7]}
                //                                    });

                var X = (MathNet.Numerics.LinearAlgebra.Matrix<double>)DenseMatrix.Create(N, 1, 0);
                for (int i = 0; i < N; i++)
                {
                    X[i, 0] = x.v[i];
                }

                X[0, 0] -= X_best[0, 0];
                X[1, 0] -= X_best[1, 0];
                var maxf = -1000.0;

                var X1 = X;

                var dX1 = X;
                var dX2 = X;
                var dX3 = X;
                var integral = X;
                var Uprav = X;

                if (group == 1) // расширенный ПИД-регулятор
                {
                    for (int j = 1; j <= K; j++)
                    {
                        //------------------------------------------------------------------------

                        //var Uprav = KP * X + KD1 * (1 / h) * X + KD2 * (1 / (h*h)) * X
                        //    - (KD1 * (1 / h) * dX1 + 5 * KD2 * (1 / (h * h)) * dX1)
                        //    + 4 * KD2 * (1 / (h * h)) * dX2
                        //    - KD2 * (1 / (h * h)) * dX3 + integral;

                        if (j >= 1)
                        {
                            Uprav = KP * X + KD1 * (1 / h) * X + KD2 * (1 / (h * h)) * X + KI * (h / 2) * integral;
                        }
                        if (j >= 2)
                        {
                            Uprav += -(KD1 * (1 / h) * dX1 + 5 * KD2 * (1 / (h * h)) * dX1);
                        }
                        if (j >= 3)
                        {
                            Uprav += 4 * KD2 * (1 / (h * h)) * dX2;
                        }
                        if (j >= 4)
                        {
                            Uprav += -KD2 * (1 / (h * h)) * dX3;
                        }

                        var X2 = A * X + B * Uprav;

                        if (j == 2) dX1 = X;
                        else if (j == 3)
                        {
                            dX2 = dX1;
                            dX1 = X;
                        }
                        else if (j >= 4)
                        {
                            dX3 = dX2;
                            dX2 = dX1;
                            dX1 = X;
                        }
                        X = X2;
                        integral += 2 * X;
                        //------------------------------------------------------------------------

                        fish.x.v[0] = X_best[0, 0] + X[0, 0];
                        if (fish.x.v[0] < oblast[0, 0]) fish.x.v[0] = oblast[0, 0];
                        if (fish.x.v[0] > oblast[0, 1]) fish.x.v[0] = oblast[0, 1];

                        fish.x.v[1] = X_best[1, 0] + X[1, 0];
                        if (fish.x.v[1] < oblast[1, 0]) fish.x.v[1] = oblast[1, 0];
                        if (fish.x.v[1] > oblast[1, 1]) fish.x.v[1] = oblast[1, 1];


                        var f = Func(fish.x.v, fish);

                        if (f > maxf)
                        {
                            maxf = f;
                            X1 = X;
                        }
                    }

                }
                else if (group == 2) // классический ПИД-регулятор
                {
                    for (int j = 1; j <= K; j++)
                    {
                        //------------------------------------------------------------------------

                        if (j >= 1)
                        {
                            Uprav = KP * X + KD1 * (1 / h) * X + KI * (h / 2) * integral;
                        }
                        if (j >= 2)
                        {
                            Uprav += -(KD1 * (1 / h) * dX1);
                        }

                        var X2 = A * X + B * Uprav;

                        if (j >= 2)
                        {
                            dX1 = X;
                        }
                        X = X2;
                        integral += 2 * X;
                        //------------------------------------------------------------------------

                        fish.x.v[0] = X_best[0, 0] + X[0, 0];
                        if (fish.x.v[0] < oblast[0, 0]) fish.x.v[0] = oblast[0, 0];
                        if (fish.x.v[0] > oblast[0, 1]) fish.x.v[0] = oblast[0, 1];

                        fish.x.v[1] = X_best[1, 0] + X[1, 0];
                        if (fish.x.v[1] < oblast[1, 0]) fish.x.v[1] = oblast[1, 0];
                        if (fish.x.v[1] > oblast[1, 1]) fish.x.v[1] = oblast[1, 1];


                        var f = Func(fish.x.v, fish);

                        if (f > maxf)
                        {
                            maxf = f;
                            X1 = X;
                        }
                    }

                }
                else if (group == 3) //  ПИ-регулятор
                {
                    for (int j = 1; j <= K; j++)
                    {
                        //------------------------------------------------------------------------

                        var X2 = A * X + B * (KP * X + KI * (h / 2) * integral);
                        X = X2;
                        integral += 2 * X;

                        //------------------------------------------------------------------------

                        fish.x.v[0] = X_best[0, 0] + X[0, 0];
                        if (fish.x.v[0] < oblast[0, 0]) fish.x.v[0] = oblast[0, 0];
                        if (fish.x.v[0] > oblast[0, 1]) fish.x.v[0] = oblast[0, 1];

                        fish.x.v[1] = X_best[1, 0] + X[1, 0];
                        if (fish.x.v[1] < oblast[1, 0]) fish.x.v[1] = oblast[1, 0];
                        if (fish.x.v[1] > oblast[1, 1]) fish.x.v[1] = oblast[1, 1];


                        var f = Func(fish.x.v, fish);

                        if (f > maxf)
                        {
                            maxf = f;
                            X1 = X;
                        }
                    }

                }
                else if (group == 4) //  ПД-регулятор
                {
                    for (int j = 1; j <= K; j++)
                    {
                        //------------------------------------------------------------------------

                        if (j >= 1)
                        {
                            Uprav = KP * X + KD1 * (1 / h) * X;
                        }
                        if (j >= 2)
                        {
                            Uprav += -(KD1 * (1 / h) * dX1);
                        }

                        var X2 = A * X + B * Uprav;

                        if (j >= 2)
                        {
                            dX1 = X;
                        }
                        X = X2;
                        integral += 2 * X;
                        //------------------------------------------------------------------------

                        fish.x.v[0] = X_best[0, 0] + X[0, 0];
                        if (fish.x.v[0] < oblast[0, 0]) fish.x.v[0] = oblast[0, 0];
                        if (fish.x.v[0] > oblast[0, 1]) fish.x.v[0] = oblast[0, 1];

                        fish.x.v[1] = X_best[1, 0] + X[1, 0];
                        if (fish.x.v[1] < oblast[1, 0]) fish.x.v[1] = oblast[1, 0];
                        if (fish.x.v[1] > oblast[1, 1]) fish.x.v[1] = oblast[1, 1];


                        var f = Func(fish.x.v, fish);

                        if (f > maxf)
                        {
                            maxf = f;
                            X1 = X;
                        }
                    }
                }

                X = X1;

                //var X2 = (A - B * Q.Inverse() * B.Transpose() * P) * X;

                X[0, 0] += X_best[0, 0];
                X[1, 0] += X_best[1, 0];

                for (int i = 0; i < N; i++)
                {
                    fish.x.v[i] = X[i, 0];
                    if (fish.x.v[i] < oblast[i, 0]) fish.x.v[i] = oblast[i, 0];
                    if (fish.x.v[i] > oblast[i, 1]) fish.x.v[i] = oblast[i, 1];
                }

                fish.UpdateF();

                //X_new = X_best + X2;
            }
            RecalcMaxf();

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

        //шаг 7. Сокращение популяции
        public void reduction_pop()
        {
            b2 = (int)Math.Round(NP * 0.1);
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

                if (double.IsNaN(member.x.v[0]) && Population.Count == 30)
                    Console.Write("");

                if (i == 0) fMin = fMax = member.f;
                else
                {
                    if (member.f < fMin) { fMin = member.f; jMin = i; }
                    if (member.f > fMax) { fMax = member.f; jMax = i; }
                }

                if (double.IsNaN(member.x.v[0]) && Population.Count == 30)
                    Console.Write("");
                Population.Add(member);
                // Mf += member.f;     // ???
            }
        }

        public void last_step_sort()
        {
            Population.Sort(new FishComparer());
            Analyze(true);
        }

        public string toprotocol = "";
        public void Work()//работа алгоритма
        {
            for (k = 0; k < Pmax; k++)
            {
                move_1group();

                //Population.Sort(new FishComparer());
                //reduction_pop();

                //if (k < Pmax - 1)
                //{
                //    refill_pop();
                //}


            }
            Population.Sort(new FishComparer());
            Analyze(true);

            if ((File.Exists("protocol.dt")))
            {
                FileStream fs = new FileStream("protocol.dt", FileMode.Append, FileAccess.Write);
                StreamWriter r1 = new StreamWriter(fs);
                //Population[0].x.v[0]

                string pr = "";
                pr += "\r\n" + "\r\n";
                pr += "Мультиагентный алгоритм на основе ПИД-регуляторов" + "\r\n" + "\r\n";
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


        public double[] x_opt = new double[2];
        public double I_opt;
        public double[] traj;
        public int _numb_tr;

        public double[] ta;

        public void Work_krill()//работа алгоритма
        {
            Analyze(true);
        }

        void Analyze(bool final = false)
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

            }
            RecalcMaxf();
            if (!final) return;

            Mfend = Math.Round(Mf, 6);
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
            //BestF.Add(maxDf);
            Mf = S / Population.Count;
            //AverF.Add(Mf);
        }

        public void RecalcMaxf()
        {
            maxDf = 0;
            var S = 0.0;
            for (int j = 0; j < Population.Count; j++)
            {
                Fish member = Population[j];
                if (member.f > maxDf) maxDf = member.f;
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
