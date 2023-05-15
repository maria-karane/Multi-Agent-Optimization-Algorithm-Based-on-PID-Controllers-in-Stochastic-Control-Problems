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
    class Linear_regulator
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
        int Pmax;
        double ks;
        double kl;
        double h;
        //


        public void Init_NEW(int ip, int k, int n, ITasks task_,
            /* параметры настройки алгоритма: */ int Pmax_, double ks_, double kl_, double h_
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
            //z = z_;

            //~~~~~~~~~
            Pmax = Pmax_;//число проходов
            ks = ks_;
            kl = kl_;
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
            public Linear_regulator LR;

            public double[] rx;

            //public Algoritm A;

            private Vector _x;
            public double[] dx;
            public double W;
            public double f;
            public double f2;
            public double normf;
            //public List<Fish> Kloni;

            public double df { get { return f2 - f; } }

            public Fish(Linear_regulator lr)
            {
                LR = lr;
                dx = new double[LR.N];
                x2 = new double[LR.N];
                W = 0.5 * LR.WScale;
            }

            public Vector x
            {
                get { return _x; }
                set
                {
                    _x = value;
                    for (int i = 0; i < LR.N; i++) x2[i] = _x.v[i];
                    f = LR.Func(_x.v, this);
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
                    for (int i = 0; i < LR.N; i++)
                        if (x2[i] < LR.oblast[i, 0] || x2[i] > LR.oblast[i, 1]) return false;
                    return true;
                }
            }

            public void CalcNextValueOfF()
            {
                f2 = LR.Func(x2, this);
            }

            public void UpdateF()
            {
                f = LR.Func(x.v, this);
            }

            public void DropNexValueOfF()
            {
                for (int i = 0; i < LR.N; i++)
                {
                    dx[i] = 0;
                    x2[i] = x.v[i];
                }
                f2 = f;
            }

            public void Move()
            {
                for (int i = 0; i < LR.N; i++) x.v[i] = x.v[i] + dx[i];
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
            //                                        {fish1.x.v[0]}, 
            //                                        {fish1.x.v[1]},
            //                                        {fish1.x.v[2]},
            //                                        {fish1.x.v[3]}
            //                                        //{fish1.x.v[4]},
            //                                        //{fish1.x.v[5]},
            //                                        //{fish1.x.v[6]},
            //                                        //{fish1.x.v[7]}
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


                //2N*2N
                //var A = DenseMatrix.OfArray(new[,] { { 0.0, 0.0, 1.0, 0.0 }, 
                //                                     { 0.0, 0.0, 0.0, 1.0 }, 
                //                                     { 0.0, 0.0, 0.0, 0.0 }, 
                //                                     { 0.0, 0.0, 0.0, 0.0 }});
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
                //2N*N
                //var B = DenseMatrix.OfArray(new[,] { { 0.0, 0.0 }, 
                //                                 { 0.0, 0.0}, 
                //                                 { 1.0, 0.0 }, 
                //                                 { 0.0, 1.0 }});
                //------------------------------------------------------------
                var B = DenseMatrix.Create(N, N / 2, 0.0);

                for (int i = N / 2; i < N; i++)
                {
                    for (int j = 0; j < N / 2; j++)
                    {
                        if (j == (i - N / 2)) { B[i, j] = 1; }
                        else { B[i, j] = 0; }
                    }
                }

                //------------------------------------------------------------
                //N*N
                //var Q = DenseMatrix.OfArray(new[,] { { 1.0, 0.0 }, 
                //                                 { 0.0, 1.0}});
                //------------------------------------------------------------
                var Q = DenseMatrix.Create(N / 2, N / 2, 0.0);
                for (int i = 0; i < N / 2; i++)
                {
                    for (int j = 0; j < N / 2; j++)
                    {
                        if (i == j) { Q[i, j] = 1; }
                        else { Q[i, j] = 0; }
                    }
                }
                //------------------------------------------------------------
                //2N*2N
                //var S = DenseMatrix.OfArray(new[,] {{ ks*1.0, 0.0, 0.0, 0.0 }, 
                //                                    { 0.0, ks*1.0, 0.0, 0.0 }, 
                //                                    { 0.0, 0.0, 0.0, 0.0 }, 
                //                                    { 0.0, 0.0, 0.0, 0.0 }});
                //------------------------------------------------------------
                var S = DenseMatrix.Create(N, N, 0.0);
                for (int i = 0; i < (N - N / 2); i++)
                {
                    for (int j = 0; j < N / 2; j++)
                    {
                        if (i == j) { S[i, j] = ks; }
                        else { S[i, j] = 0; }
                    }
                }

                //------------------------------------------------------------
                //Lambda!!!!!!!!!
                //2N*2N
                //var P = (MathNet.Numerics.LinearAlgebra.Matrix<double>)DenseMatrix.OfArray(new[,] { { kl*1.0, 0.0,    0.0, 0.0 }, 
                //                                                                                    { 0.0,    kl*1.0, 0.0, 0.0 }, 
                //                                                                                    { 0.0,    0.0,    0.0, 0.0 }, 
                //                                                                                    { 0.0,    0.0,    0.0, 0.0 }});   
                //------------------------------------------------------------
                var P = (MathNet.Numerics.LinearAlgebra.Matrix<double>)DenseMatrix.Create(N, N, 0.0);
                for (int i = 0; i < N / 2; i++)
                {
                    for (int j = 0; j < N / 2; j++)
                    {
                        if (i == j) { P[i, j] = kl; }
                        else { P[i, j] = 0; }
                    }
                }
                //------------------------------------------------------------
                ////var P = (MathNet.Numerics.LinearAlgebra.Matrix<double>)DenseMatrix.OfArray(new[,] { { 0.0, 0.0, 0.0, 0.0 }, 
                ////                                                                                { 0.0, 0.0, 0.0, 0.0 }, 
                ////                                                                                { 0.0, 0.0, 0.0, 0.0 }, 
                ////                                                                                { 0.0, 0.0, 0.0, 0.0 }});
                if (group == 1)
                {
                    //NMAX это K!!!!!
                    for (int j = K; j > 0; j--)
                    {
                        var P2 = -A.Transpose() * P * h - P * A * h + P * B * Q.Inverse() * B.Transpose() * P * h - S * h;
                        P = P2;
                    }
                }
                else if (group == 2)
                {
                    //алг ур риккати
                    //-A.Transpose() * P - P * A + P * B * Q.Inverse() * B.Transpose() * P  - S = 0   |*(-1)
                    //A.Transpose() * P + P * A - P * B * Q.Inverse() * B.Transpose() * P  + S = 0   
                    var S1 = B * Q.Inverse() * B.Transpose();



                    //var J = DenseMatrix.Create(8, 8, 0.0);
                    //for (int i = 0; i < 4; i++)
                    //{
                    //    for (int j = 4; j < 8; j++)
                    //    {
                    //        if (i == (j - 4)) { J[i, j] = 1; }
                    //        else { J[i, j] = 0; }
                    //    }
                    //}
                    //for (int i = 4; i < 8; i++)
                    //{
                    //    for (int j = 0; j < 4; j++)
                    //    {
                    //        if ((i-4) == j) { J[i, j] = -1; }
                    //        else { J[i, j] = 0; }
                    //    }
                    //}
                    // 0  I
                    //-I  0

                    var J = DenseMatrix.Create(2 * N, 2 * N, 0.0);
                    for (int i = 0; i < N; i++)
                    {
                        for (int j = N; j < 2 * N; j++)
                        {
                            if (i == (j - N)) { J[i, j] = 1; }
                        }
                    }
                    for (int i = N; i < 2 * N; i++)
                    {
                        for (int j = 0; j < N; j++)
                        {
                            if ((i - N) == j) { J[i, j] = -1; }
                        }
                    }


                    //var W = -S   -A.Transpose()
                    //        -A   S1
                    //var W = (MathNet.Numerics.LinearAlgebra.Matrix<double>)DenseMatrix.Create(8, 8, 0.0);
                    //for (int i = 0; i < 4; i++)
                    //{
                    //    for (int j = 0; j < 4; j++)
                    //    {
                    //        W[i,j] = -S[i,j];
                    //    }
                    //    for (int j = 4; j < 8; j++)
                    //    {
                    //        W[i, j] = -A.Transpose()[i,j-4];
                    //    }
                    //}

                    //for (int i = 4; i < 8; i++)
                    //{
                    //    for (int j = 0; j < 4; j++)
                    //    {
                    //        W[i, j] = -A[i-4, j];
                    //    }
                    //    for (int j = 4; j < 8; j++)
                    //    {
                    //        W[i, j] = S1[i-4, j-4];
                    //    }
                    //}

                    var W = (MathNet.Numerics.LinearAlgebra.Matrix<double>)DenseMatrix.Create(2 * N, 2 * N, 0.0);
                    for (int i = 0; i < N; i++)
                    {
                        for (int j = 0; j < N; j++)
                        {
                            W[i, j] = -S[i, j];
                        }
                        for (int j = N; j < 2 * N; j++)
                        {
                            W[i, j] = -A.Transpose()[i, j - N];
                        }
                    }

                    for (int i = N; i < 2 * N; i++)
                    {
                        for (int j = 0; j < N; j++)
                        {
                            W[i, j] = -A[i - N, j];
                        }
                        for (int j = N; j < 2 * N; j++)
                        {
                            W[i, j] = S1[i - N, j - N];
                        }
                    }




                    int n_eps = 100;
                    //for (int k = 1; k < n_eps; k++)
                    //{
                    //    double c = Math.Pow(W.Determinant(), 1/8.0);
                    //    W = (1 / (2 * c)) * (W + c * c * J * W.Inverse() * J); 
                    //}

                    for (int k = 1; k < n_eps; k++)
                    {
                        double c = Math.Pow(W.Determinant(), 1 / (2 * N));
                        W = (1 / (2 * c)) * (W + c * c * J * W.Inverse() * J);
                    }

                    //сформировать матрицы М и N
                    //var M = (MathNet.Numerics.LinearAlgebra.Matrix<double>)DenseMatrix.Create(8, 4, 0.0);
                    //var NN = (MathNet.Numerics.LinearAlgebra.Matrix<double>)DenseMatrix.Create(8, 4, 0.0);
                    //for (int i = 0; i < 4; i++)
                    //{
                    //    for (int j = 0; j < 4; j++)
                    //    {
                    //        M[i, j] = W[i+4, j+4];
                    //        NN[i, j] = J[i + 4, j + 4] - W[i + 4, j];
                    //    }
                    //}
                    //for (int i = 4; i < 8; i++)
                    //{
                    //    for (int j = 0; j < 4; j++)
                    //    {
                    //        M[i, j] = W[i - 4, j + 4];
                    //        if (i - 4 == j) { M[i, j]++; }
                    //        NN[i, j] = -W[i - 4, j];
                    //    }
                    //}


                    var M = (MathNet.Numerics.LinearAlgebra.Matrix<double>)DenseMatrix.Create(4 * (N - N / 2), 2 * (N - N / 2), 0.0);
                    var NN = (MathNet.Numerics.LinearAlgebra.Matrix<double>)DenseMatrix.Create(4 * (N - N / 2), 2 * (N - N / 2), 0.0);
                    for (int i = 0; i < 2 * (N - N / 2); i++)
                    {
                        for (int j = 0; j < 2 * (N - N / 2); j++)
                        {
                            M[i, j] = W[i + 2 * (N - N / 2), j + 2 * (N - N / 2)];
                            NN[i, j] = J[i + 2 * (N - N / 2), j + 2 * (N - N / 2)] - W[i + 2 * (N - N / 2), j];
                        }
                    }
                    for (int i = 2 * (N - N / 2); i < 4 * (N - N / 2); i++)
                    {
                        for (int j = 0; j < 2 * (N - N / 2); j++)
                        {
                            M[i, j] = W[i - 2 * (N - N / 2), j + 2 * (N - N / 2)];
                            if ((i - 2 * (N - N / 2)) == j) { M[i, j]++; }
                            NN[i, j] = -W[i - 2 * (N - N / 2), j];
                        }
                    }




                    //// M * P = NN;

                    M.Svd().Solve(NN, P);


                    //var W11 = (MathNet.Numerics.LinearAlgebra.Matrix<double>)DenseMatrix.Create(4, 4, 0.0);
                    //var W12 = (MathNet.Numerics.LinearAlgebra.Matrix<double>)DenseMatrix.Create(4, 4, 0.0);
                    //var W21 = (MathNet.Numerics.LinearAlgebra.Matrix<double>)DenseMatrix.Create(4, 4, 0.0);
                    //var W22 = (MathNet.Numerics.LinearAlgebra.Matrix<double>)DenseMatrix.Create(4, 4, 0.0);

                    //for (int i = 0; i < 4; i++)
                    //{
                    //    for (int j = 0; j < 4; j++)
                    //    {
                    //        W11[i, j] = W[i, j];
                    //        W12[i, j] = W[i, j + 4];
                    //        W21[i, j] = W[i + 4, j];
                    //        W22[i, j] = W[i + 4, j + 4];
                    //    }
                    //}

                    var W11 = (MathNet.Numerics.LinearAlgebra.Matrix<double>)DenseMatrix.Create(2 * (N - N / 2), 2 * (N - N / 2), 0.0);
                    var W12 = (MathNet.Numerics.LinearAlgebra.Matrix<double>)DenseMatrix.Create(2 * (N - N / 2), 2 * (N - N / 2), 0.0);
                    var W21 = (MathNet.Numerics.LinearAlgebra.Matrix<double>)DenseMatrix.Create(2 * (N - N / 2), 2 * (N - N / 2), 0.0);
                    var W22 = (MathNet.Numerics.LinearAlgebra.Matrix<double>)DenseMatrix.Create(2 * (N - N / 2), 2 * (N - N / 2), 0.0);

                    for (int i = 0; i < 2 * (N - N / 2); i++)
                    {
                        for (int j = 0; j < 2 * (N - N / 2); j++)
                        {
                            W11[i, j] = W[i, j];
                            W12[i, j] = W[i, j + 2 * (N - N / 2)];
                            W21[i, j] = W[i + 2 * (N - N / 2), j];
                            W22[i, j] = W[i + 2 * (N - N / 2), j + 2 * (N - N / 2)];
                        }
                    }



                }
                else if (group == 3)
                {
                    for (int j = K; j > 0; j--)
                    {
                        var P2 = -A.Transpose() * P * h - P * A * h - S * h;
                        P = P2;
                    }

                }
                else if (group == 4)
                {
                    //P = Lambda

                }


                //Решение ДУ dX/dt = (A - B * Q^-1 * B^T * P) * X => X

                //!!!!!! для n>2
                //var X = (MathNet.Numerics.LinearAlgebra.Matrix<double>)DenseMatrix.OfArray(new[,] { 
                //                                    {x.v[0]}, 
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
                for (int j = 1; j <= K; j++)
                {
                    var X2 = (A - B * Q.Inverse() * B.Transpose() * P) * X;
                    X = X2;


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


        }





        //шаг 7. Агенты второй группы.

        // A, B, Q, S те же
        //решение ур риккати
        //-A.Transpose() * P - P * A + P * B * Q.Inverse() * B.Transpose() * P - S = 0
        // var X2 = (A - B * Q.Inverse() * B.Transpose() * P) * X;
        //X_new = X_best + X2;


        //шаг 8. Агенты третьей группы.

        // A, B, Q, S те же
        //var P2 = -A.Transpose() * P * h - P * A * h - S * h;
        //P = P2;
        //var X2 = (A - B * Q.Inverse() * B.Transpose() * P) * X;
        //X_new = X_best + X2;






        //решение матр ур Риккати => P


        //решение матр ду Риккати => P
        //P[k+1] = -Atr * P[k] * dt - P[k] * A * dt + P[k] * B * Q-1 * Btr * P[k] * dt - S * dt


        //Решение ДУ dX/dt = (A - Q^-1 * B^T * P) * X => X


















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

        //public void updateNeibours()
        // {
        //     int NP = Population.Count;
        //     for (int j = 0; j < NP; j++)
        //     {
        //         int n1 = 0;
        //         Vector dxBar = Vector.Zero(N);
        //         double dfBar = 0;
        //         alphaLoc[j] = Vector.Zero(N);

        //         for (int k = 0; k < NP; k++)
        //         {
        //             if (k != j)
        //             {
        //                 if (dist[j, k] < ds[j])
        //                 {
        //                     n1++;
        //                     dxBar = (Population[k].x - Population[j].x) / (mu + dist[j, k]);
        //                     dfBar += (Population[j].f - Population[k].f) / (fMin - fMax);
        //                     alphaLoc[j] = alphaLoc[j] + dfBar * dxBar;
        //                 }
        //             }
        //         }
        //         neib[j] = n1;
        //     }
        // }

        // public void updateAlphaTarget(int k)
        // {
        //     double cbest = 2*(R.NextDouble() + 1.0*k/K);
        //     //
        //     int NP = Population.Count;
        //     for (int i = 0; i < NP; i++)
        //     {
        //         double dfBar = (Population[i].f - fMax) / (fMin - fMax);
        //         fbest[i] = dfBar;
        //         Vector dxBar = (fishMax.x - Population[i].x) / (mu + dist[i, jMax]);
        //         alphaTarget[i] = cbest * dfBar * dxBar;
        //         alpha[i] = alphaLoc[i] + alphaTarget[i];
        //     }
        // }

        // public void calcNewN(int k)
        // {
        //     // ? // double omega = 0.8 * double(k)/K + 0.2 * R.NextDouble();
        //     //
        //     int NP = Population.Count;
        //     for (int i = 0; i < NP; i++)
        //     {
        //         double omega = 0.8 - (0.8 * k/K) + 0.2 * R.NextDouble();
        //         Nnew[i] = Nmax * alpha[i] + omega * Nold[i];
        //     }
        // }

        // private string xfoodLock = "";

        // public void calcXFood()
        // {
        //     xfood = Vector.Zero(N);
        //     double D = 0;
        //     //
        //     lock (xfoodLock)
        //     {
        //         int NP = Population.Count;
        //         for (int i = 0; i < NP; i++)
        //         {
        //             double fi = 1.0 / Population[i].f;
        //             D += fi;
        //             xfood = xfood + fi * Population[i].x;
        //         }
        //         xfood = xfood / D;
        //     }
        // }

        // public void calcBetaFood(int k)
        // {
        //     double ff = Func(xfood.v);
        //     //
        //     int NP = Population.Count;
        //     for (int i = 0; i < NP; i++)
        //     {
        //         Vector xj = Population[i].x;
        //         double dfBar = (Population[i].f - ff) / (fMin - ff);
        //         Vector dxBar = (xfood - xj) / (mu + xj.distance(xfood));
        //         double cfood = 2*(1 - 1.0*k/K);
        //         betafood[i] = cfood * dfBar * dxBar;
        //     }
        // }

        // public void calcBetaBest()
        // {
        //     int NP = Population.Count;
        //     for (int i = 0; i < NP; i++)
        //     {
        //         Vector xj = Population[i].x;
        //         double dfBar = (Population[i].f - fjbest[i]) / (fMin - fMax); /// dF^{j,jbest}         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! // fbest[i];
        //         Vector dxBar = (xjbest[i] - Population[i].x) / (mu + xjbest[i].distance(xj));
        //         betabest[i] = dfBar * dxBar;
        //         beta[i] = betafood[i] + betabest[i];
        //     }
        // }

        // public void callsNewF(int k)
        // {
        //     int NP = Population.Count;
        //     for (int j = 0; j < NP; j++)
        //     {
        //         double omegaf = 0.8 - (0.8 * k / K) + 0.2 * R.NextDouble();
        //         Fnew[j] = Vf * beta[j] + omegaf * Fold[j];
        //     }
        // }

        // public void steps567(int k)
        // {
        //     int NP = Population.Count;
        //     for (int j = 0; j < NP; j++)
        //     {
        //         Vector D = new Vector(N);
        //         for (int i = 0; i < N; i++) D.v[i] = Dmax * (1 - 1.0 * k / K) * (-1 + 2 * R.NextDouble());
        //         //
        //         Vector V = D + Fnew[j] + Nnew[j];
        //         xNew[j] = Population[j].x + dt * V;
        //         //
        //         for (int i = 0; i < N; i++)
        //         {
        //             if (xNew[j].v[i] < oblast[i,0] || xNew[j].v[i] >= oblast[i,1])
        //                 xNew[j].v[i] = oblast[i,0] + (oblast[i,1] - oblast[i,0]) * R.NextDouble();
        //         }
        //     }
        // }




        // public void crossover()
        // {
        //     int NP = Population.Count;





        //     for (int j = 0; j < NP; j++)
        //     {
        //         double Cr = 0.2 * (Population[j].f - fMax) / (fMin - fMax);
        //         //
        //         xCr[j] = Vector.Zero(N);
        //         if ((fMin - fMax) == 0)
        //             continue;
        //         for (int i = 0; i < N; i++)
        //         {
        //             if (R.NextDouble() < Cr)
        //             {
        //                 int r = j;
        //                 while (r == j) r = R.Next() % NP;
        //                 xCr[j].v[i] = xNew[r].v[i];
        //             }
        //             else
        //                 xCr[j].v[i] = xNew[j].v[i];
        //         }
        //         int stop = 0;
        //         if (N > 1) 
        //         if (double.IsNaN(xCr[j].v[0]) || double.IsNaN(xCr[j].v[1]))
        //             stop = 1;
        //     }
        // }

        // public void mutation()
        // {
        // }
        public string toprotocol = "";
        public void Work()//работа алгоритма
        {
            //initAlgorithm();

            for (k = 0; k < Pmax; k++)
            {
                move_1group();
                //migration();
                //frontal_search();
                //reduction_pop();
                //refill_pop();
                //Population.Sort(new FishComparer());
                //if (k > 0)
                //{
                //    Analyze();
                //    Console.WriteLine("k=" + k + " I=" + I_opt);
                //    if (I_opt < -2.77) break;
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
                pr += "Мультиагентный алгоритм на основе линейных регуляторов управления движением агентов" + "\r\n" + "\r\n";
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
            //initAlgorithm();

            //for (k = 0; k < K; k++)
            //{
            //    updateNeibours();
            //    updateAlphaTarget(k);
            //    calcNewN(k);
            //    calcXFood();
            //    calcBetaFood(k);
            //    calcBetaBest();
            //    callsNewF(k);
            //    steps567(k);
            //    crossover();
            //    mutation();
            //    updateDistances();



            //    // var fish = Population[0];
            //    /*
            //    LogIteration();
            //    if (k > 0) Food();//питание 
            //    Swiming();//плавание
            //    Reproduction();//размножение
            //    Delete();//удаление
            //    */
            //}
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
            RecalcMaxDf();
            if (!final) return;

            Mfend = Math.Round(Mf, 6);
        }

        //операция плавания

        //public void Swiming() 
        //{
        //    calcBetaBest();
        //    callsNewF(k);
        //    steps567(k);
        //}
        ////операция питания
        //public void Food() 
        //{
        //    updateNeibours();
        //    updateAlphaTarget(k);
        //    calcNewN(k);
        //    calcXFood();
        //    calcBetaFood(k);
        //}
        ////операция разножения
        //public void Reproduction() 
        //{
        //    crossover();
        //}
        ////удаление рыб с наим весом
        //public void Delete() 
        //{
        //    mutation();
        //    updateDistances();
        //    RecalcMaxDf();
        //}

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
