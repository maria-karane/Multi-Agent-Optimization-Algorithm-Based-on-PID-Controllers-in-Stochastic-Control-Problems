using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace Управление
{
    class Imperialist
    {
        public int Ip;
        public int K;
        public int N;
        public double[,] oblast;

        public ITasks task;

        public int z;
        public double thr;
        public double WScale;
        //
        public int N_imp;
        public double Beta;

        public double Ksi;  // Xi
        public double Gamma;
        //
        public double Mfend = 0;

        public List<double> AverF = new List<double>();
        public List<double> BestF = new List<double>();

        Random R = new Random();

        public int Nm;
        int m;

        public int Nf = 0;

        public Vector optimalStart;
        public Vector optimalD;


        public void Optimize(int ip, int k, int N, ITasks task_, int n_imp, double beta, double gamma, double ksi)
        {
            int n = task_.n();

            double delta = task_.x_delta(); //ШАГ DELTA X
            int[] M = new int[n];
            double[] step = new double[n];

            for (int i = 0; i < n; i++)
            {
                //double z1 = 
                M[i] = (int)Math.Round((task_.x1_a()[i] - task_.x0_a()[i]) / delta);
                step[i] = (task_.x1_a()[i] - task_.x0_a()[i]) / M[i];
            }

            if (n == 1)
            {
                int MM = M[0];
                double[] I_opt_arr = new double[MM];
                Trajectory[] traj_arr = new Trajectory[MM];

                Trajectory[] x_opt_arr = new Trajectory[MM];
                Trajectory[] coeff_arr = new Trajectory[MM];

                double max_I = 0;
                int X1 = 0, Y1 = 0;

                for (int X = 0; X < M[0]; X++)
                {
                    double[] x = new double[n];

                    x[0] = task_.x0_a()[0] + (0.5 + X) * step[0];


                    task_.setPosition(x);
                    var alg = new Imperialist();
                    alg.Init(ip, k, N, task_, n_imp, beta, gamma, ksi);
                    alg.Create_Pop();
                    alg.Work();
                    I_opt_arr[X] = alg.I_opt;
                    traj_arr[X] = new Trajectory(alg.traj);

                    x_opt_arr[X] = new Trajectory(alg.x_opt);
                    coeff_arr[X] = new Trajectory(alg.ta);


                    if (X == 0) max_I = alg.I_opt;
                    else
                    {
                        if (alg.I_opt < max_I)
                        {
                            max_I = alg.I_opt;
                            X1 = X;

                        }
                    }

                }

                this.I_opt = max_I;
                this.traj = traj_arr[X1].data;
                this.x_opt = x_opt_arr[X1].data;
                this.ta = coeff_arr[X1].data;
            }

            else if (n == 2)
            {
                int MM = M[0] * M[1];
                double[] I_opt_arr = new double[MM];
                Trajectory[] traj_arr = new Trajectory[MM];

                Trajectory[] x_opt_arr = new Trajectory[MM];
                Trajectory[] coeff_arr = new Trajectory[MM];

                double max_I = 0;
                int X1 = 0, Y1 = 0;

                for (int X = 0; X < M[0]; X++)
                {
                    double[] x = new double[n];

                    x[0] = task_.x0_a()[0] + (0.5 + X) * step[0];

                    for (int Y = 0; Y < M[1]; Y++)
                    {
                        x[1] = task_.x0_a()[1] + (0.5 + Y) * step[1];

                        task_.setPosition(x);
                        var alg = new Imperialist();
                        alg.Init(ip, k, N, task_, n_imp, beta, gamma, ksi);
                        alg.Create_Pop();
                        alg.Work();
                        I_opt_arr[M[0] * Y + X] = alg.I_opt;
                        traj_arr[M[0] * Y + X] = new Trajectory(alg.traj);

                        x_opt_arr[M[0] * Y + X] = new Trajectory(alg.x_opt);
                        coeff_arr[M[0] * Y + X] = new Trajectory(alg.ta);


                        if (X == 0 && Y == 0) max_I = alg.I_opt;
                        else
                        {
                            if (alg.I_opt < max_I)
                            {
                                max_I = alg.I_opt;
                                X1 = X;
                                Y1 = Y;
                            }
                        }
                    }
                }

                this.I_opt = max_I;
                this.traj = traj_arr[M[0] * Y1 + X1].data;
                this.x_opt = x_opt_arr[M[0] * Y1 + X1].data;
                this.ta = coeff_arr[M[0] * Y1 + X1].data;
            }

            else
            { throw new Exception("Unsupported dimension: " + n.ToString()); }

        }


        public void Init(int ip, int k, int n, ITasks task_,
            /* параметры настройки алгоритма: */ int n_imp, double beta, double gamma, double ksi)
        {
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
            //z = z_;
            //
            N_imp = n_imp;
            Beta = beta;

            Ksi = ksi;
            Gamma = gamma;

            //
            stepVolInitial = 0.001;
            stepVolFinal = 0.0001;
            //
            optimalStart = new Vector(n);
            optimalD = new Vector(n);
            for (int j = 0; j < n; j++) optimalStart.v[j] = -1;
        }

        //double mu = 3;
        //double Nmax = 0.01;
        //double NP = 40;
        //double Vf = 0.02;
        //double Dmax = 0.005;
        //double ci = 0.2;





        private double Func(double[] ta, Fish fish)
        {
            Array.Sort(ta);
            // var tp = ta[0];
            var upr = new U(ta, task.get_x_L(), task.U1(), task.U2(), task.getUprType(), task.getViewCos(), task.getChebyshev(),
                task.t0(), task.t1(), task.x0(), task.x1());

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

            public double norm2()
            {
                var S = 0.0;
                for (int i = 0; i < v.Length; i++)
                {
                    var d = v[i];
                    S += d * d;
                }
                return S;
            }

            public double norm()
            {
                return Math.Sqrt(norm2());
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

            public static double operator *(Vector x, Vector y)
            {
                double s = 0;
                for (int i = 0; i < x.v.Length; i++) s += x.v[i] * x.v[i];
                return s;
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

            static Random R0 = new Random();

            private static double N01()
            {
                double S = 0;
                for (int i = 0; i < 9; i++) S += -1 + 2 * R0.NextDouble();
                return S / 3.0;
            }

            public static Vector RandomGaussian(int n)
            {
                Vector v = new Vector(n);
                for (int k = 0; k < v.v.Length; k++)
                {
                    v.v[k] = N01();
                }
                return v;
            }

            public static Vector RandomSphere(int n)
            {
                var u = RandomGaussian(n);
                u.normalize();
                return u;
            }

            public static Vector RandomOrthogonal(Vector h)
            {
                var u = RandomGaussian(h.v.Length);
                u = u - ((u * h) / (h * h)) * h;
                u.normalize();
                return u;
            }
        }

        public class Fish
        {
            public Imperialist I;

            public double[] rx;

            public int impIndex = 0;

            private Vector _x;
            public double[] dx;
            public double W;
            public double f;
            public double f2;
            public double normf;
            //public List<Fish> Kloni;

            public double df { get { return f2 - f; } }

            public Fish(Imperialist i)
            {
                I = i;
                dx = new double[I.N];
                x2 = new double[I.N];
                W = 0.5 * I.WScale;
            }

            public Vector x
            {
                get { return _x; }
                set
                {
                    _x = value;
                    for (int i = 0; i < I.N; i++) x2[i] = _x.v[i];
                    f = I.Func(_x.v, this);
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
                    for (int i = 0; i < I.N; i++)
                        if (x2[i] < I.oblast[i, 0] || x2[i] > I.oblast[i, 1]) return false;
                    return true;
                }
            }

            public void CalcNextValueOfF()
            {
                f2 = I.Func(x2, this);
            }

            public void DropNexValueOfF()
            {
                for (int i = 0; i < I.N; i++)
                {
                    dx[i] = 0;
                    x2[i] = x.v[i];
                }
                f2 = f;
            }

            public void Move()
            {
                for (int i = 0; i < I.N; i++) x.v[i] = x.v[i] + dx[i];
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
                return Population.Count;
            }
            set
            {
                //throw new NotImplementedException();
                k = value;
            }
        }

        
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
            Analyze_Imp();
        }

        public double[] x_opt = new double[2];
        public double I_opt;

        void Analyze_Imp(bool final = false)
        {
            double i0 = 1000;
            Fish f_opt = null;

            foreach (Fish member in Population)
            {
                if (task.I(member.rx) < i0)
                {
                    I_opt = i0 = task.I(member.rx);
                    x_opt = member.rx;
                    f_opt = member;
                }

            }

            if (f_opt != null)
            {
                var ta = f_opt.x.v;
                Array.Sort(ta);

                var upr = new U(ta, task.get_x_L(), task.U1(), task.U2(), task.getUprType(), task.getViewCos(), task.getChebyshev(),
                    task.t0(), task.t1(), task.x0(), task.x1());
                traj = task.SolveT(upr);
               // Console.WriteLine(ta);
            }

            RecalcMaxDf_Imp();
            if (!final) return;

            Mfend = Math.Round(Mf, 6);
        }

        private void RecalcMaxDf_Imp()
        {
            maxDf = 0;
            var S = 0.0;
            for (int j = 0; j < N_imp; j++)
            {
                var member = Empire(j);
                if (member.df > maxDf) maxDf = member.df;
                S += member.f;
            }
            BestF.Add(maxDf);
            Mf = S / N_imp;
            AverF.Add(Mf);
        }




        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // основной цикл

        int N_col { get { return Population.Count - N_imp; } }

        // Vector[] xImp;
        // double[] fxImp;
        // double[] tfxImp;
        // double[] P;
        int[] Nc;

        void initAlgorithm()
        {
            int NP = Population.Count;
            int NP0 = NP;
            //
            // P = new double[N_imp];
            // xImp = new Vector[N_imp];
            // fxImp = new double[N_imp];
            // tfxImp = new double[N_imp];
            Nc = new int[N_imp];
        }
        public string toprotocol = "";
        public void Work()//работа алгоритма
        {
            initAlgorithm();

            Step4();

            // for (k = 0; k < K; k++)
            while (true)
            {
                Step5();
                Step67();
                printColonies();
                if (N_imp == 1) break;
            }
            if (N_imp > 1)
                throw new Exception("Не хватило итераций, N_imp = " + N_imp);
            Analyze(true);

            if ((File.Exists("protocol.dt")))
            {
                FileStream fs = new FileStream("protocol.dt", FileMode.Append, FileAccess.Write);
                StreamWriter r1 = new StreamWriter(fs);
                //Population[0].x.v[0]

                string pr = "";
                pr += "\r\n" + "\r\n";
                pr += "Метод, имитирующий империалистическую конкуренцию" + "\r\n" + "\r\n";
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

        class FishComparer : IComparer<Fish>
        {
            public int Compare(Fish x, Fish y)
            {
                if (x.f < y.f) return 1;
                if (x.f > y.f) return -1;
                return 0;
            }
        }

        void Step3()
        {
            Population.Sort(new FishComparer());
        }

        public void Step4()
        {
            int NP = Population.Count;

            Step3();

            var tfxImp = new double[N_imp];
            for (int j = 0; j < N_imp; j++) tfxImp[j] = Empire(j).f - MainEmpire.f;

            double S1 = 0;
            for (int j = 0; j < N_imp; j++) S1 += tfxImp[j];

            if (Math.Abs(S1) < 0.000000001) S1 = S1 < 0 ? -1 : 1;

            var P = new double[N_imp];
            for (int j = 0; j < N_imp; j++) P[j] = Math.Abs(tfxImp[j] / S1);

            int N_col_last = N_col;
            int N_col_R = N_col - 2;
            for (int j = 1; j < N_imp; j++)
            {
                Nc[j] = (int)Math.Floor(P[j] * N_col_R);
                if (Nc[j] == 0) Nc[j] = 1;
                if (Nc[j] == 0)
                    throw new Exception("Империя " + j + " без колоний!");
                N_col_last -= Nc[j];
            }
            Nc[0] = N_col_last;
            if (N_col_last <= 0)
                throw new Exception("Империя " + (N_imp - 1) + " без колоний! " + N_col_last);

            // 4.4

            for (int k = 0; k < NP; k++) Population[k].impIndex = -1;

            int m = N_col;
            int[] H = new int[N_imp];
            for (int j = 0; j < N_imp; j++) H[j] = 0;
            while (m > 0)
            {
                int j0 = R.Next() % N_imp;
                if (H[j0] == Nc[j0]) continue;
                int k1 = (R.Next() % (NP - N_imp)) + N_imp;
                var country = Population[k1];
                if (country.impIndex != -1) continue;
                // oK:
                country.impIndex = j0;
                H[j0]++;
                m--;
            }

            printColonies();
            //Console.WriteLine("----");
        }

        void testColonies()
        {
            int NP = Population.Count;

            for (int j = 1; j < N_imp; j++)
            {
                var m = Nc[j];
                var x = 0;
                for (int k = 0; k < NP; k++)
                {
                    var colony = Population[k];
                    if (colony.impIndex == j) x++;
                }
                if (x != m)
                    throw new Exception("Empire " + j + ": Nc = " + m + " (" + x + ")");
            }
        }

        int getColony(int j)//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        {
            int NP = Population.Count;

            testColonies();

            if (Nc[j] <= 0)
                throw new Exception("Empire " + j + " Nc = " + Nc[j]);

            if (Nc[j] == 1)
            {
                //Console.Write("/one colony left/");
                for (int k2 = 0; k2 < NP; k2++)
                    if (Population[k2].impIndex == j) return k2;
            }

            while (true)
            {
                int k1 = (R.Next() % (NP - N_imp)) + N_imp;
                if (Population[k1].impIndex != j) continue;
                return k1;
            }
        }

        string colonyText(int j)
        {
            return "{" + j + ":" + Nc[j] + "}";
        }

        public void Step5()
        {
           // Console.Write("[step5: ");
            for (int j = 0; j < N_imp; j++)
            {
            //    Console.Write(colonyText(j));
                // var I = Population[...];

                int m = 0; // сколько сдвигов к данному моменту
                while (true) // начало 5.1
                {
                    int stop = 0;

                    int k1 = getColony(j);
                    var colony = Population[k1]; // *** 
                    var V1 = Empire(j).x - colony.x; // *** ЗДЕСЬ ОСТАНОВИЛИСЬ (1)
                    if (V1.norm() < 0.000001)
                        stop = 1;
                    V1.normalize();              // ***
                    if (V1.norm() < 1)
                        V1 = Vector.RandomSphere(N);
                    var V2 = Vector.RandomOrthogonal(V1);
                    var dc = (xEmpire(j) - colony.x).norm();
                    //var gamma = Math.PI / 4; // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    var xNew = colony.x + dc * Beta * (R.NextDouble() * V1 + Math.Tan(-Gamma * (-1 + 2 * R.NextDouble())) * V2);

                  //  if (Double.IsNaN(colony.x[0]) || Double.IsNaN(colony.x[1]) || Double.IsNaN(xNew[0]) || Double.IsNaN(xNew[1]))
                   //     stop = 1;
                    colony.x = xNew;
                    var fNew = Func(xNew.v, null);

                    // !!! // colony.x = xNew;
                    m++;

                    if (fNew <= colony.f)
                    {
                    }
                    else
                    {
                        Population[k1] = Empire(j);
                        Population[j] = colony;
                        Population[k1].impIndex = j;
                        Empire(j).impIndex = -1;
                    }

                    if (m >= Nc[j])
                        break;
                }
            }

            var s0 = 0;
            for (int k = 0; k < NP; k++)
                if (Population[k].impIndex != -1) s0++;

         //   Console.WriteLine(" I=" + N_imp + " C=" + s0 + " NP=" + NP + " (" + (NP - N_imp - s0) + ")]");
        }

        Fish Empire(int j)
        {
            return Population[j];
        }

        Fish MainEmpire { get { return Empire(0); } }

        Vector xEmpire(int j)
        {
            return Empire(j).x;
        }

        void Step67()
        {
           // Console.Write("[step67: ");

            int j0, j1;
            Step6(out j0, out j1);
            Step7(j0, j1);
        }

        public void Step6(out int out_j0, out int out_j1)
        {
            int NP = Population.Count;

            int stop = 0;

            var TC = new double[N_imp];
            var NTC = new double[N_imp];
            var P = new double[N_imp];
            double mTC = 0;
            for (int j = 0; j < N_imp; j++)
            {
                if (Nc[j] <= 0)
                    stop = 1;

                var s = 0.0;
                for (int k1 = N_imp; k1 < NP; k1++)
                {
                    if (Population[k1].impIndex != j) continue;
                    s += Population[k1].f;
                }

                if (s <= 0.00000001)
                    stop = 1;

                TC[j] = Empire(j).f + Ksi * s / Nc[j];

                if (j == 0) mTC = TC[j];
                else
                {
                    if (TC[j] < mTC) mTC = TC[j];
                }
            }

            double sNTC = 0;
            for (int j = 0; j < N_imp; j++)
            {
                NTC[j] = TC[j] - mTC;
                sNTC += Math.Abs(NTC[j]);
            }
            if (sNTC < 0.0000001) sNTC = 1;

            var j0 = 0;
            double minP = 0;
            for (int j = 0; j < N_imp; j++)
            {
                P[j] = Math.Abs(NTC[j]) / sNTC;
                if (j == 0)
                {
                    minP = P[j];
                    j0 = j;
                }
                else
                {
                    if (P[j] < minP)    // min OR max
                    {
                        minP = P[j];
                        j0 = j;
                    }
                }
            }

          //  Console.Write("W");

            double _f = 0;
            bool found1 = false;
            var weakColI = -1;
            for (int k1 = N_imp; k1 < NP; k1++)
            {
                var colony = Population[k1];
                if (colony.impIndex != j0) continue;
                if (!found1)
                {
                    _f = colony.f;
                    weakColI = k1;
                    found1 = true;
                }
                if (colony.f > _f)
                {
                    _f = colony.f;
                    weakColI = k1;
                }
            }
            if (weakColI == -1)
                throw new Exception("weakColI == -1");

           // Console.Write("/weak " + weakColI + "/");

            int j1 = -1;
            double _maxD = -10000;
            for (int j = 0; j < N_imp; j++)
            {
                var D = P[j] + R.NextDouble();
                if (D > _maxD)
                {
                    _maxD = D;
                    j1 = j;
                }
            }

            Nc[j0]--;
            Population[weakColI].impIndex = j1;
            Nc[j1]++;

            out_j0 = j0;
            out_j1 = j1;
        }

        public void Step7(int j0, int j1)
        {

            int nColJ0 = Nc[j0];

            // (страны).................................|012345678      0..8
            //                                          | ^     *
            // (страны).................................1|02345678      0..7
            //                                            01234567

           // Console.Write("R");

            if (nColJ0 == 0) // распад империи j0
            {
                Empire(j0).impIndex = j1;
                Nc[j1] += 1;
                //
                // swap империй 0 <--> j0
                //
                if (j0 < N_imp - 1)
                {
                    var tmp = Population[N_imp - 1];
                    Population[N_imp - 1] = Population[j0];
                    Population[j0] = tmp;
                    //
                    for (int k1 = 0; k1 < NP; k1++)
                    {
                        var colony = Population[k1];
                        if (colony.impIndex == N_imp - 1) colony.impIndex = j0;
                    }
                    Nc[j0] = Nc[N_imp - 1];
                }
                //
                //for (int k1 = N_imp; k1 < NP; k1++)
                //{
                //    var colony = Population[k1];
                //    colony.impIndex--;
                //}
                //
                // переформатирование списков (Nc)
                //
                //for (int j = 0; j < N_imp - 1; j++)
                //{
                //    Nc[j] = Nc[j + 1];
                //}
                //
                N_imp--;
            }
         //   Console.WriteLine("]");

            //printColonies();

            RecalcMaxDf();
        }

        void printColonies()
        {
            for (int k2 = 0; k2 < NP; k2++)
            {
                var colony = Population[k2];
                //if (colony.impIndex == -1) Console.Write("I");
                //else Console.Write("" + colony.impIndex);
            }
           // Console.WriteLine("");
        }


















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

        public void initAlgorithm_1()
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

        void updateDistances()
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
        }

        void updateNeibours()
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

        public void Work_1()//работа алгоритма
        {
            initAlgorithm_1();

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
                var upr = new U(ta, task.get_x_L(), task.U1(), task.U2(), task.getUprType(), task.getViewCos(), task.getChebyshev(),
                    task.t0(), task.t1(), task.x0(), task.x1());
                traj = task.SolveT(upr);
                _numb_tr = task.numb_tr();
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
                //Console.Write(fish.W + "  ");
                //Console.WriteLine();
            }
        }
    }
}
