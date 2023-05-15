using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace Управление
{
    public class FishSchoolSearch
    {
        public int Ip;
        public int K;
        public int N;
        public double[,] oblast;

        
        
        //public int z;
        public ITasks task;
        public double thr;
        public double WScale;

        public double Mfend = 0;

        public List<double> AverF = new List<double>();
        public List<double> BestF = new List<double>();

        Random R = new Random();

        public int Nm;
        int m;

        public int Nf = 0;

        public Vector optimalStart;
        public Vector optimalD;

        public void Optimize(int ip, int k, int N, ITasks task_, double thr_, double wScale, double stepInd0_ = 0.1, double stepInd1_ = 0.01)
        {
            int n = task_.n();

            double delta = task_.x_delta(); //ШАГ DELTA X
            int[] M = new int[n];
            double[] step = new double[n];

            for (int i = 0; i < n; i++)
            {
                //double z1 = 
                M[i] = (int)Math.Round((task_.x1_a()[i] - task_.x0_a()[i])/delta);
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
                    var alg = new FishSchoolSearch();
                    alg.Init(ip, k, N, task_, thr_, wScale, stepInd0_, stepInd1_);
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
                        var alg = new FishSchoolSearch();
                        alg.Init(ip, k, N, task_, thr_, wScale, stepInd0_, stepInd1_);
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
        public void Init(int ip, int k, int n, ITasks task_, double thr_, double wScale, double stepInd0_ = 0.1, double stepInd1_ = 0.01)
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
                    if (i < n / 2)
                    {
                        oblast[i, 0] = task.c_down();
                        oblast[i, 1] = task.c_up();
                    }
                    else
                    {
                        oblast[i, 0] = task.c_down();
                        oblast[i, 1] = task.c_up();
                    }
                }
                else
                {
                    oblast[i, 0] = task.t0();
                    oblast[i, 1] = task.t1();
                }
            }

            //----
            Ip = ip;
            
            thr = thr_;
            thr = 0.8 * wScale; // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            WScale = wScale;
            //
            stepIndInitial = stepInd0_;
            stepIndFinal = stepInd1_;
            stepVolInitial = 0.001;
            stepVolFinal = 0.0001;
            //
            optimalStart = new Vector(n);
            optimalD = new Vector(n);
            for (int j = 0; j < n; j++) optimalStart.v[j] = -1;
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

            public double this[int i]
            {
                get
                {
                    return v[i];
                }
            }
        }

        private double Func(double[] ta, Fish fish)
        {
            Array.Sort(ta);
           // var tp = ta[0];
            //var upr = new U(ta, task.U1(), task.U2(), task.getUprType(), task.getViewCos(), task.getFST(), 
            //    task.t0(), task.t1(), task.x0_a(), task.x1_a(), task.LL());
            var _x0 = task.x0();
            var _x1 = task.x1();

            var upr = new U(ta, task.get_x_L(), task.U1(), task.U2(), task.getUprType(), task.getViewCos(), task.getChebyshev(), task.t0(), task.t1(), _x0, _x1);


            var rx = task.Solve(upr);
            if (fish != null) fish.rx = rx;
            return -task.I(rx);

            //double func = 0;
            //if (z == 0)
            //{
            //    func = x[0] * Math.Sin(Math.Sqrt(Math.Abs(x[0]))) + x[1] * Math.Sin(Math.Sqrt(Math.Abs(x[1])));
            //}
            //else { }
            //Nf++;
            //return func;
        }

    
        public class Fish
        {
            public FishSchoolSearch FSS;

            public double[] rx;
            
            private Vector _x;
            public double[] dx;
            public double W;
            public double f;
            public double f2;
            public double normf;
            //public List<Fish> Kloni;

            public double df { get { return f2 - f; } }

            public Fish(FishSchoolSearch fss)
            {
                FSS = fss;
                dx = new double[FSS.N];
                x2 = new double[FSS.N];
                W = 0.5 * FSS.WScale;
            }

            public Vector x
            {
                get { return _x; }
                set
                {
                    _x = value;
                    for (int i = 0; i <FSS.N; i++) x2[i] = _x.v[i];
                    f = FSS.Func(_x.v,this);
                    f2 = f;
                }
            }

            public double[] X { get { return x.v; } }

            public double[] x2;

            public bool X2inOblast
            {
                get
                {
                    //var x2 = X2;
                    for (int i = 0; i < FSS.N; i++)
                        if (x2[i] < FSS.oblast[i, 0] || x2[i] > FSS.oblast[i, 1]) return false;
                    return true;
                }
            }

            public void CalcNextValueOfF()
            {
                f2 = FSS.Func(x2, this);
            }

            public void DropNexValueOfF()
            {
                for (int i = 0; i < FSS.N; i++)
                {
                    dx[i] = 0;
                    x2[i] = x.v[i];
                }
                f2 = f;
            }

            public void Move()
            {
                for (int i = 0; i < FSS.N; i++) x.v[i] = x.v[i] + dx[i];
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

       
        //создание начальной популяции

        public void Create_Pop()
        {
            //создаем рыбок заданное количество раз
            //заносим их в список
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

                Population.Add(member);
                // Mf += member.f;     // ???
            }

            if ((File.Exists("protocol.dt")))
            {
                FileStream fs = new FileStream("protocol.dt", FileMode.Append, FileAccess.Write);
                StreamWriter r1 = new StreamWriter(fs);


                string pr = "";
                pr += "\r\n" + "\r\n";
                pr += "Метод, имитирующий поведение стаи рыб" + "\r\n" + "\r\n";
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



        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // основной цикл

        /* F SRDF SRDF SRDF ... SRD   (F)
        */
        public string toprotocol = "";
        public void Work()//работа алгоритма
        {
            for (k = 0; k < K; k++)
            {
                //LogIteration();
                if (k > 0) Food();//питание 
                Swiming();//плавание
                Reproduction();//размножение
                Delete();//удаление
            }

            if ((File.Exists("protocol.dt")))
            {
                FileStream fs = new FileStream("protocol.dt", FileMode.Append, FileAccess.Write);
                StreamWriter r1 = new StreamWriter(fs);


                string pr = "";
                pr += "\r\n" + "\r\n";
                pr += "Метод, имитирующий поведение стаи рыб" + "\r\n" + "\r\n";
                pr += "КОНЕЧНАЯ ПОПУЛЯЦИЯ" + "\r\n" + "\r\n";
                pr += "|-----------------|--------------------------------------------------|----------------------------|\r\n";
                pr += "|   Номер агента  |                 Коэффициенты                     |         Критерий           |\r\n";
                pr += "|-----------------|--------------------------------------------------|----------------------------|\r\n";

                toprotocol = pr;
                r1.Write(toprotocol);

                r1.Close();
                fs.Close();
            }

            Analyze(true);

            //if ((File.Exists("protocol.dt")))
            //{              
            //    FileStream fs = new FileStream("protocol.dt", FileMode.Append, FileAccess.Write);
            //    StreamWriter r1 = new StreamWriter(fs);
            //    //Population[0].x.v[0]

            //    string pr = "";
            //    pr += "\r\n" + "\r\n";
            //    pr += "Метод, имитирующий поведение стаи рыб" + "\r\n" + "\r\n";
            //    pr += "КОНЕЧНАЯ ПОПУЛЯЦИЯ" + "\r\n" + "\r\n";
            //    pr += "|-----------------|--------------------------------------------------|----------------------------|\r\n";
            //    pr += "|   Номер агента  |                 Коэффициенты                     |         Критерий           |\r\n";
            //    pr += "|-----------------|--------------------------------------------------|----------------------------|\r\n";


            //    for (int j = 0; j < Ip; j++)
            //    {
            //        pr += "|        " + (j + 1) + "        |";
            //        for (int i = 0; i < N; i++)
            //        {

            //            pr += Math.Round(Population[j].x.v[i], 2) + " ";
            //        }
            //        pr += "|    " + (-Math.Round(Population[j].f, 4)) + "     |" + "\r\n";
                   
            //    }

            //    toprotocol = pr;
            //    r1.Write(toprotocol);
                
            //    r1.Close();
            //    fs.Close();
            //}
        }

        public double[] x_opt = new double[2];
        public double I_opt;
        public double[] traj;
        public int _numb_tr;

        public double[] ta;


        void Analyze(bool final = false)
        {
            double i0 = 1000000;
            Fish f_opt = null;

           
            int j = 1;
            foreach (Fish member in Population)
            {
                
                double v = task.I(member.rx);
                if ( v < i0)
                {
                    I_opt = i0 = v;
                    x_opt = member.rx;
                    f_opt = member;
                    
                }
                if ((File.Exists("protocol.dt")))
                {
                    FileStream fs = new FileStream("protocol.dt", FileMode.Append, FileAccess.Write);
                    StreamWriter r1 = new StreamWriter(fs);
                   

                    string pr = "";
                    //pr += "\r\n" + "\r\n";
                    //pr += "Метод, имитирующий поведение стаи рыб" + "\r\n" + "\r\n";
                    //pr += "КОНЕЧНАЯ ПОПУЛЯЦИЯ" + "\r\n" + "\r\n";
                    //pr += "|-----------------|--------------------------------------------------|----------------------------|\r\n";
                    //pr += "|   Номер агента  |                 Коэффициенты                     |         Критерий           |\r\n";
                    //pr += "|-----------------|--------------------------------------------------|----------------------------|\r\n";


                    
                        pr += "|        " + j + "        |";

                        for (int i = 0; i < N; i++)
                        {                           
                            pr += Math.Round(f_opt.x.v[i], 2) + " ";
                        }
                        pr += "         |    " + Math.Round(I_opt, 4) + "             |" + "\r\n";

                    

                    toprotocol = pr;
                    r1.Write(toprotocol);

                    j++;
                    r1.Close();
                    fs.Close();
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
            var sumDf = 0.0;
            var sumDx = new double[N]; // все нули

            //  а)
            foreach (Fish member in Population)
            {
                // z)
                var __x2 = member.x2;
                for (int i = 0; i < N; i++)
                    member.x.v[i] = __x2[i];

                // ----------------------------------------------------------------------------------
                for (int i = 0; i < N; i++)
                {
                    double delta = -0.5 + R.NextDouble();
                    member.dx[i] = delta * stepInd;
                    if (i >= 2) member.dx[i] *= 0.3;
                    member.x2[i] = member.x.v[i] + member.dx[i];
                }
                if (!member.X2inOblast)
                {
                    for (int i = 0; i < N; i++)
                    {
                        member.dx[i] = 0;
                        member.x2[i] = member.x[i];
                    }
                }
                //
                member.CalcNextValueOfF();
                if (member.df < 0)
                    member.DropNexValueOfF();
                //
                if (member.df < 0)
                    Console.WriteLine("!");
                //
                if (member.df > 0)
                {
                    sumDf += member.df;
                    for (int i = 0; i < N; i++)
                        sumDx[i] += member.dx[i] * member.df;
                }
            }
            // b)
            if (sumDf > 0)
            {
                var dy = new double[N];
                for (int i = 0; i < N; i++)
                    dy[i] = sumDx[i] / sumDf;
                //
                foreach (Fish member in Population)
                {
                    for (int i = 0; i < N; i++)
                    {
                        ////// ОТКЛЮЧАЕМ! // 
                        member.dx[i] += dy[i];
                        member.x2[i] = member.x.v[i] + member.dx[i];
                    }


                    if (!member.X2inOblast)
                    {
                        for (int i = 0; i < N; i++)
                        {
                            member.dx[i] = 0;
                            member.x2[i] = member.x[i];
                        }
                    }


                    member.CalcNextValueOfF();
                    //if (member.df < 0)
                    //    member.DropNexValueOfF();
                    //if (member.df < 0)
                    //{
                    //    Console.WriteLine("!");
                    //}
                }
            }
            // c)
            var Bari = new double[N]; // нули?
            var sumW = 0.0;
            foreach (Fish member in Population)
            {
                for (int i = 0; i < N; i++)
                {
                    Bari[i] += member.x2[i] * member.W;
                }
                sumW += member.W;
            }
            if (sumW > 0)
                for (int i = 0; i < N; i++) Bari[i] /= sumW;
            var s = dW > 0 ? -1 : 1;
            // завершаем: 
            // if (false)
            foreach (Fish member in Population)
            {
                var u = new Vector(N);
                for (int i = 0; i < N; i++)
                {
                    u.v[i] = member.x2[i] - Bari[i];
                }
                u.normalize();
                var _r = R.NextDouble();
                for (int i = 0; i < N; i++)
                {
                    // ОТКЛЮЧАЕМ! // 
                    member.dx[i] += u.v[i] * stepVol * _r;
                    member.x2[i] = member.x.v[i] + member.dx[i];
                }

                if (!member.X2inOblast)
                {
                    for (int i = 0; i < N; i++)
                    {
                        member.dx[i] = 0;
                        member.x2[i] = member.x[i];
                    }
                }



                member.CalcNextValueOfF();
                if (member.df < 0)
                    member.DropNexValueOfF();
                if (member.df < 0)
                    Console.WriteLine("!");
            }

            
            
            
            RecalcMaxDf();

        }
        //операция питания
        public void Food()
        {
            if (maxDf <= 0)
                return;

            dW = 0.0;
            foreach (Fish member in Population)
            {
                //if (member.df < 0)
                //    Console.WriteLine("!");

                int stop = 0;
                if (maxDf > 0.5)
                    stop = 1;

                double a0 = 0.1 * 0.5 * (maxDf + 0.33 * 0.5 * WScale / K);
                a0 = 2 * (WScale - Population[0].W) / WScale;
                if (a0 <= 0)
                    a0 = 1;
                //
                var newW = member.df > 0 ? member.W + member.df / maxDf : member.W;
                if (newW > WScale) newW = WScale;
                dW += newW - member.W;
                member.W = newW;
            }
        }
        //операция разножения
        public void Reproduction()
        {
            var L = new List<Fish>();
            foreach (Fish member in Population)
            {
                if (member.W >= thr)
                {
                    member.pairIndex = -1;
                    L.Add(member);
                }
            }
            var n = L.Count;
            while (n > 1)
            {
                var k1 = R.Next(L.Count);
                var fish = L[k1];
                if (fish.pairIndex >= 0) continue;
                //
                double m = 0;
                int k2 = -1;
                for (int k = 0; k < L.Count; k++)
                {
                    var fish2_ = L[k];
                    if (fish2_.pairIndex >= 0) continue;
                    var m_ = fish2_.W / fish2_.x.distance(fish.x);
                    if (m_ > m)
                    {
                        m = m_;
                        k2 = k;
                    }
                }
                var fish2 = L[k2];
                //
                fish.pairIndex = k2;
                fish2.pairIndex = k1;
                //
                // размножение
                var child = new Fish(this);
                var y = new Vector(N);
                for (int i = 0; i < N; i++) y.v[i] = 0.5 * (fish.x.v[i] + fish2.x.v[i]);
                child.x = y;
                child.W = 0.5 * (fish.W + fish2.W);
                Population.Add(child);
                //
                n -= 2;
            }
        }
        //удаление рыб с наим весом
        public void Delete()
        {
            Population.Sort(delegate(Fish p1, Fish p2) { return p2.W.CompareTo(p1.W); });
            //
            int r = -1;
            for (int j = 0; j < Population.Count; j++)
            {
                if (Population[j].W == WScale) r = j;
            }
            if (r >= 1)
            {
                var fish0 = Population[0];
                var j1 = 1 + (R.Next() % r);
                Population[0] = Population[j1];
                Population[j1] = fish0;
            }
            //
            int stop = 0;
            if (Population[0].W >= this.WScale)
                stop = 1;
            //
            if (Population.Count > Ip)
                Population.RemoveRange(Ip, Population.Count - Ip);
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
