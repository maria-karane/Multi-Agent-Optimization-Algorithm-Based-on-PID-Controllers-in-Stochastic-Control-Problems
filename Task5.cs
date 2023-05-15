using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Управление
{
    public class Task5 : BaseTask
    {

        double _t0, _t1, _x10, _x20, _dt, _c_down, _c_up, _x_delta, param_cheb;

        public override double t0() {return _t0;}

        public override double t1() { return _t1; }
        public override double x_delta() { return _x_delta; }

        Random R = new Random();

        //для разложения С
        public override double c_down() { return _c_down; }

        public override double c_up() { return _c_up; }
        //
        //для разложения полиномы Чебышева
        //public override double u_cheb_min() { return _u_cheb_min; }
        //public override double u_cheb_max() { return _u_cheb_max; }
        //public override double x0(int i) { return i == 0 ? _x10 : _x20; }

        public override double Dt() { return _dt; }

        public Task5(bool chebyshev) 
        {

            _x_delta = 0.02;

           
            if (chebyshev)
            {
                //
                _c_down = -2.2;
                _c_up = 1.1;
                //
                _t0 = -1;
                _t1 = 1;
                param_cheb = 0.8;
            }
            else
            {
                _c_down = 0.0;
                _c_up = 350.0;//было 50
                _t0 = 0;
                _t1 = 1.6;
                param_cheb = 1;
            }

            //_x10 = 1;
            //_x20 = 0;
            _dt = 0.05;
            _u_neg = -2;
            _u_pos = 1;


            _x0 = new double[2];
            _x1 = new double[2];


            //-------------------------------------------------------------------------------------
            //НАЧАЛЬНЫЕ УСЛОВИЯ ЗАДАНЫ ОПРЕДЕЛЕННО
            _x1[0] = 4.5;
            _x0[0] = 1.0;
            _x1[1] = 13.0;
            _x0[1] = 0.0;


            //-------------------------------------------------------------------------------------                                     
            //НАЧАЛЬНЫЕ УСЛОВИЯ ЗАДАЮТСЯ В ВИДЕ РАЗБИЕНИЯ МНОЖЕСТВА
            //var _x0_b = new double[2];
            //var _x1_b = new double[2];
            //_x0_b[0] = 0.95;
            //_x0_b[1] = 1.05;

            //_x1_b[0] = -0.05;
            //_x1_b[1] = 0.05;

            //int n0 = (int)Math.Round((_x0_b[1] - _x0_b[0]) / _x_delta);
            //int n1 = (int)Math.Round((_x1_b[1] - _x1_b[0]) / _x_delta);
            //double dx0 = (_x0_b[1] - _x0_b[0]) / n0;
            //double dx1 = (_x1_b[1] - _x1_b[0]) / n1;
            //int d = (n0) * (n1);
            //_x0_a = new double[d];
            //_x1_a = new double[d];

            //int s = 0;
            //for (int i = 0; i < n0; i++)
            //{
            //    double x0 = (i + 0.5) * dx0 + _x0_b[0];
            //    for (int j = 0; j < n1; j++)
            //    {
            //        double x1 = (j + 0.5) * dx1 + _x1_b[0];
            //        _x0_a[s] = x0;
            //        _x1_a[s] = x1;
            //        s++;
            //    }
            //}

            //-------------------------------------------------------------------------------------
            //НАЧАЛЬНЫЕ УСЛОВИЯ ЗАДАЮТСЯ СЛУЧАЙНО С ПОМОЩЬЮ ЗАКОНА РАСПРЕДЕЛЕНИЯ

            var _x0_b = new double[2];
            var _x1_b = new double[2];
            _x0_b[0] = 0.95;
            _x0_b[1] = 1.05;

            _x1_b[0] = -0.05;
            _x1_b[1] = 0.05;

            int n0 = (int)Math.Round((_x0_b[1] - _x0_b[0]) / _x_delta);
            int n1 = (int)Math.Round((_x1_b[1] - _x1_b[0]) / _x_delta);

            int d = (n0) * (n1);
            _x0_a = new double[d];
            _x1_a = new double[d];
            int s = 0;

            for (int i = 0; i < n0; i++)
            {
                double x0 = R.NextDouble() * (_x0_b[1] - _x0_b[0]) + _x0_b[0];
                for (int j = 0; j < n1; j++)
                {
                    double x1 = R.NextDouble() * (_x1_b[1] - _x1_b[0]) + _x1_b[0];
                    _x0_a[s] = x0;
                    _x1_a[s] = x1;
                    s++;
                }
            }

            //-------------------------------------------------------------------------------------
        
        }

        public override int n()
        {
            return 2;
        }

        public override double v(int i, double u, double[] x)
        {


            if (i == 0)
            {
                return param_cheb*(1/(Math.Cos(x[0])+2) + 3*Math.Sin(x[1]) + u); //Math.Sin(x[0]) + x[1] * Math.Cos(x[1]) - u * u;
            }
            else 
            {
                return param_cheb*(x[0] + x[1] + u);//x[0] - 6 * Math.Sin(x[0]) + u;
            }
        }

        public override double I(double[] x) 
        {
            //среднее
            if (I_count)
            {
                double w = 0;
                int d = _x0_a.Length;
                for (int s = 0; s < d; s++)
                {
                    w += x[2 * s] - 0.5 * x[2 * s + 1];

                }
                return w / d;
            }
            //return x[0] - 0.5 * x[1]; 

            //гарантирующее
            if (I_count2)
            {
                int d = _x0_a.Length;
                double max_I;
                max_I = x[0] - 0.5 * x[1];
                for (int s = 1; s < d; s++)
                {
                    double m = x[2 * s] - 0.5 * x[2 * s + 1];
                    if (m > max_I) max_I = m;
                }
                //
                return max_I;
            }

            return 0;
            
        }
    }
}
