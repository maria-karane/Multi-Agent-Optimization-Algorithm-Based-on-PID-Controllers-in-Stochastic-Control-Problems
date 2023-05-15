using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Управление
{
    public class Task2 : BaseTask
    {
        double _t0, _t1, _x10, _x20, _dt, _c_down, _c_up, _x_delta;

        public override double t0() {return _t0;}

        public override double t1() { return _t1; }
        public override double x_delta() { return _x_delta; }


        //для разложения С
        public override double c_down() { return _c_down; }

        public override double c_up() { return _c_up; }
        //

        public override double x0(int i) { return i == 0 ? _x10 : _x20; }

        public override double Dt() { return _dt; }

        public Task2() 
        {
            //
            _c_down = 0.0;
            _c_up = 50.0;
            //
            _x_delta = 0.05;
            _t0 = 0;
            _t1 = 1.6;
            _x10 = 1;
            _x20 = 0;
            _dt = 0.01;
            _u_neg = -2;
            _u_pos = 1;
        }

        public override int n()
        {
            return 2;
        }

        public override double v(int i, double u, double[] x)
        {
            if (i == 0)
            {
                return 1/(Math.Cos(x[0])+2) + 3*Math.Sin(x[1]) + u; //Math.Sin(x[0]) + x[1] * Math.Cos(x[1]) - u * u;
            }
            else 
            {
                return x[0] + x[1] + u;//x[0] - 6 * Math.Sin(x[0]) + u;
            }
        }

        public override double I(double[] x) { return x[0] - 0.5 * x[1]; }
    }
}
