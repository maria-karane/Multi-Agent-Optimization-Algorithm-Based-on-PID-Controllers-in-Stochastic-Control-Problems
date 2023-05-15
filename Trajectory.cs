using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Управление
{
    class Trajectory
    {
        public double[] data;
        public Trajectory(double[] a)
        {
            data = new double[a.Length];

            for (int j = 0; j < a.Length; j++)
            {
                data[j] = a[j];
            }
        }
    }
}
