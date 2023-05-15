using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace Управление
{
    public partial class Graph : Form
    {
        public Graph()
        {
            InitializeComponent();
        }

        public double[] traj;
        public int d;


        private void pictureBox1_Paint(object sender, PaintEventArgs e)
        {
            var g = e.Graphics;
            var b0 = Brushes.AliceBlue;
            var w = pictureBox1.Width;
            var h = pictureBox1.Height;
            g.FillRectangle(b0, new Rectangle(0,0,w , h));
            var pen = new Pen(Color.Blue, 2);
            var pen1 = new Pen(Color.Red, 2);
            var pen2 = new Pen(Color.Green, 2);
            var pen3 = new Pen(Color.Black, 2);
            Font font1 = new Font("TimesNewRoman", 12, FontStyle.Bold);
            Font font2 = new Font("TimesNewRoman", 8);
            float a = 30;

            var x0 = 20;
            var y0 = 20;

            if (traj == null) 
            {
                //Ох
                g.DrawLine(pen, new Point(0, h - 3*y0), new Point(w - x0, h - 3*y0));

                //Оу
                g.DrawLine(pen, new Point(x0, h), new Point(x0, y0));

                //стрелка у
                g.DrawLine(pen, new Point((int)(1.4*x0), 2*y0), new Point(x0, y0));
                g.DrawLine(pen, new Point((int)(0.6 * x0), 2 * y0), new Point(x0, y0));

                //стрелка х
                g.DrawLine(pen, new Point(w - x0, h - 3 * y0), new Point(w - (int)(1.8 * x0), h - (int)(3.4 * y0)));
                g.DrawLine(pen, new Point(w - x0, h - 3 * y0), new Point(w - (int)(1.8 * x0), h - (int)(2.6 * y0)));


                g.DrawString("t", font1, Brushes.Black, w - 50, h - 55);
                e.Graphics.DrawString("x1, x2, u", font1, Brushes.Black, 2, 40);

                return;
            }

            int m = traj.Length/d/3 - 1;
            int dx = w / m;

            var m0 = traj[0];
            var m1 = traj[0];

            for (int s = 0; s < d; s++)
            {


               

                for (int j = 1; j <= m; j++)
                {
                    if (traj[j + s * 3 * (m + 1)] < m0) m0 = traj[j + s * 3 * (m + 1)];

                    if (traj[m + 1 + j + s * 3 * (m + 1)] < m0) m0 = traj[m + 1 + j + s * 3 * (m + 1)];
                    if (traj[2 * (m + 1) + j + s * 3 * (m + 1)] < m0) m0 = traj[2 * (m + 1) + j + s * 3 * (m + 1)];

                    if (traj[j + s * 3 * (m + 1)] > m1) m1 = traj[j + s * 3 * (m + 1)];

                    if (traj[m + 1 + j + s * 3 * (m + 1)] > m1) m1 = traj[m + 1 + j + s * 3 * (m + 1)];
                    if (traj[2 * (m + 1) + j + s * 3 * (m + 1)] > m1) m1 = traj[2 * (m + 1) + j + s * 3 * (m + 1)];

                }
            }


                if (m0 > 0) { m0 = -1.2; }
                if (m1 < 0) { m1 = 1.2; }
                if (m1 < m0 + 0.1) { m1 = m0 + 0.1; }



                m0 = m0 * 1.25;
                m1 = m1 * 1.25;

                var q = (h - 2 * y0) / (m1 - m0);


                var dy = (int)(q * (-m0));


                //Ох
                g.DrawLine(pen3, new Point(0, h - y0 - dy), new Point(w - x0, h - y0 - dy));

                //Оу
                g.DrawLine(pen3, new Point(x0, h - dy - (int)(q * m0)), new Point(x0, h - y0 - dy - (int)(q * m1)));

                //стрелка х
                g.DrawLine(pen3, new Point(w - x0, h - y0 - dy), new Point(w - (int)(1.8 * x0), h - y0 - dy - (int)(0.2 * y0)));
                g.DrawLine(pen3, new Point(w - x0, h - y0 - dy), new Point(w - (int)(1.8 * x0), h - y0 - dy - (int)(-0.2 * y0)));

                //стрелка у
                g.DrawLine(pen3, new Point((int)(1.2 * x0), h - dy - y0 + (int)(0.8 * y0) - (int)(q * m1)), new Point(x0, h - dy - y0 - (int)(q * m1)));
                g.DrawLine(pen3, new Point((int)(0.8 * x0), h - dy - y0 + (int)(0.8 * y0) - (int)(q * m1)), new Point(x0, h - dy - y0 - (int)(q * m1)));

                g.DrawString("t", font1, Brushes.Black, w - 50, h - 20 - dy);
                g.DrawString("x1, x2, u", font1, Brushes.Black, 20, 40);

                //int count = 0;
                for (int j = x0; j <= (w - x0 - 15); j++)
                {
                    if (j % 50 == 0)
                    {
                        g.DrawLine(pen3, new Point(j, h - y0 - dy + 3), new Point(j, h - y0 - dy - 3));
                        //g.DrawString("0", font1, Brushes.Black, j, h - y0 - dy + 5);
                        //count++;
                    }


                }
            for (int s = 0; s < d; s++)
            {
                for (int j = 1; j <= m; j++)
                {
                    if (s == 0 && j == 1)
                    {
                        g.DrawString("x1", font1, Brushes.Blue, x0 + m * dx - dx, h - y0 - dy - (int)(q * traj[m + s * 3 * (m + 1)]));
                        //g.DrawString("-1", font1, Brushes.Blue,x0 + j * dx - dx, h - y0 - dy - (int)(q * traj[j - 1 + s * 3 * (m + 1)]));
                    }
                    g.DrawLine(pen, new Point(x0 + j * dx - dx, h - y0 - dy - (int)(q * traj[j - 1 + s * 3 * (m + 1)])),
                          new Point(x0 + j * dx, h - y0 - dy - (int)(q * traj[j + s * 3 * (m + 1)])));
                    if (s == 0 && j == 1)
                    {
                        g.DrawString("x2", font1, Brushes.Red, x0 + m * dx - dx, h - y0 - dy - (int)(q * traj[m + 1 + m + s * 3 * (m + 1)]));
                    }
                    g.DrawLine(pen1, new Point(x0 + j * dx - dx, h - y0 - dy - (int)(q * traj[m + 1 + j - 1 + s * 3 * (m + 1)])),
                           new Point(x0 + j * dx, h - y0 - dy - (int)(q * traj[m + 1 + j + s * 3 * (m + 1)])));
                    if (s == 0 && j == 1)
                    {
                        g.DrawString("u", font1, Brushes.Green, x0 + m * dx - dx, h - 20 - dy - (int)(q * traj[2 * (m + 1) + m + s * 3 * (m + 1)]));
                    }
                    g.DrawLine(pen2, new Point(x0 + j * dx - dx, h - y0 - dy - (int)(q * traj[2 * (m + 1) + j - 1 + s * 3 * (m + 1)])),
                             new Point(x0 + j * dx, h - y0 - dy - (int)(q * traj[2 * (m + 1) + j + s * 3 * (m + 1)])));
                }
            }

        }
    }
}
