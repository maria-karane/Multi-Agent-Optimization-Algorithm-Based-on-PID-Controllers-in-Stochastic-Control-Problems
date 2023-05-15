using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.IO;
using System.Diagnostics;



namespace Управление
{
    public partial class Form1 : Form
    {
        
        ITasks task;

        public Form1()
        {
            
            InitializeComponent();
            dataGridView1.RowCount = 7;
            dataGridView2.RowCount = 3;
            pictureBox1.Refresh();

            if ((File.Exists("protocol.dt")))
            {
                File.Delete("protocol.dt");
            }
            FileStream vipoln_test = new FileStream("protocol.dt", FileMode.OpenOrCreate);
            vipoln_test.Close();
            FileStream fs = new FileStream("protocol.dt", FileMode.Append, FileAccess.Write);

            StreamWriter r = new StreamWriter(fs);
            r.Write("***************************************************************************************************************" + "\r\n"
                   + "**                                                                                                          **" + "\r\n"
                   + "**         Протокол работы алгоритма поиска оптимального управления на основе мультиагентных методов        **" + "\r\n"
                   + "**                                           " + DateTime.Now.ToString() + "                                **" + "\r\n"
                   + "**                                                                                                          **" + "\r\n"
                   + "**************************************************************************************************************" + "\r\n");
            r.Close();
            fs.Close();
        }

        private void Form1_Load(object sender, EventArgs e)
        {
        }

        private void comboBox1_SelectedIndexChanged_1(object sender, EventArgs e)
        {
            if (comboBox1.SelectedIndex == 0)
            {
               // task = new Task1();
                pictureBox1.Image = Properties.Resources.task1_bundle;
            }

            else if (comboBox1.SelectedIndex == 1)
            {
                pictureBox1.Image = Properties.Resources.task2_bundle;
            }

            else if (comboBox1.SelectedIndex == 2)
            {
                pictureBox1.Image = Properties.Resources.task3_bundle;
            }

            pictureBox1.Refresh();
        }

        private void radioButton1_CheckedChanged(object sender, EventArgs e)
        {
            dataGridView1.Rows[0].Cells[0].Value = "Размер начальной популяции";
            dataGridView1.Rows[0].Cells[1].Value = "5";
            dataGridView1.Rows[1].Cells[0].Value = "Количество итераций";
            dataGridView1.Rows[1].Cells[1].Value = "50";
            dataGridView1.Rows[2].Cells[0].Value = "Максимальный вес";
            dataGridView1.Rows[2].Cells[1].Value = "150";
            dataGridView1.Rows[3].Cells[0].Value = "Пороговый вес";
            dataGridView1.Rows[3].Cells[1].Value = "100";
            dataGridView1.Rows[4].Cells[0].Value = "Шаг stepInd0";
            dataGridView1.Rows[4].Cells[1].Value = "0,1";
            dataGridView1.Rows[5].Cells[0].Value = "Шаг stepInd1";
            dataGridView1.Rows[5].Cells[1].Value = "0,01";
            dataGridView1.Rows[6].Cells[0].Value = " ";
            dataGridView1.Rows[6].Cells[1].Value = " ";

            dataGridView1.Refresh();
        }

        private void radioButton2_CheckedChanged(object sender, EventArgs e)
        {
            dataGridView1.Rows[0].Cells[0].Value = "Размер начальной популяции";
            dataGridView1.Rows[0].Cells[1].Value = "5";
            dataGridView1.Rows[1].Cells[0].Value = "Количество итераций";
            dataGridView1.Rows[1].Cells[1].Value = "10";
            dataGridView1.Rows[2].Cells[0].Value = "Макс. скорость движения криля";
            dataGridView1.Rows[2].Cells[1].Value = "0,01";
            dataGridView1.Rows[3].Cells[0].Value = "Малое положительное число";
            dataGridView1.Rows[3].Cells[1].Value = "0,5";
            dataGridView1.Rows[4].Cells[0].Value = "Макс. скорость передвижения к источнику пищи";
            dataGridView1.Rows[4].Cells[1].Value = "0,02";
            dataGridView1.Rows[5].Cells[0].Value = "Максимальная скорость диффузии криля";
            dataGridView1.Rows[5].Cells[1].Value = "0,005";
            dataGridView1.Rows[6].Cells[0].Value = "Параметр Ci";
            dataGridView1.Rows[6].Cells[1].Value = "0,2";
          
            dataGridView1.Refresh();
        }

        private void radioButton3_CheckedChanged(object sender, EventArgs e)
        {
            dataGridView1.Rows[0].Cells[0].Value = "Число стран";
            dataGridView1.Rows[0].Cells[1].Value = "10";
            dataGridView1.Rows[1].Cells[0].Value = "Число империалистических стран";
            dataGridView1.Rows[1].Cells[1].Value = "3";
            dataGridView1.Rows[2].Cells[0].Value = "Количество итераций";
            dataGridView1.Rows[2].Cells[1].Value = "50";
            dataGridView1.Rows[3].Cells[0].Value = "Параметр сдвига колоний beta";
            dataGridView1.Rows[3].Cells[1].Value = "0,2";
            dataGridView1.Rows[4].Cells[0].Value = "Параметр сдвига колоний gamma";
            dataGridView1.Rows[4].Cells[1].Value = "0,02";
            dataGridView1.Rows[5].Cells[0].Value = "Параметр учета влияния колоний ksi";
            dataGridView1.Rows[5].Cells[1].Value = "0,01";
            dataGridView1.Rows[6].Cells[0].Value = " ";
            dataGridView1.Rows[6].Cells[1].Value = " ";

            dataGridView1.Refresh();
        }


        private void radioButton8_CheckedChanged(object sender, EventArgs e)
        {
            dataGridView1.Rows[0].Cells[0].Value = "Размер начальной популяции";
            dataGridView1.Rows[0].Cells[1].Value = "20";
            dataGridView1.Rows[1].Cells[0].Value = "Количество итераций";
            dataGridView1.Rows[1].Cells[1].Value = "5";
            dataGridView1.Rows[2].Cells[0].Value = "Кол-во точек (Безье)";
            dataGridView1.Rows[2].Cells[1].Value = "1";
            dataGridView1.Rows[3].Cells[0].Value = "Кол-во точек (B-сплайн)";
            dataGridView1.Rows[3].Cells[1].Value = "7";
            dataGridView1.Rows[4].Cells[0].Value = "Активность поиска по координате";
            dataGridView1.Rows[4].Cells[1].Value = "0,9";
            dataGridView1.Rows[5].Cells[0].Value = "Кол-во возм. положений членов поп.";
            dataGridView1.Rows[5].Cells[1].Value = "5";
            dataGridView1.Rows[6].Cells[0].Value = "Кол-во особей для сокр. популяции";
            dataGridView1.Rows[6].Cells[1].Value = "8";

            dataGridView1.Refresh();
        }

        private void radioButton9_CheckedChanged(object sender, EventArgs e)
        {
            dataGridView1.Rows[0].Cells[0].Value = "Размер начальной популяции";
            dataGridView1.Rows[0].Cells[1].Value = "9";
            dataGridView1.Rows[1].Cells[0].Value = "Число итераций за время прохода";
            dataGridView1.Rows[1].Cells[1].Value = "2";
            dataGridView1.Rows[2].Cells[0].Value = "Максимальное число проходов";
            dataGridView1.Rows[2].Cells[1].Value = "2";
            dataGridView1.Rows[3].Cells[0].Value = "Параметр Ks";
            dataGridView1.Rows[3].Cells[1].Value = "0,1";
            dataGridView1.Rows[4].Cells[0].Value = "Параметр Kl";
            dataGridView1.Rows[4].Cells[1].Value = "5";
            dataGridView1.Rows[5].Cells[0].Value = "Шаг интегрирования";
            dataGridView1.Rows[5].Cells[1].Value = "0,0001";
            dataGridView1.Rows[6].Cells[0].Value = " ";
            dataGridView1.Rows[6].Cells[1].Value = " ";


            dataGridView1.Refresh();
        }
        private void radioButton11_CheckedChanged(object sender, EventArgs e)
        {
            dataGridView1.Rows[0].Cells[0].Value = "Размер начальной популяции";
            dataGridView1.Rows[0].Cells[1].Value = "9";
            dataGridView1.Rows[1].Cells[0].Value = "Число итераций за время прохода";
            dataGridView1.Rows[1].Cells[1].Value = "2";
            dataGridView1.Rows[2].Cells[0].Value = "Максимальное число проходов";
            dataGridView1.Rows[2].Cells[1].Value = "2";
            dataGridView1.Rows[3].Cells[0].Value = "Параметр Kp";
            dataGridView1.Rows[3].Cells[1].Value = "0,1";
            dataGridView1.Rows[4].Cells[0].Value = "Параметр Kd1";
            dataGridView1.Rows[4].Cells[1].Value = "10";
            dataGridView1.Rows[5].Cells[0].Value = "Параметр Kd2";
            dataGridView1.Rows[5].Cells[1].Value = "0,01";
            dataGridView1.Rows[6].Cells[0].Value = "Параметр Ki";
            dataGridView1.Rows[6].Cells[1].Value = "10";


            dataGridView1.Refresh();
        }
        
        private void radioButton4_CheckedChanged(object sender, EventArgs e)
        {
            numericUpDown1.Visible = true;
            numericUpDown2.Visible = false;
            textBox2.Visible = false;
            textBox3.Visible = false;

            label3.Text = "Количество точек переключения: ";
            label3.Visible = true;
        }
        private void radioButton5_CheckedChanged(object sender, EventArgs e)
        {
            numericUpDown1.Visible = false;
            numericUpDown2.Visible = true;
            textBox2.Visible = true;
            textBox3.Visible = true;

            label3.Text = "Количество коэффициентов \n в разложении: ";
            label3.Visible = true;
        }
        private void radioButton10_CheckedChanged(object sender, EventArgs e)
        {
            numericUpDown1.Visible = false;
            numericUpDown2.Visible = true;
            textBox2.Visible = false;
            textBox3.Visible = false;

            label3.Text = "Количество коэффициентов \n в разложении: ";
            label3.Visible = true;
        }
        

        double[] traj = null;
        int _numb_tr = 1;
        double[] TA = null;

        int reduce = 0;

        private void button1_Click(object sender, EventArgs e)
        {
            Stopwatch stopwatch = new Stopwatch();
            stopwatch.Start(); //считаем затраченное время
            
            dataGridView1.Refresh();

            double[] rx = new double[2];
            double I_opt = 0;
            double[] U_opt = new double[2];
            double rx1_max = 0, rx1_min = 0, rx2_max = 0, rx2_min = 0;
            int error = 0;

            double u0 = Convert.ToDouble(textBox1.Text);
              
            //в релейном виде
            bool f_basis = radioButton4.Checked;

            //с разложением по косинусам
            bool view_cos = radioButton5.Checked;

            //с разложением по полиномам Чебышева
            bool chebyshev = radioButton10.Checked;
                        
            //критерий средний или гарантирующий
            bool I_count = radioButton6.Checked;
            bool I_count2 = radioButton7.Checked;

            //int fst = 2;

            //-----------------------------------------------------------
            int L0 = 0;
            int L1 = 0;
            int L2 = 0;

            try
            {
                L0 = (int)numericUpDown2.Value;
                L1 = Convert.ToInt32(textBox2.Text);
                L2 = Convert.ToInt32(textBox3.Text);
            }
            catch { error = 1; }
            if ((L1 < 0) || (L2 < 0)) error = 1;
            if (L0 == 0) error = 1;
            

            int[] x_L = { L1, L2 };
            //-----------------------------------------------------------
            int z = comboBox1.SelectedIndex;

            if (comboBox1.SelectedIndex == -1) MessageBox.Show("Выберите задачу!");
            else
            {
                if (radioButton1.Checked)//решение задачи 1 2 3 методом рыб
                {
                   
                    
                    //параметры
                    int pop = 0;
                    int iter = 0;
                    int scale = 0;
                    int thr = 0;
                    double ind0 = 0;
                    double ind1 = 0;

                    try
                    {
                        pop = Convert.ToInt32(dataGridView1.Rows[0].Cells[1].Value);
                        iter = Convert.ToInt32(dataGridView1.Rows[1].Cells[1].Value);
                        scale = Convert.ToInt32(dataGridView1.Rows[2].Cells[1].Value);
                        thr = Convert.ToInt32(dataGridView1.Rows[3].Cells[1].Value);
                        ind0 = Convert.ToDouble(dataGridView1.Rows[4].Cells[1].Value);
                        ind1 = Convert.ToDouble(dataGridView1.Rows[5].Cells[1].Value);
                    }
                    catch { error = 2; }
                    if ((pop <= 0) || (iter <= 0)) error = 2;
                    //-----------------------------------------------------------

                    if (error == 0)
                    {
                        //int z = comboBox1.SelectedIndex;
                        ITasks task = null;
                        if (z == 0)
                        {
                            task = new Task4(chebyshev);
                            task.setUprType(u0, f_basis, view_cos, chebyshev, L0, x_L);
                            task.setIcount(I_count, I_count2);

                        }
                        else if (z == 1)
                        {
                            task = new Task5(chebyshev);
                            task.setUprType(u0, f_basis, view_cos, chebyshev, L0, x_L);
                            task.setIcount(I_count, I_count2);
                        }
                        else if (z == 2)
                        {
                            task = new Task6(chebyshev);
                            task.setUprType(u0, f_basis, view_cos, chebyshev, L0, x_L);
                            task.setIcount(I_count, I_count2);
                        }

                        //запуск алгоритма
                        var alg = new FishSchoolSearch();

                        int N = (int)numericUpDown1.Value;
                        if (view_cos)
                        {
                            if (L1 == 0 && L2 != 0) N = L0 * L2;
                            else if (L2 == 0 && L1 != 0) N = L0 * L1;
                            else if (L1 == 0 && L2 == 0) N = L0;
                            else N = L0 * L1 * L2;

                        }
                        else if (chebyshev) N = L0;
                        else if (f_basis) { }
                        else
                        {
                            MessageBox.Show("Выберите систему базисных функций!");
                            return;
                        }

                        alg.Init(pop, iter, N, task, scale, thr, ind0, ind1);
                        alg.Create_Pop();


                        //for (int i = 0; i < 1; i++)
                        //{
                            alg.Work();
                            I_opt = alg.I_opt;
                        //    if (I_opt < -0.128)
                        //    {
                        //        Console.Write(i);
                        //        break;
                        //    }

                        //}

                        //I_opt = alg.I_opt;
                        rx = alg.x_opt;
                        //------
                        rx1_max = rx[0];
                        rx1_min = rx[0];
                        rx2_max = rx[1];
                        rx2_min = rx[1];
                        for (int i = 1; i < rx.Length; i++)
                        {
                            if (i % 2 == 0)
                            {
                                if (rx[i] > rx1_max) rx1_max = rx[i];
                                if (rx[i] < rx1_min) rx1_min = rx[i];
                            }
                            else
                            {
                                if (rx[i] > rx2_max) rx2_max = rx[i];
                                if (rx[i] < rx2_min) rx2_min = rx[i];
                            }
                        }
                        //-----------
                        traj = alg.traj;
                        TA = alg.ta;

                        _numb_tr = alg._numb_tr;
                        reduce = TA.Length;
                    }
                    else switch (error)
                        {
                            case 1:
                                MessageBox.Show("Неверно введено количество коэффициентов!");
                                break;
                            case 2:
                                MessageBox.Show("Неверно введены параметры алгоритма!");
                                break;
                        }

                }

                else if (radioButton2.Checked) //решение задачи 1 2 3 методом криля
                {
                    //int L0 = (int)numericUpDown2.Value;
                    //int L1 = Convert.ToInt32(textBox2.Text);
                    //int L2 = Convert.ToInt32(textBox3.Text);
                    //int[] x_L = { L1, L2 };
                    //int z = comboBox1.SelectedIndex;
                    
                    //-----------------------------------------------------------
                     //параметры
                    int pop = 0;
                    int iter = 0;
                    double Nmax = 0;
                    double mu = 0;
                    double Vf = 0;
                    double Dmax = 0;
                    double ci = 0;
                    try
                    {
                        pop = Convert.ToInt32(dataGridView1.Rows[0].Cells[1].Value);
                        iter = Convert.ToInt32(dataGridView1.Rows[1].Cells[1].Value);
                        Nmax = Convert.ToDouble(dataGridView1.Rows[2].Cells[1].Value);
                        mu = Convert.ToDouble(dataGridView1.Rows[3].Cells[1].Value);
                        Vf = Convert.ToDouble(dataGridView1.Rows[4].Cells[1].Value);
                        Dmax = Convert.ToDouble(dataGridView1.Rows[5].Cells[1].Value);
                        ci = Convert.ToDouble(dataGridView1.Rows[6].Cells[1].Value);
                    }
                    catch { error = 2; }
                    if ((pop <= 0) || (iter <= 0)) error = 2;
                    //-----------------------------------------------------------

                    if (error == 0)
                    {
                        ITasks task = null;

                        if (z == 0)
                        {
                            task = new Task4(chebyshev);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                            task.setUprType(u0, f_basis, view_cos, chebyshev, L0, x_L);
                            task.setIcount(I_count, I_count2);
                        }
                        else if (z == 1)
                        {
                            task = new Task5(chebyshev);
                            task.setUprType(u0, f_basis, view_cos, chebyshev, L0, x_L);
                            task.setIcount(I_count, I_count2);
                        }
                        else if (z == 2)
                        {
                            task = new Task6(chebyshev);
                            task.setUprType(u0, f_basis, view_cos, chebyshev, L0, x_L);
                            task.setIcount(I_count, I_count2);
                        }
                        

                        //запуск алгоритма
                        var alg = new Krill_Herd();

                        int N = (int)numericUpDown1.Value;

                        if (view_cos)
                        {
                            if (L1 == 0 && L2 != 0) N = L0 * L2;
                            else if (L2 == 0 && L1 != 0) N = L0 * L1;
                            else if (L1 == 0 && L2 == 0) N = L0;
                            else N = L0 * L1 * L2;

                        }
                        else if (chebyshev) N = L0;
                        else if (f_basis) { }
                        else
                        {
                            MessageBox.Show("Выберите систему базисных функций!");
                            return;
                        }
                     
                        alg.Init(pop, iter, N, task, mu, Nmax, Vf, Dmax, ci);
                        alg.Create_Pop();
                        alg.Work();
                        I_opt = alg.I_opt;

                        rx = alg.x_opt;
                        //------
                        rx1_max = rx[0];
                        rx1_min = rx[0];
                        rx2_max = rx[1];
                        rx2_min = rx[1];
                        for (int i = 1; i < rx.Length; i++)
                        {
                            if (i % 2 == 0)
                            {
                                if (rx[i] > rx1_max) rx1_max = rx[i];
                                if (rx[i] < rx1_min) rx1_min = rx[i];
                            }
                            else
                            {
                                if (rx[i] > rx2_max) rx2_max = rx[i];
                                if (rx[i] < rx2_min) rx2_min = rx[i];
                            }
                        }
                        //-----------
                        traj = alg.traj;

                        TA = alg.ta;
                        _numb_tr = alg._numb_tr;
                        reduce = TA.Length;
                    }
                    else switch (error)
                        {
                            case 1:
                                MessageBox.Show("Неверно введено количество коэффициентов!");
                                break;
                            case 2:
                                MessageBox.Show("Неверно введены параметры алгоритма!");
                                break;
                        }
                }

                else if (radioButton3.Checked) //решение задачи 1 2 3 методом империалистов
                {
                    //int L0 = (int)numericUpDown2.Value;
                    //int L1 = Convert.ToInt32(textBox2.Text);
                    //int L2 = Convert.ToInt32(textBox3.Text);
                    //int[] x_L = { L1, L2 };
                    //int z = comboBox1.SelectedIndex;

                    //-----------------------------------------------------------
                    //параметры
                    int pop = 0;
                    int imper = 0;
                    int iter = 0;
                    double betta = 0;
                    double gamma = 0;
                    double ksi = 0;
                    try 
                    {
                        pop = Convert.ToInt32(dataGridView1.Rows[0].Cells[1].Value);
                        imper = Convert.ToInt32(dataGridView1.Rows[1].Cells[1].Value);
                        iter = Convert.ToInt32(dataGridView1.Rows[2].Cells[1].Value);
                        betta = Convert.ToDouble(dataGridView1.Rows[3].Cells[1].Value);
                        gamma = Convert.ToDouble(dataGridView1.Rows[4].Cells[1].Value);
                        ksi = Convert.ToDouble(dataGridView1.Rows[5].Cells[1].Value);
                    }
                    catch { error = 2; }
                    if ((pop <= 0) || (iter <= 0)) error = 2;
                    //-----------------------------------------------------------

                    if (error == 0)
                    {
                        ITasks task = null;

                        if (z == 0)
                        {
                            task = new Task4(chebyshev);
                            task.setUprType(u0, f_basis, view_cos, chebyshev, L0, x_L);
                            task.setIcount(I_count, I_count2);
                        }
                        else if (z == 1)
                        {
                            task = new Task5(chebyshev);
                            task.setUprType(u0, f_basis, view_cos, chebyshev, L0, x_L);
                            task.setIcount(I_count, I_count2);
                        }
                        else if (z == 2)
                        {
                            task = new Task6(chebyshev);
                            task.setUprType(u0, f_basis, view_cos, chebyshev, L0, x_L);
                            task.setIcount(I_count, I_count2);
                        }
                        else { MessageBox.Show("Выберите задачу!"); }

                        //запуск алгоритма
                        var alg = new Imperialist();

                        int N = (int)numericUpDown1.Value;
                        if (view_cos)
                        {
                            if (L1 == 0 && L2 != 0) N = L0 * L2;
                            else if (L2 == 0 && L1 != 0) N = L0 * L1;
                            else if (L1 == 0 && L2 == 0) N = L0;
                            else N = L0 * L1 * L2;

                        }
                        else if (chebyshev) N = L0;
                        else if (f_basis) { }
                        else
                        {
                            MessageBox.Show("Выберите систему базисных функций!");
                            return;
                        }
                        alg.Init(pop, iter, N, task, imper, betta, gamma, ksi);
                        alg.Create_Pop();
                        alg.Work();
                        I_opt = alg.I_opt;
                        rx = alg.x_opt;
                        //------
                        rx1_max = rx[0];
                        rx1_min = rx[0];
                        rx2_max = rx[1];
                        rx2_min = rx[1];
                        for (int i = 1; i < rx.Length; i++)
                        {
                            if (i % 2 == 0)
                            {
                                if (rx[i] > rx1_max) rx1_max = rx[i];
                                if (rx[i] < rx1_min) rx1_min = rx[i];
                            }
                            else
                            {
                                if (rx[i] > rx2_max) rx2_max = rx[i];
                                if (rx[i] < rx2_min) rx2_min = rx[i];
                            }
                        }
                        //-----------
                        traj = alg.traj;

                        TA = alg.ta;
                        _numb_tr = alg._numb_tr;
                        reduce = TA.Length;
                    }
                    else switch (error)
                        {
                            case 1:
                                MessageBox.Show("Неверно введено количество коэффициентов!");
                                break;
                            case 2:
                                MessageBox.Show("Неверно введены параметры алгоритма!");
                                break;
                        }
                }
                else if (radioButton8.Checked) //решение задачи 1 2 3 интерполяционным методом
                {
                    //int L0 = (int)numericUpDown2.Value;
                    //int L1 = Convert.ToInt32(textBox2.Text);
                    //int L2 = Convert.ToInt32(textBox3.Text);
                    //int[] x_L = { L1, L2 };

                    //if (L1 == 0 || L2 == 0)
                    //{
                    //    L1 = L1 + 1;
                    //    L2 = L2 + 1;

                    //}

                    //int z = comboBox1.SelectedIndex;
                    //параметры
                    int pop = 0;
                    int iter = 0;
                    int M1 = 0;
                    int M2 = 0;
                    double PRT = 0;
                    int nstep = 0;
                    int b2 = 0;
                    try 
                    {
                        //параметры
                        pop = Convert.ToInt32(dataGridView1.Rows[0].Cells[1].Value);
                        iter = Convert.ToInt32(dataGridView1.Rows[1].Cells[1].Value);
                        M1 = Convert.ToInt32(dataGridView1.Rows[2].Cells[1].Value);
                        M2 = Convert.ToInt32(dataGridView1.Rows[3].Cells[1].Value);
                        PRT = Convert.ToDouble(dataGridView1.Rows[4].Cells[1].Value);
                        nstep = Convert.ToInt32(dataGridView1.Rows[5].Cells[1].Value);
                        b2 = Convert.ToInt32(dataGridView1.Rows[6].Cells[1].Value);
                    }
                    catch { error = 2; }
                    if ((pop <= 0) || (iter <= 0)) error = 2;
                    //-----------------------------------------------------------

                    if (error == 0)
                    {
                        ITasks task = null;

                        if (z == 0)
                        {
                            task = new Task4(chebyshev);
                            task.setUprType(u0, f_basis, view_cos, chebyshev, L0, x_L);
                            task.setIcount(I_count, I_count2);
                        }
                        else if (z == 1)
                        {
                            task = new Task5(chebyshev);
                            task.setUprType(u0, f_basis, view_cos, chebyshev, L0, x_L);
                            task.setIcount(I_count, I_count2);
                        }
                        else if (z == 2)
                        {
                            task = new Task6(chebyshev);
                            task.setUprType(u0, f_basis, view_cos, chebyshev, L0, x_L);
                            task.setIcount(I_count, I_count2);
                        }
                        
                       
                        //запуск алгоритма
                        var alg = new Interpolation_Search();

                        int N = (int)numericUpDown1.Value;
                        if (view_cos)
                        {
                            if (L1 == 0 && L2 != 0) N = L0 * L2;
                            else if (L2 == 0 && L1 != 0) N = L0 * L1;
                            else if (L1 == 0 && L2 == 0) N = L0;
                            else N = L0 * L1 * L2;
                        }
                        else if (chebyshev) N = L0;
                        else if (f_basis) { }
                        else
                        {
                            MessageBox.Show("Выберите систему базисных функций!");
                            return;
                        }
                        alg.Init(pop, iter, N, task, M1, M2, PRT, nstep, b2);

                        //Вывод
                        //-------------------------------------------------------------------
                        //for (int i = 0; i < 1; i++)
                        //{
                        //    alg.Create_Pop();
                        //    alg.Work();
                        //    I_opt = alg.I_opt;
                        //    Console.WriteLine(" " + Math.Round(I_opt, 4));
                        //    if (I_opt < -0.133)
                        //    {
                        //        Console.WriteLine(" " + i);
                        //        break;
                        //    }
                        //}

                        //Stream statisticsDataFile = null;
                        //var path1 = "statistics-" + z + ".txt";
                        //if (File.Exists(path1)) File.Delete(path1);
                        //statisticsDataFile = File.OpenWrite(path1);
                        //-------------------------------------------------------------------


                        alg.Create_Pop();
                        alg.Work();
                        I_opt = alg.I_opt;

                        rx = alg.x_opt;
                        //------
                        rx1_max = rx[0];
                        rx1_min = rx[0];
                        rx2_max = rx[1];
                        rx2_min = rx[1];
                        for (int i = 1; i < rx.Length; i++)
                        {
                            if (i % 2 == 0)
                            {
                                if (rx[i] > rx1_max) rx1_max = rx[i];
                                if (rx[i] < rx1_min) rx1_min = rx[i];
                            }
                            else
                            {
                                if (rx[i] > rx2_max) rx2_max = rx[i];
                                if (rx[i] < rx2_min) rx2_min = rx[i];
                            }
                        }
                        //-----------
                        traj = alg.traj;

                        TA = alg.ta;
                        _numb_tr = alg._numb_tr;
                        reduce = TA.Length;
                        //Вывод продолжение
                        //-------------------------------------------------------------------
                        //string line = "";
                        //for (int i = 0; i < 41; i++)
                        //{
                        //    double t_coord = task.t0() + i * task.Dt();
                        //    line += Math.Round(t_coord, 2).ToString();

                        //    for (int j = 0; j < traj.Length / 41; j++)
                        //    {
                        //        line += "\t" + Math.Round(traj[i + j * 41], 4).ToString();
                        //    }
                        //    line += "\r\n";
                        //}

                        //var lineBuf = Encoding.ASCII.GetBytes(line);
                        //statisticsDataFile.Write(lineBuf, 0, lineBuf.Length);

                        //if (statisticsDataFile != null)
                        //{
                        //    statisticsDataFile.Close();
                        //}
                        //-------------------------------------------------------------------
                    }
                    else switch (error)
                        {
                            case 1:
                                MessageBox.Show("Неверно введено количество коэффициентов!");
                                break;
                            case 2:
                                MessageBox.Show("Неверно введены параметры алгоритма!");
                                break;
                        }
                }

                else if (radioButton9.Checked) //решение задачи 1 2 3  методом с линейными регуляторами
                {
                    //int L0 = (int)numericUpDown2.Value;
                    //int L1 = Convert.ToInt32(textBox2.Text);
                    //int L2 = Convert.ToInt32(textBox3.Text);
                    //int[] x_L = { L1, L2 };
                    //int z = comboBox1.SelectedIndex;
                    
                    //-----------------------------------------------------------
                    //параметры
                    int pop = 0;
                    int K = 0;
                    int Pmax = 0;
                    double ks = 0;
                    double kl = 0;
                    double h = 0;
                    try 
                    {                   
                        pop = Convert.ToInt32(dataGridView1.Rows[0].Cells[1].Value);
                        K = Convert.ToInt32(dataGridView1.Rows[1].Cells[1].Value);
                        Pmax = Convert.ToInt32(dataGridView1.Rows[2].Cells[1].Value);
                        ks = Convert.ToDouble(dataGridView1.Rows[3].Cells[1].Value);
                        kl = Convert.ToDouble(dataGridView1.Rows[4].Cells[1].Value);
                        h = Convert.ToDouble(dataGridView1.Rows[5].Cells[1].Value);
                    }
                    catch { error = 2; }
                    if ((pop <= 0) || ((pop - 1) % 4 != 0) || (K <= 0)) error = 2;
                    //-----------------------------------------------------------

                    

                    if (error == 0)
                    {
                        ITasks task = null;

                        if (z == 0)
                        {
                            task = new Task4(chebyshev);
                            task.setUprType(u0, f_basis, view_cos, chebyshev, L0, x_L);
                            task.setIcount(I_count, I_count2);
                        }
                        else if (z == 1)
                        {
                            task = new Task5(chebyshev);
                            task.setUprType(u0, f_basis, view_cos, chebyshev, L0, x_L);
                            task.setIcount(I_count, I_count2);
                        }
                        else if (z == 2)
                        {
                            task = new Task6(chebyshev);
                            task.setUprType(u0, f_basis, view_cos, chebyshev, L0, x_L);
                            task.setIcount(I_count, I_count2);
                        }
                        
                        //запуск алгоритма   
                        var alg = new Linear_regulator();
                        int N = (int)numericUpDown1.Value;
                        if (view_cos)
                        {
                            if (L1 == 0 && L2 != 0) N = L0 * L2;
                            else if (L2 == 0 && L1 != 0) N = L0 * L1;
                            else if (L1 == 0 && L2 == 0) N = L0;
                            else N = L0 * L1 * L2;
                        }
                        else if (chebyshev) N = L0;
                        else if (f_basis) { }
                        else
                        {
                            MessageBox.Show("Выберите систему базисных функций!");
                            return;
                        }

                        alg.Init_NEW(pop, K, 2 * N, task, Pmax, ks, kl, h);
                        alg.Create_Pop();
                        alg.Work();
                        I_opt = alg.I_opt;

                        rx = alg.x_opt;

                        //------
                        rx1_max = rx[0];
                        rx1_min = rx[0];
                        rx2_max = rx[1];
                        rx2_min = rx[1];
                        for (int i = 1; i < rx.Length; i++)
                        {
                            if (i % 2 == 0)
                            {
                                if (rx[i] > rx1_max) rx1_max = rx[i];
                                if (rx[i] < rx1_min) rx1_min = rx[i];
                            }
                            else
                            {
                                if (rx[i] > rx2_max) rx2_max = rx[i];
                                if (rx[i] < rx2_min) rx2_min = rx[i];
                            }
                        }
                        //-----------

                        traj = alg.traj;

                        TA = alg.ta;
                        _numb_tr = alg._numb_tr;
                        reduce = TA.Length / 2;
                    }
                    else switch (error)
                        {
                            case 1:
                                MessageBox.Show("Неверно введено количество коэффициентов!");
                                break;
                            case 2:
                                MessageBox.Show("Неверно введены параметры алгоритма!");
                                break;
                        }
                }
                else if (radioButton11.Checked) //решение задачи 1 2 3  методом с ПИД-регуляторами
                {
                    //-----------------------------------------------------------
                    //параметры
                    int pop = 0;
                    int K = 0;
                    int Pmax = 0;
                    double kp = 0;
                    double kd1 = 0;
                    double kd2 = 0;
                    double ki = 0;
                    double h = 0.0001;
                    try
                    {
                        pop = Convert.ToInt32(dataGridView1.Rows[0].Cells[1].Value);
                        K = Convert.ToInt32(dataGridView1.Rows[1].Cells[1].Value);
                        Pmax = Convert.ToInt32(dataGridView1.Rows[2].Cells[1].Value);
                        kp = Convert.ToDouble(dataGridView1.Rows[3].Cells[1].Value);
                        kd1 = Convert.ToDouble(dataGridView1.Rows[4].Cells[1].Value);
                        kd2 = Convert.ToDouble(dataGridView1.Rows[5].Cells[1].Value);
                        ki = Convert.ToDouble(dataGridView1.Rows[6].Cells[1].Value);                        
                    }
                    catch { error = 2; }
                    if ((pop <= 0) || ((pop - 1) % 4 != 0) || (K <= 0)) error = 2;
                    //-----------------------------------------------------------
                    if (error == 0)
                    {
                        ITasks task = null;

                        if (z == 0)
                        {
                            task = new Task4(chebyshev);
                            task.setUprType(u0, f_basis, view_cos, chebyshev, L0, x_L);
                            task.setIcount(I_count, I_count2);
                        }
                        else if (z == 1)
                        {
                            task = new Task5(chebyshev);
                            task.setUprType(u0, f_basis, view_cos, chebyshev, L0, x_L);
                            task.setIcount(I_count, I_count2);
                        }
                        else if (z == 2)
                        {
                            task = new Task6(chebyshev);
                            task.setUprType(u0, f_basis, view_cos, chebyshev, L0, x_L);
                            task.setIcount(I_count, I_count2);
                        }

                        //запуск алгоритма   
                        var alg = new PID_regulators();
                        int N = (int)numericUpDown1.Value;
                        if (view_cos)
                        {
                            if (L1 == 0 && L2 != 0) N = L0 * L2;
                            else if (L2 == 0 && L1 != 0) N = L0 * L1;
                            else if (L1 == 0 && L2 == 0) N = L0;
                            else N = L0 * L1 * L2;
                        }
                        else if (chebyshev) N = L0;
                        else if (f_basis) { }
                        else
                        {
                            MessageBox.Show("Выберите систему базисных функций!");
                            return;
                        }

                        alg.Init_NEW(pop, K, 2 * N, task, Pmax, kp, kd1, kd2, ki, h);

                        //Вывод
                        //-------------------------------------------------------------------
                        for (int i = 0; i < 1; i++)
                        {
                            alg.Create_Pop();
                            alg.Work();
                            I_opt = alg.I_opt;
                            Console.WriteLine(" " + Math.Round(I_opt, 4));
                            if (I_opt < -2.55)
                            {
                                Console.WriteLine(" " + i);
                                break;
                            }
                        }

                        Stream statisticsDataFile = null;
                        var path1 = "statistics-" + z + ".txt";
                        if (File.Exists(path1)) File.Delete(path1);
                        statisticsDataFile = File.OpenWrite(path1);
                        //-------------------------------------------------------------------
                        
                        //alg.Create_Pop();
                        //alg.Work();
                        //I_opt = alg.I_opt;

                        rx = alg.x_opt;

                        //------
                        rx1_max = rx[0];
                        rx1_min = rx[0];
                        rx2_max = rx[1];
                        rx2_min = rx[1];
                        for (int i = 1; i < rx.Length; i++)
                        {
                            if (i % 2 == 0)
                            {
                                if (rx[i] > rx1_max) rx1_max = rx[i];
                                if (rx[i] < rx1_min) rx1_min = rx[i];
                            }
                            else
                            {
                                if (rx[i] > rx2_max) rx2_max = rx[i];
                                if (rx[i] < rx2_min) rx2_min = rx[i];
                            }
                        }
                        //-----------

                        traj = alg.traj;

                        TA = alg.ta;
                        _numb_tr = alg._numb_tr;
                        reduce = TA.Length;
                        //reduce = TA.Length / 2;

                        //Вывод продолжение
                        //-------------------------------------------------------------------
                        string line = "";
                        int points = (int)((task.t1() - task.t0()) / task.Dt()) + 1;
                        for (int i = 0; i < points; i++)
                        {
                            double t_coord = task.t0() + i * task.Dt();
                            line += Math.Round(t_coord, 2).ToString();

                            for (int j = 0; j < traj.Length / points; j++)
                            {
                                line += "\t" + Math.Round(traj[i + j * points], 4).ToString();
                            }
                            line += "\r\n";
                        }

                        var lineBuf = Encoding.ASCII.GetBytes(line);
                        statisticsDataFile.Write(lineBuf, 0, lineBuf.Length);

                        if (statisticsDataFile != null)
                        {
                            statisticsDataFile.Close();
                        }
                        //-------------------------------------------------------------------

                    }
                    else switch (error)
                        {
                            case 1:
                                MessageBox.Show("Неверно введено количество коэффициентов!");
                                break;
                            case 2:
                                MessageBox.Show("Неверно введены параметры алгоритма!");
                                break;
                        }
                }
                else
                {
                    MessageBox.Show("Выберите метод оптимизации!");
                    error = 1;
                }

                if (error == 0)
                {

                    dataGridView2.Rows[0].Cells[0].Value = "x1; x2";
                    //dataGridView2.Rows[0].Cells[1].Value = rx[0] + "; " + rx[1];

                    dataGridView2.Rows[0].Cells[1].Value = "[" + Math.Round(rx1_min, 4) + ";" + Math.Round(rx1_max, 4) + "]"
                                                 + ", " + "[" + Math.Round(rx2_min, 4) + ";" + Math.Round(rx2_max, 4) + "]";

                    dataGridView2.Rows[1].Cells[1].Value = I_opt;

                    if (f_basis)
                    {
                        dataGridView2.Rows[2].Cells[0].Value = "Точки переключения";
                    }
                    else
                    {
                        dataGridView2.Rows[2].Cells[0].Value = "Коэффициенты в разложении";
                    }


                    if (I_count)
                    {
                        dataGridView2.Rows[1].Cells[0].Value = "I(среднее)";
                    }

                    else if (I_count2)
                    {
                        dataGridView2.Rows[1].Cells[0].Value = "I(гарантирующее)";
                    }

                    dataGridView2.Rows[2].Cells[1].Value = array_str(TA, reduce);

                    stopwatch.Stop();
                    TimeSpan ts = stopwatch.Elapsed;
                    string elapsedTime = String.Format("{0:00}:{1:00}:{2:00}.{3:00}", ts.Hours, ts.Minutes, ts.Seconds, ts.Milliseconds / 10);
                    Console.WriteLine("RunTime " + elapsedTime);
                    label7.Text = "Время выполнения программы: " + elapsedTime;
                    label7.Visible = true;                          
                }
            }
            
        }

        static string array_str(double[] a, int reduce)
        {
            
            var s = "";
            for (int j = 0; j < reduce; j++)
            {
                if (j > 0) s += "; ";
                s += Math.Round(a[j],2);

            }
            return s;
        }

        Graph gr;
        private void button2_Click(object sender, EventArgs e)
        {
           // if(gr==null) 
            if (traj == null) MessageBox.Show("Решите задачу. Нажмите на кнопку \"Решить\". ");
            else
            {
                gr = new Graph();
                gr.d = _numb_tr;
                gr.traj = traj;
                gr.Show();
            }
        }

        private void label4_Click(object sender, EventArgs e)
        {

        }

        private void dataGridView1_CellContentClick(object sender, DataGridViewCellEventArgs e)
        {

        }

        private void button4_Click(object sender, EventArgs e)
        {         

        }

        private void button6_Click(object sender, EventArgs e) //выход
        {
            Close();
        }

        private void button3_Click(object sender, EventArgs e) //протокол
        {
            Process.Start("protocol.dt");
        }

        private void button5_Click(object sender, EventArgs e) //Справка
        {
            Process.Start("control.pdf");
        }

       
    }
}
