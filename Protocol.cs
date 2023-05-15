using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Управление
{
    class Protocol
    {

        public string toprotocol = "";

        public void Prot2(List<FishSchoolSearch.Fish> pop, double I, int N)
        {
            string pr = "";
            pr += "\r\n" + "\r\n";
            pr += "КОНЕЧНАЯ ПОПУЛЯЦИЯ" + "\r\n" + "\r\n";

            pr += "|-----------------|-----------------|-----------------|----------------------------|\r\n";
            pr += "|   Номер агента  |        x1       |        x2       |      Приспособленность     |\r\n";
            pr += "|-----------------|-----------------|-----------------|----------------------------|\r\n";
            for (int i = 0; i < pop.Count; i++)
            {

                string g = Convert.ToString(i + 1);
                pr += "|  ";
                if (g.Length == 1) pr += "  ";
                else if (g.Length == 2) pr += " ";
                pr += Convert.ToString(i + 1) + "            |  ";
                for (int j = 0; j < N; j++)
                {
                    //string prob = "";
                    //string prob1 = "";
                    //string prob2 = "";
                    //if (pop[i].x[j] > 0) prob1 = " ";
                    //else prob2 = " ";
                    //if (j < 2)
                    //{
                    //    int length = 14 - Math.Round(pop[i].x[j], 6).ToString().Length;
                    //    for (int l = 0; l < length; l++) prob += " ";
                    //}

                    //pr += prob1 + Math.Round(pop[i].x[j], 6) + prob2 + prob + "|  ";
                    pr += Math.Round(pop[i].x[j], 4) + ", ";
                }

                //string prob0 = "";
                //string prob10 = "";
                //string prob20 = "";
                //if (pop[i].f > 0) prob10 = " ";
                //else prob20 = " ";

                //{
                //    int length2 = 25 - Math.Round(pop[i].f, 6).ToString().Length;
                //    for (int l = 0; l < length2; l++) prob0 += " ";
                //}
                //pr += prob10 + Math.Round(pop[i].f, 6) + prob20 + prob0 + "|  ";

                pr += Math.Round(pop[i].f, 6) + "|  ";

                pr += "\r\n";

            }
            pr += "|-----------------|-----------------|-----------------|----------------------------|\r\n";

            pr += "\r\n";
       
            //pr += "Количество локальных поисков: " + L.ToString() + "\r\n";
            pr += "Размер конечной популяции: " + pop.Count.ToString() + "\r\n";
            //pr += "Наилучшая приспособленность: " + Math.Round(pop[0].f,8).ToString() + "\r\n";
            pr += "Агент с наилучшей приспособленностью: x1 = " + Math.Round(pop[0].x[0], 8).ToString();
            pr += ", x2 = " + Math.Round(pop[0].x[1], 8).ToString() + "\r\n";
    
            toprotocol = pr;
        }
    }
}
