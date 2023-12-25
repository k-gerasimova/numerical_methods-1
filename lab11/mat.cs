using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace lab11
{
    internal class MAT
    {
        //класс с методами, нужными для лаб
        public void output(double[,] matrix, int n)
        {
            //вывод матрицы
            for(int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    Console.Write("{0:f4}", matrix[i, j]);
                    if (j != n - 1) Console.Write(" ");
                }
                Console.WriteLine();
            }
        }

        public void input(double[, ] A, int n)
        {
            //ввод матрицы
            string[] help = new string[n];
            for (int i = 0; i < n; i++)
            {
                help = Console.ReadLine().Split();
                for (int j = 0; j < n; j++)
                {
                    A[i, j] = Convert.ToInt32(help[j]);
                }

            }
        }

        public double[,] SwapRows(double[,] matrix, int k, int i, int n)
        {   // Свап строк
            for (int j = 0; j < n; j++)
            {
                double temp;
                temp = matrix[k, j];
                matrix[k, j] = matrix[i, j];
                matrix[i, j] = temp;
            }
            return matrix;
        }

        public double[,] E(int n)
        {
            //генерация единичной матрицы
            double[,] res = new double[n, n];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    if (i == j) res[i, j] = 1;
                    else res[i, j] = 0;
                }
            }
            return res;
        }

        public double[,] Comp(double[,] A, double[, ]B, int n)
        {
            //перемножение двух матриц n x n
            if (A.GetLength(1) != B.GetLength(0))
                throw new Exception("err");

            double[,] res = new double[n, n];
            double help = 0;
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    for (int k = 0; k < n; k++)
                    {
                        help += A[i, k] * B[k, j];
                    }
                    res[i, j] = help;
                    help = 0;
                }
            }
            return res;

        }

        public double[] Comp(double[,] A, double[]B, int n)
        {
            //перемножение матрицы n x n и столбца B n x 1
            double[] res = new double[n];
            for (int i = 0; i<n; i++)
            {
                for(int j = 0; j < n; j++)
                {
                    res[i] += A[i, j] * B[j];
                }
            }
            return res;
        }

        public double Comp(double[] A, double[] B, int n)
        {
            //перемножение двух строк/столбцов
            double res = 0;
                for (int j = 0; j < n; j++)
                {
                    res += A[j] * B[j];
                }
            return res;
        }


        public double[,] A_1(double[,] A, int n)
        {
            double[,] res = new double[n, n];

            double[,] dopoln = new double[n, n];
            double det = Det(A, n);
            //Console.WriteLine("det:{0}", det);
            double minor;
            if (det == 0)
            {
                throw new Exception("No A^-1");
            }
            for(int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    minor = Det(GetMatr(A, n, j, i), n - 1);
                  //  Console.WriteLine("Minor {0}{1}:{2}", i, j, minor);
                    dopoln[i, j] = Math.Pow(-1, i + j + 2) * minor;

                }
            }
            double [,] pre_res = Transpose(dopoln, n);
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    res[i, j] = 1 / det * pre_res[i, j];

                }
            }

            return res;
        }

        public double[,] Transpose(double[,] A, int n)
        {
            //транспонирование матрицы
            double[,] res = new double[A.GetLength(0), A.GetLength(1)];
            for (int i = 0; i < A.GetLength(0); i++)
            {
                for (int j = 0; j < A.GetLength(1); j++)
                {
                    res[i, j] = A[j, i];
                }
            }
            return res;
        }

        public double Det(double[, ] A, int n)
        {
            //поиск определителя
            double det = 0;
            int k = 1;
            if (n < 1)
            {
                Console.WriteLine("Cannot find det.");
                return 0;
            }
            else if (n == 1) return A[0, 0];
            else if (n == 2) return A[0, 0] * A[1, 1] - (A[1, 0] * A[0, 1]);
            else if (n > 2)
            {
                double[,] smal;
                for (int i = 0; i < n; i++)
                {
                    smal = GetMatr(A, n, i, 0);
                    det = det + k * A[0, i] * Det(smal, n - 1);
                    k = -k;
                }
            }
            return det;

        }

        private double[,] GetMatr(double[, ]A, int n,  int j, int i)
        {
            //вспомогательный метод для генерации матрицы без строки + столбца
            double[,] res = new double[A.GetLength(0) -1, A.GetLength(1) - 1];
            int ki, kj, di, dj;
            di = 0;
            for (ki = 0; ki < n - 1; ki++) { // проверка индекса строки
                if (ki == i) di = 1;
                dj = 0;
                for (kj = 0; kj < n - 1; kj++) { // проверка индекса столбца
                    if (kj == j) dj = 1;
                    res[ki, kj] = A[ki + di, kj + dj];
                }
            }
            return res;
        }

        public double[,] Reversed_m(double[,] m)
        {
            //обратная матрица
            int n = m.GetLength(0);

            double[,] rev = new double[n, n];
            for (int i = 0; i < n; i++)
                rev[i, i] = 1;

            double[,] big_m = new double[n, 2 * n];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                {
                    big_m[i, j] = m[i, j];
                    big_m[i, j + n] = rev[i, j];
                }

            for (int k = 0; k < n; k++) // строки
            {
 
                if (big_m[k, k] == 0)
                {
                    for (int l = k + 1; l < n; l++)
                    {
                        if (big_m[l, l] == 0)
                            continue;
                        else
                        {
                            for (int j = 0; j < n * 2; j++)
                            {
                                double temp;
                                temp = big_m[l, j];
                                big_m[l, j] = big_m[k, j];
                                big_m[k, j] = temp;
                            }
                            for (int j = 0; j < n; j++)
                            {
                                double temp;
                                temp = m[l, j];
                                m[l, j] = m[k, j];
                                m[k, j] = temp;
                            }
                            break;
                        }
                    }
                }
                for (int i = 0; i < 2 * n; i++) 
                {
                    big_m[k, i] = big_m[k, i] / m[k, k];
                }
                for (int i = k + 1; i < n; i++) 
                {
                    double K = big_m[i, k] / big_m[k, k]; //Коэффициент
                    for (int j = 0; j < 2 * n; j++) //j-номер столбца следующей строки после k
                        big_m[i, j] = big_m[i, j] - big_m[k, j] * K; //Зануление элементов матрицы ниже первого члена, преобразованного в единицу
                }
                for (int i = 0; i < n; i++) //Обновление, внесение изменений в начальную матрицу
                    for (int j = 0; j < n; j++)
                        m[i, j] = big_m[i, j];
            }

            //Обратный ход (Зануление верхнего правого угла)
            for (int k = n - 1; k > -1; k--) //k-номер строки
            {
                for (int i = 2 * n - 1; i > -1; i--) //i-номер столбца
                    big_m[k, i] = big_m[k, i] / m[k, k];
                for (int i = k - 1; i > -1; i--) //i-номер следующей строки после k
                {
                    double K = big_m[i, k] / big_m[k, k];
                    for (int j = 2 * n - 1; j > -1; j--) //j-номер столбца следующей строки после k
                        big_m[i, j] = big_m[i, j] - big_m[k, j] * K;
                }
            }

            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    rev[i, j] = big_m[i, j + n];

            return rev;
        }


        public double Normed(double[,]A, int n)
        {
            //поиск нормы для матрицы
            double[] s = new double[n];
            double norm;
            for (int i = 0; i < n; i++)
            {
                for (int j = 0, k =0; k <n; k++, j++)
                {
                    s[i] += Math.Abs(A[i, j]);
                }

            }
            norm = s.Max();
            return norm;
        }

        public double Normed(double[] X, int size)
        {
            //поиск нормы для столбца/строки
            double norm = 0;
            for (int i = 0; i < size; i++)
            {
                X[i] = Math.Abs(X[i]);
            }
            norm = X.Max();
            return norm;
        }

        public double[] Sub(double[] a, double[] b)
        {
            //поиск разности между строками/столбцами
            double[] r = new double[b.GetLength(0)];
            for (int i = 0; i < a.GetLength(0); i++)
            {
                r[i] = a[i] - b[i];
            }
            return r;
        }
        public double[,] Sub(double[,] a, double[,] b) //поиск разности между матрицами
        {
            double[,] r = new double[a.GetLength(0), a.GetLength(1)];
            for (int i = 0; i < a.GetLength(1); i++)
            {
                for (int j = 0; j < a.GetLength(0); j++)
                {
                    r[i, j] = a[i, j] - b[i, j];
                }
            }
            return r;
        }

        public double[] Add(double[] a, double[] b) // скалярное сложение векторов
        {
            double[] r = new double[a.GetLength(0)];
            for (int i = 0; i < a.GetLength(0); i++)
            {
                r[i] = a[i] + b[i];
            }
            return r;
        }

        public double[,] Add_sbstr(double[]a, double[]b, int n) //перемножение столбца и строки
        {
            double[,]r = new double[n, n];
            for(int i = 0; i < n; i++)
            {
                for (int j = 0; j< n; j++)
                {
                    r[i, j] = a[i] * b[j];
                }
            }
            return r;
        }
        public double[,] del(double[,]a, double z, int n) //деление матрицы на число
        {
            for (int i = 0;i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    a[i, j] /= z;
                }
            }
            return a;
        }

        public double[,] Comp(double[, ] A, double b, int n) //умножение матрицы на число
        {
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    A[i, j] *= b;
                }
            }
            return A;
        }
        public double Add_strsb(double[] a, double[]b, int n) // перемножение строки на столбец
        {
            double r = 0;
            for (int i = 0; i< n;i++)
            {
                r += a[i] * b[i];
            }
            return r;
        }


    }




}
