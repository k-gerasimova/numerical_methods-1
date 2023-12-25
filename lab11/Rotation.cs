using lab11;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace lab14
{
    internal class rot
    {

        public double[,] A_k;
        public double[,] U;
        
        public void rotate(double[,] A, int n, double eps)
        {
            A_k = A;
            MAT m = new MAT();
            double max = double.MinValue;
            int max_i = -1;
            int max_j = -1;

            double[,] U_k = new double[n, n];

            //поиск наиб внедиаг элемента
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    if (i != j && Math.Abs(A_k[i, j]) > max)
                    {
                        max = Math.Abs(A_k[i, j]);
                        max_i = i;
                        max_j = j;
                    }
                } 
            }

            double phi = 0.5 * Math.Atan(2 * A_k[max_i, max_j] / (A_k[max_i, max_i] - A_k[max_j, max_j]));

            //создание матрицы U

            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    if ( i == max_i && j == max_j) 
                    {
                        U_k[i, j] = - Math.Sin(phi);
                    } else if(i == max_j && j == max_i)
                    {
                        U_k[i, j] = Math.Sin(phi);
                    } else if (i == max_i && j == max_i)
                    {
                        U_k[i, j] = Math.Cos(phi);
                    } else if (i == max_j && j == max_j)
                    {
                        U_k[i, j] = Math.Cos(phi);
                    } else if (i == j)
                    {
                        U_k[i, j] = 1;
                    } else
                    {
                        U_k[i, j] = 0;
                    }
                }
            }

            A_k = new double[n, n];
            A_k = m.Comp(m.Comp(m.Transpose(U_k, n), A, n), U_k, n);
            U = m.Comp(U, U_k, n);
            double krit = 0;
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    if ( i != j)
                    {
                        krit += Math.Pow(A_k[i, j], 2);
                    }
                }
            }
            krit = Math.Pow(krit, 0.5);
            if (krit > eps) rotate(A_k, n, eps);
            else
            {

                Console.WriteLine("\nСобственные значения:");
                for(int i = 0; i < n; i++)
                {
                    Console.WriteLine("lambda{0}:{1}", i, A_k[i, i]);
                }

                Console.WriteLine();
                Console.WriteLine("Собственные векторы:");
                for (int i = 0; i < n; i++)
                {
                    Console.Write("x{0}:", i);
                    for(int j = 0; j < n; j++)
                    {
                        Console.Write("{0} ", U[i, j]);
                    }
                    Console.WriteLine();

                }
                Console.WriteLine();
                Console.WriteLine("проверка ортогональности:");
                double[] x = new double[n];
                double[] x_1 = new double[n];
                for (int i = 0; i < n; i++)
                {
                    for (int j = 0; j < n; j++)
                    {
                        if (i != j && i < j)
                        {
                            for (int k = 0; k < n; k++)
                            {
                                x[k] = U[i, k];
                                x_1[k] = U[j, k];
                            }
                            Console.WriteLine("(x{0}, x{1}):{2}", i, j, m.Comp(x, x_1, n));
                        }
                    }
                }
                Console.ReadKey();
            }
            
        }


    }

    internal class QR
    {
        public double[,] Q;
        public double[,] R;
        public double[,] A_k;
        MAT m = new MAT();
        

       public void search_sz(double[,] A, int n, double eps)
        {
            
            double d;
            double[] sz = new double[n];
            do {
                QR_raz(A, n);
                A_k = m.Comp(A_k, Q, n);
                Q = m.E(n);
                d = 0;
                for (int i = 0; i < n; i++)
                {
                    if (d < Math.Abs(sz[i] - A_k[i, i]))
                    {
                        d = Math.Abs(sz[i] - A_k[i, i]);
                    }
                    sz[i] = A_k[i, i];
                }
                for(int i = 0; i < n; i++)
                {
                    for (int j = 0; j < n; j++)
                    {
                        A[i, j] = A_k[i, j];
                    }
                }

            } while (d > eps);


            Console.WriteLine("\nСобственные значения:");
            for (int i = 0; i < n; i++)
            {
                Console.WriteLine("labmda{0}:{1}", i, sz[i]);
            }
            Console.ReadKey();
        }
        public void QR_func(double[,] A, int n)
        {
            QR_raz(A, n);
            Console.WriteLine("Q");


            m.output(Q, n);
            Console.WriteLine("R");
            m.output(A_k, n);
            
            Q = m.E(n);
            R = new double[n, n];
            A_k = new double[n, n];

        }
        public void QR_raz(double[,] A, int n)
        {
            double[,] A_kk = A;
            double[] v = new double[n];
            double[,] H;
            Q = m.E(n);
            for(int i = 0; i < n - 1; i++)
            {
                double[] b = make_b(A_kk, n, i);
                for(int j = 0; j < n; j++)
                {
                    if (j < i) v[j] = 0;
                    else if (j == i) v[i] = A_kk[i, j] + Math.Sign(A_kk[i, j]) * norm_b(b);
                    else
                    {
                        v[j] = A_kk[j, i];
                    }
                }
                H = m.Sub(m.E(n), m.Comp(m.del(m.Add_sbstr(v, v, n), m.Add_strsb(v, v, n), n), 2, n));
                A_kk = m.Comp(H, A_kk, n);
                Q = m.Comp(Q, H, n);

            }
            A_k = A_kk;
        }

        public double[] make_b(double[,] a, int n, int i)
        {
            double[] b = new double[n - i];
            int l = 0;
            for (int k = 0; k < n; k++)
            {
                if (k >= i) b[l++] = a[k, i];
            }
            return b;
        }
        public double norm_b(double[] b)
        {
            double ret = 0;
            for(int i = 0; i < b.GetLength(0); i++)
            {
                ret += Math.Pow(b[i], 2);
            }
            ret = Math.Pow(ret, 0.5);
            return ret;
        }
    }
}
