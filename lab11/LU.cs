using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace lab11
{
    internal class LU
    {
        public double[,] A;
        public double[] B;
        public double[,] A_ish;
        int p = 0;
        public void Lu(double[,]A,  double[] B, int n)
        {
            MAT matrix_method = new MAT();

            //copy
            A_ish = new double[n, n];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    A_ish[i, j] = A[i, j];
                }
            }

            Console.WriteLine("\nmy method\n");

            double[,] E = matrix_method.E(n);
            for (int i = 0; i < n; i++)
            {
                double main_element = 0;
                int max_i = -1;

                //поиск главного элемента
                for (int k = i; k < n; k++)
                {
                    if (Math.Abs(A[k, i]) > main_element)
                    {
                        main_element = Math.Abs(A[k, i]);
                        max_i = k;
                    }
                }

                if (main_element == 0)
                    throw new Exception("Матрица вырождена");

                if (max_i != i) p++;
                E = matrix_method.SwapRows(E, max_i, i, n);
                A = matrix_method.SwapRows(A, max_i, i, n);

                for (int j = i + 1; j < n; j++)
                {
                    //приравниваем поддиагональные элементы 
                    //доп множителю l
                    A[j, i] /= A[i, i];
                    for (int k = i + 1; k < n; k++)
                    {
                        A[j, k] -= A[j, i] * A[i, k];
                    }

                }

            }
            double[,] L = new double[n, n];
            double[,] U = new double[n, n];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    if (i == j)
                    {
                        L[i, j] = 1;
                        U[i, j] = A[i, j];
                    }
                    else if (i > j)
                    {
                        L[i, j] = A[i, j];
                        U[i, j] = 0;
                    }
                    else
                    {
                        L[i, j] = 0;
                        U[i, j] = A[i, j];
                    }

                }
            }

            double[] res = new double[n];
            double[] y = new double[n];
            double[] B_upd = matrix_method.Comp(E, B, n);
            for (int i = 0; i < n; i++)
            {
                if (i == 0) y[0] = B_upd[0];
                else
                {
                    double help1 = 0;
                    for (int j = 0; j < i; j++)
                    {
                        help1 += y[j] * L[i, j];
                    }
                    y[i] = B_upd[i] - help1;
                }



            }
            for (int i = n - 1; i >= 0; i--)
            {

                if (i == n - 1) res[n - 1] = y[n - 1] / U[n - 1, n - 1];
                else
                {
                    double help1 = 0;
                    for (int j = i + 1; j < n; j++)
                    {
                        help1 += U[i, j] * res[j];
                    }
                    res[i] = 1 / U[i, i] * (y[i] - help1);
                }
            }
            Console.Write("lu:\n");
            double[,] lu = new double[n, n];
            lu = matrix_method.Comp(matrix_method.Transpose(E, n), matrix_method.Comp(L, U, n), n);
            matrix_method.output(A, n);
            Console.Write("l:\n");
            matrix_method.output(L, n);
            Console.Write("u:\n");
            matrix_method.output(U, n);

            // matrix_method.output(U, n);
            for (int i = 0; i < n; i++)
            {
                Console.Write("\nx{0}:{1}", i, res[i]);

            }
            double det = 1;
            for (int i = 0; i < n; i++)
            {
                det *= A[i, i];
            }
            
            Console.WriteLine("\n\nDeterminant:{0}", det);
            double det_1 = matrix_method.Det(A_ish, n) * Math.Pow(-1, p);
            Console.WriteLine("\nDeterminant from prev A:{0}", det_1);

            Console.WriteLine("\nReverted Matrix:");
            double[,] reverse = matrix_method.A_1(A_ish, n);
            matrix_method.output(reverse, n);

            Console.WriteLine("\nReverted Matrix LU:");
            double[,] reverse_1 = matrix_method.Reversed_m(A_ish);
            matrix_method.output(reverse_1, n);

            Console.ReadKey();
        }
    }

    internal class Progonka
    {
        public double[,] A;
        public double[] B;

        private int IsCorrect(double[,] matrix, int n)
        {
            for (int i = 1; i < n - 1; i++) {
                if (Math.Abs(matrix[i, i]) < Math.Abs(matrix[i, i - 1]) + Math.Abs(matrix[i, i + 1])) {
                    return 1;
                }
            }


            if (Math.Abs(matrix[0, 0]) < Math.Abs(matrix[0, 1]) || (Math.Abs(matrix[n - 1, n - 1]) < Math.Abs(matrix[n - 1, n - 2]))) {
                return 1;
            }

            for (int i = 0; i < n; i++)
            {
                if (matrix[i, i] == 0)
                {
                    return 2;
                }
            } 


            return 0;
        }
        public void Prog(double[, ]A, double[]B, int n)
        {
            MAT matrix_method = new MAT();


            double[] v = new double[n];
            double[] u = new double[n];



            double[] x = new double[n];
            if (IsCorrect(A, n) != 0)
            {
                Console.WriteLine("not valid matrix");
                Console.ReadKey();
                return;
            }
            v[0] = A[0, 1] / (-A[0, 0]);
            u[0] = (-B[0]) / (-A[0, 0]);

            for (int i = 1; i < n-1; i++)
            {
                v[i] = -A[i, i + 1] / (A[i, i] + A[i, i - 1] * v[i - 1]);
                u[i] = (-A[i, i - 1] * u[i - 1] + B[i]) / (A[i, i] + A[i, i - 1] * v[i - 1]);
            }

            v[n - 1] = 0;
            u[n - 1] = (A[n - 1, n - 2] * u[n - 2] - B[n - 1]) / (-A[n - 1, n - 1] - A[n - 1, n - 2] * v[n - 2]);
            x[n - 1] = u[n - 1];
            for (int i = n-2; i >= 0; i--)
            {
                x[i] = v[i] * x[i + 1] + u[i];
            }
            for (int i = 0; i < n; i++)
            {
                Console.Write("\nx{0}:{1}", i, x[i]);

            }
            Console.WriteLine("\nПрогоночные коэфициенты:\nu:");
            for (int i = 0; i < n; i++)
            {
                Console.Write("\nu{0}:{1}", i, u[i]);

            }
            Console.WriteLine("\nv:");
            for (int i = 0; i < n; i++)
            {
                Console.Write("\nv{0}:{1}", i, v[i]);

            }
            Console.ReadKey();


        }

    }

    internal class Zeidel
    {
        //Метод Зейделя
        public void z_m(double[,] A, double[]B, int n, double eps)
        {
            double[,] C = new double[n, n];
            double[,] D = new double[n, n];
            double[,] E= new double[n, n];
            double[,] alpha = new double[n, n];

            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    if (i == j) E[i, j] = 1;
                    else E[i, j] = 0;
                    if (j > i)
                    {
                        D[i, j] = -A[i, j]/A[i, i];
                        alpha[i, j] = D[i, j];
                        C[i, j] = 0;
                    } else
                    {
                        C[i, j] = -A[i, j]/A[i, i];
                        alpha[i, j] = C[i, j];
                        D[i, j] = 0;
                    } if (i == j)
                    {
                        C[i, j] = 0;
                        D[i, j] = 0;
                        alpha[i, j] = 0;
                    }
                }
            }
            MAT m = new MAT();
            

            double eps_1 = 1;
            double[] X = new double[n];
            double[] X_p = new double[n];

            int k = 0;

            for (int i = 0; i < n; i++)
            {
                B[i] = B[i] / A[i, i];
                X_p[i] = B[i];
            }

            while (eps_1 > eps)
            {
                X = m.Add(m.Comp(m.Comp(m.A_1(m.Sub(E, C), n), D, n), X_p, n), m.Comp(m.A_1(m.Sub(E, C), n), B, n));
                if (m.Normed(alpha, n) >= 1)
                    eps_1 = m.Normed(m.Sub(X, X_p), n);
                else
                    eps_1 = m.Normed(D, n) * (m.Normed(m.Sub(X, X_p), n)) / (1 - m.Normed(alpha, n));
                for (int i = 0; i < n; i++)
                {
                    X_p[i] = X[i];
                    X[i] = 0;
                }
                k++;
            }



            Console.WriteLine("\nZeidel:");
            for (int i = 0; i < n; i++)
            {
                Console.Write("\nx{0}:{1:f4}", i, X_p[i]);

            }
            Console.WriteLine("\niterations:{0}", k);
            Console.ReadKey();
        }


        //метод простых итераций
        public void prost(double[, ]A, double[] B, int n, double eps)
        {
            MAT m = new MAT();
            double[,] alpha = new double[n, n];
            for(int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    if (i == j) alpha[i, j] = 0;
                    else alpha[i, j] = -A[i, j] / A[i, i];
                    
                }
            }

            double eps_1 = 1;
            double[] X = new double[n];
            double[] X_p = new double[n];
            int k = 0;
            for (int i = 0; i < n; i++)
            {
                X_p[i] = B[i]/A[i, i];
            }
            while(eps_1 > eps && m.Normed(m.Sub(X, X_p), n) > eps)
            {
                for (int i = 0; i <n; i++)
                {
                    X[i] = B[i] / A[i,i] + m.Comp(alpha, X_p ,n)[i];
                }
                if (m.Normed(alpha, n) >= 1)
                    eps_1 = m.Normed(m.Sub(X, X_p), n);
                else
                    eps_1 = m.Normed(alpha, n) * (m.Normed(m.Sub(X, X_p), n)) / (1 - m.Normed(alpha, n));
                for (int i = 0; i < n; i++)
                {
                    X_p[i] = X[i];
                    X[i] = 0;
                }
                k++;
            }



            Console.WriteLine("\nEasy Iteration:");
            for (int i = 0; i < n; i++)
            {
                Console.Write("\nx{0}:{1:f4}", i, X_p[i]);

            }
            Console.WriteLine("\niterations:{0}", k);
            // Console.ReadKey();

        }
    }
}