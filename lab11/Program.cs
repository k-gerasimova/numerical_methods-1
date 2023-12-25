using lab14;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;



namespace lab11
{

    internal class Program
    {
        static void Main(string[] args)
        {
            Console.WriteLine("Choose Method:");
            Console.WriteLine("1. LU");
            Console.WriteLine("2. Progonka");
            Console.WriteLine("3. Zeidel and simple iteration");
            Console.WriteLine("4. Rotation method");
            Console.WriteLine("5. QR");
            int method;
            method = Convert.ToInt32(Console.ReadLine());
            MAT matrix_method = new MAT();

            Console.WriteLine("Введите размерность матрицы А:");
            int n = Convert.ToInt32(Console.ReadLine());
            string[] help = new string[n];
            Console.WriteLine("Введите матрицу А:");
            double [,] A = new double[n, n];
            matrix_method.input(A, n);
            double[] B = new double[n];

            if (method != 4 && method != 5)
            {
                Console.WriteLine("Введите вектор В:");
                help = Console.ReadLine().Split(' ');
                for (int i = 0; i < n; i++) B[i] = Convert.ToDouble(help[i]);
            }

            switch (method)
            {
                case 1:
                    LU method_1 = new LU();
                    method_1.Lu(A, B, n);
                    break;
                case 2:
                    Progonka method_2 = new Progonka();
                    method_2.Prog(A, B, n);
                    break;
                case 3:
                    double eps;
                    Console.WriteLine("Введите погрешность:");
                    eps = Convert.ToDouble(Console.ReadLine());
                    Zeidel z = new Zeidel();
                    z.prost(A, B, n, eps);
                    z.z_m(A, B, n, eps);
                    break;
                case 4:
                    double eps2;
                    Console.WriteLine("Введите погрешность:");
                    eps2 = Convert.ToDouble(Console.ReadLine());
                    rot rotation = new rot();
                    rotation.U = matrix_method.E(n);
                    rotation.rotate(A,n, eps2);
                    break;
                case 5:
                    double eps3;
                    Console.WriteLine("Введите погрешность:");
                    eps3 = Convert.ToDouble(Console.ReadLine());
                    QR qr = new QR();
                    qr.QR_func(A, n);
                    qr.search_sz(A, n, eps3);
                    break;

            }
            

        }
    }
}