// МВСем7Зад1.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//


#define _USE_MATH_DEFINES
#include <iostream>
#include <math.h>
#include <windows.h>
#include <iomanip>
#include <cstdlib>
using namespace std;
//строит нулевую матрицу порядка n
double** get_null_mat(int n)
{
    double** A = new double* [2 * n];
    for (int i = 0; i <= n; i++)
    {
        A[i] = new double[2 * n];
        for (int j = 0; j <= n; j++)
            A[i][j] = 0.0;
    }
    return A;
}
//строит нулевой вектор порядка n
double* get_null_vec(int n)
{
    double* A = new double[2 * n];
    for (int i = 0; i <= n; i++)
    {
            A[i] = 0.0;
    }
    return A;
}
//Копирует матрицу b в матрицу a
void cop(double** a, double** b, int n)
{
    for (int i = 0; i <= n; i++)
        for (int j = 0; j <= n; j++)
            a[i][j] = b[i][j];
}
//Выводит матрицу
void show_mat(double** A, int n)
{
    for (int i = 0; i <= n; i++)
    {
        for (int j = 0; j <= n; j++)
        {
            cout << fixed << setprecision(8)<< " \t" << A[i][j] ;
        }
        cout << endl;
    }
}
//Выводит вектор
void show_vec(double* A, int n)
{
    for (int i = 0; i <= n; i++)
    {
            cout << fixed << setprecision(4) << " vec["<<i<<"] = " << A[i] << "   ";
        cout << endl;
    }
}
//Считает норму матрицы 
double norm_mat(double** a, int n)
{
    double sum1 = 0.0;
    double sum2 = 0.0;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
            sum1 += fabs(a[i][j]);
        //cout<<sum1<<endl;
        if (sum1 > sum2)
        {
            sum2 = sum1;
        }
        sum1 = 0.0;
    }
    return sum2;
}
//Считает евклидову норму вектора
double EVnorm_vec(double* vec, int n)
{
    double norm = 0.0;
    for (int i = 0; i < n; i++)
    {
        norm += pow(vec[i], 2);
    }
    norm = sqrt(norm);
    return norm;
}
//Считает норму вектора
double norm_vec(double* vec, int n)
{
    double norm = vec[0];
    for (int i = 1; i < n; i++)
    {
        if (fabs(vec[i]) >= fabs(norm))
        {
            norm = fabs(vec[i - 1]);
        }
    }
    return norm;
}
//Переводит матрицу в вектор */
double* Trans_mat_to_vec(double** A, int n)
{
    double* vec = get_null_vec((n + 1) * (n + 1));
    int k = 0;
    for (int i = 1; i < n; i++)
    {
        for (int j = 1; j < n; j++)
        {
            vec[k] = A[i][j];
            // cout << A[i][j] << "--" << vec[k] << endl;
            k++;
        }
    }

    show_vec(vec, (n + 1) * (n + 1) - 1);
    // cout << "всего " << k << endl;
    return vec;
}
double real_solution(double x, double y)
{
    return((x * x * y * y) / 2);
}
double f(double x, double y)
{
    return((-1)*(pow(x, 2) + pow(y, 2)));
}
double mu_x(double x)
{
    return(pow(x, 2)/2);
}
double mu_y(double y)
{
    return(pow(y, 2) / 2);
}
//Приближение дифференциального оператора L разностным
double** L(double** U_0, double h_x, int n)
{
    double** L = get_null_mat(n + 1);
    cop(L, U_0, n);
    for (int i = 1; i < n; i++)
        for (int j = 1; j < n; j++)
            L[i][j] = (U_0[i + 1][j] + U_0[i - 1][j] + U_0[i][j + 1] + U_0[i][j - 1] - 4 * U_0[i][j]) / pow(h_x, 2);
    return L;
}

//Вычисляет приближение к спектральному радиусу
double Spectral_Radious(double h_x)
{
    double sigma = 8 * pow(sin(M_PI * h_x / 2), 2) / pow(h_x, 2);
    double delta = 8 * pow(cos(M_PI * h_x / 2), 2) / pow(h_x, 2);
    double spectral_radious = (delta - sigma) / (sigma + delta);
    return spectral_radious;
}
//метод простой итерации
double** Easy_iterations(double** U,double** f, double h_x, int n, double EPS)
{
    double m = 2.0 * log10(1.0 / EPS) / pow((M_PI * h_x), 2);
    cout <<"число итераций в методе Прост. итер. = " << (int)m << endl;
    double** U_next = get_null_mat(n);
    cop(U_next, U, n);
    for (int k = 1; k <= m; k++)
    {
        for (int i = 1; i < n; i++)
            for (int j = 1; j < n; j++)
            {
                U_next[i][j] = (U[i - 1][j] + U[i + 1][j] + U[i][j - 1] + U[i][j + 1] + (pow(h_x, 2) * f[i][j])) * 0.25;
            }
        cop(U, U_next, n);
    }
    return U;
}
//Метод Зейделя
double** Seidel(double** U, double** f, double h_x, int n, double EPS)
{
    double m = log10(1.0 / EPS) / pow((M_PI * h_x), 2);
    cout << "число итераций в методе Зейделя = " <<(int) m << endl;
    double** U_next = get_null_mat(n);
    cop(U_next, U, n);
    for (int k = 1; k <= m; k++)
    {
        for (int i = 1; i < n; i++)
            for (int j = 1; j < n; j++)
            {
                U_next[i][j] = (U_next[i - 1][j] + U[i + 1][j] + U_next[i][j - 1] + U[i][j + 1] + (pow(h_x, 2) * f[i][j])) * 0.25;
            }
        cop(U, U_next, n);
    }
    return U;
}
//Метод Верхней релаксации
double** Top_relaxation(double** U, double** f, double h_x, int n, double EPS)
{
    double m = 2.0 * log10(1.0 / EPS) / (M_PI * h_x);
    cout << "число итераций в методе Верхней Релаксации = " << (int)m << endl;
    double** U_next = get_null_mat(n);
    double** compare = get_null_mat(n);
    double w_opt = 2.0 / (1 + sin(M_PI * h_x));
    cop(U_next, U, n);
    for (int k = 1; k <= m; k++)
    {
        for (int i = 1; i < n; i++)
            for (int j = 1; j < n; j++)
            {
                U_next[i][j] = U[i][j] + w_opt * (pow(h_x, 2) * f[i][j] + U_next[i - 1][j] + U[i + 1][j] + U_next[i][j - 1] + U[i][j + 1] - 4 * U[i][j]) * 0.25;
            }
        cop(U, U_next, n);
    }
    return U;
}
//Строит матрицу А системы -AU=F
double** Build_mat_A(int n)
{
    int N = (n - 1) * (n - 1);
    double** A = get_null_mat((n-1)*(n-1));
    A[0][0] = 4;
    A[0][1] = -1;
    A[0][1 + (n - 1)] = -1;
    for (int i = 1; i <= N; i ++)
    {
        A[i][i] = 4;
        A[i][i - 1] = -1;
        A[i][i + 1] = -1;
        A[i][i + (n)] = -1;
        A[i][i - (n)] = -1;
    }
    show_mat(A, N);
    return A;
}
//Метод переменных треугольников
double** Triangle(double** U_k, double** f, double h_x, int n, double EPS)
{
    double m = 0.29 * log10(2.0 / EPS) / (sqrt(h_x));
    cout << "число итераций в Попеременно треугольном методе = " << (int)m << endl;
    double help = sin(M_PI * h_x * 0.5);
    double eta = pow(help, 2);
    double delta = 8 / pow(h_x, 2);
    double sigma = delta * eta;
    double omega = pow(h_x, 2) / (4 * help);
    double gamma_1 = (4 * eta) / (pow(h_x, 2) * (1 + help));
    double gamma_2 = (2 * help) / pow(h_x, 2);
    double k = omega / pow(h_x, 2);
    double F_i_j = 0.0;
    double** Lu_k = get_null_mat(n + 2);
    double** U_k_1 = get_null_mat(n+2);
    double** w_bar = get_null_mat(n + 2);
    double** w = get_null_mat(n + 2);
    double tau_k = 0.0;
    cop(U_k_1, U_k, n);
    for (int k = 0; k <= m; k++)
    {
        tau_k = 2 / (gamma_1 + gamma_2 + (gamma_2 - gamma_1) * cos((2 * k - 1) * M_PI * 0.5 * n));
        Lu_k = L(U_k_1, h_x, n);
        for(int i = 1; i <n ;i++)
            for (int j = 1; j < n; j++)
            {
                F_i_j = Lu_k[i][j] + f[i][j];
                w_bar[i][j] = (k * w_bar[i - 1][j] + k * w_bar[i][j-1] + F_i_j) / (1 + 2 * k);
            }
        for (int i = n-1; i > 0; i--)
            for (int j = n-1; j > 0; j--)
            {
                w[i][j] = (k * w_bar[i + 1][j] + k * w_bar[i][j+1] + w_bar[i][j]) / (1 + 2 * k);
            }
        for (int i = 1; i < n; i++)
            for (int j = 1; j < n; j++)
            {
                U_k[i][j] = U_k_1[i][j] + tau_k * w[i][j];
            }
        cop(U_k_1, U_k, n);
    }
    return U_k;
}

int main()
{
    system("chcp 1251");
    system("cls");
    double N, M;
    double l_x = 1;
    double l_y = 1;
    cout << "Кисиев Тимур, группа 20-Б06-мм, Вариант 9 " << endl;
    cout << "Введите число разбиений отрезка [0,1]: ";
    cin >> N;
    M = N;
    double h_x = 1 / N;
    double h_y = 1 / M;
    int n = (int)N;
    int m = (int)M;
    m = n;
    cout << "Приюлижение к спектральному радиусу матрицы H: rho = ";
    cout << Spectral_Radious(h_x) << endl;
    double EPS = 0.001;
    double** U = get_null_mat(n+1); //матрица (0...n) x (0...n) со значениями решения во всех узлах 
    double** F = get_null_mat(n+1); //матрица (0...n) x (0...n) со значениями правой части во всех узлах 
    double** real_sol = get_null_mat(n + 1); //точное решение 
    for (int i = 0; i <= n; i++)
    {
        U[i][0] = 0.0;              //значения решения на нижней границе
        U[i][m] = mu_x(h_x * i);    //значения решения на верхней границе 
        for (int j = 0; j <= n; j++)
        {
            F[i][j] = f(h_x * i, h_y * j); //значения правой части в узлах
            real_sol[i][j] = real_solution(h_x * i, h_y * j);
        }       
    }
    double** FF = get_null_mat(n);
    cop(FF, F, n);
    for (int j = 1; j < m; j++)
    {
        U[0][j] = 0.0;             //значения решения на левой грани
        U[n][j] = mu_y(h_y * j);   //значения решения на правой грани
    }
    U[0][0] = 0;
    U[0][m] = 0;
    U[n][0] = 0;
    U[n][m] = mu_x(l_x);
    double** TEST = get_null_mat(2 * n);    
    cout <<" --------------------------------------------------------------------------------------------------------------------------------------" << endl;
    cout << "Таблица значений решения в узлах:" << endl;
    show_mat(U, n);
    cop(TEST, U, n);
    cout << " " << endl;
    cout << "Таблица значений правой части в узлах:" << endl;
    show_mat(F, n);
    cout << " --------------------------------------------------------------------------------------------------------------------------------------" << endl;
    cout << "Таблица значений решения, метод простой итерации: " << endl;
    double** Easy_iter_vec = Easy_iterations(U, F, h_x, n, EPS);
    show_mat(Easy_iter_vec, n);
    cout << " --------------------------------------------------------------------------------------------------------------------------------------" << endl;
    cout << "Таблица значений решения, метод Зейделя: " << endl;
    double** Seid_vec = Seidel(U, F, h_x, n, EPS);
    show_mat(Seid_vec, n);
    cout << " --------------------------------------------------------------------------------------------------------------------------------------" << endl;
    cout << "Таблица значений решения, метод Врехней релаксации: " << endl;
    double** Upper_relax_vec = Top_relaxation(U, F, h_x, n, EPS);
    show_mat(Upper_relax_vec, n);
    cout << " --------------------------------------------------------------------------------------------------------------------------------------" << endl;
    cout << "Таблица значений решения, Попеременно-треугольный метод: " << endl;
    double** Triangle_vec = Triangle(U, F, h_x, n, EPS);
    show_mat(Triangle_vec, n);
    cout << "Значения правильного решения на сетке:" << endl;
    show_mat(real_sol, n);
    return 0;
}