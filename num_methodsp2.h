#pragma once
#include <fstream>
#include <cmath>
#include <vector>



namespace nums {
    // Норма максимума матрицы U (макс. по всем i=0..N–1, j=0..n–1 |U[i][j]|)
    // Находится в скрытой части практикума!
    double CNorm(int N, int n, double** U) {
        double max_abs = 0.0;
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < n; ++j) {
                max_abs = std::max(max_abs, std::fabs(U[i][j]));
            }
        }
        return max_abs;
    }

    void Gauss(int N, double** A, double* B, double* ANS) {
        for (int i = 0; i < N; i++) {
            if (A[i][i] == 0) {
                for (int j = i + 1; j < N; j++) {
                    if (A[j][i] != 0) {
                        for (int k = i; k < N; k++) {
                            double trans = A[i][k];
                            A[i][k] = A[j][k];
                            A[j][k] = trans;
                        }
                        double trans = B[i];
                        B[i] = B[j];
                        B[j] = trans;
                        break;
                    }
                }
            }


            double del = A[i][i];
            for (int j = i; j < N; j++) {
                A[i][j] /= del;
            }
            B[i] /= del;



            for (int j = i + 1; j < N; j++) {
                double m = A[j][i] / A[i][i];
                for (int k = i; k < N; k++) {
                    A[j][k] -= m * A[i][k];

                }
                B[j] -= m * B[i];
            }
            for (int i = 0; i < N; i++) {
                ANS[i] = B[i];
            }
        }





        for (int i = N - 1; i > -1; i--) {
            for (int j = i - 1; j > -1; j--) {
                B[j] -= A[j][i] / A[i][i] * B[i];
            }
        }
    }






    // Явная схема Рунге–Кутты 2-го порядка (метод Рунге–Хейна)
    //n     — размерность вектор - функции решения(количество функций в системе).
    //N     — количество узлов сетки.
    //t     — сетка(введённая на независимой переменной).
    //u     — численное решение(массив[n,N] из столбцов - значений численного решения в каждой точке).
    //v     — столбец, содержащий начальные условия.
    //t0    — начальное значение независимой переменной.
    //T     — конечное значение независимой переменной.
    //F     — вектор - функция правой части.
    void ERK_2(int n, int N, double* t, double** u, double* v,
        double t0, double T,
        void F(double* ans, double t, double* u))
    {
        // Дополнительные массивы для k1, промежуточного u и k2
        double* k1 = new double[n];
        double* utmp = new double[n];
        double* k2 = new double[n];

        // Шаг сетки h
        double h = (T - t0) / (N - 1);

        // Начальные условия
        t[0] = t0;
        for (int i = 0; i < n; ++i) {
            u[0][i] = v[i];
        }

        // Основной цикл по узлам сетки
        for (int step = 0; step < N - 1; ++step) {
            double tk = t[step];

            // 1) Вычисляем k1 = F(tk, u[step])
            F(k1, tk, u[step]);

            // 2) Промежуточное значение u*
            for (int i = 0; i < n; ++i) {
                utmp[i] = u[step][i] + h * k1[i];
            }

            // 3) Вычисляем k2 = F(tk + h, u*)
            F(k2, tk + h, utmp);

            // 4) Обновляем значение решения и время
            t[step + 1] = tk + h;
            for (int i = 0; i < n; ++i) {
                u[step + 1][i] = u[step][i] + (h / 2.0) * (k1[i] + k2[i]);
            }
        }

        // Освобождаем память
        delete[] k1;
        delete[] utmp;
        delete[] k2;
    }

    // ВНИМАНИЕ! Эта функция отличается от той, что в задании 1!
    // Её НЕ НУЖНО использовать для задания 1, у них разные сигнатуры
    // Она нужна только для проверки работы ERK_2_prec!
    // ПЕРЕГРУЗКА ERK_2 для использования в ERK_2_prec с адаптивным масштабированием сетки
    // n  — размерность системы уравнений
    // N  — число точек на исходной (грубой) сетке
    // u  — массив[N][n] для решения на исходной сетке
    // v  — начальные условия
    // t0 — начальное значение независимой переменной
    // T  — конечное значение независимой переменной
    // F  — правая часть системы
    // r  — масштаб сетки (количество подшагов на каждом интервале)
    void ERK_2(int n, int N, double** u, double* v,
        double t0, double T,
        void F(double* ans, double t, double* u), int r) {
        // число точек уточнённой сетки
        int Nf = (N - 1) * r + 1;
        // выделение временных массивов для уточнённой сетки
        double* tf = new double[Nf];
        double** uf = new double* [Nf];
        for (int i = 0; i < Nf; ++i) {
            uf[i] = new double[n];
        }
        // вызов базовой реализации ERK-2 на уточнённой сетке
        ERK_2(n, Nf, tf, uf, v, t0, T, F);

        // выбор значений решения на исходной сетке (каждый r-й шаг)
        for (int i = 0; i < N; ++i) {
            int idx = i * r;
            for (int j = 0; j < n; ++j) {
                u[i][j] = uf[idx][j];
            }
        }

        // освобождение временных массивов
        for (int i = 0; i < Nf; ++i) {
            delete[] uf[i];
        }
        delete[] uf;
        delete[] tf;
    }

    // ВНИМАНИЕ! Эта функция отличается от той, что в задании 1!
    // Её НЕ НУЖНО использовать для задания 1, у них разные сигнатуры
    // Она нужна только для проверки работы метода стрельбы!
    void ERK_2(int N, double* t, double* y,
        double y0, double dy0,
        double t0, double T,
        void F(double* ans, double t, double* u))
    {
        const int n = 2;  // система: u и v = u'
        // Выделяем память под полное решение [N][n]
        double** U = new double* [N];
        for (int i = 0; i < N; ++i) {
            U[i] = new double[n];
        }

        // Вектор начальных условий v = {u(a), v(a)}
        double v[2] = { y0, dy0 };

        // Вызываем исходную функцию
        ERK_2(n, N, t, U, v, t0, T, F);

        // Копируем только компоненту u (U[i][0]) в выходной массив y[i]
        for (int i = 0; i < N; ++i) {
            y[i] = U[i][0];
        }

        // Освобождаем память
        for (int i = 0; i < N; ++i) {
            delete[] U[i];
        }
        delete[] U;
    }

    // Одностадийная схема Розенброка ROS_1 (a = 0.5), порядок точности 2
    // Решает систему u' = F(t,u), используя якобиан Fu(t,u).
    // n     — размерность системы
    // N     — число узлов сетки
    // t     — массив узлов сетки размера N
    // u     — двумерный массив [N][n] значений решения в узлах
    // v     — начальные условия (размер n)
    // t0,T  — начало и конец отрезка интегрирования
    // F     — функция правой части: F(ans, t, u)
    // Fu    — функция, заполняющая якобиан: Fu(ans, t, u)
    void ROS_1_a05(int n, int N, double* t, double** u, double* v,
        double t0, double T,
        void F(double* ans, double t, double* u),
        void Fu(double** ans, double t, double* u)) {
        // Параметр схемы
        const double gamma = 0.5;

        // Дополнительные массивы
        double* w1 = new double[n];
        double* f = new double[n];
        double** M = new double* [n];
        for (int i = 0; i < n; ++i) {
            M[i] = new double[n];
        }

        // Расчет шага по времени
        double h = (T - t0) / (N - 1);

        // Начальные условия
        for (int j = 0; j < n; ++j) {
            u[0][j] = v[j];
        }
        t[0] = t0;

        // Основной цикл по шагам
        for (int i = 0; i < N - 1; ++i) {
            // 1) вычисляем f = F(t_i, u_i)
            F(f, t[i], u[i]);

            // 2) вычисляем якобиан J = Fu(t_i, u_i) в M
            Fu(M, t[i], u[i]);

            // 3) собираем матрицу (I - gamma*h*J)
            for (int r = 0; r < n; ++r) {
                for (int c = 0; c < n; ++c) {
                    M[r][c] = -gamma * h * M[r][c];
                }
                M[r][r] += 1.0;
            }

            // 4) решаем систему M * w1 = f
            Gauss(n, M, f, w1);

            // 5) обновляем u_{i+1} = u_i + h*w1 и t_{i+1}
            for (int j = 0; j < n; ++j) {
                u[i + 1][j] = u[i][j] + h * w1[j];
            }
            t[i + 1] = t[i] + h;
        }

        // Освобождаем память
        delete[] w1;
        delete[] f;
        for (int i = 0; i < n; ++i) {
            delete[] M[i];
        }
        delete[] M;
    }




    // Схема Розенброка ROS_1_a1 (L1-устойчивая)
    // n  — размерность системы
    // N  — число узлов сетки
    // t  — массив узлов по независимой переменной (размер N)
    // u  — массив решений: u[k][i] = u_i(t_k), размер [N][n]
    // v  — начальные условия: v[i] = u_i(t0)
    // t0 — начальное значение t
    // T  — конечное значение t
    // F  — правая часть системы: F(ans, t, u), ans длины n
    // Fu — матрица Якоби: Fu(ans, t, u), ans — n x n
    void ROS_1_a1(int n, int N, double* t, double** u, double* v,
        double t0, double T,
        void F(double* ans, double t, double* u),
        void Fu(double** ans, double t, double* u)) {
        // Параметр схемы (gamma = 1 для первого порядка)
        const double gamma = 1.0;

        // Дополнительные массивы
        double* w1 = new double[n];      // вектор приращений
        double* f = new double[n];      // F(t_i, u_i)
        double** M = new double* [n];     // матрица (I - gamma·h·J)
        for (int i = 0; i < n; ++i) {
            M[i] = new double[n];
        }

        // Шаг по времени
        double h = (T - t0) / (N - 1);

        // Начальные условия
        for (int j = 0; j < n; ++j) {
            u[0][j] = v[j];
        }
        t[0] = t0;

        // Основной цикл по узлам сетки
        for (int i = 0; i < N - 1; ++i) {
            // 1) вычисляем f = F(t_i, u_i)
            F(f, t[i], u[i]);

            // 2) вычисляем якобиан J = Fu(t_i, u_i) в M
            Fu(M, t[i], u[i]);

            // 3) собираем матрицу (I - gamma·h·J)
            for (int r = 0; r < n; ++r) {
                for (int c = 0; c < n; ++c) {
                    M[r][c] = -gamma * h * M[r][c];
                }
                M[r][r] += 1.0;
            }

            // 4) решаем систему M * w1 = f
            Gauss(n, M, f, w1);

            // 5) обновляем u_{i+1} = u_i + h * w1 и t_{i+1}
            for (int j = 0; j < n; ++j) {
                u[i + 1][j] = u[i][j] + h * w1[j];
            }
            t[i + 1] = t[i] + h;
        }

        // Освобождение памяти
        delete[] w1;
        delete[] f;
        for (int i = 0; i < n; ++i) {
            delete[] M[i];
        }
        delete[] M;
    }



    // Явная схема Рунге–Кутты 2-го порядка с адаптивным контролем точности.
    // n     — размер системы (количество компонент).
    // N     — число узлов исходной сетки.
    // t     — массив из N значений t[i], равномерно разбиенных на [t0,T].
    // u     — выход: массив N x n, u[i][j] ~ решение в точке t[i], j-ая компонента.
    // v     — длина-n, начальные условия u( t0 ) = v.
    // t0,T  — границы по времени.
    // F     — функция правой части: ans[0..n-1] = F( t, u[0..n-1] ).
    // eps   — требуемая точность (в той же норме).
    // r     — множитель для сгущения сетки (целое >=2).
    void ERK_2_prec(int n, int N, double** u, double* v, double t0, double T,
        void F(double* ans, double t, double* u), int r, double eps) {
        double R = 2 * eps; //Задаём больше, чем величина ошибки, чтобы выполнился while
        int r1 = 1;
        int r2 = r;

        // выделение памяти для временных массивов решений
        double** u1 = new double* [N];
        double** u2 = new double* [N];
        double** ud = new double* [N];
        for (int i = 0; i < N; ++i) {
            u1[i] = new double[n];
            u2[i] = new double[n];
            ud[i] = new double[n];
        }
        // инициализация массивов для безопасного использования (чтобы подавить предупреждения о неинициализированной памяти)
        for (int j = 0; j < n; ++j) {
            u1[0][j] = v[j];
            u2[0][j] = v[j];
        }

        // уточнение, пока оценочная погрешность R > eps
        while (R > eps) {
            // грубая сетка с масштабом r1
            ERK_2(n, N, u1, v, t0, T, F, r1);
            // уточнённая сетка с масштабом r2
            ERK_2(n, N, u2, v, t0, T, F, r2);

            // вычисление разности между грубым и уточнённым решениями
            for (int i = 0; i < N; ++i) {
                for (int j = 0; j < n; ++j) {
                    ud[i][j] = u2[i][j] - u1[i][j];
                }
            }
            double norm = CNorm(N, n, ud);
            // схема ERK-2 имеет порядок p=2, поэтому погрешность масштабируется как r2^2
            R = norm / (std::pow(r2, 2) - 1.0);

            if (R > eps) {
                r1 = r2;
                r2 *= r;
                // определить новые начальные данные для следующей итерации
                for (int j = 0; j < n; ++j) {
                    u1[0][j] = v[j];
                    u2[0][j] = v[j];
                }
            }
        }

        // копирование итогового решения в u на исходной сетке
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < n; ++j) {
                u[i][j] = u2[i][j];
            }
        }

        // освобождение памяти временных массивов
        for (int i = 0; i < N; ++i) {
            delete[] u1[i];
            delete[] u2[i];
            delete[] ud[i];
        }
        delete[] u1;
        delete[] u2;
        delete[] ud;
    }



    // Функция стрельбы для краевой задачи u'' = f(t,u,u')
    //double eps    — желаемая точность решения нелинейного уравнения.
    //double g_0    — левая граница отрезка на котором локализовано искомое значение вспомогательного параметра.
    //double g_1    — правая граница отрезка на котором локализовано искомое значение вспомогательного параметра.
    //int N         — количество узлов сетки.
    //double* t     — сетка(введённая на независимой переменной).
    //double** u    — численное решение(массив[n, N] из столбцов - значений численного решения в каждой точке).
    //double* u     — искомое численное решение задачи.
    //double a      — левая граница интервала.
    //double b      — правая граница интервала(гарантируется, что a < b).
    //double alpha  — значение искомой функции на левой границе.
    //double beta   — значение искомой функции на правой границе.
    //void F        — вектор - функция правой части вспомогательной задачи Коши, определяющая конкретную задачу.
    void ShootingMethod(double eps,
        double g0, double g1,
        int N,
        double* t, double* u,
        double a, double b,
        double alpha, double beta,
        void F(double* ans, double t, double* u))
    {
        // Массивы для промежуточных интеграций
        double* y_low = new double[N];
        double* y_high = new double[N];

        // Первый прогон для g0
        ERK_2(N, t, y_low, alpha, g0, a, b, F);
        double phi_low = y_low[N - 1] - beta;

        // Второй прогон для g1
        ERK_2(N, t, y_high, alpha, g1, a, b, F);
        double phi_high = y_high[N - 1] - beta;

        // Проверка изменения знака
        if (phi_low * phi_high > 0.0) {
            std::cerr << "ShootingMethod error: no sign change on [g0,g1]\n";
            delete[] y_low;
            delete[] y_high;
            return;
        }

        double s_mid = 0.0, phi_mid;
        // Бисекция по s
        while (std::fabs(g1 - g0) > eps) {
            s_mid = 0.5 * (g0 + g1);

            // Прогон для s_mid
            double* y_mid = new double[N];
            ERK_2(N, t, y_mid, alpha, s_mid, a, b, F);
            phi_mid = y_mid[N - 1] - beta;
            delete[] y_mid;

            // Корректировка интервала
            if (phi_low * phi_mid <= 0.0) {
                g1 = s_mid;
                phi_high = phi_mid;
            }
            else {
                g0 = s_mid;
                phi_low = phi_mid;
            }
        }

        // Финальный прогон и запись решения
        double s_final = 0.5 * (g0 + g1);
        ERK_2(N, t, u, alpha, s_final, a, b, F);

        // Очистка
        delete[] y_low;
        delete[] y_high;
    }

    // ВНИМАНИЕ! Эта функция нужна только для работы сеточного метода.
    // Предполагается, что она будет в скрытой части практикума.
    // Метод прогонки (алгоритм Томаса) для решения СЛАУ с трехдиагональной матрицей
    // N    — размер системы
    // a[i] — элементы нижней диагонали (i=0..N-2)
    // b[i] — элементы главной диагонали (i=0..N-1)
    // c[i] — элементы верхней диагонали (i=0..N-2)
    // d[i] — правая часть системы (i=0..N-1)
    // ANS  — массив для записи решения (i=0..N-1)
    void Progonka(int N, double* a, double* b, double* c, double* d, double* ANS) {
        // Массивы прогоночных коэффициентов
        double* ksi = new double[N];
        double* eta = new double[N];

        // Прямой ход
        // Для i = 0
        ksi[0] = -c[0] / b[0];
        eta[0] = d[0] / b[0];
        // Для i = 1..N-2
        for (int i = 1; i < N - 1; ++i) {
            double denom = b[i] + a[i - 1] * ksi[i - 1];
            ksi[i] = -c[i] / denom;
            eta[i] = (d[i] - a[i - 1] * eta[i - 1]) / denom;
        }
        // Для i = N-1 (используем только eta)
        double denom_last = b[N - 1] + a[N - 2] * ksi[N - 2];
        eta[N - 1] = (d[N - 1] - a[N - 2] * eta[N - 2]) / denom_last;

        // Обратный ход
        ANS[N - 1] = eta[N - 1];
        for (int i = N - 2; i >= 0; --i) {
            ANS[i] = ksi[i] * ANS[i + 1] + eta[i];
        }

        // Освобождение памяти
        delete[] ksi;
        delete[] eta;
    }


    // Реализация сеточного метода для краевой задачи u'' + p(t) u' + q(t) u = f(t)
    // N       — число узлов сетки
    // t       — массив узлов сетки (заполняется внутри функции)
    // u       — искомое решение на узлах сетки
    // a, b    — границы интервала [a, b]
    // alpha   — граничное условие u(a) = alpha
    // beta    — граничное условие u(b) = beta
    // p(t), q(t), f(t) — заданные функции коэффициентов и правой части
    void GridMethod(int N,
        double* t,
        double* u,
        double a,
        double b,
        double alpha,
        double beta,
        double p(double),
        double q(double),
        double f(double)) {
        // Шаг сетки
        double h = (b - a) / (N - 1);
        // Заполнение массива узлов сетки
        for (int i = 0; i < N; ++i) {
            t[i] = a + i * h;
        }

        // Три диагональные матрицы и правая часть
        double* A = new double[N - 1];  // поддиагональ
        double* B = new double[N];      // главная диагональ
        double* C = new double[N - 1];  // наддиагональ
        double* d = new double[N];      // правая часть

        // Левое граничное условие
        B[0] = 1.0;
        C[0] = 0.0;
        d[0] = alpha;

        // Внутренние узлы
        for (int i = 1; i < N - 1; ++i) {
            double ti = t[i];
            // Коэффициенты для дискретизации
            A[i - 1] = 1.0 / (h * h) - p(ti) / (2.0 * h);
            B[i] = -2.0 / (h * h) + q(ti);
            C[i] = 1.0 / (h * h) + p(ti) / (2.0 * h);
            d[i] = f(ti);
        }

        // Правое граничное условие
        A[N - 2] = 0.0;
        B[N - 1] = 1.0;
        d[N - 1] = beta;

        // Запуск прогонки (метод прогонки реализован отдельно)
        Progonka(N, A, B, C, d, u);

        // Освобождение памяти
        delete[] A;
        delete[] B;
        delete[] C;
        delete[] d;
    }

}