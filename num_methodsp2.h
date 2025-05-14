#pragma once
#include <fstream>
#include <cmath>
#include <vector>



namespace nums {
    // ����� ��������� ������� U (����. �� ���� i=0..N�1, j=0..n�1 |U[i][j]|)
    // ��������� � ������� ����� ����������!
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






    // ����� ����� ���������� 2-�� ������� (����� ����������)
    //n     � ����������� ������ - ������� �������(���������� ������� � �������).
    //N     � ���������� ����� �����.
    //t     � �����(�������� �� ����������� ����������).
    //u     � ��������� �������(������[n,N] �� �������� - �������� ���������� ������� � ������ �����).
    //v     � �������, ���������� ��������� �������.
    //t0    � ��������� �������� ����������� ����������.
    //T     � �������� �������� ����������� ����������.
    //F     � ������ - ������� ������ �����.
    void ERK_2(int n, int N, double* t, double** u, double* v,
        double t0, double T,
        void F(double* ans, double t, double* u))
    {
        // �������������� ������� ��� k1, �������������� u � k2
        double* k1 = new double[n];
        double* utmp = new double[n];
        double* k2 = new double[n];

        // ��� ����� h
        double h = (T - t0) / (N - 1);

        // ��������� �������
        t[0] = t0;
        for (int i = 0; i < n; ++i) {
            u[0][i] = v[i];
        }

        // �������� ���� �� ����� �����
        for (int step = 0; step < N - 1; ++step) {
            double tk = t[step];

            // 1) ��������� k1 = F(tk, u[step])
            F(k1, tk, u[step]);

            // 2) ������������� �������� u*
            for (int i = 0; i < n; ++i) {
                utmp[i] = u[step][i] + h * k1[i];
            }

            // 3) ��������� k2 = F(tk + h, u*)
            F(k2, tk + h, utmp);

            // 4) ��������� �������� ������� � �����
            t[step + 1] = tk + h;
            for (int i = 0; i < n; ++i) {
                u[step + 1][i] = u[step][i] + (h / 2.0) * (k1[i] + k2[i]);
            }
        }

        // ����������� ������
        delete[] k1;
        delete[] utmp;
        delete[] k2;
    }

    // ��������! ��� ������� ���������� �� ���, ��� � ������� 1!
    // Ÿ �� ����� ������������ ��� ������� 1, � ��� ������ ���������
    // ��� ����� ������ ��� �������� ������ ERK_2_prec!
    // ���������� ERK_2 ��� ������������� � ERK_2_prec � ���������� ���������������� �����
    // n  � ����������� ������� ���������
    // N  � ����� ����� �� �������� (������) �����
    // u  � ������[N][n] ��� ������� �� �������� �����
    // v  � ��������� �������
    // t0 � ��������� �������� ����������� ����������
    // T  � �������� �������� ����������� ����������
    // F  � ������ ����� �������
    // r  � ������� ����� (���������� �������� �� ������ ���������)
    void ERK_2(int n, int N, double** u, double* v,
        double t0, double T,
        void F(double* ans, double t, double* u), int r) {
        // ����� ����� ��������� �����
        int Nf = (N - 1) * r + 1;
        // ��������� ��������� �������� ��� ��������� �����
        double* tf = new double[Nf];
        double** uf = new double* [Nf];
        for (int i = 0; i < Nf; ++i) {
            uf[i] = new double[n];
        }
        // ����� ������� ���������� ERK-2 �� ��������� �����
        ERK_2(n, Nf, tf, uf, v, t0, T, F);

        // ����� �������� ������� �� �������� ����� (������ r-� ���)
        for (int i = 0; i < N; ++i) {
            int idx = i * r;
            for (int j = 0; j < n; ++j) {
                u[i][j] = uf[idx][j];
            }
        }

        // ������������ ��������� ��������
        for (int i = 0; i < Nf; ++i) {
            delete[] uf[i];
        }
        delete[] uf;
        delete[] tf;
    }

    // ��������! ��� ������� ���������� �� ���, ��� � ������� 1!
    // Ÿ �� ����� ������������ ��� ������� 1, � ��� ������ ���������
    // ��� ����� ������ ��� �������� ������ ������ ��������!
    void ERK_2(int N, double* t, double* y,
        double y0, double dy0,
        double t0, double T,
        void F(double* ans, double t, double* u))
    {
        const int n = 2;  // �������: u � v = u'
        // �������� ������ ��� ������ ������� [N][n]
        double** U = new double* [N];
        for (int i = 0; i < N; ++i) {
            U[i] = new double[n];
        }

        // ������ ��������� ������� v = {u(a), v(a)}
        double v[2] = { y0, dy0 };

        // �������� �������� �������
        ERK_2(n, N, t, U, v, t0, T, F);

        // �������� ������ ���������� u (U[i][0]) � �������� ������ y[i]
        for (int i = 0; i < N; ++i) {
            y[i] = U[i][0];
        }

        // ����������� ������
        for (int i = 0; i < N; ++i) {
            delete[] U[i];
        }
        delete[] U;
    }

    // ������������� ����� ���������� ROS_1 (a = 0.5), ������� �������� 2
    // ������ ������� u' = F(t,u), ��������� ������� Fu(t,u).
    // n     � ����������� �������
    // N     � ����� ����� �����
    // t     � ������ ����� ����� ������� N
    // u     � ��������� ������ [N][n] �������� ������� � �����
    // v     � ��������� ������� (������ n)
    // t0,T  � ������ � ����� ������� ��������������
    // F     � ������� ������ �����: F(ans, t, u)
    // Fu    � �������, ����������� �������: Fu(ans, t, u)
    void ROS_1_a05(int n, int N, double* t, double** u, double* v,
        double t0, double T,
        void F(double* ans, double t, double* u),
        void Fu(double** ans, double t, double* u)) {
        // �������� �����
        const double gamma = 0.5;

        // �������������� �������
        double* w1 = new double[n];
        double* f = new double[n];
        double** M = new double* [n];
        for (int i = 0; i < n; ++i) {
            M[i] = new double[n];
        }

        // ������ ���� �� �������
        double h = (T - t0) / (N - 1);

        // ��������� �������
        for (int j = 0; j < n; ++j) {
            u[0][j] = v[j];
        }
        t[0] = t0;

        // �������� ���� �� �����
        for (int i = 0; i < N - 1; ++i) {
            // 1) ��������� f = F(t_i, u_i)
            F(f, t[i], u[i]);

            // 2) ��������� ������� J = Fu(t_i, u_i) � M
            Fu(M, t[i], u[i]);

            // 3) �������� ������� (I - gamma*h*J)
            for (int r = 0; r < n; ++r) {
                for (int c = 0; c < n; ++c) {
                    M[r][c] = -gamma * h * M[r][c];
                }
                M[r][r] += 1.0;
            }

            // 4) ������ ������� M * w1 = f
            Gauss(n, M, f, w1);

            // 5) ��������� u_{i+1} = u_i + h*w1 � t_{i+1}
            for (int j = 0; j < n; ++j) {
                u[i + 1][j] = u[i][j] + h * w1[j];
            }
            t[i + 1] = t[i] + h;
        }

        // ����������� ������
        delete[] w1;
        delete[] f;
        for (int i = 0; i < n; ++i) {
            delete[] M[i];
        }
        delete[] M;
    }




    // ����� ���������� ROS_1_a1 (L1-����������)
    // n  � ����������� �������
    // N  � ����� ����� �����
    // t  � ������ ����� �� ����������� ���������� (������ N)
    // u  � ������ �������: u[k][i] = u_i(t_k), ������ [N][n]
    // v  � ��������� �������: v[i] = u_i(t0)
    // t0 � ��������� �������� t
    // T  � �������� �������� t
    // F  � ������ ����� �������: F(ans, t, u), ans ����� n
    // Fu � ������� �����: Fu(ans, t, u), ans � n x n
    void ROS_1_a1(int n, int N, double* t, double** u, double* v,
        double t0, double T,
        void F(double* ans, double t, double* u),
        void Fu(double** ans, double t, double* u)) {
        // �������� ����� (gamma = 1 ��� ������� �������)
        const double gamma = 1.0;

        // �������������� �������
        double* w1 = new double[n];      // ������ ����������
        double* f = new double[n];      // F(t_i, u_i)
        double** M = new double* [n];     // ������� (I - gamma�h�J)
        for (int i = 0; i < n; ++i) {
            M[i] = new double[n];
        }

        // ��� �� �������
        double h = (T - t0) / (N - 1);

        // ��������� �������
        for (int j = 0; j < n; ++j) {
            u[0][j] = v[j];
        }
        t[0] = t0;

        // �������� ���� �� ����� �����
        for (int i = 0; i < N - 1; ++i) {
            // 1) ��������� f = F(t_i, u_i)
            F(f, t[i], u[i]);

            // 2) ��������� ������� J = Fu(t_i, u_i) � M
            Fu(M, t[i], u[i]);

            // 3) �������� ������� (I - gamma�h�J)
            for (int r = 0; r < n; ++r) {
                for (int c = 0; c < n; ++c) {
                    M[r][c] = -gamma * h * M[r][c];
                }
                M[r][r] += 1.0;
            }

            // 4) ������ ������� M * w1 = f
            Gauss(n, M, f, w1);

            // 5) ��������� u_{i+1} = u_i + h * w1 � t_{i+1}
            for (int j = 0; j < n; ++j) {
                u[i + 1][j] = u[i][j] + h * w1[j];
            }
            t[i + 1] = t[i] + h;
        }

        // ������������ ������
        delete[] w1;
        delete[] f;
        for (int i = 0; i < n; ++i) {
            delete[] M[i];
        }
        delete[] M;
    }



    // ����� ����� ���������� 2-�� ������� � ���������� ��������� ��������.
    // n     � ������ ������� (���������� ���������).
    // N     � ����� ����� �������� �����.
    // t     � ������ �� N �������� t[i], ���������� ���������� �� [t0,T].
    // u     � �����: ������ N x n, u[i][j] ~ ������� � ����� t[i], j-�� ����������.
    // v     � �����-n, ��������� ������� u( t0 ) = v.
    // t0,T  � ������� �� �������.
    // F     � ������� ������ �����: ans[0..n-1] = F( t, u[0..n-1] ).
    // eps   � ��������� �������� (� ��� �� �����).
    // r     � ��������� ��� �������� ����� (����� >=2).
    void ERK_2_prec(int n, int N, double** u, double* v, double t0, double T,
        void F(double* ans, double t, double* u), int r, double eps) {
        double R = 2 * eps; //����� ������, ��� �������� ������, ����� ���������� while
        int r1 = 1;
        int r2 = r;

        // ��������� ������ ��� ��������� �������� �������
        double** u1 = new double* [N];
        double** u2 = new double* [N];
        double** ud = new double* [N];
        for (int i = 0; i < N; ++i) {
            u1[i] = new double[n];
            u2[i] = new double[n];
            ud[i] = new double[n];
        }
        // ������������� �������� ��� ����������� ������������� (����� �������� �������������� � �������������������� ������)
        for (int j = 0; j < n; ++j) {
            u1[0][j] = v[j];
            u2[0][j] = v[j];
        }

        // ���������, ���� ��������� ����������� R > eps
        while (R > eps) {
            // ������ ����� � ��������� r1
            ERK_2(n, N, u1, v, t0, T, F, r1);
            // ��������� ����� � ��������� r2
            ERK_2(n, N, u2, v, t0, T, F, r2);

            // ���������� �������� ����� ������ � ��������� ���������
            for (int i = 0; i < N; ++i) {
                for (int j = 0; j < n; ++j) {
                    ud[i][j] = u2[i][j] - u1[i][j];
                }
            }
            double norm = CNorm(N, n, ud);
            // ����� ERK-2 ����� ������� p=2, ������� ����������� �������������� ��� r2^2
            R = norm / (std::pow(r2, 2) - 1.0);

            if (R > eps) {
                r1 = r2;
                r2 *= r;
                // ���������� ����� ��������� ������ ��� ��������� ��������
                for (int j = 0; j < n; ++j) {
                    u1[0][j] = v[j];
                    u2[0][j] = v[j];
                }
            }
        }

        // ����������� ��������� ������� � u �� �������� �����
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < n; ++j) {
                u[i][j] = u2[i][j];
            }
        }

        // ������������ ������ ��������� ��������
        for (int i = 0; i < N; ++i) {
            delete[] u1[i];
            delete[] u2[i];
            delete[] ud[i];
        }
        delete[] u1;
        delete[] u2;
        delete[] ud;
    }



    // ������� �������� ��� ������� ������ u'' = f(t,u,u')
    //double eps    � �������� �������� ������� ����������� ���������.
    //double g_0    � ����� ������� ������� �� ������� ������������ ������� �������� ���������������� ���������.
    //double g_1    � ������ ������� ������� �� ������� ������������ ������� �������� ���������������� ���������.
    //int N         � ���������� ����� �����.
    //double* t     � �����(�������� �� ����������� ����������).
    //double** u    � ��������� �������(������[n, N] �� �������� - �������� ���������� ������� � ������ �����).
    //double* u     � ������� ��������� ������� ������.
    //double a      � ����� ������� ���������.
    //double b      � ������ ������� ���������(�������������, ��� a < b).
    //double alpha  � �������� ������� ������� �� ����� �������.
    //double beta   � �������� ������� ������� �� ������ �������.
    //void F        � ������ - ������� ������ ����� ��������������� ������ ����, ������������ ���������� ������.
    void ShootingMethod(double eps,
        double g0, double g1,
        int N,
        double* t, double* u,
        double a, double b,
        double alpha, double beta,
        void F(double* ans, double t, double* u))
    {
        // ������� ��� ������������� ����������
        double* y_low = new double[N];
        double* y_high = new double[N];

        // ������ ������ ��� g0
        ERK_2(N, t, y_low, alpha, g0, a, b, F);
        double phi_low = y_low[N - 1] - beta;

        // ������ ������ ��� g1
        ERK_2(N, t, y_high, alpha, g1, a, b, F);
        double phi_high = y_high[N - 1] - beta;

        // �������� ��������� �����
        if (phi_low * phi_high > 0.0) {
            std::cerr << "ShootingMethod error: no sign change on [g0,g1]\n";
            delete[] y_low;
            delete[] y_high;
            return;
        }

        double s_mid = 0.0, phi_mid;
        // �������� �� s
        while (std::fabs(g1 - g0) > eps) {
            s_mid = 0.5 * (g0 + g1);

            // ������ ��� s_mid
            double* y_mid = new double[N];
            ERK_2(N, t, y_mid, alpha, s_mid, a, b, F);
            phi_mid = y_mid[N - 1] - beta;
            delete[] y_mid;

            // ������������� ���������
            if (phi_low * phi_mid <= 0.0) {
                g1 = s_mid;
                phi_high = phi_mid;
            }
            else {
                g0 = s_mid;
                phi_low = phi_mid;
            }
        }

        // ��������� ������ � ������ �������
        double s_final = 0.5 * (g0 + g1);
        ERK_2(N, t, u, alpha, s_final, a, b, F);

        // �������
        delete[] y_low;
        delete[] y_high;
    }

    // ��������! ��� ������� ����� ������ ��� ������ ��������� ������.
    // ��������������, ��� ��� ����� � ������� ����� ����������.
    // ����� �������� (�������� ������) ��� ������� ���� � ���������������� ��������
    // N    � ������ �������
    // a[i] � �������� ������ ��������� (i=0..N-2)
    // b[i] � �������� ������� ��������� (i=0..N-1)
    // c[i] � �������� ������� ��������� (i=0..N-2)
    // d[i] � ������ ����� ������� (i=0..N-1)
    // ANS  � ������ ��� ������ ������� (i=0..N-1)
    void Progonka(int N, double* a, double* b, double* c, double* d, double* ANS) {
        // ������� ����������� �������������
        double* ksi = new double[N];
        double* eta = new double[N];

        // ������ ���
        // ��� i = 0
        ksi[0] = -c[0] / b[0];
        eta[0] = d[0] / b[0];
        // ��� i = 1..N-2
        for (int i = 1; i < N - 1; ++i) {
            double denom = b[i] + a[i - 1] * ksi[i - 1];
            ksi[i] = -c[i] / denom;
            eta[i] = (d[i] - a[i - 1] * eta[i - 1]) / denom;
        }
        // ��� i = N-1 (���������� ������ eta)
        double denom_last = b[N - 1] + a[N - 2] * ksi[N - 2];
        eta[N - 1] = (d[N - 1] - a[N - 2] * eta[N - 2]) / denom_last;

        // �������� ���
        ANS[N - 1] = eta[N - 1];
        for (int i = N - 2; i >= 0; --i) {
            ANS[i] = ksi[i] * ANS[i + 1] + eta[i];
        }

        // ������������ ������
        delete[] ksi;
        delete[] eta;
    }


    // ���������� ��������� ������ ��� ������� ������ u'' + p(t) u' + q(t) u = f(t)
    // N       � ����� ����� �����
    // t       � ������ ����� ����� (����������� ������ �������)
    // u       � ������� ������� �� ����� �����
    // a, b    � ������� ��������� [a, b]
    // alpha   � ��������� ������� u(a) = alpha
    // beta    � ��������� ������� u(b) = beta
    // p(t), q(t), f(t) � �������� ������� ������������� � ������ �����
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
        // ��� �����
        double h = (b - a) / (N - 1);
        // ���������� ������� ����� �����
        for (int i = 0; i < N; ++i) {
            t[i] = a + i * h;
        }

        // ��� ������������ ������� � ������ �����
        double* A = new double[N - 1];  // ������������
        double* B = new double[N];      // ������� ���������
        double* C = new double[N - 1];  // ������������
        double* d = new double[N];      // ������ �����

        // ����� ��������� �������
        B[0] = 1.0;
        C[0] = 0.0;
        d[0] = alpha;

        // ���������� ����
        for (int i = 1; i < N - 1; ++i) {
            double ti = t[i];
            // ������������ ��� �������������
            A[i - 1] = 1.0 / (h * h) - p(ti) / (2.0 * h);
            B[i] = -2.0 / (h * h) + q(ti);
            C[i] = 1.0 / (h * h) + p(ti) / (2.0 * h);
            d[i] = f(ti);
        }

        // ������ ��������� �������
        A[N - 2] = 0.0;
        B[N - 1] = 1.0;
        d[N - 1] = beta;

        // ������ �������� (����� �������� ���������� ��������)
        Progonka(N, A, B, C, d, u);

        // ������������ ������
        delete[] A;
        delete[] B;
        delete[] C;
        delete[] d;
    }

}