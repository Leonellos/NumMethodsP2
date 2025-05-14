#pragma once
#include <cmath>
#include <iomanip>
//#include "nums_alexander.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


// ��������� ������� ��� ������������
namespace nums_service {

    // ��� ERK2
        // ��������� ans[i] = f_i(t, u)
    void F_exp_decay(double* ans, double t, double* u) {
        // �������� lambda ��� �������� ������ u' = -lambda u
        const double lambda = 1.0;
        // n = 1, ���� ���������
        ans[0] = -lambda * u[0];
    }

    // ������ ����� F ��� �������� ������ u' = -u
    void F_test(double* ans, double t, double* u) {
        // ans[0] = du/dt = -u[0]
        ans[0] = -u[0];
    }

    // ������� Fu ��� �������� ������: dF/du = -1
    void Fu_test(double** ans, double t, double* u) {
        // ans ��� ������� n x n; ����� n=1
        ans[0][0] = -1.0;
    }


    // ������ ����� �������: f = -lambda * u
    void F_exp(double* ans, double t, double* u) {
        const double lambda = 100.0; // ����������� � 1/�
        ans[0] = -lambda * u[0];
    }

    // ������� F �� u: df/du = -lambda
    void Fu_exp(double** ans, double t, double* u) {
        const double lambda = 100.0; // ����������� � 1/�
        ans[0][0] = -lambda;
    }

    // ��������������� F ��� ������� ������ u'' + u = 0:
    // u(0)=0, u(pi/2)=1, ������ ������� u(t)=sin(t)
    void F2_shoot(double* ans, double t, double* u) {
        // u[0] = u, u[1] = v = u'
        ans[0] = u[1];     // u' = v
        ans[1] = -u[0];    // v' = u'' = -u
    }

    // ������ u'' = -\pi^2 sin(\pi t) �� [0,1] � u(0)=u(1)=0 (������ ������� sin(\pi t))
    static double p_test(double) { return 0.0; }
    static double q_test(double) { return 0.0; }
    static double f_test(double t) {
        return -M_PI * M_PI * std::sin(M_PI * t);
    }
}

namespace nums {

    /*
    --------------------------------------
    �������� ������� ��� ��������� �������
    --------------------------------------
    */

    void test_ERK2() {
        std::cout << "\n\t ���� ERK2...\n";
        // ��������� �����
        const int n = 2;       // ������ �������
        const int N = 101;     // ����� �����
        const double t0 = 0.0; // ��������� �����
        const double T = 5.0; // �������� �����

        // �������� ������ ��� �������
        std::vector<double> t(N);
        std::vector<double> v(n);            // ��������� ������� u(0)
        std::vector<std::vector<double>> u(N, std::vector<double>(n));
        std::vector<double*> uptr(N);
        for (int i = 0; i < N; ++i) {
            uptr[i] = u[i].data();
        }

        // ����� ��������� ������� u(0)=1
        v[0] = 1.0;
        v[1] = 0;

        // ��������� ����� ��2
        nums::ERK_2(n, N, t.data(), uptr.data(), v.data(), t0, T, nums_service::F_test);
        //nums_alex::ERK_2(n, N, t.data(), uptr.data(), v.data(), t0, T, nums_service::F_test);

        // ����� ����������� � ������� � � CSV-����
        std::ofstream fout("solution.csv");
        fout << "t,u\n";
        std::cout << "\tt\t|\tu\n";
        std::cout << "-------------------------------\n";
        for (int i = 0; i < N; ++i) {
            std::cout<< "\t" << t[i] << "\t|\t" << u[i][0] << "\n";
            fout << t[i] << "," << u[i][0] << "\n";
        }
        fout.close();

        std::cout << "\n������� ��������� � ����� solution.csv\n";

    }


    void test_ROS_1_a05() {
        const int n = 1;        // ����������� �������
        const int N = 11;       // ����� ����� �����
        const double t0 = 0.0;  // ��������� �����, �
        const double T = 1.0;  // �������� �����, �

        // �������� ������ ��� ����� � �������
        double* t = new double[N];
        double** u = new double* [N];
        for (int i = 0; i < N; ++i) {
            u[i] = new double[n];
        }
        double v[n];
        v[0] = 1.0;  // u(0)=1
        //v[1] = 0.0;

        // ������ ������
        ROS_1_a05(n, N, t, u, v, t0, T, nums_service::F_test, nums_service::Fu_test);
        //nums_alex::ROS_1_a05(n, N, t, u, v, t0, T, nums_service::F_test, nums_service::Fu_test);

        // ����� ����������� � ������ �����������
        std::cout << " t [s]    u_numeric      u_exact       error\n";
        std::cout << "---------------------------------------------\n";
        for (int k = 0; k < N; ++k) {
            double t_k = t[k];
            double u_num = u[k][0];
            double u_ex = std::exp(-t_k);
            double err = std::fabs(u_num - u_ex);
            std::cout << std::fixed
                << " " << std::setw(6) << std::setprecision(3) << t_k
                << "  " << std::setw(12) << std::setprecision(6) << u_num
                << "  " << std::setw(12) << std::setprecision(6) << u_ex
                << "  " << std::setw(12) << std::setprecision(6) << err
                << "\n";
        }

        // ����������� ������
        for (int i = 0; i < N; ++i) {
            delete[] u[i];
        }
        delete[] u;
        delete[] t;
    }


    void test_ROS_1_a1() {
        // ��������� ������
        const int n = 1;            // ����������� �������
        const double t0 = 0.0;      // ��������� �����, �
        const double T = 0.1;       // �������� �����, �
        const int N = 11;           // ����� �����
        const double lambda = 100.0;// �������� ���������, 1/�

        // �������� ������ ��� ����� � �������
        double* t = new double[N];
        double** u = new double* [N];
        for (int k = 0; k < N; ++k) {
            u[k] = new double[n];
        }
        double v_init = 1.0;        // ��������� ������� u(0)=1
        double v[1] = { v_init };

        // ��������� �����
        ROS_1_a1(n, N, t, u, v, t0, T, nums_service::F_exp, nums_service::Fu_exp);

        // ����� ����������� � ������ �����������
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "  t(s)    u_num      u_exact    abs_error\n";
        for (int k = 0; k < N; ++k) {
            double t_k = t[k];
            double u_num = u[k][0];
            double u_ex = std::exp(-lambda * t_k);
            double err = std::fabs(u_num - u_ex);
            std::cout << std::setw(7) << t_k << "  "
                << std::setw(9) << u_num << "  "
                << std::setw(9) << u_ex << "  "
                << std::setw(9) << err << "\n";
        }

        // ����������� ������
        delete[] t;
        for (int k = 0; k < N; ++k) {
            delete[] u[k];
        }
        delete[] u;
    }


    // ����-������� ��� ERK_2_prec.
    void test_ERK_2_prec() {
        const int n = 2;          // ����������� �������
        const int N = 11;         // ����� ����� �� �������� �����
        const double t0 = 0.0;    // ��������� �����, �
        const double T = 1.0;     // �������� �����, �
        const int r = 2;          // ������� ��������� (���������� ��������)
        const double eps = 1e-6;  // ��������� ��������

        // �������� ������ ��� ������� u[N][n] � ��������� ������� v[n]
        double** u = new double* [N];
        for (int i = 0; i < N; ++i) {
            u[i] = new double[n];
        }
        double* v = new double[n];
        v[0] = 1.0;  // u(0)=1
        v[1] = 0.0;

        // ��������� ���������� ERK2
        ERK_2_prec(n, N, u, v, t0, T, nums_service::F_exp, r, eps);
        //nums_alex::ERK_2_prec(n, N, u, v, t0, T, nums_service::F_exp, r, eps);

        // ������������� �������: u(T) = exp(-T)
        double u_exact = std::exp(-T);
        double u_num = u[N - 1][0];
        double err = std::fabs(u_exact - u_num);

        // ����� �����������
        std::cout << "���� ERK_2_prec �� u' = -u, u(0)=1:\n";
        std::cout << " ��������� = " << u_num << "\n";
        std::cout << " ������������� = " << u_exact << "\n";
        std::cout << " ���������� ����������� = " << err << "\n";

        // ����������� ������
        for (int i = 0; i < N; ++i) {
            delete[] u[i];
        }
        delete[] u;
        delete[] v;
    }

    // �������� ������� ��� ������ ��������
    void TestShootingMethod() {
        const int N = 100;                // ����� ����� �����
        const double a = 0.0;             // ����� �����
        const double b = M_PI / 2.0;      // ������ �����
        const double alpha = 0.0;         // u(a) = 0
        const double beta = 1.0;         // u(b) = 1
        const double eps = 1e-6;        // �������� �� ��������� s
        const double g0 = 0.0;         // ������ ������� ��� u'(a)
        const double g1 = 2.0;         // ������� ������� ��� u'(a)

        // �������� �������
        double* t = new double[N];
        double* u = new double[N];

        // ��������� ����� ��������
        ShootingMethod(eps, g0, g1, N, t, u, a, b, alpha, beta, nums_service::F2_shoot);

        // ������� ������������� ������
        std::cout << "    t\t ��������� u\t ������ u\t ������\n";
        for (int i = 0; i < N; ++i) {
            double u_exact = std::sin(t[i]);
            double err = u[i] - u_exact;
            std::cout << t[i]
                << "\t" << u[i]
                << "\t" << u_exact
                << "\t" << err
                << "\n";
        }

        // ������ ���������� ���������� ������� s = u'(a)
        double h = t[1] - t[0];
        double s_approx = (u[1] - u[0]) / h;
        std::cout << "\n����������� ��������� s = u'(0): "
            << s_approx << "\n";

        // ����� ��������� �������� �� t = b
        double u_end = u[N - 1];
        double err_end = u_end - beta;
        std::cout << "�������� ����� ��� t = " << b << ":\n"
            << "  ��������� u(b) = " << u_end << "\n"
            << "  ������   u(b) = " << beta << "\n"
            << "  ������        = " << err_end << "\n";

        delete[] t;
        delete[] u;
    }


    // �������� ������� ��� GridMethod
    void TestGridMethod() {
        const int N = 11;
        const double a = 0.0, b = 1.0;
        const double alpha = 0.0, beta = 0.0;
        double* t = new double[N];
        double* u = new double[N];
        GridMethod(N, t, u, a, b, alpha, beta, nums_service::p_test, nums_service::q_test, nums_service::f_test);
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "i\tt\t\tu\n";
        for (int i = 0; i < N; ++i) {
            std::cout << i << '\t' << t[i] << '\t' << u[i] << '\n';
        }
        delete[] t; delete[] u;
    }

}