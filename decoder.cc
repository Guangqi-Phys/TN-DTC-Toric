#include "Python/Python.h"
#include "itensor/all.h"
#include "itensor/util/print_macro.h"
#include "libs/pcg-cpp/include/pcg_random.hpp"
#include "include/sampling.h"
#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <random>

using namespace itensor;
using namespace std;
using std::vector;

int main(int argc, char *argv[])
{

    int N = 12;
    int N_sample = 500;
    double theta0; // error rate
    double theta0_prob;
    double theta1; // stabilizer operator
    double theta2; // logical x operator
    int const n_t = 5;
    int const n_simu = 50;
    double measure0;
    double measure1;
    double measure2;
    auto sites = SpinHalf(N, {"ConserveQNs=", false});
    auto state = InitState(sites);
    double measur1list[n_t] = {0};
    double measur2list[n_t] = {0};
    string filename;

    // PrintData(state);
    theta0 = 0.02 * M_PI;
    theta1 = 0.25 * M_PI;
    theta2 = 0.5 * M_PI;

    pcg_extras::seed_seq_from<random_device> seed_source;

    pcg32 engine0(seed_source);
    pcg32 engine1(seed_source);

    uniform_real_distribution<double> dist0(-0.02, 0.02);

    // define the logical z operator
    auto logzmpo = AutoMPO(sites);
    logzmpo += pow(2, (N / 4)), "Sz", 3, "Sz", 7, "Sz", 11;

    auto args = Args("Cutoff=", 1E-9, "MaxDim=", 500);

    // Make initial MPS psi to be in the all up state
    for (auto j : range1(N))
    {
        state.set(j, "Up");
    }

    ofstream outfile;
    filename = string("data/decoder_") + string("N=") + to_string(N) + string(".dat");

    outfile.open(filename);

    Py_Initialize();

    // define and apply the sttabilizer operator to state psi
    for (int simu = 0; simu < n_simu; simu++)
    {
        auto psi = MPS(state);

        // measure0 = real(innerC(psi, logzmpo, psi));
        // PrintData(real(measure0));

        for (int time = 1; time < n_t; time++)
        {

            // adding errors

            for (int i_er = 1; i_er < N; ++i_er)
            {
                auto x_er = AutoMPO(sites);
                x_er += "Sx", i_er;
                auto xer_mpo = toMPO(x_er);
                theta0_prob = theta0 + dist0(engine0);
                auto exp_xer_mpo = toExpH(x_er, theta0_prob * Cplx_i);
                psi = applyMPO(exp_xer_mpo, psi, args).noPrime("Site");
            }

            // applying stabilizer operators
            for (int b = 0; b < N; b++)
            {
                auto stbmpo = AutoMPO(sites);
                if (b % 4 == 0)
                {
                    stbmpo += sin(theta1) * Cplx_i, "Sz", ((b % 12) + 1), "Sz", (((b + 2) % 12) + 1), "Sz", (((b + 3) % 12) + 1), "Sz", (((b + 4) % 12) + 1);
                }
                if (b % 4 == 1)
                {
                    stbmpo += sin(theta1) * Cplx_i, "Sz", ((b % 12) + 1), "Sz", (((b + 1) % 12) + 1), "Sz", (((b + 2) % 12) + 1), "Sz", (((b + 4) % 12) + 1);
                }
                if (b % 4 == 2)
                {
                    stbmpo += sin(theta1) * Cplx_i, "Sx", ((b % 12) + 1), "Sx", (((b + 2) % 12) + 1), "Sx", (((b + 3) % 12) + 1), "Sx", (((b + 4) % 12) + 1);
                }
                if (b % 4 == 3)
                {
                    stbmpo += sin(theta1) * Cplx_i, "Sx", ((b % 12) + 1), "Sx", (((b + 1) % 12) + 1), "Sx", (((b + 2) % 12) + 1), "Sx", (((b + 4) % 12) + 1);
                }
                // note that need to do noPrime operation after applyMPO() function.
                auto Hstb = toMPO(stbmpo);
                psi = sum(cos(theta1) * psi, applyMPO(Hstb, psi, args).noPrime("Site"), args);
                psi.normalize();
            }

            // PrintData(psi);

            // vector<int> vec1 = sample(psi);
            // PrintData(vec1[1]);

            // measure the logical Z operator
            if (time % 2 == 0)
            {
                measure1 = -measure_average(psi, N_sample);
            }
            else
            {
                measure1 = measure_average(psi, N_sample);
            }

            measur1list[time] = measur1list[time] + measure1 / n_simu;

            // PrintData(measure1);

            // apply logical x operator to state psi

            auto logxmpo = AutoMPO(sites);
            logxmpo += sin(theta2) * Cplx_i, "Sx", 7, "Sx", 8;
            auto Hlogx = toMPO(logxmpo);
            psi = sum(cos(theta2) * psi, applyMPO(Hlogx, psi, args).noPrime("Site"), args);
            psi.normalize();

            // auto logxmpo = AutoMPO(sites);
            // logxmpo += "Sx", 3, "Sx", 4;
            // auto Hlogx = toMPO(logxmpo);
            // auto expHlogx = toExpH(logxmpo, theta2 * Cplx_i, args);
            // psi = applyMPO(expHlogx, psi, args).noPrime("Site");
            // psi.normalize();

            // for (int i = 7; i < 9; i++)
            // {
            //     auto logxmpo = AutoMPO(sites);
            //     logxmpo += "Sx", i;
            //     auto Hlogx = toMPO(logxmpo);
            //     auto expHlogx = toExpH(logxmpo, theta2 * Cplx_i);
            //     psi = applyMPO(expHlogx, psi, args).noPrime("Site");
            // }
            psi.normalize();
            // PrintData(psi);

            // measure the logical Z operation
            if (time % 2 == 0)
            {
                measure2 = -measure_average(psi, N_sample);
            }
            else
            {
                measure2 = measure_average(psi, N_sample);
            }

            measur2list[time] = measur2list[time] + measure2 / n_simu;
            // PrintData(measure2);

            // vector<int> vec = sample(psi);

            // for (int i = 1; i < N; i++)
            //     print(vec[i]);
        }
    }

    // PrintData(measure2);

    for (int time = 1; time < n_t; time++)
    {
        outfile << time << " " << measur1list[time] << " " << measur2list[time] << endl;
    }

    Py_Finalize();
    outfile.close();

    return 0;
}