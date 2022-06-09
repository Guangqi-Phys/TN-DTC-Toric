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
    int N_sample = 8000;
    double theta0; // error rate
    double theta0_prob;
    double theta1; // stabilizer operator
    double theta2; // logical x operator
    int const n_t = 10;
    int const n_simu = 50;
    double measure0;
    double measure1 = 0.0;
    double measure2 = 0.0;
    auto sites = SpinHalf(N, {"ConserveQNs=", false});
    auto state = InitState(sites);
    double measur1list[n_t] = {0};
    double measur2list[n_t] = {0};
    string filename;

    // PrintData(state);
    theta0 = 0.02 * M_PI;
    theta1 = 0.10 * M_PI;
    theta2 = 0.5 * M_PI;

    pcg_extras::seed_seq_from<random_device> seed_source;

    pcg32 engine0(seed_source);
    pcg32 engine1(seed_source);

    uniform_real_distribution<double> dist0(-0.02, 0.02);

    // define the logical z operator
    auto logzmpo = AutoMPO(sites);
    logzmpo += pow(2, (N / 4)), "Sz", 3, "Sz", 7, "Sz", 11;

    auto args = Args("Cutoff=", 1E-16, "MaxDim=", 1000);

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
                theta0_prob = theta0 + dist0(engine0);
                auto x_er = AutoMPO(sites);
                x_er += 2 * sin(theta0_prob) * Cplx_i, "Sx", i_er;
                auto xer_mpo = toMPO(x_er);
                // auto exp_xer_mpo = toExpH(x_er, theta0_prob * Cplx_i);
                // psi = applyMPO(exp_xer_mpo, psi, args).noPrime("Site");
                psi = sum(cos(theta0_prob) * psi, applyMPO(xer_mpo, psi, args).noPrime("Site"), args);
                psi.normalize();
            }

            // applying stabilizer operators
            for (int b = 0; b < N; b++)
            {
                auto stbmpo = AutoMPO(sites);
                if (b % 4 == 0)
                {
                    stbmpo += 16 * sin(theta1) * Cplx_i, "Sz", ((b % N) + 1), "Sz", (((b + 2) % N) + 1), "Sz", (((b + 3) % N) + 1), "Sz", (((b + 4) % N) + 1);
                }
                else if (b % 4 == 1)
                {
                    stbmpo += 16 * sin(theta1) * Cplx_i, "Sz", ((b % N) + 1), "Sz", (((b + 1) % N) + 1), "Sz", (((b + 2) % N) + 1), "Sz", (((b + 4) % N) + 1);
                }
                else if (b % 4 == 2)
                {
                    stbmpo += 16 * sin(theta1) * Cplx_i, "Sx", ((b % N) + 1), "Sx", (((b + 2) % N) + 1), "Sx", (((b + 3) % N) + 1), "Sx", (((b + 4) % N) + 1);
                }
                else if (b % 4 == 3)
                {
                    stbmpo += 16 * sin(theta1) * Cplx_i, "Sx", ((b % N) + 1), "Sx", (((b + 1) % N) + 1), "Sx", (((b + 2) % N) + 1), "Sx", (((b + 4) % N) + 1);
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

            auto psi_1 = psi;

            if (time % 2 == 0)
            {
                measure1 = -measure_average(psi_1, N_sample);
            }
            else
            {
                measure1 = measure_average(psi_1, N_sample);
            }

            measur1list[time] = measur1list[time] + measure1 / n_simu;

            // PrintData(measure1);

            // apply logical x operator to state psi

            auto logxmpo = AutoMPO(sites);
            logxmpo += 4 * sin(theta2) * Cplx_i, "Sx", 7, "Sx", 8;
            auto Hlogx = toMPO(logxmpo);
            psi = sum(cos(theta2) * psi, applyMPO(Hlogx, psi, args).noPrime("Site"), args);
            psi.normalize();

            // auto logxmpo = AutoMPO(sites);
            // logxmpo += 4,"Sx", 7, "Sx", 8;
            // auto Hlogx = toMPO(logxmpo);
            // auto expHlogx = toExpH(logxmpo, theta2 * Cplx_i, args);
            // psi = ApplyMPO(expHlogx, psi, args).noPrime("Site");
            // psi.normalize();

            // for (int i = 7; i < 9; i++)
            // {
            //     auto logxmpo = AutoMPO(sites);
            //     logxmpo += sin(theta2) * Cplx_i, "Sx", i;
            //     auto Hlogx = toMPO(logxmpo);
            //     psi = sum(cos(theta2) * psi, applyMPO(Hlogx, psi, args).noPrime("Site"), args);
            //     psi.normalize();
            // }

            // measure the logical Z operation
            auto psi_2 = psi;
            if (time % 2 == 0)
            {
                measure2 = -measure_average(psi_2, N_sample);
            }
            else
            {
                measure2 = measure_average(psi_2, N_sample);
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