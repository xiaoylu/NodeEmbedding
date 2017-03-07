#include "sample.h"
#include "hierar_infer.h"
#include <iostream>
# include <chrono>
using namespace std;
using seconds = chrono::seconds;
using get_time = chrono::steady_clock;

using namespace std;

int main(int argc, char *argv[])
{
    //srand(time(NULL))
    srand(0); // re-produce results

    if (argc == 1)
    {
        cout << "Wrong parameter format. Please type:" << endl;
        cout << argv[0] << endl;
        cout << "input=<input txtfile>" << endl;
        cout << "          required." << endl;
        cout << "d=<dimension: length of the output vectors>" << endl;
        cout << "          optional. default=20" << endl;
        cout << "OW=<the observation window of a cascades, infections after OW are ignored>" << endl;
        cout << "          optional. default=1.01" << endl;
        cout << "stepsize=<the stepsize of the gradient descent>" << endl;
        cout << "          optional. default=0.0001" << endl;
        cout << "negative_sampling=<the number of negative samples chosen>" << endl;
        cout << "          optional. default=20" << endl;
        cout << "max_iterations=<the maximum number of iterations>" << endl;
        cout << "          optional. default=20" << endl;
        cout << "early_stop=<learning process terminates if the loglikelihood does not increase for early_stop iterations>" << endl;
        cout << "          optional. default=2" << endl;
        cout << "epsilon=<the cutoff threshold for projected gradient descent, i.e. set v[i] = 0 if v[i] <= epsilon>" << endl;
        cout << "          optional. default=0" << endl;
        return 0;
    }

    map<string, string> kwargs;
    kwargs["input"] = "../data/cascade-synthentic-graph.txt";
    kwargs["d"] = "20"; 
    kwargs["OW"] = "1.01"; 
    kwargs["stepsize"] = "0.0001";
    kwargs["negative_sampling"] = "20"; 
    kwargs["max_iterations"] = "20"; 
    kwargs["epsilon"] = "0";
    kwargs["early_stop"] = "2";

    for (int i = 0; i < argc; ++i) {
        char name[256], val[256];
        if (sscanf(argv[i], "%[^=]=%[^\n]", name, val) == 2) {
            kwargs[name] = val;
        }
    }

    cout << "Reading Data..." << endl;
    HierarInfer algo(kwargs["input"]);
    cout << "Data Loaded" << endl;

    auto start = get_time::now(); 

    algo.filter(0.01, 2); // delete illegal cascades
    algo.normalize(); // make the average duration to be 1.0 

    algo.infer(stoi(kwargs["d"]), stof(kwargs["OW"]),
               stof(kwargs["stepsize"]), stof(kwargs["epsilon"]),
               stoi(kwargs["negative_sampling"]), stoi(kwargs["max_iterations"]),
               stoi(kwargs["early_stop"])
              );

    auto end = get_time::now();
    auto diff = end - start;
    cout << "Elapsed time is :  " << 
            chrono::duration_cast<seconds>(diff).count()
         << " seconds " << endl;

    algo.output(); // write the A, B matrix to file

    return 0;
}
