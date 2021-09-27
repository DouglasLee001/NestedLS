#include "build_free.h"
#include "update.h"
#include "checker.h"
#include "construct_solution.h"
#include "fastds-pure.h"
#include "fastds_pure_small.h"
#include "parse_cmd.h"
#include <sstream>

int main(int argc, char *argv[])
{
    Parameters *para = new Parameters(argc, argv);
    string instanceName = para->getParameterValue("--i");
    seed = atoi(para->getParameterValue("--s", to_string(seed)).c_str());
    cutoff_time = atoi(para->getParameterValue("--t", to_string(cutoff_time)).c_str());
    delete para;
    BuildInstance(instanceName); 
    srand(seed);
    start = chrono::steady_clock::now();
    enter_ls();
    cout <<fixed<<setprecision(2)<<"best weight: \n"<<bestWeight<<endl<<"found time:\n"<< best_comp_time << endl;
    FreeMemory();
    return 0;
}
