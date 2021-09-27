
#ifndef _PARSE_CMD_H
#define _PARSE_CMD_H

#include <string.h>
#include <map>
#include <iostream>

using namespace std;

class Parameters
{
private:
    //set<parameter string, argu value>
    map<string, string> paras;
    //parsed parameter pairs count
    int parsed_parameter_count;

public:
    Parameters(const int argc, char **const argv)
    {
        if (argc % 2 == 0)
        {
            cout << "Error: parameter key and value count not match!!" << endl;
        }
        int i = 0;
        while (++i < argc)
        {
            string key = argv[i];
            string value = argv[++i];
            paras.emplace(map<string, string>::value_type(key, value));
            parsed_parameter_count++;
        }
    }
    int getParsedParameterCount() { return parsed_parameter_count; }

    string getParameterValue(const string key, string defaultValue = "")
    {
        if (paras.count(key))
            return paras.at(key);
        return defaultValue;
    }

    void printParamaterPairs()
    {
        cout << "Arguments list:" << endl;
        for (auto it = paras.begin(); it != paras.end(); ++it)
        {
            cout << it->first << "\t" << it->second << endl;
        }
    }
};

#endif
