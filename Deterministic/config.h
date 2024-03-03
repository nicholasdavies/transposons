// config.h
// Module for providing strongly-typed program options with default values,
// which can be set through configuration files and/or the command line.

#ifndef CONFIG_H
#define CONFIG_H

#include <string>
#include <vector>
#include <iostream>
#include <map>
#include "luafunc.h"
using std::string;
using std::vector;
using LuaFunc::LuaFunction;

class Parameters
{
public:
    typedef std::map<string, string> NameValueMap;

    Parameters();

    void Reset();

    void Read(std::istream& in);
    void Read(string filename);
    void Read(int argc, char* argv[]);

    void Write(std::ostream& out) const;

    int Sweep() const;
    string SweepName() const;
    int NSweeps() const;
    void NextSweep();
    void GoToSweep(int s);

    // INCLUSION PASS 1: Declare all parameters as public members of this class.
    template<typename T> using Func = LuaFunction<T>;
    #define PARAMETER(Type, Name, Default) Type Name;
    #include "config_def.h"
    #undef PARAMETER

private:
    void InterpretLines(std::istream& in);
    void InterpretLine(string line);
    void SetAllToDefault();
    void AssignFromMap(NameValueMap& nvm);

    vector<NameValueMap> _nvm_sweeps;
    vector<string> _nvm_sweep_names;

    vector<NameValueMap> _templates;
    vector<string> _template_names;
    vector<vector<string>> _template_params;

    NameValueMap _nvm_override;

    int _sweep;
    bool _assignment_virgin;
    bool _template_mode;
};

#endif  // CONFIG_H