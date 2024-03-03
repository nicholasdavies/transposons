// config.cpp

#include "config.h"
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <limits>
#include <iomanip>
#include <regex>
#include <cctype>
using namespace std;

// HELPER FUNCTIONS

// ConvertFromString
//  helper template functors to convert strings into any needed types.
//  to convert string to string, no transformation is done;
//  to convert string to bool, "false" is false and all other values are true;
//  to convert string to vector, comma-delimited tokens are passed to appropriate handler for value type;
//  to convert string to LuaFunction, LuaFunction constructor is used with string as compact form;
//  for other conversions, uses istringstream.
template <typename Type>
struct ConvertFromString {
    Type operator()(string s)
    {
        Type return_value;
        istringstream iss(s);
        iss >> return_value;
        return return_value;
    }
};

template <>
struct ConvertFromString<string> {
    string operator()(string s)
    {
        return s;
    }
};

template <>
struct ConvertFromString<bool> {
    bool operator()(string s)
    {
        return s != "false";
    }
};

template <typename Subtype>
struct ConvertFromString<vector<Subtype>> {
    vector<Subtype> operator()(string s)
    {
        ConvertFromString<Subtype> subconv;
        vector<Subtype> value_vector;
        if (s.empty())
            return value_vector;
        
        string::size_type start = 0, end = 0;
        while (end != string::npos)
        {
            end = s.find(',', start);
            value_vector.push_back(subconv(s.substr(start, end - start)));
            start = end + 1;
        }
        return value_vector;
    }
};

template <typename Subtype>
struct ConvertFromString<LuaFunction<Subtype>> {
    LuaFunction<Subtype> operator()(string s)
    {
        return LuaFunction<Subtype>(s);
    }
};

// ConvertToString
//  helper template functors to convert parameter types to strings for
//  writing to an output stream, ensuring that necessary precision is kept
//  for numerical types.
template <typename Type>
struct ConvertToString {
    string operator()(Type const& v)
    {
        ostringstream oss;
        oss << setprecision(numeric_limits<Type>::digits10 + 1) << v;
        return oss.str();
    }
};

template<>
struct ConvertToString<string> {
    string operator()(string const& v)
    {
        return v;
    }
};

template<>
struct ConvertToString<bool> {
    string operator()(bool const& v)
    {
        return v ? "true" : "false";
    }
};

template<typename Subtype>
struct ConvertToString<vector<Subtype>> {
    string operator()(vector<Subtype> const& v)
    {
        ConvertToString<Subtype> subconv;
        string total;
        for (unsigned int i = 0; i < v.size(); ++i)
        {
            total += subconv(v[i]);
            if (i + 1 != v.size()) total += ',';
        }
        return total;
    }
};

template<typename Subtype>
struct ConvertToString<LuaFunction<Subtype>> {
    string operator()(LuaFunction<Subtype> const& v)
    {
        return v.CompactForm();
    }
};

// Trim
//  helper function to trim whitespace from a string.
string Trim(string s)
{
    string::size_type start = s.find_first_not_of(" \t");
    if (start == string::npos)
        return string();
    string::size_type end = s.find_last_not_of(" \t");
    return s.substr(start, end - start + 1);
}

// SetParameterFromMap
//  assigns a value to the supplied parameter.
//  if there is an entry in the name-value map with the correct parameter name,
//  that entry is used; otherwise the parameter is unchanged.
template <typename Type>
void SetParameterFromMap(Type& parameter,
                         string parameter_name,
                         Parameters::NameValueMap& nvm)
{
    auto map_entry = nvm.find(parameter_name);
    if (map_entry != nvm.end())
    {
        ConvertFromString<Type> conv;
        parameter = conv(map_entry->second);
    }
}

// PARAMETERS METHODS

// Parameters constructor
//  sets parameters to default values.
Parameters::Parameters()
{
    Reset();
}

// Reset
//  resets the parameters, restoring them to their defaults.
void Parameters::Reset()
{
    _nvm_sweeps = vector<NameValueMap>(1);
    _nvm_sweep_names = vector<string>(1, "Main");
    _nvm_override.clear();
    _sweep = 0;
    _assignment_virgin = true;
    _template_mode = false;
    SetAllToDefault();
}

// Read (istream)
//  load parameters from an input stream, ignoring //-style comments
//  and interpreting statements of the form pname = value as assignments.
//  Statements of the form [SweepName] name new parameter sweeps.
void Parameters::Read(istream& in)
{
    InterpretLines(in);
    AssignFromMap(_nvm_sweeps[_sweep]);
}

// Read (string)
//  load parameters from a named file, using the input stream method above.
void Parameters::Read(string filename)
{
    ifstream fin(filename);
    if (!fin.good())
        throw runtime_error("Could not load parameter file " + filename);
    InterpretLines(fin);
    AssignFromMap(_nvm_sweeps[_sweep]);
}

// Read (argc, argv)
//  load parameters from the command line; parameter names are specified by
//  prefixing with '-'. Values are optionally taken from the following token;
//  if not provided, the parameter is set to true. If argv[1] is not prefixed
//  with '-', it is taken as a filename from which to load parameters before
//  overriding from the command line. If a config file is specified in that
//  manner, the next argument can be a sweep range in the format N or N-M,
//  where N and M are integers (M = N by default) specifying a range of sweeps
//  to limit the run by.
void Parameters::Read(int argc, char* argv[])
{
    // First try load parameters from parameter file, if applicable.
    int first = 1;
    if (argc > first && argv[first][0] != '-')
    {
        string filename = argv[first];
        ifstream fin(filename);
        if (!fin.good())
            throw runtime_error("Could not load parameter file " + filename);
        InterpretLines(fin);
        first = 2;

        // Attempt to read sweep range in the form e.g. "12" or "12-15"
        if (argc > first && argv[first][0] != '-')
        {
            ConvertFromString<int> conv;
            int sweep_start, sweep_end;
            string range = argv[first];
            
            auto dash = range.find('-');
            if (dash != string::npos)
            {
                sweep_start = conv(range.substr(0, dash));
                sweep_end = conv(range.substr(dash + 1));
            }
            else
            {
                sweep_start = sweep_end = conv(range);
            }
            if (sweep_start < 1 || sweep_end < 1 || sweep_start > NSweeps() || sweep_end > NSweeps() || sweep_start > sweep_end)
            {
                throw runtime_error("Invalid sweep range " + range);
            }

            // combine any skipped sweeps into the first non-skipped sweep
            for (int i = sweep_start - 2; i >= 0; --i)
            {
                _nvm_sweeps[sweep_start - 1].insert(_nvm_sweeps[i].begin(), _nvm_sweeps[i].end());
            }

            _nvm_sweeps = vector<NameValueMap>(_nvm_sweeps.begin() + sweep_start - 1, _nvm_sweeps.begin() + sweep_end);
            _nvm_sweep_names = vector<string>(_nvm_sweep_names.begin() + sweep_start - 1, _nvm_sweep_names.begin() + sweep_end);

            first = 3;
        }
    }

    // Now load overriding parameters from command line.
    string parameter_name;
    for (int a = first; a < argc; ++a)
    {
        string token = argv[a];
        // If token starts with a hyphen, interpret as parameter name, unless it seems to be a negative number.
        if (token[0] == '-' && token.size() > 1 && !isdigit(token[1]))
        {
            parameter_name = token.substr(1);
            _nvm_override[parameter_name] = "true";
        }
        else
        {
            if (parameter_name.empty())
            {
                throw runtime_error("Unexpected token in command line: " + token);
            }
            _nvm_override[parameter_name] = token;
            parameter_name.clear();
        }
    }

    AssignFromMap(_nvm_sweeps[_sweep]);
    AssignFromMap(_nvm_override);   
}

// Write
//  write parameters to the specified output stream, in a format which is
//  human- and Parameters- readable.
// NOTE: Here I prefix local identifiers with an underscore
// and fully qualify global identifiers to avoid name clashes.
void Parameters::Write(ostream& _out) const
{
    const int _col1_width = 18;
    const int _col2_width = 24;

    // INCLUSION PASS 2: Write to output stream.
    #define PARAMETER(Type, Name, Default) \
      { ConvertToString<Type> _conv; \
        _out.width(_col1_width); _out << std::left << #Name << " = "; \
        _out.width(_col2_width); _out << std::left << _conv(Name); \
        _out << " // " << #Type << ", default = " << #Default; \
        _out << std::endl; }
    #include "config_def.h"
    #undef PARAMETER
}

// Sweep
//  Get the number of the current sweep.
int Parameters::Sweep() const
{
    return _sweep;
}

// SweepName
//  Get the name of the current sweep.
string Parameters::SweepName() const
{
    return _nvm_sweep_names[_sweep];
}

// NSweeps
//  Get the total number of sweeps.
int Parameters::NSweeps() const
{
    return (int)_nvm_sweeps.size();
}

// NextSweep
//  Advances to the next sweep.
void Parameters::NextSweep()
{
    ++_sweep;
    if (_sweep == NSweeps()) // silently stay on last sweep (useful for loops)
        _sweep = NSweeps() - 1;

    AssignFromMap(_nvm_sweeps[_sweep]);
    AssignFromMap(_nvm_override);
}

// GoToSweep
//  Sets parameters to those of a given sweep.
//  This is less efficient than jumping forward a certain number of times using NextSweep
//  if that is an option.
void Parameters::GoToSweep(int s)
{
    if (s > NSweeps())
    {
        throw runtime_error("Requested sweep in GoToSweep larger than maximum.");
    }

    SetAllToDefault();
    for (int S = 0; S < s; ++S)
    {
        NextSweep();
    }
}

// InterpretLines
//  call InterpretLine on sequential lines of in.
void Parameters::InterpretLines(istream& in)
{
    while (in.good())
    {
        string line;
        getline(in, line);
        InterpretLine(line);
    }
}

// InterpretLine
//  first removes //-style comments from the line.
//  then interprets the following style of lines:
//      [SweepName]
//      SweepTemplate <A,B,C>
//      [SweepTemplateInstantiation]: SweepTemplate<0, 1, 2>
//      parameter = value
void Parameters::InterpretLine(string line)
{
    // Remove comments and leading/trailing whitespace from the line.
    string::size_type comment_start = line.find("//");
    if (comment_start != string::npos)
        line.resize(comment_start);
    line = Trim(line);
    if (line.empty())
        return;

    // Match interpretable lines.
    smatch matches;

    // Sweep heading
    if (regex_match(line, matches, regex("\\[\\s*(\\w+)\\s*\\]")))
    {
        // If we haven't seen an assignment yet, set the current sweep name to the supplied name.
        if (_assignment_virgin)
        {
            _nvm_sweep_names.back() = matches[1];
        }
        // Otherwise, make a new sweep.
        else
        {
            _nvm_sweeps.push_back(NameValueMap());
            _nvm_sweep_names.push_back(matches[1]);
        }

        _template_mode = false;
    }

    // Template declaration
    else if (regex_match(line, matches, regex("(\\w+)\\s*<(\\s*\\w+\\s*(\\s*,\\s*\\w+)*\\s*)>")))
    {
        // Make a new template.
        _templates.push_back(NameValueMap());
        _template_names.push_back(matches[1]);
        ConvertFromString<vector<string>> conv;
        _template_params.push_back(conv(matches[2]));

        _template_mode = true;
    }

    // Template instantiation
    else if (regex_match(line, matches, regex("\\[\\s*(\\w+)\\s*\\]\\s*:\\s*(\\w+)\\s*<(.*)>")))
    {
        // Make sure this is a valid template
        auto templ = find(_template_names.begin(), _template_names.end(), matches[2]);
        if (templ != _template_names.end())
        {
            auto which = templ - _template_names.begin();
            
            // Get the substitutions
            ConvertFromString<vector<string>> conv;
            auto subst = conv(matches[3]);
            if (subst.size() != _template_params[which].size())
            {
                cout << "Config: incorrect number of template parameters in line [" << line << "].\n";
            }
            else
            {
                // If we haven't seen an assignment yet, set the current sweep name to the supplied name.
                if (_assignment_virgin)
                {
                    _nvm_sweep_names.back() = matches[1];
                }
                // Otherwise, make a new sweep.
                else
                {
                    _nvm_sweeps.push_back(NameValueMap());
                    _nvm_sweep_names.push_back(matches[1]);
                }

                // Fill the sweep as though each assignment in the invoked template were now run, with appropriate substitutions.
                // Include the special substitution <$Name> -> name of this instantiation.
                for (auto& entry : _templates[which])
                {
                    string value = entry.second;
                    if (value.find('<') != string::npos)    // Don't bother trying any argument replacement if there are no <s.
                    {
                        for (unsigned int p = 0; p < _template_params[which].size(); ++p)
                            value = regex_replace(value, regex("<\\s*" + Trim(_template_params[which][p]) + "\\s*>"), Trim(subst[p]));
                        value = regex_replace(value, regex("<\\s*\\$Name\\s*>"), string(matches[1]));
                    }
                    _nvm_sweeps.back()[entry.first] = value;
                }

                _assignment_virgin = false;
            }
        }

        _template_mode = false;
    }

    // Assignment
    else if (regex_match(line, matches, regex("(\\w+)\\s*=\\s*(.*)")))
    {
        // Record the value of this parameter, either into the current sweep or the current template.
        if (_template_mode)
        {
            _templates.back()[matches[1]] = matches[2];
        }
        else
        {
            _nvm_sweeps.back()[matches[1]] = matches[2];
            _assignment_virgin = false;
        }
    }

    // Uninterpretable line
    else
    {
        cout << "Config: could not interpret line [" << line << "].\n";
    }
}

// SetAllToDefault
//  sets all defined parameters to the corresponding default value.
void Parameters::SetAllToDefault()
{
    // INCLUSION PASS 3: Set default values for each defined parameter.
    #define PARAMETER(Type, Name, Default) Name = (Type)(Default);
    #include "config_def.h"
    #undef PARAMETER
}

// AssignFromMap
//  attempts to assign a value to all defined parameters from the provided
//  name-value map.
void Parameters::AssignFromMap(NameValueMap& _nvm)
{
    // INCLUSION PASS 4: Call SetParameterFromMap on each defined parameter.
    #define PARAMETER(Type, Name, Default) \
        SetParameterFromMap<Type>(Name, #Name, _nvm);
    #include "config_def.h"
    #undef PARAMETER
}
