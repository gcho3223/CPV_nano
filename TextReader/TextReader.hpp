// TextReader v2.00

#ifndef __TextReader_hpp_Inclueded__
#define __TextReader_hpp_Inclueded__

#include <string>
#include <vector>
#include <map>
#include <iterator>

typedef std::string st;
typedef std::vector<double> vec_d;
typedef vec_d::iterator vec_d_it;
typedef std::vector<st> vec_s;
typedef vec_s::iterator vec_s_it;
typedef std::map<st, double> map_d;
typedef map_d::iterator map_d_it;
typedef std::map<st, st> map_s;
typedef std::map<st, st>::iterator map_s_it;
typedef std::map<st, vec_d> map_sd;
typedef std::map<st, vec_d>::iterator map_sd_it;
typedef std::map<st, vec_s> map_ss;
typedef std::map<st, vec_s>::iterator map_ss_it;

class TextReader {

public:

    TextReader();
    ~TextReader();
    void ReadFile(st);
    void ReadVariables();
    void PrintoutVariables();

    double GetNumber(st,bool PrintError = true);
    double GetNumber(st,int,bool PrintError = true);
    int GetNumberInt(st,bool PrintError = true);
    int GetNumberInt(st,int,bool PrintError = true);
    unsigned int GetNumberUint(st,bool PrintError = true);
    unsigned int GetNumberUint(st,int,bool PrintError = true);
    float GetNumberFloat(st,bool PrintError = true);
    float GetNumberFloat(st,int,bool PrintError = true);
    double GetNumberDouble(st,bool PrintError = true);
    double GetNumberDouble(st,int,bool PrintError = true);
    st GetText(st,bool PrintError = true);
    st GetText(st,int,bool PrintError = true);
    bool GetBool(st,bool PrintError = true);
    bool GetBool(st,int,bool PrintError = true);
    bool Check(st);
    bool CheckNumber(st);
    bool CheckText(st);
    bool CheckArray(st);
    int Size(st);

    map_d ReturnMap(map_d);
    map_s ReturnMap(map_s);
    map_sd ReturnMap(map_sd);
    map_ss ReturnMap(map_ss);

private:

    FILE *configfile;
    map_d CVf; // Cut_Variable_floating-point_type
    map_s CVt; // Cut_Variable_text_type
    map_sd CVfa; //Cut_Variable_floating-point_type_array
    map_ss CVta; // Cut_Variable_text_type_array
    map_d_it CVf_it;
    map_s_it CVt_it;
    map_sd_it CVfa_it;
    map_ss_it CVta_it;

};

#endif // __TextReader_hpp_Inclueded__
