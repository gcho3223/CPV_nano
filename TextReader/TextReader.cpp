// TextReader v2.00

#include <stdio.h>  // fopen
#include <stdlib.h> // atof
#include <iomanip>  // setw
#include <iostream> // cout
#include "TextReader.hpp"

TextReader::TextReader(){

}

TextReader::~TextReader() {

    fclose(configfile);

}
void TextReader::ReadFile(st inputconfig){
    //if ((configfile = fopen(inputconfig.c_str(),"r")) == NULL) {
    if ((configfile = fopen(inputconfig.c_str(),"r")) == 0) {
	std::cout << std::endl << "Check input file (" << inputconfig << ")." << std::endl << std::endl;
    }  // need error

}

void TextReader::ReadVariables() {

    int state = 0;
    int text_state = 0;
    int value_state = 0;
    int array_state = 0;

    st Variable_name = "";
    st Variable_value = "";
    st Variable_array = "";

    vec_d Variable_floatingpoint_array;
    vec_s Variable_text_array;

    int input_check;
    st input_buffer;

    const char *array_buffer;
    st array_input;
    int array_index = 0;

    while ((input_check = fgetc(configfile)) != EOF) {

	if (input_check == '\n') {

	    if (Variable_name.length() > 0) {
		switch (array_state) {
		    case 0:
			switch (text_state) {
			    case 0:
				CVf[Variable_name] = atof(Variable_value.c_str());
				break;
			    case 1:
				CVt[Variable_name] = Variable_value;
				break;
			}
			break;
		    case 1:
			switch (text_state) {
			    case 0:
				CVfa[Variable_name] = Variable_floatingpoint_array;
				Variable_floatingpoint_array.clear();
				break;
			    case 1:
				CVta[Variable_name] = Variable_text_array;
				Variable_text_array.clear();
				break;
			}
		}
		Variable_name = "";
		Variable_value = "";
	    } // buffer to map

	    state = 0;
	    text_state = 0;
	    value_state = 0;
	    array_state = 0;
	    array_index = 0;
	    continue; // reset state for new line

	} // if - line end

	input_buffer = char(input_check);

	switch (state) {
	    case 0: // line start
		switch (input_check) {
		    case ' ':
		    case '\t':
			break;
		    case '#':
			state = 99; //comment
			break;
		    default:
			state = 1; // variable name start
			Variable_name.append(input_buffer);
			break;
		}
		break;

	    case 1: // variable name start
		switch (input_check) {
		    case ' ':
		    case '\t':
			state = 20;
			break;
		    case ':':
			state = 21;
			break;
		    default:
			Variable_name.append(input_buffer);
			break;
		}
		break;

	    case 20: // variable name end
		switch (input_check) {
		    case ' ':
		    case '\t':
			break;
		    case ':':
			state = 21;
			break;
		    default:
			state = 99;
			break;
		}
		break;

	    case 21: // after :
		switch (input_check) {
		    case ' ':
		    case '\t':
			break;
		    case '"':
			state = 4;
			text_state = 1;
			break;
		    case '{':
			state = 5;
			array_state = 1;
			break;
		    default:
			state = 3;
			Variable_value.append(input_buffer);
			break;
		}
		break;

	    case 3: // variable value start
		Variable_value.append(input_buffer);
		break;

	    case 4: // variable text type
		switch (input_check) {
		    case '"':
			state = 99; // text end
			break;
		    default:
			Variable_value.append(input_buffer);
			break;
		}
		break;

	    case 5: // variable array type
		switch (input_check) {
		    case '}':
			Variable_value.append(input_buffer);
			array_buffer = Variable_value.c_str();

			if (array_buffer[0] == '"') {
			    text_state = 1;
			    while (array_buffer[array_index] != '}') {
				switch (array_buffer[array_index]) {
				    case '"':
					if (value_state == 0) {
					    value_state = 1;
					}
					else {
					    value_state = 0;
					    Variable_text_array.push_back(Variable_array); 
					    Variable_array = "";
					}
					break; 
				    default:
					if (value_state == 1) {
					    array_input = array_buffer[array_index];
					    Variable_array.append(array_input);
					}
					break;
				}
				array_index++;
			    }
			}
			else {
			    while (array_buffer[array_index] != '}') {
				switch (array_buffer[array_index]) {
				    case ' ':
				    case '\t':
					break;
				    case ',':
					Variable_floatingpoint_array.push_back(atof(Variable_array.c_str()));
					Variable_array = "";
					break;
				    default:
					array_input = array_buffer[array_index];
					Variable_array.append(array_input);
					break;
				}
				array_index++;
			    }
			    Variable_floatingpoint_array.push_back(atof(Variable_array.c_str()));
			    Variable_array = "";

			}
			state = 99;
			break;
		    default:
			Variable_value.append(input_buffer);
			break;
		}
		break;

	    case 99: // End
		break;

	} // switch - variable input end

    } // while - EOF end

} // ReadVariable end

double TextReader::GetNumber(st CutName, bool PrintError){
    if ((CVf_it = CVf.find(CutName.c_str())) != CVf.end()) {
	return CVf_it->second;
    }
    else if ((CVfa_it = CVfa.find(CutName.c_str())) != CVfa.end()) {
	if(PrintError) std::cout << "Warning - GetNumber(\"" << CutName << "\") : \"" << CutName << "\" is array" << std::endl;
	return (CVfa_it->second).at(0);
    }
    else {
	if(PrintError) std::cout << "Error - GetNumber(\"" << CutName << "\") : Cannot find \"" << CutName << "\"" << std::endl;
	return -999.00;
    }
}

double TextReader::GetNumber(st CutName, int array_index, bool PrintError){
    if ((CVfa_it = CVfa.find(CutName.c_str())) != CVfa.end()) {
	if (array_index > 0 && (CVfa_it->second).size() >= array_index) {
	    return (CVfa_it->second).at(array_index-1);
	}
	else {
	//    return (CVfa_it->second).at(0);
	    if(PrintError) std::cout << "Error - GetNumber(\"" << CutName << "\"," << array_index << ") : Cannot find index \"" << array_index << "\"" << std::endl;
	    return -999.00;
	}

    }
    else if ((CVf_it = CVf.find(CutName.c_str())) != CVf.end()) {
	if(PrintError) std::cout << "Warning - GetNumber(\"" << CutName << "\"," << array_index << ") : \"" << CutName << "\" is not array" << std::endl;
	return CVf_it->second;
    }
    else {
	if(PrintError) std::cout << "Error - GetNumber(\"" << CutName << "\"," << array_index << ") : Cannot find \"" << CutName << "\"" << std::endl;
	return -999.00;
    }
}

int TextReader::GetNumberInt(st CutName, bool PrintError){
    if ((CVf_it = CVf.find(CutName.c_str())) != CVf.end()) {
	return int(CVf_it->second);
    }
    else if ((CVfa_it = CVfa.find(CutName.c_str())) != CVfa.end()) {
	if(PrintError) std::cout << "Warning - GetNumberInt(\"" << CutName << "\") : \"" << CutName << "\" is array" << std::endl;
	return int((CVfa_it->second).at(0));
    }
    else {
	if(PrintError) std::cout << "Error - GetNumberInt(\"" << CutName << "\") : Cannot find \"" << CutName << "\"" << std::endl;
	return -999;
    }
}

int TextReader::GetNumberInt(st CutName, int array_index, bool PrintError){
    if ((CVfa_it = CVfa.find(CutName.c_str())) != CVfa.end()) {
	if (array_index > 0 && (CVfa_it->second).size() >= array_index) {
	    return int((CVfa_it->second).at(array_index-1));
	}
	else {
	//    return int((CVfa_it->second).at(0));
	    if(PrintError) std::cout << "Error - GetNumberInt(\"" << CutName << "\"," << array_index << ") : Cannot find index\"" << array_index << "\"" << std::endl;
	    return -999;
	}

    }
    else if ((CVf_it = CVf.find(CutName.c_str())) != CVf.end()) {
	if(PrintError) std::cout << "Warning - GetNumberInt(\"" << CutName << "\"," << array_index << ") : \"" << CutName << "\" is not array" << std::endl;
	return int(CVf_it->second);
    }
    else {
	if(PrintError) std::cout << "Error - GetNumberInt(\"" << CutName << "\"," << array_index << ") : Cannot find \"" << CutName << "\"" << std::endl;
	return -999;
    }
}

unsigned int TextReader::GetNumberUint(st CutName, bool PrintError){
    if ((CVf_it = CVf.find(CutName.c_str())) != CVf.end()) {
	return (unsigned int)(CVf_it->second);
    }
    else if ((CVfa_it = CVfa.find(CutName.c_str())) != CVfa.end()) {
	if(PrintError) std::cout << "Warning - GetNumberUint(\"" << CutName << "\") : \"" << CutName << "\" is array" << std::endl;
	return (unsigned int)((CVfa_it->second).at(0));
    }
    else {
	if(PrintError) std::cout << "Error - GetNumberUint(\"" << CutName << "\") : Cannot find \"" << CutName << "\"" << std::endl;
	return 999;
    }
}

unsigned int TextReader::GetNumberUint(st CutName, int array_index, bool PrintError){
    if ((CVfa_it = CVfa.find(CutName.c_str())) != CVfa.end()) {
	if (array_index > 0 && (CVfa_it->second).size() >= array_index) {
	    return (unsigned int)((CVfa_it->second).at(array_index-1));
	}
	else {
	//    return int((CVfa_it->second).at(0));
	    if(PrintError) std::cout << "Error - GetNumberUint(\"" << CutName << "\"," << array_index << ") : Cannot find index\"" << array_index << "\"" << std::endl;
	    return 999;
	}

    }
    else if ((CVf_it = CVf.find(CutName.c_str())) != CVf.end()) {
	if(PrintError) std::cout << "Warning - GetNumberUint(\"" << CutName << "\"," << array_index << ") : \"" << CutName << "\" is not array" << std::endl;
	return (unsigned int)(CVf_it->second);
    }
    else {
	if(PrintError) std::cout << "Error - GetNumberUint(\"" << CutName << "\"," << array_index << ") : Cannot find \"" << CutName << "\"" << std::endl;
	return 999;
    }
}

float TextReader::GetNumberFloat(st CutName, bool PrintError){
    if ((CVf_it = CVf.find(CutName.c_str())) != CVf.end()) {
	return float(CVf_it->second);
    }
    else if ((CVfa_it = CVfa.find(CutName.c_str())) != CVfa.end()) {
	if(PrintError) std::cout << "Warning - GetNumberFloat(\"" << CutName << "\") : \"" << CutName << "\" is array" << std::endl;
	return float((CVfa_it->second).at(0));
    }
    else {
	if(PrintError) std::cout << "Error - GetNumberFloat(\"" << CutName << "\") : Cannot find \"" << CutName << "\"" << std::endl;
	return -999.00;
    }
}

float TextReader::GetNumberFloat(st CutName, int array_index, bool PrintError){
    if ((CVfa_it = CVfa.find(CutName.c_str())) != CVfa.end()) {
	if (array_index > 0 && (CVfa_it->second).size() >= array_index) {
	    return float((CVfa_it->second).at(array_index-1));
	}
	else {
	//    return int((CVfa_it->second).at(0));
	    if(PrintError) std::cout << "Error - GetNumberFloat(\"" << CutName << "\"," << array_index << ") : Cannot find index\"" << array_index << "\"" << std::endl;
	    return -999.00;
	}

    }
    else if ((CVf_it = CVf.find(CutName.c_str())) != CVf.end()) {
	if(PrintError) std::cout << "Warning - GetNumberFloat(\"" << CutName << "\"," << array_index << ") : \"" << CutName << "\" is not array" << std::endl;
	return float(CVf_it->second);
    }
    else {
	if(PrintError) std::cout << "Error - GetNumberFloat(\"" << CutName << "\"," << array_index << ") : Cannot find \"" << CutName << "\"" << std::endl;
	return -999.00;
    }
}

double TextReader::GetNumberDouble(st CutName, bool PrintError){
    if ((CVf_it = CVf.find(CutName.c_str())) != CVf.end()) {
	return double(CVf_it->second);
    }
    else if ((CVfa_it = CVfa.find(CutName.c_str())) != CVfa.end()) {
	if(PrintError) std::cout << "Warning - GetNumberDouble(\"" << CutName << "\") : \"" << CutName << "\" is array" << std::endl;
	return double((CVfa_it->second).at(0));
    }
    else {
	if(PrintError) std::cout << "Error - GetNumberDouble(\"" << CutName << "\") : Cannot find \"" << CutName << "\"" << std::endl;
	return -999.00;
    }
}

double TextReader::GetNumberDouble(st CutName, int array_index, bool PrintError){
    if ((CVfa_it = CVfa.find(CutName.c_str())) != CVfa.end()) {
	if (array_index > 0 && (CVfa_it->second).size() >= array_index) {
	    return double((CVfa_it->second).at(array_index-1));
	}
	else {
	//    return int((CVfa_it->second).at(0));
	    if(PrintError) std::cout << "Error - GetNumberDouble(\"" << CutName << "\"," << array_index << ") : Cannot find index\"" << array_index << "\"" << std::endl;
	    return -999.00;
	}

    }
    else if ((CVf_it = CVf.find(CutName.c_str())) != CVf.end()) {
	if(PrintError) std::cout << "Warning - GetNumberDouble(\"" << CutName << "\"," << array_index << ") : \"" << CutName << "\" is not array" << std::endl;
	return double(CVf_it->second);
    }
    else {
	if(PrintError) std::cout << "Error - GetNumberDouble(\"" << CutName << "\"," << array_index << ") : Cannot find \"" << CutName << "\"" << std::endl;
	return -999.00;
    }
}

st TextReader::GetText(st CutName, bool PrintError){
    if ((CVt_it = CVt.find(CutName.c_str())) != CVt.end()) {
	return CVt_it->second;
    }
    else if ((CVta_it = CVta.find(CutName.c_str())) != CVta.end()) {
	if(PrintError) std::cout << "Warning - GetText(\"" << CutName << "\") : \"" << CutName << "\" is array" << std::endl;
	return (CVta_it->second).at(0);
    }
    else {
	if(PrintError) std::cout << "Error - GetText(\"" << CutName << "\") : Cannot find \"" << CutName << "\"" << std::endl;
	return "DUMMY";
    }
}

st TextReader::GetText(st CutName, int array_index, bool PrintError){
    if ((CVta_it = CVta.find(CutName.c_str())) != CVta.end()) {
	if (array_index > 0 && (CVta_it->second).size() >= array_index) {
	    return (CVta_it->second).at(array_index-1);
	}
	else {
	//    return (CVta_it->second).at(0);
	    if(PrintError) std::cout << "Error - GetText(\"" << CutName << "\"," << array_index << ") : Cannot find index\"" << array_index << "\"" << std::endl;
	    return "DUMMY";
	}

    }
    else if ((CVt_it = CVt.find(CutName)) != CVt.end()) {
	if(PrintError) std::cout << "Warning - GetText(\"" << CutName << "\"," << array_index << ") : \"" << CutName << "\" is not array" << std::endl;
	return CVt_it->second;
    }
    else {
	if(PrintError) std::cout << "Error - GetText(\"" << CutName << "\"," << array_index << ") : Cannot find \"" << CutName << "\"" << std::endl;
	return "DUMMY";
    }

}

bool TextReader::GetBool(st CutName, bool PrintError){
    double BoolDouble;
    st BoolText;
    if(Check(CutName)) {
	if(CheckArray(CutName)) {
	    if(PrintError) std::cout << "Warning - GetBool(\"" << CutName << "\") : \"" << CutName << "\" is array" << std::endl;
	    if(CheckNumber(CutName)) {
		BoolDouble = double(((CVfa.find(CutName.c_str()))->second).at(0));
		if(BoolDouble == 1) {
		    return true;
		}
		else if(BoolDouble == 0) {
		    return false;
		}
		else {
		    if(PrintError) std::cout << "Error - GetBool(\"" << CutName << "\") : \"" << CutName << "\" is not bool type" << std::endl;
		    return false;
		}
	    }
	    else {
		BoolText = ((CVta.find(CutName.c_str()))->second).at(0);
		if (BoolText == "1" || BoolText == "true" || BoolText == "True" || BoolText == "TRUE") {
		    return true;
		}
		else if (BoolText == "0" || BoolText == "false" || BoolText == "False" || BoolText == "FALSE") {
		    return false;
		}
		else {
		    if(PrintError) std::cout << "Error - GetBool(\"" << CutName << "\") : \"" << CutName << "\" is not bool type" << std::endl;
		    return false;
		}
	    }
	}
	else {
	    if(CheckNumber(CutName)) {
		BoolDouble = double((CVf.find(CutName))->second);
		if(BoolDouble == 1) {
		    return true;
		}
		else if(BoolDouble == 0) {
		    return false;
		}
		else {
		    if(PrintError) std::cout << "Error - GetBool(\"" << CutName << "\") : \"" << CutName << "\" is not bool type" << std::endl;
		    return false;
		}
	    }
	    else {
		BoolText = (CVt.find(CutName))->second;
		if (BoolText == "1" || BoolText == "true" || BoolText == "True" || BoolText == "TRUE") {
		    return true;
		}
		else if (BoolText == "0" || BoolText == "false" || BoolText == "False" || BoolText == "FALSE") {
		    return false;
		}
		else {
		    if(PrintError) std::cout << "Error - GetBool(\"" << CutName << "\") : \"" << CutName << "\" is not bool type" << std::endl;
		    return false;
		}
	    }
	}
    }
    else {
	if(PrintError) std::cout << "Error - GetBool(\"" << CutName << "\") : Cannot find \"" << CutName << "\"" << std::endl;
	return false;
    }
}

bool TextReader::GetBool(st CutName, int array_index, bool PrintError){
    double BoolDouble;
    st BoolText;
    if(Check(CutName)) {
	if(CheckArray(CutName)) {
	    if(CheckNumber(CutName)) {
		BoolDouble = double(((CVfa.find(CutName.c_str()))->second).at(array_index-1));
		if(BoolDouble == 1) {
		    return true;
		}
		else if(BoolDouble == 0) {
		    return false;
		}
		else {
		    if(PrintError) std::cout << "Error - GetBool(\"" << CutName << "\"," << array_index << ") : \"" << array_index << "\"th \"" << CutName << "\" is not bool type" << std::endl;
		    return false;
		}
	    }
	    else {
		BoolText = ((CVta.find(CutName.c_str()))->second).at(array_index-1);
		if (BoolText == "1" || BoolText == "true" || BoolText == "True" || BoolText == "TRUE") {
		    return true;
		}
		else if (BoolText == "0" || BoolText == "false" || BoolText == "False" || BoolText == "FALSE") {
		    return false;
		}
		else {
		    if(PrintError) std::cout << "Error - GetBool(\"" << CutName << "\"," << array_index << ") : \"" << array_index << "\"th \"" << CutName << "\" is not bool type" << std::endl;
		    return false;
		}
	    }
	}
	else {
	    if(PrintError) std::cout << "Warning - GetBool(\"" << CutName << "\"," << array_index << ") : \"" << CutName << "\" is not array" << std::endl;
	    if(CheckNumber(CutName)) {
		BoolDouble = double((CVf.find(CutName))->second);
		if(BoolDouble == 1) {
		    return true;
		}
		else if(BoolDouble == 0) {
		    return false;
		}
		else {
		    if(PrintError) std::cout << "Error - GetBool(\"" << CutName << "\"," << array_index << ") : \"" << CutName << "\" is not bool type" << std::endl;
		    return false;
		}
	    }
	    else {
		BoolText = (CVt.find(CutName))->second;
		if (BoolText == "1" || BoolText == "true" || BoolText == "True" || BoolText == "TRUE") {
		    return true;
		}
		else if (BoolText == "0" || BoolText == "false" || BoolText == "False" || BoolText == "FALSE") {
		    return false;
		}
		else {
		    if(PrintError) std::cout << "Error - GetBool(\"" << CutName << "\"," << array_index << ") : \"" << CutName << "\" is not bool type" << std::endl;
		    return false;
		}
	    }
	}
    }
    else {
	if(PrintError) std::cout << "Error - GetBool(\"" << CutName << "\"," << array_index << ") : Cannot find \"" << CutName << "\"" << std::endl;
	return false;
    }
}

void TextReader::PrintoutVariables() {

    std::cout << std::endl;
    for (CVt_it = CVt.begin(); CVt_it != CVt.end(); ++CVt_it) {
	std::cout << std::setw(30) << CVt_it->first << std::setw(8) << " : " << '\t' << "\"" << CVt_it->second << "\"" << std::endl;
    }
    std::cout << std::endl;
    for (CVf_it = CVf.begin(); CVf_it != CVf.end(); ++CVf_it) {
	std::cout << std::setw(30) << CVf_it->first << std::setw(8) << " : " << '\t' << CVf_it->second << std::endl;
    }
    std::cout << std::endl;
    for (CVta_it = CVta.begin(); CVta_it != CVta.end(); ++CVta_it) {
	std::cout << std::setw(30) << CVta_it->first << " [" << (CVta_it->second).size() << "]" << std::setw(4) << " : " << '\t' << "{";
	for (vec_s_it s_it = (CVta_it->second).begin(); s_it != (CVta_it->second).end(); ++s_it) {
	    if(s_it != (CVta_it->second).begin()) std::cout << ", ";
	    std::cout << "\"" << *s_it << "\"";
	}
	std::cout << "}" << std::endl;
    }
    std::cout << std::endl;
    for (CVfa_it = CVfa.begin(); CVfa_it != CVfa.end(); ++CVfa_it) {
	std::cout << std::setw(30) << CVfa_it->first << " [" << (CVfa_it->second).size() <<  "]" << std::setw(4) << " : " << '\t' << "{";
	for (vec_d_it f_it = (CVfa_it->second).begin(); f_it != (CVfa_it->second).end(); ++f_it) {
	    if(f_it != (CVfa_it->second).begin()) std::cout << ", ";
	    std::cout << *f_it;
	}
	std::cout << "}" << std::endl;
    }
    std::cout << std::endl;

}

bool TextReader::Check(st CutName){
    bool Check_Cut = false;
    if (CVf.find(CutName.c_str())  != CVf.end())  { Check_Cut = true; } 
    if (CVfa.find(CutName.c_str()) != CVfa.end()) { Check_Cut = true; } 
    if (CVt.find(CutName.c_str())  != CVt.end())  { Check_Cut = true; } 
    if (CVta.find(CutName.c_str()) != CVta.end()) { Check_Cut = true; } 
    return Check_Cut;
}

bool TextReader::CheckNumber(st CutName){
    bool Check_Cut = false;
    if (CVf.find(CutName.c_str())  != CVf.end())  { Check_Cut = true; } 
    if (CVfa.find(CutName.c_str()) != CVfa.end()) { Check_Cut = true; } 
    return Check_Cut;
}

bool TextReader::CheckText(st CutName){
    bool Check_Cut = false;
    if (CVt.find(CutName.c_str())  != CVt.end())  { Check_Cut = true; } 
    if (CVta.find(CutName.c_str()) != CVta.end()) { Check_Cut = true; } 
    return Check_Cut;
}

bool TextReader::CheckArray(st CutName){
    bool Check_Cut = false;
    if (CVfa.find(CutName.c_str()) != CVfa.end()) { Check_Cut = true; } 
    if (CVta.find(CutName.c_str()) != CVta.end()) { Check_Cut = true; } 
    return Check_Cut;
}

int TextReader::Size(st CutName){
    int Variable_size = 0;
    if (CVf.find(CutName.c_str())  != CVf.end())  { Variable_size = 1; }
    if (CVt.find(CutName.c_str())  != CVt.end())  { Variable_size = 1; }
    if ((CVfa_it = CVfa.find(CutName.c_str())) != CVfa.end()) { Variable_size = (CVfa_it->second).size(); }
    if ((CVta_it = CVta.find(CutName.c_str())) != CVta.end()) { Variable_size = (CVta_it->second).size(); }
    return Variable_size;
}

map_d TextReader::ReturnMap(map_d MapType){
    return CVf;
}

map_s TextReader::ReturnMap(map_s MapType){
    return CVt;
}

map_sd TextReader::ReturnMap(map_sd MapType){
    return CVfa;
}

map_ss TextReader::ReturnMap(map_ss MapType){
    return CVta;
}

