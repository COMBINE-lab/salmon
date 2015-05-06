/* * * * *
 *  cmdline.h
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Copyright (C) 2012, Rongkai Xia, shinekaixia@gmail.com
 *
 * * * * */

#ifndef libTM_cmdline_h
#define libTM_cmdline_h
#include <string>
#include <vector>
#include <map>

static int param_offset=25;
static int line_max_char=80;
class cmdline {
    
public:
    
    bool parse_cmdline(int argc,const char * argv[]);
    
    void print_help();
    
    std::string get_value(const std::string& param);
    
    bool check_param(const std::string *param,int len);
    
    bool has_param(const std::string& param);
    
    void register_help(const std::string& param,const std::string& help);
    
    void register_help(const std::string *param,const std::string *help,int len);
    
    void clear_help();
protected:
    bool parse_name(std::string& name);
    std::map<std::string,std::string> value;
    
    std::vector<std::pair<std::string, std::string> > help;
    
};

bool cmdline::has_param(const std::string &param){
    if (value.find(param)!=value.end()) {
        return true;
    }
    return false;
}
std::string cmdline::get_value(const std::string &param){
    if (value.find(param)==value.end()) {
        return "";
    }
    return value[param];
}

bool cmdline::parse_name(std::string& name){
    if (name.length()>=1 && name[0]!='-') {
        return false;
    }
    if (name.length()>=2 && name[1]=='-') {
        name=name.substr(2);
    }else {
        name=name.substr(1);
    }
    return true;
}

bool cmdline::parse_cmdline(int argc,const char * argv[]){
    int i = 1;
    while (i < argc) {
        std::string s(argv[i]);
        if (parse_name(s)) {
            if (value.find(s) != value.end()) {
                fprintf(stderr, "Error: The parameter -%s is already specified.\n",s.c_str());
                return false;
            }
            if ((i+1) < argc) {
                std::string s_next(argv[i+1]);
                if (! parse_name(s_next)) {
                    value[s] = s_next;
                    i++;
                } else {
                    value[s] = "";
                }
            } else {
                value[s] = "";
            }
        } else {
            fprintf(stderr, "Error: Unknow parameter %s.\n",s.c_str());
            print_help();
            return false;
        }
        i++;
    }
    return true;
}

void cmdline::print_help() {
    std::vector< std::pair<std::string, std::string> >::const_iterator pv;
    
    for ( pv = help.begin(); pv != help.end(); ++pv) {
        std::cout << "-" << pv->first;
        for (int i=pv->first.size()+1; i < param_offset; i++) { std::cout << " "; } 
        std::string s_out = pv->second;
        while (s_out.size() > 0) {
            if (s_out.size() > (line_max_char-param_offset)) {
                size_t p = s_out.substr(0, line_max_char-param_offset).find_last_of(" \t");
                if (p == 0) {
                    p = line_max_char-param_offset;
                }
                std::cout << s_out.substr(0, p) << std::endl;
                s_out = s_out.substr(p+1, s_out.length()-p);            
            } else {
                std::cout << s_out << std::endl;
                s_out = "";  
            }
            if (s_out.size() > 0) {
                for (int i=0; i < param_offset; i++) { std::cout << " "; }
            }
        }
    }
}

bool cmdline::check_param(const std::string *param, int len){
    int i;
    for (std::map<std::string,std::string>::const_iterator iter=value.begin(); iter!=value.end(); iter++) {
        for (i=0; i<len; i++) {
            if (iter->first==param[i]) {
                break;
            }
        }
        if (i==len) {
            fprintf(stderr, "Error: Unknow parameter %s.\n",iter->first.c_str());
//            print_help();
            return false;
        }
    }
    return true;
}

void cmdline::register_help(const std::string &param, const std::string &help){
    this->help.push_back(make_pair(param, help));
}

void cmdline::register_help(const std::string *param,const std::string *help,int len){
    clear_help();
    for (int i=0; i<len; i++) {
        register_help(param[i], help[i]);
    }
}

void cmdline::clear_help(){
    this->help.clear();
}
#endif
