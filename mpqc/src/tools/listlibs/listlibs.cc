
#include <stdio.h>
#include <stdlib.h>
#include <iostream.h>
#include <fstream.h>
#include <string>
#include <list>
#include <map>
#include <set>

static const bool debug = false;

bool
include_to_filename(string &filename)
{
  if (filename.substr(0,8) == "#include") {
    filename.remove(0,filename.find("<")+1);
    filename.remove(filename.find(">"),
                    filename.size() - filename.find(">") + 1);
    return true;
    }
  return false;
}

void
process_file(const string &passedfilename,
             map<string,list<string>,less<string> > &read_files,
             const list<string> &includes)
{
  string filename(passedfilename);

  if (debug) cout << "process_file: filename: " << filename << endl;

  // find and open the file
  string ifile;
  ifstream file;
  for (list<string>::const_iterator i = includes.begin();
       i != includes.end();
       i++) {
    ifile = filename;
    if (filename[0] != '/') {
      ifile = string(*i) + "/" + ifile;
      }
    if (debug) cout << "process_file: trying: " << ifile << endl;
    file.open(ifile.c_str());
    if (file.good()) break;
    }

  if (!file.good()) {
    cerr << "listlibs: couldn't find " << filename << endl;
    exit(1);
    }

  if (debug) cout << "process_file: ifile: " << ifile << endl;

  // read the contents of the file and insert into the file contents map
  list<string> filecontents;
  const int max_line_length = 100;
  char c_line[max_line_length];
  while (file.good()) {
    file.getline(c_line, max_line_length);
    string line(c_line);
    if (line == "" ) continue;
    filecontents.push_back(line);
    if (debug) cout << "process_file: contents: " << line << endl;
    }
  file.close();
  read_files[filename] = filecontents;

  // read in any other referenced files that have not yet been read
  for (list<string>::iterator i = filecontents.begin();
       i != filecontents.end();
       i++) {
    string line(*i);
    if (include_to_filename(line)
        && read_files.find(line) == read_files.end()) {
      process_file(line, read_files, includes);
      }
    }
}

void
find_libraries(list<string> &libraries,
               const string &passedfilename,
               const map<string,list<string>,less<string> > &read_files,
               set<string,less<string> > &current_includes)
{
  if (current_includes.find(passedfilename) != current_includes.end()) {
    cerr << "listlibs: recursive include detected for "
         << passedfilename << ":" << endl;
    for (set<string,less<string> >::iterator i = current_includes.begin();
         i != current_includes.end();
         i++) {
      cerr << "  " << (*i) << endl;
      }
    exit(1);
    }
#ifndef sgi
  current_includes.insert(passedfilename);
#endif
  string filename(passedfilename);
  const list<string> &filecontents = (*(read_files.find(filename))).second;
  for (list<string>::const_iterator i = filecontents.begin();
       i != filecontents.end();
       i++) {
    string line(*i);
    if (include_to_filename(line)) {
      find_libraries(libraries, line, read_files, current_includes);
      }
    else {
      if (debug) cout << "listlibs: find_libraries: found: " << line << endl;
      libraries.push_back(line);
      }
    }
#ifndef sgi
  current_includes.erase(passedfilename);
#endif
}

void
eliminate_redundant(list<string> &libraries)
{
  list<string> nonredund_libraries;
  set<string,less<string> > known_libs;
  for (list<string>::reverse_iterator inext, i = libraries.rbegin();
       i != libraries.rend();
       i++) {
    if (known_libs.find(*i) == known_libs.end()) {
      known_libs.insert(*i);
      nonredund_libraries.push_front(*i);
      }
    }
  libraries = nonredund_libraries;
}

void
substibute_defines(list<string> &libraries,
                   map<string,string,less<string> > &defines)
{
  for (map<string,string,less<string> >::iterator i = defines.begin();
       i != defines.end();
       i++) {
    const string &macro = (*i).first;
    const string &value = (*i).second;
    int macrosize = macro.size();
    for (list<string>::iterator j = libraries.begin();
         j != libraries.end();
         j++) {
      string &lib = (*j);
      for (int npass = 0, index = lib.find(macro);
           index != lib.npos;
           npass++,index = lib.find(macro)) {
        lib.replace(index,macrosize,value);
        if (npass > 20) {
          cerr << "lib = " << lib << endl;
          cerr << "listlibs: recursive macro replacement "
               << macro << "=" << value
               << endl;
          exit(1);
          }
        }
      }
    }
}

int
main(int argc, char *argv[])
{
  list<string> libraries;
  list<string> includes;
  map<string,string,less<string> > defines;
  map<string,list<string>,less<string> > read_files;
  string filename;
  includes.push_back(".");

  // process the arguments
  for (int i=1; i<argc; i++) {
    string arg(argv[i]);
    if (arg.substr(0,2) == "-D") {
      string def(&argv[i][2]);
      string symbol(def);
      string value(def);
      symbol.remove(symbol.find("="),symbol.size()-symbol.find("=")+1);
      def.remove(0, def.find("=") + 1);
      if (debug) cout << "Defining " << symbol << " to be " << def << endl;
      defines[symbol] = def;
      }
    else if (arg.substr(0,2) == "-I") {
      string incdir(arg.substr(2,arg.size()-2));
      if (debug) cout << "main: incdir: " << incdir << endl;
      includes.push_back(incdir);
      }
    else {
      filename = arg;
      }
    }

  if (filename.size() == 0) {
    cerr << "listlibs: requires a filename" << endl;
    exit(1);
    }

  process_file(filename, read_files, includes);

  set<string,less<string> > current_includes;
  find_libraries(libraries, filename, read_files, current_includes);

  eliminate_redundant(libraries);

  substibute_defines(libraries, defines);

  for (list<string>::iterator i = libraries.begin();
       i != libraries.end();
       i++) {
    if (i != libraries.begin()) cout << " ";
    cout << (*i);
    }
  cout << endl;

  return 0;
}
