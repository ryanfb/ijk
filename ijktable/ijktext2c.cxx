/// \file ijktext2c.cxx
/// convert text files to C #define and output to standard output

#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

char * outfilename = NULL;
vector<string> infilename;

// routines
void convert_text_to_c_define(ostream & out);
void parse_command_line(int argc, char **argv);
void usage_error();

int main(int argc, char **argv)
{
  parse_command_line(argc, argv);

  if (outfilename == NULL) {
    convert_text_to_c_define(cout);
  }
  else {
    ofstream outfile(outfilename, ios::out);

    if (!outfile) {
      cerr << "Error opening file " << outfilename << "." << endl;
      cerr << "Exiting." << endl;
      exit(15);
    }

    convert_text_to_c_define(outfile);
    outfile.close();
  }
}

void convert_text_to_c_define(ostream & out)
{
  for (int k = 0; k < infilename.size(); k++) {
    ifstream infile(infilename[k].c_str(), ios::in);
    if (!infile) {
      cerr << "Error opening file " << infilename[k] << "." << endl;
      cerr << "Exiting." << endl;
      exit(10);
    }

    string def_name = infilename[k];
    for (int i = 0; i < def_name.length(); i++) {
      def_name[i] = toupper(def_name[i]);
      if (!isalnum(def_name[i])) {
	def_name[i] = '_';
      }
    }

    // add an '_' prefix
    def_name = string("_") + def_name;

    out << "#define " << def_name << " = \\" << endl;

    char c;
    while (infile.get(c)) {
      if (c == '\"') {
	out << "\\\"";
      }
      else if (c == '\n') {
	out << "\\" << endl;
      }
      else {
	out << c;
      }
    }

    out << endl << endl;
    infile.close();
  }

}

void parse_command_line(int argc, char **argv)
{
  int iarg = 1;
  while (iarg < argc) {
    if (strcmp(argv[iarg], "-o") == 0) {
      iarg++;
      if (iarg >= argc)
	usage_error();
      outfilename = argv[iarg];
    }
    else {
      break;
    }
    iarg++;
  };

  while (iarg < argc) {
    infilename.push_back(argv[iarg]);
    iarg++;
  }

  if (infilename.size() < 1) { usage_error(); }
}

void usage_error()
{
  cerr << "Usage: ijktext2c {-o output-file} file1 {file2} {file3} ... "
       << endl;
  exit(40);
}
