#ifdef __GNUG__
#pragma implementation
#endif

#ifndef _chemistry_qc_psi_input_timpl_h
#define _chemistry_qc_psi_input_timpl_h

using namespace std;
using namespace sc;

#include<chemistry/qc/psi/psiinput.h>

template <typename T>
void
PsiInput::write_keyword_array(const char *keyword, int num, const std::vector<T>& values)
{
  write_indent();
  file_ << scprintf("%s = (", keyword);
  file_ << setw(24) << setprecision(12);
  for (int i=0; i<num; i++) {
    file_ << values[i] << " ";
  }
  file_ << ")" << endl;
}

#endif /* header guard */
