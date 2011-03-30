#ifndef _chemistry_qc_psi_input_timpl_h
#define _chemistry_qc_psi_input_timpl_h

#include<chemistry/qc/psi/psiinput.h>

template <typename T>
void
sc::PsiInput::write_keyword_array(const char *keyword, const std::vector<T>& values)
{
  if (!can_run_on_me()) return;
  write_indent();
  file_ << scprintf("%s = (", keyword);
  file_ << std::setw(24) << std::setprecision(12);
  const size_t size = values.size();
  for (int i=0; i<size; i++) {
    file_ << values[i] << " ";
  }
  file_ << ")" << std::endl;
}

#endif /* header guard */
