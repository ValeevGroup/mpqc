
#include "scostream.h"

main()
{
  SCostream& os = SCostream::cout;

  {
    os++;
    os.indent() << "|" << endl;
    {
      os++;
      os.indent() << "|"; os << os.get_column() << endl;
      os.indent() << "|" << os.get_column() << endl;
      os--;
    }
    os.indent() << "|" << endl;
    
    os--;
  }

  os << "x" << endl;
  os << 'y' << endl;

  
}
