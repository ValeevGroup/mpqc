
#include <util/keyval/keyval.h>
#include "oogl.h"
#include "object.h"
#include "sphere.h"

int
main()
{
  RefKeyVal keyval = new ParsedKeyVal(SRCDIR "/rentest.in");
  printf("getting render\n"); fflush(stdout);
  RefRender render = keyval->describedclassvalue("render");
  printf("getting object\n"); fflush(stdout);
  RefRenderedObject object = keyval->describedclassvalue("object");

  printf("rendering object\n"); fflush(stdout);
  render->render(object);
  printf("rendered object\n"); fflush(stdout);

  printf("getting rid of keyval\n"); fflush(stdout);
  keyval = 0;
  printf("getting rid of render\n"); fflush(stdout);
  render = 0;
  printf("getting rid of object\n"); fflush(stdout);
  object = 0;

  printf("main done\n"); fflush(stdout);

  return 0;
}
