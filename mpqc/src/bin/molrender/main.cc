
#include <stdio.h>
#include <util/keyval/keyval.h>
#include <util/render/render.h>
#include <util/render/oogl.h>
#include <util/render/object.h>
#include <math/isosurf/surf.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/molecule/molrender.h>
#include <chemistry/molecule/molshape.h>

int
main(int argc, char** argv)
{
  int i;
  const char* model = "stick";
  const char* keyword = "molecule";
  const char* inputfile = "molrender.in";
  int level = 3;
  enum InputType { PDB, KEYVAL } input = KEYVAL;
  const char* render = 0;
  int quiet = 0;

  for (i=1; i<argc; i++) {
      if (!strcmp(argv[i], "-model")) {
          i++;
          model = argv[i];
        }
      else if (!strcmp(argv[i], "-quiet")) {
          quiet = 1;
        }
      else if (!strcmp(argv[i], "-render")) {
          i++;
          render = argv[i];
        }
      else if (!strcmp(argv[i], "-keyword")) {
          i++;
          keyword = argv[i];
        }
      else if (!strcmp(argv[i], "-pdb")) {
          i++;
          input = PDB;
          inputfile = argv[i];
        }
      else if (!strcmp(argv[i], "-keyval")) {
          i++;
          input = KEYVAL;
          inputfile = argv[i];
        }
      else if (!strcmp(argv[i], "-level")) {
          i++;
          level = atoi(argv[i]);
        }
      else {
          fprintf(stderr,"%s: unknown option: \"%s\"\n", argv[0], argv[i]);
          abort();
        }
    }

  // Find the molecule.
  RefMolecule mol;
  if (input == PDB) {
      RefAssignedKeyVal keyval = new AssignedKeyVal();
      keyval->assign("pdb_file", inputfile);
      mol = new Molecule(keyval.pointer());
    }
  else {
      RefKeyVal keyval = new ParsedKeyVal(inputfile);
      mol = new Molecule(new PrefixKeyVal(keyword, keyval));
    }

  // Set up the rendered molecule object.
  RefAssignedKeyVal tmpkv = new AssignedKeyVal();
  RefAssignedKeyVal keyval = new AssignedKeyVal();

  keyval->assign("molecule", mol.pointer());
  keyval->assign("model", model);

  RefDescribedClass atominfo = new AtomInfo(tmpkv.pointer());
  keyval->assign("atominfo", atominfo);
  tmpkv->clear();

  RefRenderedObject molobject;
  if (!strcmp(model,"stick")) {
      molobject = new RenderedStickMolecule(keyval.pointer());
    }
  else if (!strcmp(model,"ball")) {
      molobject = new RenderedBallMolecule(keyval.pointer());
    }
  else if (!strcmp(model,"connolly")) {
      tmpkv->assign("molecule", mol);
      tmpkv->assign("atominfo", atominfo);
      RefDescribedClass volume = new ConnollyShape(tmpkv.pointer());
      tmpkv->clear();
      tmpkv->assignboolean("verbose", !quiet);
      RefDescribedClass trisurf = new TriangulatedSurface(tmpkv.pointer());
      tmpkv->clear();
      tmpkv->assign("surface", trisurf);
      tmpkv->assign("volume", volume);
      tmpkv->assign("resolution", 1.0);
      tmpkv->assignboolean("remove_short_edges", 0);
      tmpkv->assignboolean("remove_slender_edges", 0);
      RefDescribedClass surface
          = new TriangulatedImplicitSurface(tmpkv.pointer());
      tmpkv->clear();
      keyval->assign("surface", surface);
      molobject = new RenderedMolecularSurface(keyval.pointer());
    }
  else {
      fprintf(stderr,"%s: unknown model \"%s\"\n", argv[0], model);
      abort();
    }

  RefRenderedObjectSet object;

  if (render) {
      object = new RenderedObjectSet;
      object->add(molobject);
      RefAppearance appearance = new Appearance;
      appearance->level().set(level);
      object->appearance(appearance);
      if (object.null()) {
          fprintf(stderr,"%s: got a null object to render\n",argv[0]);
          abort();
        }

      // Set up the renderer.
      RefRender renderer;
      if (!strcmp("oogl", render)) {
          renderer = new OOGLRender;
        }
      else {
          fprintf(stderr,"%s: unknown renderer: \"%s\"\n", argv[0], render);
        }

      // Render the object.
      renderer->render(object.pointer());
    }

  if (!quiet) {
      ConnollyShape::print_counts();
      CS2Sphere::print_counts();
    }

  fflush(stdout);
  fflush(stderr);
  return 0;
}
