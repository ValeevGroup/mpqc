
#ifdef __GNUC__
#pragma implementation
#endif

#include <math/scmat/matrix.h>
#include <math/isosurf/surf.h>

/////////////////////////////////////////////////////////////////////////
// TriangulatedSurface

TriangulatedSurface::TriangulatedSurface():
  _triangle_vertex(0),
  _triangle_edge(0),
  _edge_vertex(0),
  _integrator(new GaussTriangleIntegrator(1))
{
  clear();
}

TriangulatedSurface::~TriangulatedSurface()
{
  clear();
  delete _integrator;
}

void
TriangulatedSurface::set_integrator(TriangleIntegrator*i)
{
  delete _integrator;
  _integrator = i;
}

TriangleIntegrator*
TriangulatedSurface::integrator(int itri)
{
  // currently itri is ignored
  return _integrator;
}

void
TriangulatedSurface::clear()
{
  _completed_surface = 0;

  if (_triangle_vertex) {
      for (int i=0; i<_triangles.length(); i++) {
          delete[] _triangle_vertex[i];
        }
      delete[] _triangle_vertex;
    }
  _triangle_vertex = 0;

  if (_triangle_edge) {
      for (int i=0; i<_triangles.length(); i++) {
          delete[] _triangle_edge[i];
        }
      delete[] _triangle_edge;
    }
  _triangle_edge = 0;

  if (_edge_vertex) {
      for (int i=0; i<_edges.length(); i++) {
          delete[] _edge_vertex[i];
        }
      delete[] _edge_vertex;
    }
  _edge_vertex = 0;

  _have_values = 0;
  _values.clear();

  _vertices.clear();
  _edges.clear();
  _triangles.clear();
}

void
TriangulatedSurface::remove_short_edges(double length_cutoff)
{
  if (!_completed_surface) {
      complete_surface();
    }
  
  int nvertex = _vertices.length();
  int nedge = _edges.length();
  int ntriangle = _triangles.length();

  printf("TriangulatedSurface::remove_short_edges:\n");
  printf("  nvertex nedge ntriangle\n");
  printf("    %3d    %3d     %3d     (%s)\n",
         nvertex,nedge,ntriangle,"initial");

  int* vertex_map = new int[nvertex];
  int* edge_map = new int[nedge];
  int* triangle_map = new int[ntriangle];
  
  for (int ie=0; ie<nedge; ie++) {
      // the edge with potentially two deleted points
      RefEdge E2D(_edges[ie]);
      if (E2D->length() < length_cutoff) {
          int i;
          
          // find the deleted vertices
          ArraysetRefVertex D;
          D.add(_edges[ie]->vertex(0));
          D.add(_edges[ie]->vertex(1));

          // the edges with one deleted point which border triangles with
          // one deleted point.
          ArraysetRefEdge E1D_T1D;
          // the edges with one deleted point which border the triangles
          // with two deleted points.
          ArraysetRefEdge E1D_T2DA;
          ArraysetRefEdge E1D_T2DB;

          // the triangles with two deleted points
          RefTriangle T2DA;
          RefTriangle T2DB;
          // the triangles with one deleted point
          ArraysetRefTriangle T1D;

          // find the triangles T2DA, T2DB, and the triangle set T1D
          // also find the edge arrays E1D_T2DA and E1D_T2DB
          for (i=0; i<ntriangle; i++) {
              if (_triangles[i]->edge(0) == E2D) {
                  if (T2DA.null()) {
                      T2DA = _triangles[i];
                      E1D_T2DA.add(_triangles[i]->edge(1));
                      E1D_T2DA.add(_triangles[i]->edge(2));
                    }
                  else {
                      T2DB = _triangles[i];
                      E1D_T2DB.add(_triangles[i]->edge(1));
                      E1D_T2DB.add(_triangles[i]->edge(2));
                    }
                }
              else if (_triangles[i]->edge(1) == E2D) {
                  if (T2DA.null()) {
                      T2DA = _triangles[i];
                      E1D_T2DA.add(_triangles[i]->edge(0));
                      E1D_T2DA.add(_triangles[i]->edge(2));
                    }
                  else {
                      T2DB = _triangles[i];
                      E1D_T2DB.add(_triangles[i]->edge(0));
                      E1D_T2DB.add(_triangles[i]->edge(2));
                    }
                }
              else if (_triangles[i]->edge(2) == E2D) {
                  if (T2DA.null()) {
                      T2DA = _triangles[i];
                      E1D_T2DA.add(_triangles[i]->edge(0));
                      E1D_T2DA.add(_triangles[i]->edge(1));
                    }
                  else {
                      T2DB = _triangles[i];
                      E1D_T2DB.add(_triangles[i]->edge(0));
                      E1D_T2DB.add(_triangles[i]->edge(1));
                    }
                }
              else if (_triangles[i]->vertex(0) == D[0]
                       || _triangles[i]->vertex(0) == D[1]
                       || _triangles[i]->vertex(1) == D[0]
                       || _triangles[i]->vertex(1) == D[1]
                       || _triangles[i]->vertex(2) == D[0]
                       || _triangles[i]->vertex(2) == D[1]) {
                  T1D.add(_triangles[i]);
                }
            }

          // find the edge array E1D_T1D
          for (i=0; i<T1D.length(); i++) {
              for (int j=0; j<3; j++) {
                  if (T1D[i]->edge(j)->vertex(0) == D[0]
                      ||T1D[i]->edge(j)->vertex(0) == D[1]
                      ||T1D[i]->edge(j)->vertex(1) == D[0]
                      ||T1D[i]->edge(j)->vertex(1) == D[1]) {
                      E1D_T1D.add(T1D[i]->edge(j));
                    }
                }
            }

          // find the vertex replacement, D_rep
          RefSCVector newp(E2D->vertex(0)->dimension());
          RefSCVector newv(E2D->vertex(0)->dimension());
          RefVertex D_rep(new Vertex(newp,newv));
          for (i=0; i<3; i++) {
              D_rep->point()[i] =
                0.5*(D[0]->point()[i]
                     + D[1]->point()[i]);
              D_rep->gradient().operator[](i) =
                0.5*(D[0]->gradient().operator[](i)
                     + D[1]->gradient().operator[](i));
            }

          // find the edge replacements, E1D_T2DA_rep
          RefEdge E1D_T2DA_rep;
          if (D[0] != E1D_T2DA[0]->vertex(0)
              && D[1] != E1D_T2DA[0]->vertex(0)) {
              E1D_T2DA_rep = new Edge(E1D_T2DA[0]->vertex(0),D_rep);
            }
          else {
              E1D_T2DA_rep = new Edge(E1D_T2DA[0]->vertex(1),D_rep);
            }

          // find the edge replacements, E1D_T2DB_rep
          RefEdge E1D_T2DB_rep;
          if (D[0] != E1D_T2DB[0]->vertex(0)
              && D[1] != E1D_T2DB[0]->vertex(0)) {
              E1D_T2DB_rep = new Edge(E1D_T2DB[0]->vertex(0),D_rep);
            }
          else {
              E1D_T2DB_rep = new Edge(E1D_T2DB[0]->vertex(1),D_rep);
            }

          // find the edge replacements, E1D_T1D_rep
          ArrayRefEdge E1D_T1D_rep(E1D_T1D.length());
          for (i=0; i<E1D_T1D.length(); i++) {
              if (E1D_T1D[i]->vertex(0) == D[0]
                  || E1D_T1D[i]->vertex(0) == D[1]) {
                  E1D_T1D_rep[i] = new Edge(D_rep,E1D_T1D[i]->vertex(1));
                }
              else {
                  E1D_T1D_rep[i] = new Edge(D_rep,E1D_T1D[i]->vertex(0));
                }
            }

          // find the triangle replacments, T1D_rep
          ArrayRefTriangle T1D_rep(T1D.length());
          for (i=0; i<T1D.length(); i++) {
              ArrayRefEdge edges(3);
              for (int j=0; j<3; j++) {
                  edges[j] = T1D[i]->edge(j);
                  int iedge = E1D_T1D.iseek(edges[j]);
                  if (iedge >= 0) {
                      edges[j] = E1D_T1D_rep[iedge];
                      continue;
                    }
                  iedge = E1D_T2DA.iseek(edges[j]);
                  if (iedge >= 0) {
                      edges[j] = E1D_T2DA_rep;
                      continue;
                    }
                  iedge = E1D_T2DB.iseek(edges[j]);
                  if (iedge >= 0) {
                      edges[j] = E1D_T2DB_rep;
                      continue;
                    }
                }
              T1D_rep[i] = new Triangle(edges[0],edges[1],edges[2]);
            }

          for (i=0; i<ntriangle; i++) triangle_map[i] = i;
          for (i=0; i<nvertex; i++) vertex_map[i] = i;
          for (i=0; i<nedge; i++) edge_map[i] = i;

          // compute the mapping of old vertex indices to new vertex indices
          for (i=0; i<2; i++) {
              int j;
              int jstart = _vertices.iseek(D[i]);
              vertex_map[jstart] = nvertex-2; // delete 2 add 1
              for (j=jstart+1; j<nvertex; j++) {
                  vertex_map[j]--;
                }
            }

          // compute the mapping of old edge indices to new edge indices
          int istart = _edges.iseek(E2D);
          edge_map[istart] = -1;
          for (i=istart+1; i<nedge; i++) edge_map[i]--;
          for (i=0; i<2; i++) {
              int j;
              int jstart = _edges.iseek(E1D_T2DA[i]);
              edge_map[jstart] = nedge - 5; // delete 5 add 1
              for (j=jstart+1; j<nedge; j++) edge_map[j]--;
              jstart = _edges.iseek(E1D_T2DB[i]);
              edge_map[jstart] = nedge - 4; // delete 5 add 2
              for (j=jstart+1; j<nedge; j++) edge_map[j]--;
            }

          // compute the mapping of old triangle indices to new
          istart = _triangles.iseek(T2DA);
          triangle_map[istart] = -1;
          for (i=istart+1; i<ntriangle; i++) triangle_map[i]--;
          istart = _triangles.iseek(T2DB);
          triangle_map[istart] = -1;
          for (i=istart+1; i<ntriangle; i++) triangle_map[i]--;
          
          // replace the E1D_T1D edges
          for (i=0; i<E1D_T1D.length(); i++) {
              _edges[_edges.iseek(E1D_T1D[i])] = E1D_T1D_rep[i];
            }

          // replace the T1D triangles
          for (i=0; i<T1D.length(); i++) {
              _triangles[_triangles.iseek(T1D[i])] = T1D_rep[i];
            }

          // delete the unused objects
          _triangles.del(T2DA);
          _triangles.del(T2DB);
          _edges.del(E2D);
          _vertices.del(D[0]);
          _vertices.del(D[1]);
          for (i=0; i<2; i++) {
              _edges.del(E1D_T2DA[i]);
              _edges.del(E1D_T2DB[i]);
            }

          // put the new vertex on the end
          _vertices.add(D_rep);

          // put the new edges on the end
          _edges.add(E1D_T2DA_rep);
          _edges.add(E1D_T2DB_rep);

          if (_triangle_vertex) {
              int ii;
              for (i=ii=0; ii<ntriangle; ii++) {
                  if (triangle_map[ii] < 0) {
                      delete[] _triangle_vertex[ii];
                      continue;
                    }
                  _triangle_vertex[i] = _triangle_vertex[ii];
                  for (int j=0; j<3; j++) {
                      _triangle_vertex[i][j]
                        = vertex_map[_triangle_vertex[i][j]];
                    }
                  i++;
                }
            }
          if (_triangle_edge) {
              int ii;
              for (i=ii=0; ii<ntriangle; ii++) {
                  if (triangle_map[ii] < 0) {
                      delete[] _triangle_edge[ii];
                      continue;
                    }
                  _triangle_edge[i] = _triangle_edge[ii];
                  for (int j=0; j<3; j++) {
                      _triangle_edge[i][j]
                        = edge_map[_triangle_edge[i][j]];
                      if (_triangle_edge[i][j] < 0
                          || _triangle_edge[i][j] > 1000) {
                          printf("funny edge\n");
                          sqrt(-1.0);
                          abort();
                        }
                    }
                  i++;
                }
            }
          if (_edge_vertex) {
              int ii;
              for (i=ii=0; ii<nedge; ii++) {
                  if (edge_map[ii] < 0) {
                      delete[] _edge_vertex[ii];
                      continue;
                    }
                  _edge_vertex[i] = _edge_vertex[ii];
                  for (int j=0; j<2; j++) {
                      _edge_vertex[i][j]
                        = vertex_map[_edge_vertex[i][j]];
                    }
                  i++;
                }
            }
          ntriangle -= 2;
          nedge -= 3;
          nvertex -= 1;
          
        } /* if below cutoff */
    }
  delete[] edge_map;
  delete[] vertex_map;

  printf("    %3d    %3d     %3d     (%s)\n",
         nvertex,nedge,ntriangle,"final");

  printf("    %3d    %3d     %3d     (%s)\n",
         _vertices.length(),_edges.length(),_triangles.length(),"actual");
}

void
TriangulatedSurface::complete_surface()
{
  int i;
  int ntri = ntriangle();
  for (i=0; i<ntri; i++) {
      RefTriangle tri = triangle(i);
      _edges.add(tri->edge(0));
      _edges.add(tri->edge(1));
      _edges.add(tri->edge(2));
    }
  int ne = nedge();
  for (i=0; i<ne; i++) {
      RefEdge e = edge(i);
      _vertices.add(e->vertex(0));
      _vertices.add(e->vertex(1));
    }

  // construct the array that converts the triangle number and vertex
  // number within the triangle to the overall vertex number
  _triangle_vertex = new int*[ntri];
  for (i=0; i<ntri; i++) {
      _triangle_vertex[i] = new int[3];
      for (int j=0; j<3; j++) _triangle_vertex[i][j] =
                                 _vertices.iseek(triangle(i)->vertex(j));
    }

  // construct the array that converts the triangle number and edge number
  // within the triangle to the overall edge number
  _triangle_edge = new int*[ntri];
  for (i=0; i<ntri; i++) {
      _triangle_edge[i] = new int[3];
      for (int j=0; j<3; j++) _triangle_edge[i][j] =
        _edges.iseek(triangle(i)->edge(j));
    }

  // construct the array that converts the edge number and vertex number
  // within the edge to the overall vertex number
  _edge_vertex = new int*[ne];
  for (i=0; i<ne; i++) {
      _edge_vertex[i] = new int[2];
      for (int j=0; j<2; j++)
        _edge_vertex[i][j] = _vertices.iseek(_edges[i]->vertex(j));
    }

  _completed_surface = 1;
}

void
TriangulatedSurface::compute_values(RefVolume&vol)
{
  int n = _vertices.length();
  _values.set_length(n);

  for (int i=0; i<n; i++) {
      vol->set_x(_vertices[i]->point());
      _values[i] = vol->value();
    }
  _have_values = 1;
}

double
TriangulatedSurface::area()
{
  double result = 0.0;
  for (int i=0; i<_triangles.length(); i++) {
      result += _triangles[i]->area();
    }
  return result;
}

double
TriangulatedSurface::volume()
{
  double result = 0.0;
  for (int i=0; i<_triangles.length(); i++) {

      // get the vertices of the triangle
      RefSCVector A(_vertices[triangle_vertex(i,0)]->point());
      RefSCVector B(_vertices[triangle_vertex(i,1)]->point());
      RefSCVector C(_vertices[triangle_vertex(i,2)]->point());

      // project the vertices onto the xy plane
      RefSCVector Axy(A.dim()); Axy.assign(A); Axy[2] = 0.0;
      RefSCVector Bxy(B.dim()); Bxy.assign(B); Bxy[2] = 0.0;
      RefSCVector Cxy(C.dim()); Cxy.assign(C); Cxy[2] = 0.0;

      // construct the legs of the triangle in the xy plane
      RefSCVector BAxy = Bxy - Axy;
      RefSCVector CAxy = Cxy - Axy;

      // find the lengths of the legs of the triangle in the xy plane
      double baxy = sqrt(BAxy.dot(BAxy));
      double caxy = sqrt(CAxy.dot(CAxy));

      // if one of the legs is of length zero, then there is
      // no contribution from this triangle
      if (baxy < 1.e-16 || caxy < 1.e-16) continue;

      // find the sine of the angle between the legs of the triangle
      // in the xy plane
      double costheta = BAxy.dot(CAxy)/(baxy*caxy);
      double sintheta = sqrt(1.0 - costheta*costheta);

      // the area of the triangle in the xy plane
      double areaxy = 0.5 * baxy * caxy * sintheta;

      // the height of the three corners of the triangle
      // (relative to the z plane)
      double hA = A[2];
      double hB = B[2];
      double hC = C[2];

      // the volume of the space under the triangle
      double volume = areaxy * (hA + (hB + hC - 2.0*hA)/3.0);

      // the sum of the gradients of the triangle's three vertices
      // along the projection axis (z)
      double zgrad
        = _vertices[triangle_vertex(i,0)]->gradient()[2]
        + _vertices[triangle_vertex(i,1)]->gradient()[2]
        + _vertices[triangle_vertex(i,2)]->gradient()[2];

      if (zgrad > 0.0) {
          result += volume;
        }
      else {
          result -= volume;
        }

    }

  // If the volume is negative, then the surface gradients were
  // opposite in sign to the direction assumed.  Flip the sign
  // to fix.
  return fabs(result);
}

void
TriangulatedSurface::add_triangle(const RefTriangle&t)
{
  if (_completed_surface) clear();
  _triangles.add(t);
}

// If a user isn't keeping track of edges while add_triangle is being
// used to build the surface, then this can be called to see if an edge
// already exists (at a great performance cost).
RefEdge
TriangulatedSurface::find_edge(const RefVertex& v1, const RefVertex& v2)
{
  int i;

  for (i=0; i<_triangles.length(); i++) {
      RefTriangle t = _triangles[i];
      RefEdge e1 = t->edge(0);
      RefEdge e2 = t->edge(1);
      RefEdge e3 = t->edge(2);
      if (e1->vertex(0) == v1 && e1->vertex(1) == v2) return e1;
      if (e1->vertex(1) == v1 && e1->vertex(0) == v2) return e1;
      if (e2->vertex(0) == v1 && e2->vertex(1) == v2) return e2;
      if (e2->vertex(1) == v1 && e2->vertex(0) == v2) return e2;
      if (e3->vertex(0) == v1 && e3->vertex(1) == v2) return e3;
      if (e3->vertex(1) == v1 && e3->vertex(0) == v2) return e3;
    }

  return 0;
}

void
TriangulatedSurface::print(FILE*fp)
{
  fprintf(fp,"TriangulatedSurface:\n");
  int i;

  int np = nvertex();
  fprintf(fp," %3d Vertices:\n",np);
  for (i=0; i<np; i++) {
      RefVertex p = vertex(i);
      int dim = p->dimension();
      fprintf(fp,"  %3d:",i);
      for (int j=0; j<dim; j++) {
          fprintf(fp," % 15.10f", (double)p->point()[j]);
        }
      fprintf(fp,"\n");
    }

  int ne = nedge();
  fprintf(fp," %3d Edges:\n",ne);
  for (i=0; i<ne; i++) {
      RefEdge e = edge(i);
      fprintf(fp,"  %3d: %3d %3d\n",i,
              _vertices.iseek(e->vertex(0)),
              _vertices.iseek(e->vertex(1)));
    }

  int nt = ntriangle();
  fprintf(fp," %3d Triangles:\n",nt);
  for (i=0; i<nt; i++) {
      RefTriangle tri = triangle(i);
      fprintf(fp,"  %3d: %3d %3d %3d\n",i,
              _edges.iseek(tri->edge(0)),
              _edges.iseek(tri->edge(1)),
              _edges.iseek(tri->edge(2)));
    }
}

void
TriangulatedSurface::print_vertices_and_triangles(FILE*fp)
{
  fprintf(fp,"TriangulatedSurface:\n");
  int i;

  int np = nvertex();
  fprintf(fp," %3d Vertices:\n",np);
  for (i=0; i<np; i++) {
      RefVertex p = vertex(i);
      int dim = p->dimension();
      fprintf(fp,"  %3d:",i);
      for (int j=0; j<dim; j++) {
          fprintf(fp," % 15.10f", (double)p->point()[j]);
        }
      fprintf(fp,"\n");
    }

  int nt = ntriangle();
  fprintf(fp," %3d Triangles:\n",nt);
  for (i=0; i<nt; i++) {
      RefTriangle tri = triangle(i);
      fprintf(fp,"  %3d: %3d %3d %3d\n",i,
              _triangle_vertex[i][0],
              _triangle_vertex[i][1],
              _triangle_vertex[i][2]);
    }
}

void
TriangulatedSurface::print_geomview_format(FILE*fp)
{
  fprintf(fp,"OFF\n");

  fprintf(fp,"%d %d %d\n",nvertex(),ntriangle(),nedge());
  int i;

  int np = nvertex();
  for (i=0; i<np; i++) {
      RefVertex p = vertex(i);
      int dim = p->dimension();
      for (int j=0; j<dim; j++) {
          fprintf(fp," % 15.10f", (double)p->point()[j]);
        }
      fprintf(fp,"\n");
    }

  int nt = ntriangle();
  for (i=0; i<nt; i++) {
      RefTriangle tri = triangle(i);
      fprintf(fp," 3 %3d %3d %3d\n",
              _triangle_vertex[i][0],
              _triangle_vertex[i][1],
              _triangle_vertex[i][2]);
    }
}

//////////////////////////////////////////////////////////////////////
// TriangulatedSurface10

TriangulatedSurface10::TriangulatedSurface10(RefVolume&vol,double isovalue):
  _vol(vol),
  _isovalue(isovalue)
{
  // improve the integrator
  set_integrator(new GaussTriangleIntegrator(7));
}

TriangulatedSurface10::~TriangulatedSurface10()
{
}

void
TriangulatedSurface10::complete_surface()
{
  TriangulatedSurface::complete_surface();

  _edges4.set_length(nedge());
  _triangles10.set_length(ntriangle());

  int i;
  for (i=0; i<nedge(); i++) {
      _edges4[i] = new Edge4(vertex(edge_vertex(i,0)),
                             vertex(edge_vertex(i,1)),
                             _vol,_isovalue);
      _edges[i] = _edges4[i].pointer();
    }

  for (i=0; i<ntriangle(); i++) {
      _triangles10[i] = new Triangle10(_edges4[triangle_edge(i,0)],
                                       _edges4[triangle_edge(i,1)],
                                       _edges4[triangle_edge(i,2)],
                                       _vol,_isovalue);
      _triangles[i] = _triangles10[i].pointer();
    }
}

double
TriangulatedSurface10::area()
{
  TriangulatedSurfaceIntegrator tsi(*this);

  double area = 0.0;
  for (tsi = 0; tsi; tsi++) {
      area += tsi.w();
    }
  
  return area;
}

double
TriangulatedSurface10::volume()
{
  TriangulatedSurfaceIntegrator tsi(*this);

  double volume = 0.0;
  RefVertex current = tsi.current();
  for (tsi = 0; tsi; tsi++) {
      RefSCVector norm(current->gradient());
      norm.normalize();
      volume += tsi.w() * norm[0] * current->point()[0];
    }
  
  return volume;
}

void
TriangulatedSurface10::remove_short_edges(double length_cutoff)
{
  fprintf(stderr,"TriangulatedSurface10::remove_short_edges(): "
          "not implemented\n");
  abort();
}

//////////////////////////////////////////////////////////////////////
// TriangulatedSurfaceIntegrator

TriangulatedSurfaceIntegrator::
  TriangulatedSurfaceIntegrator(TriangulatedSurface&ts)
{
  _ts = &ts;

  _itri = 0;
  _irs = 0;

  if (_ts->nvertex()) {
      RefSCDimension n = _ts->vertex(0)->point().dim();
      RefSCVector p(n);
      RefSCVector g(n);
      _current = new Vertex(p,g);
    }
}

TriangulatedSurfaceIntegrator::
  ~TriangulatedSurfaceIntegrator()
{
}

int
TriangulatedSurfaceIntegrator::
  vertex_number(int i)
{
  return _ts->triangle_vertex(_itri,i);
}

RefVertex
TriangulatedSurfaceIntegrator::
  current()
{
  return _current;
}

TriangulatedSurfaceIntegrator::
  operator int()
{
  if (_itri < 0 || _itri >= _ts->ntriangle()) return 0;

  TriangleIntegrator* i = _ts->integrator(_itri);
  _s = i->s(_irs);
  _r = i->r(_irs);
  _weight = i->w(_irs);
  _surface_element = _ts->triangle(_itri)->interpolate(_r,_s,_current);

  //printf("%3d: r=%5.3f s=%5.3f w=%8.5f dA=%8.5f ",
  //       _itri, _r, _s, _weight, _surface_element);
  //_current->print(stdout);

  return (int) 1;
}

void
TriangulatedSurfaceIntegrator::
  operator ++()
{
  int n = _ts->integrator(_itri)->n();
  if (_irs == n-1) {
      _irs = 0;
      _itri++;
    }
  else {
      _irs++;
    }
}

int
TriangulatedSurfaceIntegrator::
  operator = (int i)
{
  _itri = i;
  _irs = 0;
  return i;
}

