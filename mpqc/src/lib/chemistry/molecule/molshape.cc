
#ifdef __GNUC__
#pragma implementation
#endif

extern "C" {
# include <math.h>
  }

#include <chemistry/molecule/molshape.h>
#include <chemistry/molecule/molecule.h>
#include <math/scmat/matrix3.h>

////////////////////////////////////////////////////////////////////////
// VDWShape

#define CLASSNAME VDWShape
#define PARENTS public UnionShape
#define HAVE_KEYVAL_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
VDWShape::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = UnionShape::_castdown(cd);
  return do_castdowns(casts,cd);
}

VDWShape::VDWShape(const RefMolecule&mol)
{
  initialize(mol);
}

VDWShape::VDWShape(const RefKeyVal&keyval)
{
  RefMolecule mol = keyval->describedclassvalue("molecule");
  initialize(mol);
}

void
VDWShape::initialize(const RefMolecule&mol)
{
  _shapes.clear();
  for (int i=0; i<mol->natom(); i++) {
      SCVector3 r;
      for (int j=0; j<3; j++) r[j] = mol->operator[](i)[j];
      add_shape(
          new SphereShape(r,mol->operator[](i).element().vdw_radius())
          );
    }
}

VDWShape::~VDWShape()
{
}

////////////////////////////////////////////////////////////////////////
// static functions for ConnollyShape and ConnollyShape2

static double
find_atom_size(const RefAtomInfo& a, ChemicalElement&element)
{
  return a->radius(element);
}

////////////////////////////////////////////////////////////////////////
// ConnollyShape

#define CLASSNAME ConnollyShape
#define PARENTS public UnionShape
#define HAVE_KEYVAL_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
ConnollyShape::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = UnionShape::_castdown(cd);
  return do_castdowns(casts,cd);
}

ConnollyShape::ConnollyShape(const RefKeyVal&keyval)
{
  RefMolecule mol = keyval->describedclassvalue("molecule");
  double probe_radius = keyval->doublevalue("probe_radius");
  if (keyval->error() != KeyVal::OK) {
      probe_radius = 2.6456173;
    }
  atominfo_ = keyval->describedclassvalue("atominfo");
  initialize(mol,probe_radius);
}

void
ConnollyShape::initialize(const RefMolecule&mol,double probe_radius)
{
  _shapes.clear();
  ArraysetRefSphereShape spheres;
  for (int i=0; i<mol->natom(); i++) {
      SCVector3 r;
      for (int j=0; j<3; j++) r[j] = mol->operator[](i)[j];
      RefSphereShape
        sphere(
            new SphereShape(r,find_atom_size(atominfo_,
                                             mol->operator[](i).element()))
            );
      add_shape(sphere.pointer());
      spheres.add(sphere);
    }

  ////////////////////// Leave out the other shapes
  //return;

  for (i=0; i<spheres.length(); i++) {
      for (int j=0; j<i; j++) {
          RefShape th =
            UncappedTorusHoleShape::newUncappedTorusHoleShape(probe_radius,
                                              *(spheres[i].pointer()),
                                              *(spheres[j].pointer()));
          if (th.null()) continue;
          add_shape(th);

          ////////////////////// Leave out the three sphere shapes
          //continue;
          
          // now check for excluding volume for groups of three spheres
          for (int k=0; k<j; k++) {
              RefShape e =
                Uncapped5SphereExclusionShape::
              newUncapped5SphereExclusionShape(probe_radius,
                                               *(spheres[i].pointer()),
                                               *(spheres[j].pointer()),
                                               *(spheres[k].pointer()));
              if (e.nonnull()) add_shape(e);
            }
        }
    }
}

ConnollyShape::~ConnollyShape()
{
}

////////////////////////////////////////////////////////////////////////
// ConnollyShape2

#define CLASSNAME ConnollyShape2
#define PARENTS public Shape
#define HAVE_KEYVAL_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
ConnollyShape2::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = Shape::_castdown(cd);
  return do_castdowns(casts,cd);
}

ConnollyShape2::ConnollyShape2(const RefKeyVal&keyval)
{
  sphere = 0;
  RefMolecule mol = keyval->describedclassvalue("molecule");
  probe_r = keyval->doublevalue("probe_radius");
  if (keyval->error() != KeyVal::OK) {
      probe_r = 2.6456173;
    }
  atominfo_ = keyval->describedclassvalue("atominfo");
  initialize(mol,probe_r);
}

#if COUNT_CONNOLLY2
int ConnollyShape2::n_total_ = 0;
int ConnollyShape2::n_inside_vdw_ = 0;
int ConnollyShape2::n_with_nsphere_[CONNOLLYSHAPE2_N_WITH_NSPHERE_DIM];
#endif

void
ConnollyShape2::print_counts(FILE*fp)
{
  fprintf(fp,"ConnollyShape2::print_counts():\n");
#if COUNT_CONNOLLY2
  fprintf(fp,"  n_total = %d\n",
          n_total_);
  fprintf(fp,"  n_inside_vdw = %d\n",
          n_inside_vdw_);
  for (int i=0; i<CONNOLLYSHAPE2_N_WITH_NSPHERE_DIM-1; i++) {
      fprintf(fp,"  n with nsphere = %2d: %d\n",
              i, n_with_nsphere_[i]);
    }
  fprintf(fp,"  n with nsphere >= %d: %d\n",
          CONNOLLYSHAPE2_N_WITH_NSPHERE_DIM-1, n_with_nsphere_[CONNOLLYSHAPE2_N_WITH_NSPHERE_DIM-1]);
#else
  fprintf(fp,"  No count information is available.\n");
#endif
}

void
ConnollyShape2::initialize(const RefMolecule&mol,double probe_radius)
{
  if (sphere) delete[] sphere;
  n_spheres = mol->natom();
  sphere = new CS2Sphere[n_spheres];

  for (int i=0; i<n_spheres; i++) {
      SCVector3 r;
      for (int j=0; j<3; j++) r[j] = mol->operator[](i)[j];
      sphere[i].initialize(r,find_atom_size(atominfo_,
                                            mol->operator[](i).element())
                              + probe_r);
    }
}

ConnollyShape2::~ConnollyShape2()
{
  if (sphere) delete[] sphere;
}

double
ConnollyShape2::distance_to_surface(const SCVector3&r, SCVector3*grad) const
{
#if COUNT_CONNOLLY2
  n_total_++;
#endif
  
  // can't compute grad so zero it if it is requested
  if (grad) {
      *grad = 0.0;
    }

  CS2Sphere probe_centers(r,probe_r);

  const int max_local_spheres = 60;
  CS2Sphere local_sphere[max_local_spheres];

  const double outside = 1.0;
  const double inside = -1.0;

  // find out which spheres are near the probe_centers sphere
  int n_local_spheres = 0;
  for (int i=0; i<n_spheres; i++) {
      double distance = sphere[i].distance(probe_centers);
      double r_i = sphere[i].radius();
      if (distance < r_i - probe_r) {
#if COUNT_CONNOLLY2
          n_inside_vdw_++;
#endif
          return inside;
        }
      else if (distance < r_i + probe_r) {
          if (n_local_spheres == max_local_spheres) {
              fprintf(stderr,"ConnollyShape2::distance_to_surface:"
                      " max_local_spheres exceeded\n");
              abort();
            }
          local_sphere[n_local_spheres] = sphere[i];
          n_local_spheres++;
        }
    }

#if COUNT_CONNOLLY2
  if (n_local_spheres >= CONNOLLYSHAPE2_N_WITH_NSPHERE_DIM) {
      n_with_nsphere_[CONNOLLYSHAPE2_N_WITH_NSPHERE_DIM-1]++;
    }
  else {
      n_with_nsphere_[n_local_spheres]++;
    }
#endif

  if (probe_centers.intersect(local_sphere,n_local_spheres)
      == 1) return inside;
  return outside;
}

void
ConnollyShape2::boundingbox(double valuemin,
                            double valuemax,
                            SCVector3& p1, SCVector3& p2)
{
  int i,j;
  if (valuemin < -1.0 || valuemax > 1.0) {
      fprintf(stderr,"ConnollyShape2::boundingbox: value out of range\n");
      abort();
    }

  if (n_spheres == 0) {
      for (i=0; i<3; i++) {
          p1[i] = 0.0;
          p2[i] = 0.0;
        }
      return;
    }

  double r = sphere[0].radius() - probe_r;
  SCVector3 v1(sphere[0].x() - r, sphere[0].y() - r, sphere[0].z() - r);
  SCVector3 v2(sphere[0].x() + r, sphere[0].y() + r, sphere[0].z() + r);

  for (i=1; i<n_spheres; i++) {
      double r = sphere[i].radius() - probe_r;
      for (j=0; j<3; j++) {
          if (v1[j] > sphere[i].center()[j] - r) {
              v1[j] = sphere[i].center()[j] - r;
            }
          if (v2[j] < sphere[i].center()[j] + r) {
              v2[j] = sphere[i].center()[j] + r;
            }
        }
    }
  
  for (i=0; i<3; i++) {
      p1[i] = v1[i] - 0.01;
      p2[i] = v2[i] + 0.01;
    }
}

////////////////////////////////////////////////////////////////////////
// interval class needed by CS2Sphere

// Simple class to keep track of regions along an interval
class interval
{
    int _nsegs;                // # disjoint segments in interval
    int _max_segs;              // # segments currently allocated

    double *_min, *_max;        // arrays of ranges for segments

  private:
    // internal member function to compact interval list--this 
    // assumes that new segment is located in last element of 
    // _min and _max
    void compact(void)
    {

        if (_nsegs==1) return;

        // case 0 new segment is disjoint and below all other segments
        if (_max[_nsegs-1] < _min[0]) 
        {
            double mintmp=_min[_nsegs-1];
            double maxtmp=_max[_nsegs-1];
            for (int i=_nsegs-2; i>=0 ; i--)
            {
                _min[i+1]=_min[i];
                _max[i+1]=_max[i];
            }
            _min[0]=mintmp;
            _max[0]=maxtmp;
            return;
        }

        // case 1: new segment is disjoint and above all other segments
        if (_min[_nsegs-1] > _max[_nsegs-2]) return;

        // Fast forward to where this interval belongs
        int icount=0;
        while (_min[_nsegs-1] > _max[icount]) icount++;

        // case 2: new segment is disjoint and between two segments
        if (_max[_nsegs-1] < _min[icount]) 
        {
            double mintmp=_min[_nsegs-1];
            double maxtmp=_max[_nsegs-1];
            for (int i=_nsegs-2; i >= icount; i--)
            {
                _min[i+1]=_min[i];
                _max[i+1]=_max[i];
            }
            _min[icount]=mintmp;
            _max[icount]=maxtmp;
            return;
        }

        // new segment must overlap lower part of segment icount,
        // so redefine icount's lower boundary
        _min[icount] = (_min[_nsegs-1] < _min[icount])?
            _min[_nsegs-1]:_min[icount];

        // Now figure how far up this new segment extends
        // case 3: if it doesn't extend beyond this segment, just exit
        if (_max[_nsegs-1] < _max[icount]) { _nsegs--; return;}

        // Search forward till we find its end
        int jcount=icount;
        while (_max[_nsegs-1] > _max[jcount]) jcount++;
        
        // Case 4
        // The new segment goes to the end of all the other segments
        if (jcount == _nsegs-1)
        {
            _max[icount]=_max[_nsegs-1];
            _nsegs=icount+1;
            return;
        }

        // Case 5 
        // The new segment ends between segments
        if (_max[_nsegs-1] < _min[jcount])
        {
            _max[icount]=_max[_nsegs-1];
            // Now clobber all the segments covered by the new one
            int kcount=icount+1;
            for (int i=jcount; i<_nsegs; i++)
            {
                _min[kcount]=_min[i];
                _max[kcount]=_max[i];
                kcount++;
            }
            _nsegs=kcount-1;
            return;
        }
        
        // Case 6 
        // The new segment ends inside a segment
        if (_max[_nsegs-1] >= _min[jcount])
        {
            _max[icount]=_max[jcount];
            // Now clobber all the segments covered by the new one
            int kcount=icount+1;
            for (int i=jcount+1; i<_nsegs; i++)
            {
                _min[kcount]=_min[i];
                _max[kcount]=_max[i];
                kcount++;
            }
            _nsegs=kcount-1;
            return;
        }

        // Shouldn't get here!
        printf(" Found no matching cases in interval::compact()\n");
        print();
        exit(1);
    }

  public:
    interval(void):_nsegs(0),_max_segs(10) 
   { _min = (double*) malloc(_max_segs*sizeof(double));   // Use malloc so
      _max = (double*) malloc(_max_segs*sizeof(double));} //we can use realloc

    ~interval() { free(_min); free(_max); }
    
    // add a new segment to interval
    void add(double min, double max)
    {
        if (min > max) {double tmp=min; min=max; max=tmp;}
        if (_nsegs == _max_segs)
        {
            _max_segs *= 2;
            _min=(double *)realloc(_min, _max_segs*sizeof(double));
            _max=(double *)realloc(_max, _max_segs*sizeof(double));
        }
        
        _min[_nsegs]=min;
        _max[_nsegs]=max;
        _nsegs++;
        compact();
    }
    
    // Test to see if the interval is complete over {min, max}
    int test_interval(double min, double max)
    {
        if (_nsegs == 0) return 0;
  
        if (min > max) {double tmp=min; min=max; max=tmp;}
        
        if (min < _min[0] || max > _max[_nsegs-1]) return 0;
        for (int i=0; i < _nsegs; i++)
        {
            if (min > _min[i] && max < _max[i]) return 1;
            if (max < _min[i]) return 0;
        }
        return 0;
    }
    
    // Print out the currect state of the interval
    void print()
    {
        printf(" _nsegs=%d; _max_segs=%d\n",_nsegs, _max_segs);
        for (int i=0; i<_nsegs; i++)
            printf("min[%d]=%7.4lf, max[%d]=%7.4lf\n",i,_min[i],i,_max[i]); 
    }

    void clear() { _nsegs = 0; }
};

////////////////////////////////////////////////////////////////////////
// CS2Sphere

#if COUNT_CONNOLLY2
int CS2Sphere::n_no_spheres_ = 0;
int CS2Sphere::n_probe_enclosed_by_a_sphere_ = 0;
int CS2Sphere::n_probe_center_not_enclosed_ = 0;
int CS2Sphere::n_surface_of_s0_not_covered_ = 0;
int CS2Sphere::n_plane_totally_covered_ = 0;
int CS2Sphere::n_internal_edge_not_covered_ = 0;
int CS2Sphere::n_totally_covered_ = 0;
#endif

void
CS2Sphere::print_counts(FILE*fp)
{
  fprintf(fp,"CS2Sphere::print_counts():\n");
#if COUNT_CONNOLLY2
  fprintf(fp,"  n_no_spheres = %d\n",
          n_no_spheres_);
  fprintf(fp,"  n_probe_enclosed_by_a_sphere = %d\n",
          n_probe_enclosed_by_a_sphere_);
  fprintf(fp,"  n_probe_center_not_enclosed = %d\n",
          n_probe_center_not_enclosed_);
  fprintf(fp,"  n_surface_of_s0_not_covered = %d\n",
          n_surface_of_s0_not_covered_);
  fprintf(fp,"  n_plane_totally_covered_ = %d\n",
          n_plane_totally_covered_);
  fprintf(fp,"  n_internal_edge_not_covered = %d\n",
          n_internal_edge_not_covered_);
  fprintf(fp,"  n_totally_covered = %d\n",
          n_totally_covered_);
#else
  fprintf(fp,"  No count information is available.\n");
#endif
}

// Function to determine if the centers of a bunch of spheres are separated
// by a plane from the center of another plane

// s0 is assumed to be at the origin.

// Return 1 if all of the points can be placed on the same side of a
// plane passing through s0's center.
static int
same_side(const CS2Sphere& s0, CS2Sphere *s, int n_spheres)
{
  if (n_spheres <= 3) return 1;

  SCVector3 perp;
  int sign;

  for (int i=0; i<n_spheres; i++)
    {
      for (int j=0; j<i; j++)
        {
          perp = s[i].center().perp_unit(s[j].center());
          int old_sign=0;
          for (int k=0; k < n_spheres; k++)
            {
              if (i != k && j != k)
                {
                  sign=(perp.dot(s[k].center()) < 0)? -1:1;
                  if (old_sign && old_sign != sign)
                      goto next_plane;
                  old_sign=sign;
                }
            }
          // We found a  plane with all centers on one side
          return 1;
          next_plane:
          continue;
        }
    }
  // All of the planes had points on both sides.
  return 0;
}

double
CS2Sphere::common_radius(CS2Sphere &asphere)
{
  double d=distance(asphere);
  double s=0.5*(d+_radius+asphere._radius);
  double p = s*(s-d)*(s-_radius)*(s-asphere._radius);
  //printf("common_radius: p = %5.3f\n", p);
  if (p <= 0.0) return 0.0;
  return 2.*sqrt(p)/d;
}

#define PRINT_SPECIAL_CASES 0
#if PRINT_SPECIAL_CASES
static void
print_spheres(const CS2Sphere& s0, CS2Sphere* s, int n_spheres)
{
  static int output_number;
  char filename[80];
  sprintf(filename,"spherelist_%d.oogl",output_number);
  FILE* fp = fopen(filename,"w");
  fprintf(fp,"LIST\n");
  fprintf(fp,"{\n");
  fprintf(fp,"  appearance {\n");
  fprintf(fp,"      material {\n");
  fprintf(fp,"         ambient 0.5 0.1 0.1\n");
  fprintf(fp,"         diffuse 1.0 0.2 0.2\n");
  fprintf(fp,"       }\n");
  fprintf(fp,"    }\n");
  fprintf(fp," = SPHERE\n");
  fprintf(fp," %15.8f %15.8f %15.8f %15.8f\n",
          s0.radius(), s0.x(), s0.y(), s0.z());
  fprintf(fp,"}\n");
  for (int i=0; i<n_spheres; i++) {
      fprintf(fp,"{ = SPHERE\n");
      fprintf(fp," %15.8f %15.8f %15.8f %15.8f\n",
              s[i].radius(), s[i].x(), s[i].y(), s[i].z());
      fprintf(fp,"}\n");
    }
  fclose(fp);
  output_number++;
}
#endif

// Function to determine if there is any portion of s0 that 
// is not inside one or more of the spheres in s[]
int
CS2Sphere::intersect(CS2Sphere *s, int n_spheres) const
{
    if (n_spheres == 0) {
        n_no_spheres_++;
        return 0;
      }
    CS2Sphere s0;
    s0 = *this;
    // Declare an interval object to manage overlap information
    // it is static so it will only call malloc twice
    static interval intvl;
    // First make sure that at least one sphere in s[] contains
    // the center of s0 and that s0 is not contained inside
    // one of the spheres
    int center_is_contained = 0;
    for (int i=0; i<n_spheres; i++)
    {
        double d=s0.distance(s[i]);
        if (d+s0.radius() < s[i].radius()) {
            n_probe_enclosed_by_a_sphere_++;
            return 1;
          }
        if (d < s[i].radius()) center_is_contained = 1;
    }
    if (!center_is_contained) {
        n_probe_center_not_enclosed_++;
        return 0;
      }
    
    // Let's first put s0 at the origin
    for (i=0; i<n_spheres; i++)
        s[i].recenter(s0.center());
    s0.recenter(s0.center());
            
    // Now check to make sure that the surface of s0 is completely
    // included in spheres in s[], by making sure that all the
    // circles describing the intersections of every sphere with
    // s0 are included in at least one other sphere.
    double epsilon=1.e-8;
    for (i=0; i<n_spheres; i++)
    {
        // calculate radius of the intersection of s0 and s[i]
        double cr = s0.common_radius(s[i]);
        if (cr == 0.0) {
            continue;
        }
        
        // We're chosing that the intersection of s[i] and s0 
        // occurs parallel to the x-y plane, so we'll need to rotate the
        // center of s[j] appropriately.  
        // Create a rotation matrix that take the vector from 
        // the centers of s0 to s[i] and puts it on the z axis
        static const SCVector3 Zaxis(0.0, 0.0, 1.0);
        SCMatrix3 rot = rotation_mat(s0.center_vec(s[i]),Zaxis);
        
        // Now calculate the Z position of the intersection of 
        // s0 and s[i]
        double d=s0.distance(s[i]);
        double z_plane;
        if (s[i].radius()*s[i].radius() < d*d+s0.radius()*s0.radius())
            z_plane=sqrt(s0.radius()*s0.radius()-cr*cr);
        else
            z_plane=-sqrt(s0.radius()*s0.radius()-cr*cr);
        
        // Initialize the interval object
        intvl.clear();
        
        // Loop over the other spheres
        for (int j=0; j<n_spheres; j++)
            if (i != j)
            {
                // Rotate the center of s[j] to appropriate refence frame
                SCVector3 rcent = rot*s0.center_vec(s[j]);
                
                double x0=rcent.x();
                double y0=rcent.y();
                double z0=rcent.z();
                
                // Does this sphere even reach the plane where
                // the intersection of s0 and s[i] occurs?
                // If not, let's go to the next sphere
                double z_dist=s[j].radius()*s[j].radius()-
                    (z0-z_plane)*(z0-z_plane);
                if (z_dist < 0.0)
                    continue;
                
                // Calculate radius of circular projection of s[j]
                // onto s0-s[i] intersection plane
                double r_2=z_dist;
                
                // Precalculate a bunch of factors 
                double cr_2=cr*cr;
                double x0_2=x0*x0; double y0_2=y0*y0;
                double dist=sqrt(x0_2+y0_2);
                
                // If the projection of s[j] on x-y doesn't reach the
                // intersection of s[i] and s0, continue.
                if (r_2 < (dist-cr)*(dist-cr))
                    continue;
                
                // If the projection of s[j] on x-y engulfs the intersection
                // of s[i] and s0, cover interval and continue
                if (r_2 > (dist+cr)*(dist+cr))
                {
                    intvl.add(0, 2.*M_PI);
                    continue;
                }
                
                // Calculation the radical in the quadratic equation
                // determining the overlap of the two circles
                double radical=x0_2*(-cr_2*cr_2 + 2*cr_2*r_2 - 
                                     r_2*r_2 + 2*cr_2*x0_2 + 
                                     2*r_2*x0_2 - x0_2*x0_2 + 
                                     2*cr_2*y0_2 + 2*r_2*y0_2 - 
                                     2*x0_2*y0_2 - y0_2*y0_2);
                
                // Check to see if there's any intersection at all
                // I.e. if one circle is inside the other  (Note that
                // we've already checked to see if s[j] engulfs
                // the intersection of s0 and s[i])
                if (radical <= 0.0) continue;
                
                // Okay, go ahead and calculate the intersection points
                double x_numer = cr_2*x0_2 - r_2*x0_2 + x0_2*x0_2 + x0_2*y0_2;
                double x_denom = 2*x0*x0_2 + 2*x0*y0_2;
                double y_numer = cr_2*y0 - r_2*y0 + x0_2*y0 + y0*y0_2;
                double y_denom = 2*(x0_2 + y0_2);
                
                double sqrt_radical = sqrt(radical);
                
                double x_0=(x_numer - y0*sqrt_radical)/x_denom;
                double y_0=(y_numer + sqrt_radical)/y_denom;
                double x_1=(x_numer + y0*sqrt_radical)/x_denom;
                double y_1=(y_numer - sqrt_radical)/y_denom;
                
                // Now calculate the angular range of these ordered
                // points and place them on the first Riemann sheet.
                // and sort their order
                double theta1=atan2(y_0, x_0);
                double theta2=atan2(y_1, x_1);
                if (theta1 < 0.0) theta1+=2.*M_PI;
                if (theta2 < 0.0) theta2+=2.*M_PI;
                if (theta1 > theta2)
                {
                    double tmptheta=theta1;
                    theta1=theta2;
                    theta2=tmptheta;
                }
                
                // Determine which of the two possible chords 
                // is inside s[j]
                double dor=(x0-cr)*(x0-cr)+y0*y0;
                if (dor < r_2)
                {
                    intvl.add(0, theta1);
                    intvl.add(theta2, 2.*M_PI);
                }
                else            
                {
                    intvl.add(theta1, theta2);
                }
                
                // Now test to see if the range is covered
                if (intvl.test_interval(epsilon, 2.*M_PI-epsilon))
                {
                    // No need to keep testing, move on to next i
                    break;
                }
                
            }
        // If the intersection wasn't totally covered, the sphere
        // intersection is incomplete
        if (!intvl.test_interval(epsilon, 2.*M_PI-epsilon)) {
            n_surface_of_s0_not_covered_++;
            // goto next_test;
            return 0;
        }
    }

    // for the special case of all sphere's centers on one side of
    // a plane passing through s0's center we are done; the probe
    // must be completely intersected.
    if (same_side(s0,s,n_spheres)) {
        n_plane_totally_covered_++;
        return 1;
      }
    
    // As a final test of the surface coverage, make sure that all
    // of the intersection surfaces between s0 and s[] are included 
    // inside more than one sphere.
    int angle_segs;
    double max_angle[2], min_angle[2];
    for (i=0; i<n_spheres; i++)
    {
        // For my own sanity, let's put s[i] at the origin 
        for (int k=0; k<n_spheres; k++)
            if (k != i)
                s[k].recenter(s[i].center());
        s0.recenter(s[i].center());
        s[i].recenter(s[i].center());
        
        for (int j=0; j<i; j++)
        {

            // calculate radius of the intersection of s[i] and s[j]
            double cr = s[i].common_radius(s[j]);
            if (cr == 0.0) {
                continue;                   // s[i] and s[j] don't intersect
            }

            // We're chosing that the intersection of s[i] and s[j]
            // occurs parallel to the x-y plane, so we'll need to rotate the
            // center of all s[]'s and s0 appropriately.  
            // Create a rotation matrix that take the vector from 
            // the centers of s0 to s[i] and puts it on the z axis
            static const SCVector3 Zaxis(0.0, 0.0, 1.0);
            SCMatrix3 rot = rotation_mat(s[i].center_vec(s[j]),Zaxis);
            
            // Now calculate the Z position of the intersection of 
            // s[i] and s[j]
            double d=s[i].distance(s[j]);
            double z_plane;
            if (s[j].radius()*s[j].radius() < s[i].radius()*s[i].radius()+d*d)
                z_plane=sqrt(s[i].radius()*s[i].radius()-cr*cr);
            else
                z_plane=-sqrt(s[i].radius()*s[i].radius()-cr*cr);
            
            // Determine which part of the this intersection
            // occurs within s0
            // Rotate the center of s0 to appropriate refence frame
            SCVector3 rcent = rot*s[i].center_vec(s0);
            
            double x0=rcent.x();
            double y0=rcent.y();
            double z0=rcent.z();
            
            // Does this s0 even reach the plane where
            // the intersection of s[i] and s[j] occurs?
            // If not, let's go to the next sphere j
            double z_dist=s0.radius()*s0.radius()-
                (z0-z_plane)*(z0-z_plane);
            if (z_dist < 0.0)
                continue;
            
            // Calculate radius of circular projection of s0
            // onto s[i]-s[j] intersection plane
            double r_2=z_dist;
            
            // Precalculate a bunch of factors 
            double cr_2=cr*cr;
            double x0_2=x0*x0; double y0_2=y0*y0;
            double dist=sqrt(x0_2+y0_2);
            
            // If the projection of s[j] on x-y doesn't reach the
            // intersection of s[i] and s0, continue.
            if (r_2 < (dist-cr)*(dist-cr))
                continue;
            
            // If the projection of s0 on x-y engulfs the intersection
            // of s[i] and s[j], the intersection interval is 0 to 2pi
            if (r_2 > (dist+cr)*(dist+cr))
            {
                angle_segs=1;
                min_angle[0]=0.0;
                max_angle[0]=2.*M_PI;
            }
            
            // Calculation the radical in the quadratic equation
            // determining the overlap of the two circles
            double radical=x0_2*(-cr_2*cr_2 + 2*cr_2*r_2 - 
                                 r_2*r_2 + 2*cr_2*x0_2 + 
                                 2*r_2*x0_2 - x0_2*x0_2 + 
                                 2*cr_2*y0_2 + 2*r_2*y0_2 - 
                                 2*x0_2*y0_2 - y0_2*y0_2);
            
            // Check to see if there's any intersection at all
            // I.e. if one circle is inside the other  (Note that
            // we've already checked to see if s0 engulfs
            // the intersection of s[i] and s[j]), so this
            // must mean that the intersection of s[i] and s[j]
            // occurs outside s0
            if (radical <= 0.0) continue;

            // Okay, go ahead and calculate the intersection points
            double x_numer = cr_2*x0_2 - r_2*x0_2 + x0_2*x0_2 + x0_2*y0_2;
            double x_denom = 2*x0*x0_2 + 2*x0*y0_2;
            double y_numer = cr_2*y0 - r_2*y0 + x0_2*y0 + y0*y0_2;
            double y_denom = 2*(x0_2 + y0_2);
            
            double sqrt_radical = sqrt(radical);
            
            double x_0=(x_numer - y0*sqrt_radical)/x_denom;
            double y_0=(y_numer + sqrt_radical)/y_denom;
            double x_1=(x_numer + y0*sqrt_radical)/x_denom;
            double y_1=(y_numer - sqrt_radical)/y_denom;
            
            // Now calculate the angular range of these ordered
            // points and place them on the first Riemann sheet.
            // and sort their order
            double theta1=atan2(y_0, x_0);
            double theta2=atan2(y_1, x_1);
            if (theta1 < 0.0) theta1+=2.*M_PI;
            if (theta2 < 0.0) theta2+=2.*M_PI;
            if (theta1 > theta2)
            {
                double tmptheta=theta1;
                theta1=theta2;
                theta2=tmptheta;
            }
            //printf("theta1=%lf, theta2=%lf\n",theta1,theta2);

            // Determine which of the two possible chords 
            // is inside s0

            // But first see if s0 is inside this intersection:
            double origin_dist=((x0-cr)*(x0-cr)+(y0*y0));
            if (origin_dist < r_2) // it's the angle containing
                                       // the origin
            {
                angle_segs=2;
                min_angle[0]=0.0;
                max_angle[0]=theta1;
                min_angle[1]=theta2;
                max_angle[1]=2.*M_PI;
            }
            else            // it's the angle not including the origin
            {
                angle_segs=1;
                min_angle[0]=theta1;
                max_angle[0]=theta2;
            }
            
            // Initialize the interval object
            intvl.clear();
            
            // Loop over the other spheres
            for (k=0; k<n_spheres; k++)
            {
                if (k != i && k != j)
                {
                    // Rotate the center of s[k] to appropriate reference frame
                    rcent = rot*s[i].center_vec(s[k]);
                    
                    double x0=rcent.x();
                    double y0=rcent.y();
                    double z0=rcent.z();
                    
                    // Does this sphere even reach the plane where
                    // the intersection of s[i] and s[j] occurs?
                    // If not, let's go to the next sphere
                    double z_dist=s[k].radius()*s[k].radius()-
                        (z0-z_plane)*(z0-z_plane);
                    if (z_dist < 0.0)
                        continue;
                    
                    // Calculate radius of circular projection of s[k]
                    // onto s[i]-s[j] intersection plane
                    double r_2=z_dist;
                    
                    // Precalculate a bunch of factors 
                    double cr_2=cr*cr;
                    double x0_2=x0*x0; double y0_2=y0*y0;
                    double dist=sqrt(x0_2+y0_2);
                    
                    // If the projection of s[k] on x-y doesn't reach the
                    // intersection of s[i] and s[j], continue.
                    if (r_2 < (dist-cr)*(dist-cr))
                        continue;
                    
                    // If the projection of s[k] on x-y engulfs the intersection
                    // of s[i] and s0, cover interval and continue
                    if (r_2 > (dist+cr)*(dist+cr))
                    {
                        intvl.add(0, 2.*M_PI);
                        continue;
                    }
                    
                    // Calculation the radical in the quadratic equation
                    // determining the overlap of the two circles
                    radical=x0_2*(-cr_2*cr_2 + 2*cr_2*r_2 - 
                                  r_2*r_2 + 2*cr_2*x0_2 + 
                                  2*r_2*x0_2 - x0_2*x0_2 + 
                                  2*cr_2*y0_2 + 2*r_2*y0_2 - 
                                  2*x0_2*y0_2 - y0_2*y0_2);
                    
                    // Check to see if there's any intersection at all
                    // I.e. if one circle is inside the other  (Note that
                    // we've already checked to see if s[k] engulfs
                    // the intersection of s[i] and s[j])
                    if (radical <= 0.0) continue;
                    
                    // Okay, go ahead and calculate the intersection points
                    x_numer = cr_2*x0_2 - r_2*x0_2 + x0_2*x0_2 + x0_2*y0_2;
                    x_denom = 2*x0*x0_2 + 2*x0*y0_2;
                    y_numer = cr_2*y0 - r_2*y0 + x0_2*y0 + y0*y0_2;
                    y_denom = 2*(x0_2 + y0_2);
                    
                    sqrt_radical = sqrt(radical);
                    
                    double x_0=(x_numer - y0*sqrt_radical)/x_denom;
                    double y_0=(y_numer + sqrt_radical)/y_denom;
                    double x_1=(x_numer + y0*sqrt_radical)/x_denom;
                    double y_1=(y_numer - sqrt_radical)/y_denom;

                    // Now calculate the angular range of these ordered
                    // points and place them on the first Riemann sheet.
                    // and sort their order
                    theta1=atan2(y_0, x_0);
                    theta2=atan2(y_1, x_1);
                    if (theta1 < 0.0) theta1+=2.*M_PI;
                    if (theta2 < 0.0) theta2+=2.*M_PI;
                    if (theta1 > theta2)
                    {
                        double tmptheta=theta1;
                        theta1=theta2;
                        theta2=tmptheta;
                    }
                    //printf("In k loop, k=%d, theta1=%lf, theta2=%lf\n",
                    //       k,theta1, theta2);
                    // Determine which of the two possible chords 
                    // is inside s[k]
                    double origin_dist=((x0-cr)*(x0-cr)+(y0*y0));
                    if (origin_dist < r_2) // it's got the origin
                    {
                        intvl.add(0, theta1);
                        intvl.add(theta2, 2.*M_PI);
                    }
                    else        // it doesn't have the origin
                    {
                        intvl.add(theta1, theta2);
                    }

                    // Now test to see if the range is covered
                    if (intvl.test_interval(min_angle[0]+epsilon,
                                            max_angle[0]-epsilon) &&
                        (angle_segs!=2 ||
                         intvl.test_interval(min_angle[1]+epsilon,
                                             max_angle[1]-epsilon)))
                    {
                        goto next_j;
                    }
                }
            }
            if (!intvl.test_interval(min_angle[0]+epsilon,
                                     max_angle[0]-epsilon))
            {
                // No need to keep testing, return 0
                n_internal_edge_not_covered_++;
                return 0;
                //printf(" Non-internal coverage(1)\n");
                //goto next_test;
            }
            if (angle_segs==2)
            {
                if  (!intvl.test_interval(min_angle[1]+epsilon,
                                          max_angle[1]-epsilon))
                {
                    n_internal_edge_not_covered_++;
                    return 0;
                    //printf(" Non-internal coverage(2)\n");
                    //goto next_test;
                }
                else
                {
                    goto next_j;
                }
            }
          next_j:
            continue;
        }
    }
    
    // Since we made it past all of the sphere intersections, the
    // surface is totally covered
    n_totally_covered_++;
    return 1;
}


