
#ifndef _chemistry_qc_scf_fockdist_h
#define _chemistry_qc_scf_fockdist_h

#include <vector>
#include <util/misc/scint.h>
#include <util/group/message.h>
#include <util/group/thread.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/basis/petite.h>

namespace sc {

class FockBlocks: public sc::RefCount {
    int nblock_;
    std::vector<int> begin_;
    std::vector<int> end_;
    std::vector<int> shell_to_block_;
  public:
    typedef enum { Shell, Atom } Method;

    FockBlocks(const sc::Ref<sc::GaussianBasisSet> &, Method);

    int nblock() const { return nblock_; }
    int begin(int iblock) const {return begin_[iblock];}
    int end(int iblock) const {return end_[iblock];}
    int size(int iblock) const {return end_[iblock]-begin_[iblock];}
    int shell_to_block(int ishell) const { return shell_to_block_[ishell]; }
};

class FockDist: public sc::RefCount {
  protected:
    sc::Ref<FockBlocks> blocks_;

    sc::Ref<sc::MessageGrp> msg_;
    sc::Ref<sc::PetiteList> pl_;
    int nthread_, mythread_;
    int nproc_, myproc_;
    // This returns true if any of the shells in block i are in p1.
    bool in_p1(int i);
  public:
    FockDist(const sc::Ref<sc::GaussianBasisSet> &,
             const sc::Ref<FockBlocks> &,
             const sc::Ref<sc::PetiteList> &pl,
             const sc::Ref<sc::MessageGrp> &,
             int nthread, int mythread);
    virtual ~FockDist();
    virtual void init() = 0;
    // This only returns canonical block indices.  When symmetry is used,
    // this returns only blocks in_p1.  This can be done, since only blocks
    // comprised of shells on the same atom are used.  Checks for p2 and p4
    // must be done on individual shell indices, since the results of those
    // checks depend on the particular shell number, even for shells on
    // the same atom.
    virtual bool get_blocks(int &i, int &j, int &k, int &l) = 0;
    int nblock() const { return blocks_->nblock(); }
    int begin(int iblock) const {return blocks_->begin(iblock);}
    int end(int iblock) const {return blocks_->end(iblock);}
    int size(int iblock) const {return blocks_->size(iblock);}
    const sc::Ref<FockBlocks> &fockblocks() const { return blocks_; }
    // Returns true if a given shell quartet will always end up on the same
    // node/thread.
    virtual bool fixed_integral_map() = 0;
};

class FockDistStatic: public FockDist {
  protected:
    int i_, j_, k_, l_;
    int nglobalthread_, myglobalthread_;
  public:
    FockDistStatic(const sc::Ref<sc::GaussianBasisSet> &,
                   const sc::Ref<FockBlocks> &,
                   const sc::Ref<sc::PetiteList> &pl,
                   const sc::Ref<sc::MessageGrp> &,
                   int nthread, int mythread);
    ~FockDistStatic();
    void init();
    bool fixed_integral_map();
};

class FockDistStatic4: public FockDistStatic {
  private:
    bool next_block();
    bool next_block(int n);
  public:
    FockDistStatic4(const sc::Ref<sc::GaussianBasisSet> &,
                    const sc::Ref<FockBlocks> &,
                    const sc::Ref<sc::PetiteList> &pl,
                    const sc::Ref<sc::MessageGrp> &,
                    int nthread, int mythread);
    ~FockDistStatic4();
    bool get_blocks(int&,int&,int&,int&);
};

class FockDistStatic2: public FockDistStatic {
  private:
    bool next_block();
  public:
    FockDistStatic2(const sc::Ref<sc::GaussianBasisSet> &,
                    const sc::Ref<FockBlocks> &,
                    const sc::Ref<sc::PetiteList> &pl,
                    const sc::Ref<sc::MessageGrp> &,
                    int nthread, int mythread);
    ~FockDistStatic2();
    bool get_blocks(int&,int&,int&,int&);
};

class FockDistDynamic: public FockDist {
  protected:
    typedef std::multimap<std::pair<int,int>, std::pair<int,int>, std::greater<std::pair<int,int> > > pairmap_t;

    pairmap_t ijmap_, klmap_;
    pairmap_t::const_iterator ij_iter_, kl_iter_;

    // If there are not enough processors to use dynamic load balancing,
    // then static load balancing must be used.  These members support
    // static load balancing.
    FockDist *static_;
    bool static_distribution() { return static_ != 0; }
  public:
    FockDistDynamic(const sc::Ref<sc::GaussianBasisSet> &,
                    const sc::Ref<FockBlocks> &,
                    const sc::Ref<sc::PetiteList> &pl,
                    const sc::Ref<sc::MessageGrp> &,
                    int nthread, int mythread,
                    const signed char *pmax,
                    sc::Ref<sc::TwoBodyInt> eri,
                    int l2tol);
    ~FockDistDynamic();
    bool fixed_integral_map();
};

class FockDistDynamic2: public FockDistDynamic {
  protected:
    int i_, j_;

    bool next_block();
    void run_server();
    bool server_get_blocks(int&,int&);
  public:
    FockDistDynamic2(const sc::Ref<sc::GaussianBasisSet> &,
                     const sc::Ref<FockBlocks> &,
                     const sc::Ref<sc::PetiteList> &pl,
                     const sc::Ref<sc::MessageGrp> &,
                     int nthread, int mythread,
                     const signed char *pmax,
                     sc::Ref<sc::TwoBodyInt> eri,
                     int l2tol);
    ~FockDistDynamic2();
    bool get_blocks(int&,int&,int&,int&);
    void init();
};

class FockDistDynamic4: public FockDistDynamic {
  protected:
    typedef std::vector<std::pair<int,int> > pairvec_t;
    pairvec_t ijvec_, klvec_;
    int a_, b_; // this are the iterators used to compute ij and kl
    int a_end_, b_end_;
    int nij_, nkl_;

    sc::MessageGrp::MessageHandle recv_handle_;
    sc::MessageGrp::MessageHandle send_handle_;
    bool requested_work_;
    int work_[4];

    bool next_block();
    void run_server();
    bool server_get_blocks(int&,int&,int&,int&);
    bool ijkl(int &i, int &j, int &k, int &l);
  public:
    FockDistDynamic4(const sc::Ref<sc::GaussianBasisSet> &,
                     const sc::Ref<FockBlocks> &,
                     const sc::Ref<sc::PetiteList> &pl,
                     const sc::Ref<sc::MessageGrp> &,
                     int nthread, int mythread,
                     const signed char *pmax,
                     sc::Ref<sc::TwoBodyInt> eri,
                     int l2tol);
    ~FockDistDynamic4();
    bool get_blocks(int&,int&,int&,int&);
    void init();
};

/** FockDistribution is a factory for constructing the
    desired FockDist specialization. */
class FockDistribution: virtual public sc::SavableState {
    /// If true, do dynamic shell distribution.
    int dynamic_;
    /// Gives the number of indices to distribute.  Must be 2 or 4.
    int nindex_;
    /// If true, do shell blocking.  Otherwise do atom blocking.
    int shell_;
    /// If true, cache the integrals, if the integral map is fixed.
    int cache_integrals_;
  public:
    FockDistribution(bool dynamic = 0, bool shell = 1, int nindex = 4,
                     bool cache_integrals = 1);
    FockDistribution(sc::StateIn&);
    FockDistribution(const sc::Ref<sc::KeyVal> &);
    ~FockDistribution();
    void save_data_state(sc::StateOut&);
    sc::Ref<FockBlocks> fockblocks(const sc::Ref<sc::GaussianBasisSet> &);
    sc::Ref<FockDist> fockdist(const sc::Ref<sc::GaussianBasisSet> &,
                               const sc::Ref<FockBlocks> &,
                               const sc::Ref<sc::PetiteList> &pl,
                               const sc::Ref<sc::MessageGrp> &,
                               int nthread, int mythread,
                               const signed char *pmax,
                               sc::Ref<sc::TwoBodyInt> eri,
                               int l2tol);
    void print(std::ostream&o=sc::ExEnv::out0()) const;
    bool cache_integrals() const { return cache_integrals_; }
};

}

#endif
