
#include <chemistry/qc/scf/fockdist.h>

#include <util/class/scexception.h>
#include <util/state/statein.h>
#include <util/state/stateout.h>

#include <algorithm>

#undef FOCKDIST_SERVER
#define FOCKDIST_SERVER 0

using namespace sc;

using sc::Ref;
using sc::GaussianBasisSet;
using sc::MessageGrp;
using sc::ThreadGrp;
using sc::PetiteList;

///////////////////////////////////////////////////////////////////


static sc::ClassDesc FockDistribution_cd(
  typeid(FockDistribution),"FockDistribution",1,"virtual public SavableState",
  0, sc::create<FockDistribution>, sc::create<FockDistribution>);


FockDistribution::FockDistribution(sc::StateIn&si):
  SavableState(si)
{
  si.get(dynamic_);
  si.get(nindex_);
  si.get(shell_);
  si.get(cache_integrals_);
}

FockDistribution::FockDistribution(const sc::Ref<sc::KeyVal> &keyval)
{
  sc::KeyValValueboolean def_dynamic(0);
  dynamic_ = keyval->booleanvalue("dynamic", def_dynamic);

  sc::KeyValValueint def_nindex(4);
  nindex_ = keyval->intvalue("nindex", def_nindex);

  sc::KeyValValueboolean def_shell(1);
  shell_ = keyval->booleanvalue("shell", def_shell);

  sc::KeyValValueboolean def_cache_integrals(1);
  cache_integrals_ = keyval->booleanvalue("cache_integrals", def_cache_integrals);
}

FockDistribution::FockDistribution(bool dynamic, bool shell, int nindex,
                                   bool cache_integrals):
  dynamic_(dynamic),
  shell_(shell),
  nindex_(nindex),
  cache_integrals_(cache_integrals)
{
}

FockDistribution::~FockDistribution()
{
}

void
FockDistribution::save_data_state(sc::StateOut&so)
{
  so.put(dynamic_);
  so.put(nindex_);
  so.put(shell_);
  so.put(cache_integrals_);
}

sc::Ref<FockBlocks>
FockDistribution::fockblocks(const sc::Ref<sc::GaussianBasisSet> &gbs)
{
  FockBlocks::Method blocking;
  if (shell_) {
      blocking = FockBlocks::Shell;
    }
  else {
      blocking = FockBlocks::Atom;
    }
  return new FockBlocks(gbs,blocking);
}

sc::Ref<FockDist>
FockDistribution::fockdist(const sc::Ref<sc::GaussianBasisSet> &gbs,
                           const sc::Ref<FockBlocks> &blocks,
                           const sc::Ref<sc::PetiteList> &pl,
                           const sc::Ref<sc::MessageGrp> &msg,
                           int nthread, int mythread,
                           const signed char *pmax,
                           sc::Ref<sc::TwoBodyInt> eri,
                           int l2tol)
{
  sc::Ref<FockDist> r;
  if (nindex_ == 2) {
      if (dynamic_)
          r = new FockDistDynamic2(gbs, blocks, pl, msg,
                                   nthread, mythread,
                                   pmax, eri, l2tol);
      else
          r = new FockDistStatic2(gbs, blocks, pl, msg,
                                  nthread, mythread);
    }
  else if (nindex_ == 4) {
      if (dynamic_)
          r = new FockDistDynamic4(gbs, blocks, pl, msg,
                                   nthread, mythread,
                                   pmax, eri, l2tol);
      else
          r = new FockDistStatic4(gbs, blocks, pl, msg,
                                  nthread, mythread);
    }
  else {
      throw sc::InputError("nindex must be 2 or 4",
                           __FILE__,
                           __LINE__,
                           "nindex",
                           0,
                           class_desc());
    }

  return r;
}

void
FockDistribution::print(std::ostream&o) const
{
  o << sc::indent << "FockDistribution:" << std::endl
    << sc::incindent
    << sc::indent << "dynamic         = " << dynamic_ << std::endl
    << sc::indent << "nindex          = " << nindex_ << std::endl
    << sc::indent << "shell           = " << shell_ << std::endl
    << sc::indent << "cache_integrals = " << cache_integrals_ << std::endl;
}

///////////////////////////////////////////////////////////////////

FockBlocks::FockBlocks(const Ref<GaussianBasisSet> &gbs,
                       Method blocking)
{

  if (blocking == Shell) {
      nblock_ = gbs->nshell();
      begin_.resize(nblock_);
      end_.resize(nblock_);
      for (int i=0; i<nblock_; i++) {
          begin_[i] = i;
          end_[i] = i + 1;
        }
    }
  else {
      begin_.reserve(gbs->ncenter());
      end_.reserve(gbs->ncenter());
      for (int i=0; i<gbs->ncenter(); i++) {
          begin_.push_back(gbs->shell_on_center(i,0));
          end_.push_back(begin_[i] + gbs->nshell_on_center(i));
        }
      nblock_ = begin_.size();
    }
  shell_to_block_.resize(gbs->nshell());
  for (int i=0,ishell=0; i<nblock(); i++) {
      for (int j=0; j<size(i); j++,ishell++) {
          shell_to_block_[begin(i) + j] = i;
        }
    }
}

///////////////////////////////////////////////////////////////////

FockDist::FockDist(const Ref<GaussianBasisSet> &gbs,
                   const Ref<FockBlocks> &blocks,
                   const Ref<PetiteList> &pl,
                   const Ref<MessageGrp> &msg,
                   int nthread, int mythread)
{
  msg_ = msg;

  nproc_ = msg_->n();
  myproc_ = msg_->me();
  nthread_ = nthread;
  mythread_ = mythread;

  pl_ = pl;

  blocks_ = blocks;
}

FockDist::~FockDist()
{
}

bool
FockDist::in_p1(int iblock)
{
  for (int i=blocks_->begin(iblock); i<blocks_->end(iblock); i++) {
      if (pl_->in_p1(i)) return true;
    }
  return false;
}

///////////////////////////////////////////////////////////////////

FockDistStatic::FockDistStatic(const Ref<GaussianBasisSet> &gbs,
                               const Ref<FockBlocks> &blocks,
                               const Ref<PetiteList> &pl,
                               const Ref<MessageGrp> &msg,
                               int nthread, int mythread):
  FockDist(gbs, blocks, pl, msg, nthread, mythread)
{
  nglobalthread_ = nproc_ * nthread_;
  myglobalthread_ = myproc_ * nthread_ + mythread_;
}

FockDistStatic::~FockDistStatic()
{
}

void
FockDistStatic::init()
{
  i_ = j_ = k_ = l_ = 0;
}

bool
FockDistStatic::fixed_integral_map()
{
  return true;
}

///////////////////////////////////////////////////////////////////

FockDistStatic4::FockDistStatic4(const Ref<GaussianBasisSet> &gbs,
                                 const Ref<FockBlocks> &blocks,
                                 const Ref<PetiteList> &pl,
                                 const Ref<MessageGrp> &msg,
                                 int nthread, int mythread):
  FockDistStatic(gbs,blocks,pl,msg,nthread,mythread)
{
  init();

  // advance to the first contributing block
  while (i_<blocks_->nblock() && !in_p1(i_)) {
      i_++;
    }

  // put the iterator at the correct initial position
  for (int i=0; i<myglobalthread_; i++) {
      if (!next_block()) break;
    }
}

FockDistStatic4::~FockDistStatic4()
{
}

bool
FockDistStatic4::next_block()
{
  l_++;

  // The following commented lines would be included for fully
  // canonical iteration.  However, since we are iterating over blocks
  // with potentially more than one shell in a block and are
  // subsequently testing for canonical shell indices, we must use
  // only the l_ <= k_ constraint, unless there is only a single
  // shell in the k_ == i_ block.
  // if (k_ == i_ && l_ <= j_) return true;
  // if (k_ <  i_ && l_ <= k_) return true;

  if (k_ == i_) {
      if (size(k_) == 1 && l_ <= j_) return true;
    }

  if (k_ < i_ && l_ <= k_) return true;

  l_ = 0;
  k_++;
  if (k_ <= i_) return true;
  k_ = 0;
  j_++;
  if (j_ <= i_) return true;
  j_ = 0;
  do {
      i_++;
    } while (i_<blocks_->nblock() && !in_p1(i_));
  if (i_ < blocks_->nblock()) return true;

  return false;
}

bool
FockDistStatic4::next_block(int n)
{

  while (n) {
      int lend = blocks_->nblock();
      if (k_ == i_ && size(k_) == 1) lend = j_ + 1;
      if (k_ + 1 < lend) lend = k_ + 1;

      if (lend - l_ > n) {
          l_ += n;
          return true;
        }

      n -= lend - l_;

      l_ = 0;
      k_++;
      if (k_ <= i_) continue;
      k_ = 0;
      j_++;
      if (j_ <= i_) continue;
      j_ = 0;
      do {
          i_++;
        } while (i_<blocks_->nblock() && !in_p1(i_));
      if (i_ < blocks_->nblock()) continue;

      return false;
    }

  return true;
}

bool
FockDistStatic4::get_blocks(int &i, int &j, int &k, int &l)
{
  i = i_; j = j_; k = k_; l = l_;

  if (i_ >= blocks_->nblock()) return false;

  // advance to the next block quartet

//   for (int ithread=0; ithread<nglobalthread_; ithread++) {
//       if (!next_block()) break;
//     }

  next_block(nglobalthread_);

  return true;
}

///////////////////////////////////////////////////////////////////

FockDistStatic2::FockDistStatic2(const Ref<GaussianBasisSet> &gbs,
                                 const Ref<FockBlocks> &blocks,
                                 const Ref<PetiteList> &pl,
                                 const Ref<MessageGrp> &msg,
                                 int nthread, int mythread):
  FockDistStatic(gbs,blocks,pl,msg,nthread,mythread)
{
  init();

  // advance to the first contributing block
  while (i_<blocks_->nblock() && !in_p1(i_)) {
      i_++;
    }

  // put the iterator at the correct initial position
  for (int i=0; i<myglobalthread_; i++) {
      if (!next_block()) break;
    }
}

FockDistStatic2::~FockDistStatic2()
{
}

bool
FockDistStatic2::next_block()
{
  j_++;
  if (j_ <= i_) return true;
  j_ = 0;
  do {
      i_++;
    } while (i_<blocks_->nblock() && !in_p1(i_));
  if (i_ < blocks_->nblock()) return true;

  return false;
}

bool
FockDistStatic2::get_blocks(int &i, int &j, int &k, int &l)
{
  i = i_; j = j_; k = k_; l = l_;

  if (i_ >= blocks_->nblock()) return false;

  l_++;

  // The following commented lines would be included for fully
  // canonical iteration.  However, since we are iterating over blocks
  // with potentially more than one shell in a block and are
  // subsequently testing for canonical shell indices, we must use
  // only the l_ <= k_ constraint, unless there is only a single
  // shell in the k_ == i_ block.
  // if (k_ == i_ && l_ <= j_) return true;
  // if (k_ <  i_ && l_ <= k_) return true;

  if (k_ == i_) {
      if (size(k_) == 1 && l_ <= j_) return true;
    }

  if (l_ <= k_) return true;

  l_ = 0;
  k_++;
  if (k_ <= i_) return true;
  k_ = 0;

  // advance to the next block quartet
  for (int ithread=0; ithread<nglobalthread_; ithread++) {
      if (!next_block()) break;
    }

  return true;
}

///////////////////////////////////////////////////////////////////

FockDistDynamic::FockDistDynamic(const sc::Ref<sc::GaussianBasisSet> &gbs,
                                 const Ref<FockBlocks> &blocks,
                                 const sc::Ref<sc::PetiteList> &pl,
                                 const sc::Ref<sc::MessageGrp> &msg,
                                 int nthread, int mythread,
                                 const signed char *pmax,
                                 sc::Ref<sc::TwoBodyInt> eri,
                                 int l2tol):
  FockDist(gbs,blocks,pl,msg,nthread,mythread),
  static_(0)
{
  const bool skip_bound_check = false;

  // Find the maximum of all of pmax:
  int nshell = gbs->nshell();
  int nshell2 = (nshell*(nshell+1))/2;
  int maxpmax;
  if (nshell>0) maxpmax = pmax[0];
  for (int i=0; i<nshell2; i++) {
      if (maxpmax < pmax[i]) maxpmax = pmax[i];
    }

  std::vector<int> nfunc(blocks_->nblock());
  for (int i=0; i<blocks_->nblock(); i++) {
      int tmp = 0;
      for (int j=begin(i); j<end(i); j++) {
          tmp += gbs->shell(j).nfunction();
        }
      nfunc[i] = tmp;
    }
  int ijtasks = 0, kltasks = 0;
  for (int i=0; i<blocks_->nblock(); i++) {
      for (int j=0; j<=i; j++) {
          int integral_bound = INT_MIN;
          for (int ii=begin(i); ii<end(i); ii++) {
              for (int jj=begin(j); jj<end(j); jj++) {
                  int integral_bound_iijj
                      = eri->log2_shell_bound(ii,jj,-1,-1);
                  if (integral_bound < integral_bound_iijj)
                      integral_bound = integral_bound_iijj;
                }
            }
          if (in_p1(i)) {
              if (skip_bound_check
                  || integral_bound + maxpmax > l2tol) {
                  // Using rand() in the key randomizes the order of the
                  // shell pairs having the same size, reducing the chance
                  // of many nodes needing the same block at the same time.
                  ijmap_.insert(std::make_pair(std::make_pair(
                                                   nfunc[i]*nfunc[j],
                                                   rand()),
                                               std::make_pair(i,j)));
                  ijtasks++;
                }
            }
          if (skip_bound_check
              || integral_bound + maxpmax > l2tol) {
              klmap_.insert(std::make_pair(std::make_pair(
                                               nfunc[i]*nfunc[j],
                                               rand()),
                                           std::make_pair(i,j)));
              kltasks++;
            }
        }
    }
}

FockDistDynamic::~FockDistDynamic()
{
  delete static_;
}

bool
FockDistDynamic::fixed_integral_map()
{
  return false;
}

///////////////////////////////////////////////////////////////////

FockDistDynamic2::FockDistDynamic2(const sc::Ref<sc::GaussianBasisSet> &gbs,
                                   const Ref<FockBlocks> &blocks,
                                   const sc::Ref<sc::PetiteList> &pl,
                                   const sc::Ref<sc::MessageGrp> &msg,
                                   int nthread, int mythread,
                                   const signed char *pmax,
                                   sc::Ref<sc::TwoBodyInt> eri,
                                   int l2tol):
    FockDistDynamic(gbs,blocks,pl,msg,nthread,mythread,pmax,eri,l2tol)
{
  if (nproc_ < 2) {
      static_ = new FockDistStatic2(gbs,blocks,pl,msg,nthread,mythread);
    }
  ij_iter_ = ijmap_.begin();
  kl_iter_ = klmap_.end();
}

FockDistDynamic2::~FockDistDynamic2()
{
}

bool
FockDistDynamic2::next_block()
{
  ij_iter_++;
  if (ij_iter_ != ijmap_.end()) return true;

  return false;
}

bool
FockDistDynamic2::server_get_blocks(int&i,int&j)
{
  if (ij_iter_ == ijmap_.end()) return false;

  i = ij_iter_->second.first;
  j = ij_iter_->second.second;

  next_block();

  return true;
}

void
FockDistDynamic2::run_server()
{
  int idat[2];
  int rproc;

  int i,j,k,l;
  while(server_get_blocks(i,j)) {
      idat[0] = i;
      idat[1] = j;
//       std::cout << "server waiting" << std::endl;
      msg_->recv(-1,rproc);
//       std::cout << "server sending " << i << " " << j << std::endl;
      msg_->send(rproc,idat,2);
    }

  // Stay around to send out the messages indicating we are done
  idat[0] = idat[1] = blocks_->nblock();
//   std::cout << "nproc = " << nproc_ << " nthread = " << nthread_ << std::endl;
  for (int i=1; i<nproc_; i++) {
      for (int j=0; j<nthread_; j++) {
//           std::cout << "in shutdown: waiting for request" << std::endl;
          msg_->recv(-1,rproc);
//           std::cout << "in shutdown: sending shutdown to " << rproc << std::endl;
          msg_->send(rproc,idat,2);
        }
    }
}

bool
FockDistDynamic2::get_blocks(int&i,int&j,int&k,int&l)
{
  const int server = FOCKDIST_SERVER;

  if (static_distribution()) {
      bool ret = static_->get_blocks(i,j,k,l);
//       if (ret) {
//           std::cout << "FockDistDynamic4::get_blocks: statically processing "
//                     << i << " " << j << " " << k << " " << l
//                     << std::endl;
//         }
//       else {
//           std::cout << "FockDistDynamic4::get_blocks: "
//                     << "no more blocks in static processing"
//                     << std::endl;
//         }
      return ret;
    }

  if (myproc_ == server) {
//       std::cout << "Running server"
//                 << std::endl;
      // Only one thread runs the server.
      if (mythread_ == 0) run_server();
      return false;
    }

  if (kl_iter_ == klmap_.end()) {
      int idat[2];
      msg_->send(server,myproc_);
      msg_->recv(server,idat,2);
      i_ = idat[0];
      j_ = idat[1];
//       std::cout << "FockDistDynamic2::get_blocks:"
//                 << " node " << myproc_
//                 << " thread " << mythread_
//                 << " : received "
//                 << i_ << " " << j_
//                 << " size = " << size(i_)*size(j_)
//                 << std::endl;
      if (i_ >= blocks_->nblock()) {
          return false;
        }
      kl_iter_ = klmap_.begin();
    }

  i = i_;
  j = j_;
  k = kl_iter_->second.first;
  l = kl_iter_->second.second;
  kl_iter_++;

//   std::cout << "FockDistDynamic2::get_blocks:"
//             << " node " << myproc_
//             << " thread " << mythread_
//             << " : processing "
//             << i << " " << j << " " << k << " " << l
//             << " size = " << size(i)*size(j) << " * " << size(k)*size(l)
//             << std::endl;

  return true;
}

void
FockDistDynamic2::init()
{
  ij_iter_ = ijmap_.begin();
  kl_iter_ = klmap_.end();

  if (static_distribution()) static_->init();
}

///////////////////////////////////////////////////////////////////

FockDistDynamic4::FockDistDynamic4(const sc::Ref<sc::GaussianBasisSet> &gbs,
                                   const Ref<FockBlocks> &blocks,
                                   const sc::Ref<sc::PetiteList> &pl,
                                   const sc::Ref<sc::MessageGrp> &msg,
                                   int nthread, int mythread,
                                   const signed char *pmax,
                                   sc::Ref<sc::TwoBodyInt> eri,
                                   int l2tol):
    FockDistDynamic(gbs,blocks,pl,msg,nthread,mythread,pmax,eri,l2tol)
{
  requested_work_ = false;

  if (nproc_ < 2) {
      static_ = new FockDistStatic4(gbs,blocks,pl,msg,nthread,mythread);
    }

  // place the ij and kl maps into a vector for random access
  nij_ = ijmap_.size();
  nkl_ = klmap_.size();
  ijvec_.resize(nij_);
  klvec_.resize(nkl_);
  {
    int i = 0;
    for (pairmap_t::iterator iter=ijmap_.begin();
         iter != ijmap_.end();
         iter++, i++) {
        ijvec_[i] = iter->second;
      }
  }
  {
    int i = 0;
    for (pairmap_t::iterator iter=klmap_.begin();
         iter != klmap_.end();
         iter++, i++) {
        klvec_[i] = iter->second;
      }
  }

  if (nij_ != 0 && nkl_ != 0) {
      a_end_ = 1;
      while (a_end_ < nij_) a_end_ *= 2;
      b_end_ = nkl_;
    }
  else {
      a_end_ = 0;
      b_end_ = 0;
    }

  a_ = 0;
  b_ = 0;
}

FockDistDynamic4::~FockDistDynamic4()
{
}

bool
FockDistDynamic4::next_block()
{
  int i,j,k,l;
  do {
      b_++;
      if (b_ >= b_end_) {
          b_ = 0;
          a_++;
          if (a_ >= a_end_) return false;
        }
    } while (!ijkl(i,j,k,l));

  return true;
}

bool
FockDistDynamic4::server_get_blocks(int&i,int&j,int&k,int&l)
{
  if (a_ >= a_end_) return false;

  ijkl(i,j,k,l);

  next_block();

  return true;
}

bool
FockDistDynamic4::ijkl(int &i, int &j, int &k, int &l)
{
  int ij = a_ ^ b_;
  int kl = b_;

  if (ij >= ijvec_.size()) return false;
  if (kl >= klvec_.size()) return false;

  i = ijvec_[ij].first;
  j = ijvec_[ij].second;
  k = klvec_[kl].first;
  l = klvec_[kl].second;

// These first two are satisfied by construction.
//   if (i<j) return false;
//   if (k<l) return false;
  if (i<k) return false;
  if (i==k && size(i) == 1 && j<l) return false;
  return true;
}

void
FockDistDynamic4::run_server()
{
  int idat[4];
  int rproc;

  int i,j,k,l;
  while(server_get_blocks(i,j,k,l)) {
      idat[0] = i;
      idat[1] = j;
      idat[2] = k;
      idat[3] = l;
//       std::cout << "sending out work unit "
//                 << i << " " << j << " " << k << " " << l
//                 << std::endl;
      msg_->recvt(-1,0,rproc);
      msg_->sendt(rproc,0,idat,4,true);
    }

  // Stay around to send out the messages indicating we are done
  idat[0] = idat[1] = idat[2] = idat[3] = blocks_->nblock();
//   std::cout << "nproc = " << nproc_ << " nthread = " << nthread_ << std::endl;
  for (int i=1; i<nproc_; i++) {
      for (int j=0; j<nthread_; j++) {
//           std::cout << "in shutdown: waiting for request" << std::endl;
          msg_->recvt(-1,0,rproc);
//           std::cout << "in shutdown: sending shutdown to " << rproc << std::endl;
          msg_->sendt(rproc,0,idat,4,true);
        }
    }
}

bool
FockDistDynamic4::get_blocks(int&i,int&j,int&k,int&l)
{
  const int server = FOCKDIST_SERVER;

  if (static_distribution()) {
      bool ret = static_->get_blocks(i,j,k,l);
//       if (ret) {
//           std::cout << "FockDistDynamic4::get_blocks: statically processing "
//                     << i << " " << j << " " << k << " " << l
//                     << std::endl;
//         }
//       else {
//           std::cout << "FockDistDynamic4::get_blocks: "
//                     << "no more blocks in static processing"
//                     << std::endl;
//         }
      return ret;
    }

  if (myproc_ == server) {
//       std::cout << "Running server"
//                 << std::endl;
      // Only one thread runs the server.
      if (mythread_ == 0) run_server();
      return false;
    }

  if (!requested_work_) {
      // if we haven't previously requested work,
      // then issue a request here
      msg_->nb_recvt(server,0,work_,4,recv_handle_);
      msg_->nb_sendt(server,0,myproc_,send_handle_);
    }

  // wait on the previously requested work
  msg_->wait(send_handle_);
  msg_->wait(recv_handle_);

  if (work_[0] >= blocks_->nblock()) {
//       std::cout << "FockDistDynamic4::get_blocks:"
//                 << " node " << myproc_
//                 << " thread " << mythread_
//                 << " nnode " << nproc_
//                 << " nthread " << nthread_
//                 << " : no more blocks"
//                 << std::endl;
      requested_work_ = false;
      return false;
    }


  i = work_[0];
  j = work_[1];
  k = work_[2];
  l = work_[3];

  msg_->nb_recvt(server,0,work_,4,recv_handle_);
  msg_->nb_sendt(server,0,myproc_,send_handle_);
  requested_work_ = true;

  // issue the request for the next piece of work

//   std::cout << "FockDistDynamic4::get_blocks:"
//             << " node " << myproc_
//             << " thread " << mythread_
//             << " : processing "
//             << i << " " << j << " " << k << " " << l
//             << " size = " << size(i)*size(j) << " * " << size(k)*size(l)
//             << std::endl;

  return true;
}

void
FockDistDynamic4::init()
{
  a_ = 0;
  b_ = 0;

  if (static_distribution()) static_->init();
}

///////////////////////////////////////////////////////////////////
