//
// thread_wrap.h
//
// Copyright (C) 2014 David Hollman
//
// Author: David Hollman
// Maintainer: DSH
// Created: Feb 24, 2014
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#ifndef _util_container_thread_wrap_h
#define _util_container_thread_wrap_h

#include <vector>
#include <algorithm>
#include <iterator>

namespace sc {

template <typename STLContainer>
STLContainer merge(const std::vector<STLContainer>& containers) {
  STLContainer rv;
  for(auto&& c : containers) {
    std::copy(c.begin(), c.end(), std::back_inserter(rv));
  }
  return rv;
}

namespace detail {

  template <typename Container>
  struct _Merger {
      Container operator()(const std::vector<Container>& containers) const {
        return merge(containers);
      }
  };

}


template <typename Container, typename Merger=detail::_Merger<Container>>
class ThreadReplicated
{

    int nthread_;
    std::vector<Container> thread_containers_;
    Merger merger_;

  public:

    ThreadReplicated()
      : nthread_(0),
        thread_containers_(0)
    { }

    ThreadReplicated(int nthread)
      : nthread_(nthread),
        thread_containers_(nthread)
    { }

    void set_nthread(int nthr) {
      nthread_ = nthr;
      thread_containers_.resize(nthr);
    }

    Container& mine(int ithr) {
      return thread_containers_[ithr];
    }

    operator Container() const {
      return merger_(thread_containers_);
    }

    Container merged() const {
      return merger_(thread_containers_);
    }


};





}

#endif /* _util_container_thread_wrap_h */
