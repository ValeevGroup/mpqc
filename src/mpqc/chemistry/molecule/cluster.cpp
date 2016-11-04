#include "mpqc/chemistry/molecule/cluster.h"

namespace mpqc {

void Cluster::update_center(){
    center_ = Vector3d{0,0,0};

    for(auto const &elem : elements_){
        center_ += elem.center();
    }

    center_ /= elements_.size();
}

double Cluster::sum_distances_from_center() const {
    assert(false); // I am not sure this is correct.
    return 0.0;
    // auto reduce_r = [&](double d, const Clusterable &c) {
    //     return d + std::sqrt(diff_squaredNorm(c.center(), center_));
    // };

    // return std::accumulate(elements_.begin(), elements_.end(), 0.0, reduce_r);
}

std::ostream & operator<<(std::ostream &os, Cluster const &c){
    const auto end = c.end();
    const auto last = end - 1;
    os << "Cluster: {";
    os << "Center: " << c.center().transpose() << ", elements: {";
    for(auto i = c.begin(); i != end; ++i){
        if(i != last){
            i->print(os) << ", ";
        } else {
            i->print(os) << "}}";
        }
    }

    return os;
}

} // namespace mpqc
