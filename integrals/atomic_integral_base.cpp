//
// Created by Chong Peng on 3/2/16.
//

#include "atomic_integral_base.h"
namespace mpqc{
namespace integrals{

libint2::OneBodyEngine  AtomicIntegralBase::get_one_body_engine(const Operation &operation, int64_t max_nprim, int64_t max_am){

    libint2::OneBodyEngine::operator_type itype;
    std::vector<std::pair<double, std::array<double, 3>>> q;
    if (operation.oper() == Operation::Operations::Overlap) {
        itype = libint2::OneBodyEngine::overlap;
    } else if (operation.oper() == Operation::Operations::Kinetic) {
        itype = libint2::OneBodyEngine::kinetic;
    } else if (operation.oper() == Operation::Operations::Nuclear) {
        itype = libint2::OneBodyEngine::nuclear;
        q = make_q(*(this->mol_));
    } else {
        throw std::runtime_error("Invalid One Body Operation");
    }

    libint2::OneBodyEngine engine(itype, max_nprim, static_cast<int>(max_am),0);


    if(itype == libint2::OneBodyEngine::nuclear){
        engine.set_params(std::move(q));
    }

    return engine;
}

libint2::MultiplicativeSphericalTwoBodyKernel AtomicIntegralBase::get_two_body_engine_kernel(
        const Operation &operation)
{
    libint2::MultiplicativeSphericalTwoBodyKernel kernel;
    if (operation.oper() == Operation::Operations::Coulomb) {
        kernel = libint2::Coulomb;
    }
    else if(operation.oper() == Operation::Operations::cGTGCoulomb) {
        kernel = libint2::cGTG_times_Coulomb;
    }
    else if(operation.oper() == Operation::Operations::cGTG){
        kernel = libint2::cGTG;
    }
    else if(operation.oper() == Operation::Operations::cGTG2){
        kernel = libint2::cGTG;
    }
    else if(operation.oper() == Operation::Operations::DelcGTG2){
        kernel = libint2::DelcGTG_square;
    }
    else {
        throw std::runtime_error("Invalid Two Body Operation");
    }
    return kernel;
}


void AtomicIntegralBase::parse_one_body(const Formula& formula, std::shared_ptr<EnginePool<libint2::OneBodyEngine>>& engine_pool, Barray<2>& bases){

    auto bra_indexs = formula.left_index();
    auto ket_indexs = formula.right_index();

    TA_ASSERT(bra_indexs.size() == 1);
    TA_ASSERT(ket_indexs.size() == 1);

    auto bra_index = bra_indexs[0];
    auto ket_index = ket_indexs[0];

    TA_ASSERT(bra_index.is_ao());
    TA_ASSERT(ket_index.is_ao());

    auto bra_basis = this->index_to_basis(bra_index);
    auto ket_basis = this->index_to_basis(ket_index);

    TA_ASSERT(bra_basis != nullptr);
    TA_ASSERT(ket_basis != nullptr);

    auto max_nprim = std::max(bra_basis->max_nprim(), ket_basis->max_nprim());
    auto max_am = std::max(bra_basis->max_am(), ket_basis->max_am());
    bases = mpqc::utility::make_array(*bra_basis, *ket_basis);

    // convert operation to libint operator
    auto operation = formula.operation();
    engine_pool = make_pool(get_one_body_engine(operation,max_nprim,max_am));
}

void AtomicIntegralBase::parse_two_body_two_center(const Formula& formula, libint2::MultiplicativeSphericalTwoBodyKernel& kernel, Barray<2>& bases, int64_t& max_nprim, int64_t& max_am){

    TA_USER_ASSERT(formula.notation() == Formula::Notation::Chemical, "Two Body Two Center Integral Must Use Chemical Notation");

    auto bra_indexs = formula.left_index();
    auto ket_indexs = formula.right_index();

    TA_ASSERT(bra_indexs.size() == 1);
    TA_ASSERT(ket_indexs.size() == 1);

    auto bra_index0 = bra_indexs[0];
    auto ket_index0 = ket_indexs[0];

    TA_ASSERT(ket_index0.is_ao());
    TA_ASSERT(bra_index0.is_ao());

    auto bra_basis0 = this->index_to_basis(bra_index0);
    auto ket_basis0 = this->index_to_basis(ket_index0);

    TA_ASSERT(bra_basis0 != nullptr);
    TA_ASSERT(ket_basis0 != nullptr);

    max_nprim = std::max(bra_basis0->max_nprim(), ket_basis0->max_nprim());
    max_am = std::max(bra_basis0->max_am(), ket_basis0->max_am());
    bases = mpqc::utility::make_array(*bra_basis0, *ket_basis0);

    auto operation = formula.operation();
    kernel = get_two_body_engine_kernel(operation);

}

void AtomicIntegralBase::parse_two_body_three_center(const Formula& formula, libint2::MultiplicativeSphericalTwoBodyKernel& kernel, Barray<3>& bases, int64_t& max_nprim, int64_t& max_am){

    TA_USER_ASSERT(formula.notation() == Formula::Notation::Chemical, "Three Center Integral Must Use Chemical Notation");

    auto bra_indexs = formula.left_index();
    auto ket_indexs = formula.right_index();

    TA_ASSERT(bra_indexs.size() == 1);
    TA_ASSERT(ket_indexs.size() == 2);


    TA_ASSERT(bra_indexs[0].is_ao());
    TA_ASSERT(ket_indexs[0].is_ao());
    TA_ASSERT(ket_indexs[1].is_ao());

    auto bra_basis0 = this->index_to_basis(bra_indexs[0]);
    auto ket_basis0 = this->index_to_basis(ket_indexs[0]);
    auto ket_basis1 = this->index_to_basis(ket_indexs[1]);

    TA_ASSERT(bra_basis0 != nullptr);
    TA_ASSERT(ket_basis0 != nullptr);
    TA_ASSERT(ket_basis1 != nullptr);

    max_nprim = std::max({bra_basis0->max_nprim(), ket_basis0->max_nprim()
                                 ,ket_basis1->max_nprim()});
    max_am = std::max({bra_basis0->max_am(), ket_basis0->max_am(),
                       ket_basis1->max_am()});

    bases = mpqc::utility::make_array(*bra_basis0, *ket_basis0, *ket_basis1);

    // convert operation to libint operator
    auto operation = formula.operation();
    kernel = get_two_body_engine_kernel(operation);

}

void AtomicIntegralBase::parse_two_body_four_center(const Formula& formula, libint2::MultiplicativeSphericalTwoBodyKernel& kernel, Barray<4>& bases, int64_t& max_nprim, int64_t& max_am){
    auto bra_indexs = formula.left_index();
    auto ket_indexs = formula.right_index();

    TA_ASSERT(bra_indexs.size() == 2);
    TA_ASSERT(ket_indexs.size() == 2);


    TA_ASSERT(bra_indexs[0].is_ao());
    TA_ASSERT(ket_indexs[0].is_ao());
    TA_ASSERT(bra_indexs[1].is_ao());
    TA_ASSERT(ket_indexs[1].is_ao());

    auto bra_basis0 = this->index_to_basis(bra_indexs[0]);
    auto ket_basis0 = this->index_to_basis(ket_indexs[0]);
    auto bra_basis1 = this->index_to_basis(bra_indexs[1]);
    auto ket_basis1 = this->index_to_basis(ket_indexs[1]);

    TA_ASSERT(bra_basis0 != nullptr);
    TA_ASSERT(ket_basis0 != nullptr);
    TA_ASSERT(bra_basis1 != nullptr);
    TA_ASSERT(ket_basis1 != nullptr);

    max_nprim = std::max({bra_basis0->max_nprim(), bra_basis1->max_nprim()
                                 ,ket_basis0->max_nprim(), ket_basis1->max_nprim()});
    max_am = std::max({bra_basis0->max_am(), bra_basis1->max_am(),
                       ket_basis0->max_am(),ket_basis1->max_am()});

    if (formula.notation() == Formula::Notation::Chemical){
        bases = {*bra_basis0, *bra_basis1, *ket_basis0, *ket_basis1};
    }
    else if(formula.notation() == Formula::Notation::Physical){
        bases = {*bra_basis0, *ket_basis0, *bra_basis1, *ket_basis1};
    }

    // convert operation to libint operator
    auto operation = formula.operation();

    kernel = get_two_body_engine_kernel(operation);
}


std::array<std::wstring, 3> AtomicIntegralBase::get_df_formula(const Formula &formula) {
    std::array<std::wstring, 3> result;

    //chemical notation
    if (formula.notation() == Formula::Notation::Chemical) {

        std::wstring left =
                L"( Κ |" + formula.operation().oper_string() + L"| " + formula.left_index()[1].name() + L" " +
                formula.left_index()[0].name() + L" )";
        std::wstring right =
                L"( Κ |" + formula.operation().oper_string() + L"| " + formula.right_index()[0].name() + L" " +
                formula.right_index()[1].name() + L" )";
        std::wstring center = L"( Κ |" + formula.operation().oper_string() + L"| Λ)[inv]";
        result[0] = left;
        result[1] = center;
        result[2] = right;
    }
    //physical notation
    else {
        std::wstring left =
                L"( Κ |" + formula.operation().oper_string() + L"| " + formula.right_index()[0].name() + L" " +
                formula.left_index()[0].name() + L" )";
        std::wstring right =
                L"( Κ |" + formula.operation().oper_string() + L"| " + formula.left_index()[1].name() + L" " +
                formula.right_index()[1].name() + L" )";
        std::wstring center = L"( Κ |" + formula.operation().oper_string() + L"| Λ)[inv]";
        result[0] = left;
        result[1] = center;
        result[2] = right;
    }

    return result;
}


Formula AtomicIntegralBase::get_jk_formula(const Formula &formula) {

    Formula result;

    result.set_notation(Formula::Notation::Chemical);
    Operation oper(L"G");
    result.set_operation(oper);
    if(formula.operation().oper() == Operation::Operations::J){

        result.left_index().push_back(formula.left_index()[0]);
        result.left_index().push_back(formula.right_index()[0]);
        result.right_index().push_back(OrbitalIndex(L"μ"));
        result.right_index().push_back(OrbitalIndex(L"ν"));

    }
    else{

        result.left_index().push_back(formula.left_index()[0]);
        result.left_index().push_back(OrbitalIndex(L"μ"));
        result.right_index().push_back(formula.right_index()[0]);
        result.right_index().push_back(OrbitalIndex(L"ν"));

    }
    return result;
}

std::array<Formula,3> AtomicIntegralBase::get_jk_df_formula(const Formula &formula) {
    std::array<Formula,3> result;

    if(formula.operation().oper() == Operation::Operations::J){
        std::wstring left =  L"( Κ |G| " + formula.right_index()[0].name() + L" " + formula.left_index()[0].name() + L" )";
        std::wstring right = L"( Κ |G| μ ν  )";

        result[0] = Formula(left);
        result[2] = Formula(right);
    }
    else{
        std::wstring left =  L"( Κ |G|μ " + formula.left_index()[0].name() + L" )";
        std::wstring right = L"( Κ |G| " + formula.right_index()[0].name() + L" ν )";

        result[0] = Formula(left);
        result[2] = Formula(right);
    }
    std::wstring center = L"( Κ |G| Λ)[inv]";
    result[1] = Formula(center);

    return result;
}

OrbitalIndex AtomicIntegralBase::get_jk_orbital_space(const Operation &operation) {

    if(operation.oper() == Operation::Operations::J || operation.oper() == Operation::Operations::K){
        return OrbitalIndex(L"m");
    }
    else if(operation.oper() == Operation::Operations::KAlpha){
        return OrbitalIndex(L"m_α");
    }
    else if(operation.oper() == Operation::Operations::KBeta){
        return OrbitalIndex(L"m_β");
    }
}

std::array<Formula, 4> AtomicIntegralBase::get_fock_formula(const Formula &formula) {

    std::array<Formula, 4> result;
    auto v = formula;
    v.set_operation(Operation(L"V"));
    auto t = formula;
    t.set_operation(Operation(L"T"));
    auto j = formula;
    j.operation().set_oper(Operation::Operations::J);
    auto k = formula;
    if(formula.operation().oper() == Operation::Operations::Fock){
        k.operation().set_oper(Operation::Operations::K);
    }
    else if(formula.operation().oper() == Operation::Operations::FockAlpha){
        k.operation().set_oper(Operation::Operations::KAlpha);
    }
    else if(formula.operation().oper() == Operation::Operations::FockBeta){
        k.operation().set_oper( Operation::Operations::KBeta);
    }

    result[0] = v;
    result[1] = t;
    result[2] = j;
    result[3] = k;
    return result;
}
}// end of namespace integral
}// end of namespace mpqc
