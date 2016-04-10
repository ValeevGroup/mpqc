//
// Created by Chong Peng on 10/14/15.
//

#ifndef TILECLUSTERCHEM_ATOMIC_INTEGRAL_H
#define TILECLUSTERCHEM_ATOMIC_INTEGRAL_H


#include "../f12/f12_utility.h"
#include "../expression/formula_registry.h"
#include "../expression/orbital_registry.h"
#include "atomic_integral_base.h"
#include "../utility/wcout_utf8.h"
#include "../ta_routines/array_to_eigen.h"
#include "../utility/parallel_print.h"
#include "../utility/parallel_break_point.h"
#include "../utility/time.h"

namespace mpqc{
namespace integrals{

    /// Atomic Integral Class
    //// Op is a function
    /// Op will take TA::TensorD as argument and return Tile

    template<typename Tile, typename Policy>
    class AtomicIntegral : public AtomicIntegralBase{

    public:
        using TArray = TA::DistArray<Tile,Policy>;

        typedef Tile (*Op)(TA::TensorD &&);

        template<unsigned int N, typename E>
        using IntegralBuilder = integrals::IntegralBuilder<N,E,Op>;

        AtomicIntegral() = default;

        AtomicIntegral(madness::World& world,
                Op op,
                const std::shared_ptr<molecule::Molecule>& mol,
                const std::shared_ptr <OrbitalBasisRegistry>& obs,
                const std::vector<std::pair<double,double>>& gtg_params = std::vector<std::pair<double,double>>()
        ) : AtomicIntegralBase(world,mol,obs,gtg_params), op_(op), ao_formula_registry_(), orbital_space_registry_()
        {}

        AtomicIntegral(AtomicIntegral&& ) = default;
        AtomicIntegral& operator=(AtomicIntegral&& ) = default;

        virtual ~AtomicIntegral() = default;


        TArray compute(const std::wstring& );
        TArray compute(const Formula& );

        TA::expressions::TsrExpr<TArray,true> operator() (const std::wstring& str){
            auto formula = Formula(str);
            auto array = compute(formula);
            auto& result = ao_formula_registry_.retrieve(formula);
            return result(formula.to_ta_expression());
        };

        const FormulaRegistry<TArray> &registry() const {
            return ao_formula_registry_;
        }

        FormulaRegistry<TArray> &registry(){
            return ao_formula_registry_;
        }

        void set_orbital_space_registry(const std::shared_ptr<OrbitalSpaceRegistry<TArray>> regitsry) {
            orbital_space_registry_ = regitsry;
        }


    protected:

        // compute two body integral
        template <typename Basis>
        TArray compute_two_body_integral(const libint2::MultiplicativeSphericalTwoBodyKernel& kernel, const Basis& basis, int64_t max_nprim, int64_t max_am, const Operation& operation);

        // compute integrals that has two dimension
        TArray compute2(const Formula& formula_string);
        // compute integrals that has three dimension
        TArray compute3(const Formula& formula_stirng);
        // compute integrals that has four dimension
        TArray compute4(const Formula& formula_string);


    private:

        //TODO better inverse of two center
        //TODO direct integral
        //TODO Screener for different type of integral
        // compute integral for sparse policy
        template <typename E, unsigned long N, typename U = Policy>
        TA::Array<double, N,Tile,typename std::enable_if<std::is_same<U,TA::SparsePolicy>::value, TA::SparsePolicy>::type> compute_integrals(
                madness::World& world, E const &engine, Barray<N> const &bases)
        {
            auto sreener = std::make_shared<Screener>(Screener{});
            auto result = mpqc::integrals::sparse_integrals(world,engine,bases,sreener,op_);
            return result;
        }

        // compute integral for dense policy
        template <typename E, unsigned long N, typename U = Policy>
        TA::Array<double, N,Tile,typename std::enable_if<std::is_same<U,TA::DensePolicy>::value, TA::DensePolicy>::type> compute_integrals(
                madness::World& world, E const &engine, Barray<N> const &bases)
        {
            auto sreener = std::make_shared<Screener>(Screener{});
            auto result = mpqc::integrals::dense_integrals(world,engine,bases,sreener,op_);
            return result;
        }

    private:
        FormulaRegistry<TArray> ao_formula_registry_;
        std::shared_ptr<OrbitalSpaceRegistry<TArray>> orbital_space_registry_;
        Op op_;

    };


    template <typename Tile, typename Policy>
    template<typename Basis>
    typename AtomicIntegral<Tile, Policy>::TArray AtomicIntegral<Tile, Policy>::compute_two_body_integral(
            const libint2::MultiplicativeSphericalTwoBodyKernel& kernel, const Basis& bs_array, int64_t max_nprim,
            int64_t max_am, const Operation& operation) {

        typename AtomicIntegral<Tile,Policy>::TArray result;

        if(kernel == libint2::Coulomb){
            libint2::TwoBodyEngine<libint2::Coulomb> engine(max_nprim, static_cast<int>(max_am));
            auto engine_pool = make_pool(engine);
            result = compute_integrals(this->world_,engine_pool,bs_array);
        }
        else{
            if(this->gtg_params_.empty()){
                throw std::runtime_error("Gaussian Type Genminal Parameters are empty!");
            }

            if(kernel == libint2::cGTG){
                if(operation.oper() == Operation::Operations::cGTG2){

                    auto squared_pragmas = f12::gtg_params_squared(this->gtg_params_);
                    libint2::TwoBodyEngine<libint2::cGTG> engine(max_nprim, max_am,0,std::numeric_limits<double>::epsilon(),squared_pragmas);
                    auto engine_pool = make_pool(engine);
                    result = compute_integrals(this->world_,engine_pool,bs_array);

                }else{
                    libint2::TwoBodyEngine<libint2::cGTG> engine(max_nprim, max_am,0,std::numeric_limits<double>::epsilon(),this->gtg_params_);
                    auto engine_pool = make_pool(engine);
                    result = compute_integrals(this->world_,engine_pool,bs_array);
                }
            }
            else if(kernel == libint2::cGTG_times_Coulomb){
                libint2::TwoBodyEngine<libint2::cGTG_times_Coulomb> engine(max_nprim,max_am,0,std::numeric_limits<double>::epsilon(),this->gtg_params_);
                auto engine_pool = make_pool(engine);
                result = compute_integrals(this->world_,engine_pool,bs_array);
            }
            else if(kernel == libint2::DelcGTG_square){
                libint2::TwoBodyEngine<libint2::DelcGTG_square> engine(max_nprim, max_am,0,std::numeric_limits<double>::epsilon(),this->gtg_params_);
                auto engine_pool = make_pool(engine);
                result = compute_integrals(this->world_,engine_pool,bs_array);
            }

        }
        TA_ASSERT(result.is_initialized());
        return result;
    }

    template <typename Tile, typename Policy>
    typename AtomicIntegral<Tile,Policy>::TArray AtomicIntegral<Tile,Policy>::compute(const std::wstring& formula_string) {
        auto formula = Formula(formula_string);
        return compute(formula);
    }

    template <typename Tile, typename Policy>
    typename AtomicIntegral<Tile,Policy>::TArray AtomicIntegral<Tile,Policy>::compute(const Formula& formula) {

        TArray result;

        auto iter = ao_formula_registry_.find(formula);

        if(iter != ao_formula_registry_.end()){
            result = iter->second;
            utility::print_par(world_,"Retrived AO Integral: ");
            utility::wprint_par(world_, formula.formula_string());
            double size = utility::array_size(result);
            utility::print_par(world_," Size: ", size, " GB\n");
        }else{

            if(formula.rank() == 2){
                result =  compute2(formula);
                ao_formula_registry_.insert(formula, result);
            }
            else if(formula.rank() == 3){
                result =  compute3(formula);
                ao_formula_registry_.insert(formula, result);
            }
            else if(formula.rank() == 4){
                result =  compute4(formula);
                ao_formula_registry_.insert(formula, result);
            }
        }

        // wait all process to obtain and insert result
//        world_.gop.fence();
        return result;

    }

    template <typename Tile, typename Policy>
    typename AtomicIntegral<Tile,Policy>::TArray AtomicIntegral<Tile,Policy>::compute2(const Formula& formula) {

        Barray<2> bs_array;
        double time = 0.0;
        TArray result;

        // use one body engine
        if(formula.operation().is_onebody()){

            // H = V + T
            if(formula.operation().oper() == Operation::Operations::Core){


                auto v_formula = formula;
                v_formula.operation().set_oper(Operation::Operations::Nuclear);

                auto t_formula = formula;
                t_formula.operation().set_oper(Operation::Operations::Kinetic);

                auto v = this->compute(v_formula);
                auto t = this->compute(t_formula);

                auto time0 = mpqc_time::fenced_now(world_);

                result("i,j") = v("i,j") + t("i,j");

                auto time1 = mpqc_time::fenced_now(world_);
                time+= mpqc_time::duration_in_s(time0,time1);
            }
            else {
                auto time0 = mpqc_time::fenced_now(world_);

                std::shared_ptr<EnginePool<libint2::OneBodyEngine>> engine_pool;
                parse_one_body(formula,engine_pool,bs_array);
                result = compute_integrals(this->world_,engine_pool,bs_array);

                auto time1 = mpqc_time::fenced_now(world_);
                time+= mpqc_time::duration_in_s(time0,time1);
            }
            utility::print_par(world_,"Computed One Body Integral: ");
            utility::wprint_par(world_,formula.formula_string());
            double size = utility::array_size(result);
            utility::print_par(world_," Size: ", size, " GB");
            utility::print_par(world_," Time: ", time, " s\n");
        }
        // use two body engine
        else if(formula.operation().is_twobody()){
            auto time0 = mpqc_time::fenced_now(world_);

            int64_t max_nprim, max_am;
            libint2::MultiplicativeSphericalTwoBodyKernel kernel;

            parse_two_body_two_center(formula,kernel,bs_array,max_nprim,max_am);
            result = compute_two_body_integral( kernel, bs_array, max_nprim, max_am, formula.operation());

            if(formula.operation().has_option(Operation::Options::Inverse)){

                if(formula.operation().oper() == Operation::Operations::cGTG
                    ||formula.operation().oper() == Operation::Operations::cGTGCoulomb){
                    result("i,j") = -result("i,j");
                }

//                utility::parallel_break_point(world_,0);
//                std::cout << "Before Array To Eigen" << std::endl;
//                std::cout << result << std::endl;
                MatrixD result_eig = array_ops::array_to_eigen(result);

                // compute cholesky decomposition
                auto llt_solver = Eig::LLT<MatrixD>(result_eig);

                // check success
                Eigen::ComputationInfo info = llt_solver.info();
                if(info == Eigen::ComputationInfo::NumericalIssue){
                    utility::print_par(world_,"NumericalIssue","\n");
                }
                else if(info == Eigen::ComputationInfo::NoConvergence){
                    utility::print_par(world_,"NoConvergence","\n");
                }
                else if (info == Eigen::ComputationInfo::InvalidInput){
                    utility::print_par(world_,"InvalidInput","\n");
                }


                // print the eigen value if LLT failed
                if(info != Eigen::ComputationInfo::Success) {
                    if(world_.rank() == 0){
                        Eigen::SelfAdjointEigenSolver<MatrixD> ei_solver(result_eig);
                        std::cout << ei_solver.eigenvalues() << std::endl;
                    }
                }
                TA_ASSERT( info == Eigen::ComputationInfo::Success);

                MatrixD L = MatrixD(llt_solver.matrixL());
                MatrixD L_inv_eig = L.inverse();
                result_eig = L_inv_eig.transpose() * L_inv_eig;

                auto tr_result = result.trange().data()[0];
                result = array_ops::eigen_to_array<TA::TensorD>(result.get_world(), result_eig, tr_result, tr_result);

                if(formula.operation().oper() == Operation::Operations::cGTG
                    ||formula.operation().oper() == Operation::Operations::cGTGCoulomb){
                    result("i,j") = -result("i,j");
                }

            }

            if(formula.operation().has_option(Operation::Options::InverseSquareRoot)){
                if(formula.operation().oper() == Operation::Operations::cGTG ||
                        formula.operation().oper() == Operation::Operations::cGTGCoulomb){
                    result("i,j") = -result("i,j");
                }

                auto result_eig = array_ops::array_to_eigen(result);
                MatrixD L_inv_eig = MatrixD(Eig::LLT<MatrixD>(result_eig).matrixL()).inverse();
                auto tr_result = result.trange().data()[0];
                result = array_ops::eigen_to_array<TA::TensorD>(result.get_world(), L_inv_eig, tr_result, tr_result);

                if(formula.operation().oper() == Operation::Operations::cGTG ||
                        formula.operation().oper() == Operation::Operations::cGTGCoulomb){
                    result("i,j") = -result("i,j");
                }
            }

            auto time1 = mpqc_time::fenced_now(world_);
            time+= mpqc_time::duration_in_s(time0,time1);

            utility::print_par(world_,"Computed Twobody Two Center Integral: ");
            utility::wprint_par(world_, formula.formula_string());
            double size = utility::array_size(result);
            utility::print_par(world_," Size: ", size, " GB");
            utility::print_par(world_," Time: ", time, " s\n");
        }
            //compute JK, requires orbital space registry
        else if(formula.operation().is_jk()){

            // density fitting case
            if(formula.operation().has_option(Operation::Options::DensityFitting)){
                auto three_center_formula = get_jk_df_formula(formula);

                auto left = compute(three_center_formula[0]);
                auto center = compute(three_center_formula[1]);
                auto right = compute(three_center_formula[2]);

                auto time0 = mpqc_time::fenced_now(world_);
                // find the density
                auto space_index = get_jk_orbital_space(formula.operation());
                auto space = orbital_space_registry_->retrieve(space_index);

                // J case
                if(formula.operation().oper() == Operation::Operations::J){
                    result("i,j") = center("K,Q")*right("Q,k,l")*(space("k,a")*space("l,a"))*left("K,i,j");
                }
                // K case
                else{
                    result("i,j") = (left("K,i,k")*space("k,a"))*center("K,Q")*(right("Q,j,l")*space("l,a"));
                }
                auto time1 = mpqc_time::fenced_now(world_);
                time+= mpqc_time::duration_in_s(time0,time1);
            }
                // four center case
            else{
                // convert to ao formula
                auto four_center_formula = get_jk_formula(formula);
                auto four_center = this->compute(four_center_formula);

                auto time0 = mpqc_time::fenced_now(world_);
                // find the density
                auto space_index = get_jk_orbital_space(formula.operation());
                auto space = orbital_space_registry_->retrieve(space_index);

                if(formula.operation().oper() == Operation::Operations::J){
                    result("rho,sigma") = four_center("rho,sigma,mu,nu")*(space("mu,i")*space("nu,i"));
                }
                else{
                    result("rho,sigma") = four_center("rho,mu,sigma,nu")*(space("mu,i")*space("nu,i"));
                }
                auto time1 = mpqc_time::fenced_now(world_);
                time+= mpqc_time::duration_in_s(time0,time1);

            }
            utility::print_par(world_,"Computed Coulumb/Exchange Integral: ");
            utility::wprint_par(world_, formula.formula_string());
            double size = utility::array_size(result);
            utility::print_par(world_," Size: ", size, " GB");
            utility::print_par(world_," Time: ", time, " s\n");

        }
            // hJ = H + J
        else if(formula.operation().oper() == Operation::Operations::hJ){
            auto h_formula = formula;
            h_formula.set_operation(Operation(L"H"));


            auto j_formula = formula;
            j_formula.operation().set_oper(Operation::Operations::J);

            auto h = this->compute(h_formula);
            auto j = this->compute(j_formula);

            auto time0 = mpqc_time::fenced_now(world_);

            result("i,j") = h("i,j") + 2*j("i,j");

            auto time1 = mpqc_time::fenced_now(world_);
            time+= mpqc_time::duration_in_s(time0,time1);

            utility::print_par(world_,"Computed Coulumb/Exchange Integral: ");
            utility::wprint_par(world_, formula.formula_string());
            double size = utility::array_size(result);
            utility::print_par(world_," Size: ", size, " GB");
            utility::print_par(world_," Time: ", time, " s\n");
        }
        //compute Fock, requires orbital space registry
        else if(formula.operation().is_fock()){

            auto formulas = get_fock_formula(formula);

            auto h = compute(formulas[0]);
            auto j = compute(formulas[1]);
            auto k = compute(formulas[2]);

            auto time0 = mpqc_time::fenced_now(world_);
            // if closed shell
            if(formula.operation().oper() == Operation::Operations::Fock){
                result("rho,sigma") = h("rho,sigma") + 2*j("rho,sigma") - k("rho,sigma");
            }
            // else if spin orbital
            else{
                result("rho,sigma") = h("rho,sigma") + j("rho,sigma") - k("rho,sigma");
            }
            auto time1 = mpqc_time::fenced_now(world_);
            time+= mpqc_time::duration_in_s(time0,time1);

            utility::print_par(world_,"Computed Fock Integral: ");
            utility::wprint_par(world_, formula.formula_string());
            double size = utility::array_size(result);
            utility::print_par(world_," Size: ", size, " GB");
            utility::print_par(world_," Time: ", time, " s\n");
        }


        return result;

    }

    template <typename Tile, typename Policy>
    typename AtomicIntegral<Tile,Policy>::TArray AtomicIntegral<Tile,Policy>::compute3(const Formula& formula) {

        double time;
        TArray result;

        libint2::MultiplicativeSphericalTwoBodyKernel kernel;
        Barray<3> bs_array;

        auto time0 = mpqc_time::fenced_now(world_);
        int64_t max_nprim, max_am;
        parse_two_body_three_center(formula,kernel,bs_array,max_nprim,max_am);

        result = compute_two_body_integral( kernel, bs_array, max_nprim, max_am, formula.operation());

        auto time1 = mpqc_time::fenced_now(world_);
        time+= mpqc_time::duration_in_s(time0,time1);

        utility::print_par(world_,"Computed Twobody Three Center Integral: ");
        utility::wprint_par(world_, formula.formula_string());
        double size = utility::array_size(result);
        utility::print_par(world_," Size: ", size, " GB");
        utility::print_par(world_," Time: ", time, " s\n");

        return result;
    }

    template <typename Tile, typename Policy>
    typename AtomicIntegral<Tile,Policy>::TArray AtomicIntegral<Tile,Policy>::compute4(const Formula& formula) {
        double time = 0.0;
        TArray result;

        if(formula.operation().has_option(Operation::Options::DensityFitting)){

            // convert formula to df formula
            auto formula_strings = get_df_formula(formula);

            // compute integral
            auto left = compute(formula_strings[0]);
            auto center = compute(formula_strings[1]);
            auto right = compute(formula_strings[2]);

            auto time0 = mpqc_time::fenced_now(world_);

            if(formula.notation() == Formula::Notation::Chemical){
                result("i,j,k,l") = left("q,i,j")*center("q,p")*right("p,k,l");
            }else{
                result("i,j,k,l") = left("q,i,k")*center("q,p")*right("p,j,l");
            }
            auto time1 = mpqc_time::fenced_now(world_);
            time+= mpqc_time::duration_in_s(time0,time1);

            utility::print_par(world_,"Computed Twobody Four Center Density-Fitting Integral: ");
            utility::wprint_par(world_, formula.formula_string());
            double size = utility::array_size(result);
            utility::print_par(world_," Size: ", size, " GB");
            utility::print_par(world_," Time: ", time, " s\n");


        }
        else{
            auto time0 = mpqc_time::fenced_now(world_);

            int64_t max_nprim, max_am;
            libint2::MultiplicativeSphericalTwoBodyKernel kernel;
            Barray<4> bs_array;

            parse_two_body_four_center(formula,kernel,bs_array,max_nprim,max_am);

            result = compute_two_body_integral( kernel, bs_array, max_nprim, max_am, formula.operation());

            //TODO handle permutation better

            if(formula.notation() == Formula::Notation::Physical){
                result("i,j,k,l") = result("i,k,j,l");
            }

            auto time1 = mpqc_time::fenced_now(world_);
            time+= mpqc_time::duration_in_s(time0,time1);

            utility::print_par(world_,"Computed Twobody Four Center Integral: ");
            utility::wprint_par(world_, formula.formula_string());
            double size = utility::array_size(result);
            utility::print_par(world_," Size: ", size, " GB");
            utility::print_par(world_," Time: ", time, " s\n");
        }
        return result;

    }

}
}


#endif //TILECLUSTERCHEM_ATOMIC_INTEGRAL_H
