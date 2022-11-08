// Generated by rstantools.  Do not edit by hand.

/*
    menbayes is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    menbayes is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with menbayes.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.21.0
#include <stan/model/model_header.hpp>
namespace model_f1_g4_namespace {
using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;
static int current_statement_begin__;
stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "model_f1_g4");
    reader.add_event(83, 81, "end", "model_f1_g4");
    return reader;
}
template <typename T0__>
std::vector<typename boost::math::tools::promote_args<T0__>::type>
segfreq4(const T0__& alpha,
             const int& g, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T0__>::type local_scalar_t__;
    typedef local_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
    int current_statement_begin__ = -1;
    try {
        {
        current_statement_begin__ = 9;
        validate_non_negative_index("p", "3", 3);
        std::vector<local_scalar_t__  > p(3, local_scalar_t__(DUMMY_VAR__));
        stan::math::initialize(p, DUMMY_VAR__);
        stan::math::fill(p, DUMMY_VAR__);
        current_statement_begin__ = 10;
        if (as_bool(logical_eq(g, 0))) {
            current_statement_begin__ = 11;
            stan::model::assign(p, 
                        stan::model::cons_list(stan::model::index_uni(1), stan::model::nil_index_list()), 
                        1.0, 
                        "assigning variable p");
            current_statement_begin__ = 12;
            stan::model::assign(p, 
                        stan::model::cons_list(stan::model::index_uni(2), stan::model::nil_index_list()), 
                        0.0, 
                        "assigning variable p");
            current_statement_begin__ = 13;
            stan::model::assign(p, 
                        stan::model::cons_list(stan::model::index_uni(3), stan::model::nil_index_list()), 
                        0.0, 
                        "assigning variable p");
        } else if (as_bool(logical_eq(g, 1))) {
            current_statement_begin__ = 15;
            stan::model::assign(p, 
                        stan::model::cons_list(stan::model::index_uni(1), stan::model::nil_index_list()), 
                        ((2.0 + alpha) / 4.0), 
                        "assigning variable p");
            current_statement_begin__ = 16;
            stan::model::assign(p, 
                        stan::model::cons_list(stan::model::index_uni(2), stan::model::nil_index_list()), 
                        ((2.0 * (1.0 - alpha)) / 4.0), 
                        "assigning variable p");
            current_statement_begin__ = 17;
            stan::model::assign(p, 
                        stan::model::cons_list(stan::model::index_uni(3), stan::model::nil_index_list()), 
                        (alpha / 4.0), 
                        "assigning variable p");
        } else if (as_bool(logical_eq(g, 2))) {
            current_statement_begin__ = 19;
            stan::model::assign(p, 
                        stan::model::cons_list(stan::model::index_uni(1), stan::model::nil_index_list()), 
                        ((1.0 + (2.0 * alpha)) / 6.0), 
                        "assigning variable p");
            current_statement_begin__ = 20;
            stan::model::assign(p, 
                        stan::model::cons_list(stan::model::index_uni(2), stan::model::nil_index_list()), 
                        ((4.0 * (1.0 - alpha)) / 6.0), 
                        "assigning variable p");
            current_statement_begin__ = 21;
            stan::model::assign(p, 
                        stan::model::cons_list(stan::model::index_uni(3), stan::model::nil_index_list()), 
                        ((1.0 + (2.0 * alpha)) / 6.0), 
                        "assigning variable p");
        } else if (as_bool(logical_eq(g, 3))) {
            current_statement_begin__ = 23;
            stan::model::assign(p, 
                        stan::model::cons_list(stan::model::index_uni(1), stan::model::nil_index_list()), 
                        (alpha / 4.0), 
                        "assigning variable p");
            current_statement_begin__ = 24;
            stan::model::assign(p, 
                        stan::model::cons_list(stan::model::index_uni(2), stan::model::nil_index_list()), 
                        ((2.0 * (1.0 - alpha)) / 4.0), 
                        "assigning variable p");
            current_statement_begin__ = 25;
            stan::model::assign(p, 
                        stan::model::cons_list(stan::model::index_uni(3), stan::model::nil_index_list()), 
                        ((2.0 + alpha) / 4.0), 
                        "assigning variable p");
        } else if (as_bool(logical_eq(g, 4))) {
            current_statement_begin__ = 27;
            stan::model::assign(p, 
                        stan::model::cons_list(stan::model::index_uni(1), stan::model::nil_index_list()), 
                        0.0, 
                        "assigning variable p");
            current_statement_begin__ = 28;
            stan::model::assign(p, 
                        stan::model::cons_list(stan::model::index_uni(2), stan::model::nil_index_list()), 
                        0.0, 
                        "assigning variable p");
            current_statement_begin__ = 29;
            stan::model::assign(p, 
                        stan::model::cons_list(stan::model::index_uni(3), stan::model::nil_index_list()), 
                        1.0, 
                        "assigning variable p");
        } else {
        }
        current_statement_begin__ = 33;
        return stan::math::promote_scalar<fun_return_scalar_t__>(p);
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}
struct segfreq4_functor__ {
    template <typename T0__>
        std::vector<typename boost::math::tools::promote_args<T0__>::type>
    operator()(const T0__& alpha,
             const int& g, std::ostream* pstream__) const {
        return segfreq4(alpha, g, pstream__);
    }
};
template <typename T0__, typename T1__>
Eigen::Matrix<typename boost::math::tools::promote_args<T0__, T1__>::type, Eigen::Dynamic, 1>
convolve(const std::vector<T0__>& p1,
             const std::vector<T1__>& p2,
             const int& K,
             const int& khalf, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T0__, T1__>::type local_scalar_t__;
    typedef local_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
    int current_statement_begin__ = -1;
    try {
        {
        current_statement_begin__ = 41;
        validate_non_negative_index("q", "(K + 1)", (K + 1));
        Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> q((K + 1));
        stan::math::initialize(q, DUMMY_VAR__);
        stan::math::fill(q, DUMMY_VAR__);
        current_statement_begin__ = 42;
        for (int k = 1; k <= (K + 1); ++k) {
            {
            current_statement_begin__ = 43;
            int iup(0);
            (void) iup;  // dummy to suppress unused var warning
            stan::math::fill(iup, std::numeric_limits<int>::min());
            stan::math::assign(iup,std::min((k - 1), (khalf - 1)));
            current_statement_begin__ = 44;
            int ilo(0);
            (void) ilo;  // dummy to suppress unused var warning
            stan::math::fill(ilo, std::numeric_limits<int>::min());
            stan::math::assign(ilo,std::max(0, (k - khalf)));
            current_statement_begin__ = 45;
            stan::model::assign(q, 
                        stan::model::cons_list(stan::model::index_uni(k), stan::model::nil_index_list()), 
                        0.0, 
                        "assigning variable q");
            current_statement_begin__ = 46;
            for (int i = ilo; i <= iup; ++i) {
                current_statement_begin__ = 47;
                stan::model::assign(q, 
                            stan::model::cons_list(stan::model::index_uni(k), stan::model::nil_index_list()), 
                            (stan::model::rvalue(q, stan::model::cons_list(stan::model::index_uni(k), stan::model::nil_index_list()), "q") + (get_base1(p1, (i + 1), "p1", 1) * get_base1(p2, (k - i), "p2", 1))), 
                            "assigning variable q");
            }
            }
        }
        current_statement_begin__ = 50;
        return stan::math::promote_scalar<fun_return_scalar_t__>(q);
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}
struct convolve_functor__ {
    template <typename T0__, typename T1__>
        Eigen::Matrix<typename boost::math::tools::promote_args<T0__, T1__>::type, Eigen::Dynamic, 1>
    operator()(const std::vector<T0__>& p1,
             const std::vector<T1__>& p2,
             const int& K,
             const int& khalf, std::ostream* pstream__) const {
        return convolve(p1, p2, K, khalf, pstream__);
    }
};
#include <stan_meta_header.hpp>
class model_f1_g4
  : public stan::model::model_base_crtp<model_f1_g4> {
private:
        std::vector<double> p1_gl;
        std::vector<double> p2_gl;
        std::vector<int> x;
        double drbound;
public:
    model_f1_g4(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_f1_g4(stan::io::var_context& context__,
        unsigned int random_seed__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, random_seed__, pstream__);
    }
    void ctor_body(stan::io::var_context& context__,
                   unsigned int random_seed__,
                   std::ostream* pstream__) {
        typedef double local_scalar_t__;
        boost::ecuyer1988 base_rng__ =
          stan::services::util::create_rng(random_seed__, 0);
        (void) base_rng__;  // suppress unused var warning
        current_statement_begin__ = -1;
        static const char* function__ = "model_f1_g4_namespace::model_f1_g4";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            // initialize data block variables from context__
            current_statement_begin__ = 55;
            validate_non_negative_index("p1_gl", "5", 5);
            context__.validate_dims("data initialization", "p1_gl", "double", context__.to_vec(5));
            p1_gl = std::vector<double>(5, double(0));
            vals_r__ = context__.vals_r("p1_gl");
            pos__ = 0;
            size_t p1_gl_k_0_max__ = 5;
            for (size_t k_0__ = 0; k_0__ < p1_gl_k_0_max__; ++k_0__) {
                p1_gl[k_0__] = vals_r__[pos__++];
            }
            current_statement_begin__ = 56;
            validate_non_negative_index("p2_gl", "5", 5);
            context__.validate_dims("data initialization", "p2_gl", "double", context__.to_vec(5));
            p2_gl = std::vector<double>(5, double(0));
            vals_r__ = context__.vals_r("p2_gl");
            pos__ = 0;
            size_t p2_gl_k_0_max__ = 5;
            for (size_t k_0__ = 0; k_0__ < p2_gl_k_0_max__; ++k_0__) {
                p2_gl[k_0__] = vals_r__[pos__++];
            }
            current_statement_begin__ = 57;
            validate_non_negative_index("x", "5", 5);
            context__.validate_dims("data initialization", "x", "int", context__.to_vec(5));
            x = std::vector<int>(5, int(0));
            vals_i__ = context__.vals_i("x");
            pos__ = 0;
            size_t x_k_0_max__ = 5;
            for (size_t k_0__ = 0; k_0__ < x_k_0_max__; ++k_0__) {
                x[k_0__] = vals_i__[pos__++];
            }
            current_statement_begin__ = 58;
            context__.validate_dims("data initialization", "drbound", "double", context__.to_vec());
            drbound = double(0);
            vals_r__ = context__.vals_r("drbound");
            pos__ = 0;
            drbound = vals_r__[pos__++];
            // initialize transformed data variables
            // execute transformed data statements
            // validate transformed data
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 62;
            num_params_r__ += 1;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_f1_g4() { }
    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        typedef double local_scalar_t__;
        stan::io::writer<double> writer__(params_r__, params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;
        current_statement_begin__ = 62;
        if (!(context__.contains_r("alpha")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable alpha missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("alpha");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "alpha", "double", context__.to_vec());
        double alpha(0);
        alpha = vals_r__[pos__++];
        try {
            writer__.scalar_lub_unconstrain(0.0, drbound, alpha);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable alpha: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }
    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }
    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(std::vector<T__>& params_r__,
                 std::vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {
        typedef T__ local_scalar_t__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // dummy to suppress unused var warning
        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;
        try {
            stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
            // model parameters
            current_statement_begin__ = 62;
            local_scalar_t__ alpha;
            (void) alpha;  // dummy to suppress unused var warning
            if (jacobian__)
                alpha = in__.scalar_lub_constrain(0.0, drbound, lp__);
            else
                alpha = in__.scalar_lub_constrain(0.0, drbound);
            // transformed parameters
            current_statement_begin__ = 66;
            validate_non_negative_index("glmat", "5", 5);
            validate_non_negative_index("glmat", "5", 5);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, Eigen::Dynamic> glmat(5, 5);
            stan::math::initialize(glmat, DUMMY_VAR__);
            stan::math::fill(glmat, DUMMY_VAR__);
            // transformed parameters block statements
            current_statement_begin__ = 67;
            for (int i = 1; i <= 5; ++i) {
                {
                current_statement_begin__ = 68;
                validate_non_negative_index("p1", "3", 3);
                std::vector<local_scalar_t__  > p1(3, local_scalar_t__(DUMMY_VAR__));
                stan::math::initialize(p1, DUMMY_VAR__);
                stan::math::fill(p1, DUMMY_VAR__);
                stan::math::assign(p1,segfreq4(alpha, i, pstream__));
                current_statement_begin__ = 69;
                for (int j = 1; j <= 5; ++j) {
                    {
                    current_statement_begin__ = 70;
                    validate_non_negative_index("p2", "3", 3);
                    std::vector<local_scalar_t__  > p2(3, local_scalar_t__(DUMMY_VAR__));
                    stan::math::initialize(p2, DUMMY_VAR__);
                    stan::math::fill(p2, DUMMY_VAR__);
                    stan::math::assign(p2,segfreq4(alpha, i, pstream__));
                    current_statement_begin__ = 71;
                    validate_non_negative_index("q", "5", 5);
                    Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> q(5);
                    stan::math::initialize(q, DUMMY_VAR__);
                    stan::math::fill(q, DUMMY_VAR__);
                    stan::math::assign(q,convolve(p1, p2, 4, 3, pstream__));
                    current_statement_begin__ = 72;
                    stan::model::assign(glmat, 
                                stan::model::cons_list(stan::model::index_uni(i), stan::model::cons_list(stan::model::index_uni(j), stan::model::nil_index_list())), 
                                ((multinomial_log(x, q) + get_base1(p1_gl, i, "p1_gl", 1)) + get_base1(p2_gl, j, "p2_gl", 1)), 
                                "assigning variable glmat");
                    }
                }
                }
            }
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            current_statement_begin__ = 66;
            size_t glmat_j_1_max__ = 5;
            size_t glmat_j_2_max__ = 5;
            for (size_t j_1__ = 0; j_1__ < glmat_j_1_max__; ++j_1__) {
                for (size_t j_2__ = 0; j_2__ < glmat_j_2_max__; ++j_2__) {
                    if (stan::math::is_uninitialized(glmat(j_1__, j_2__))) {
                        std::stringstream msg__;
                        msg__ << "Undefined transformed parameter: glmat" << "(" << j_1__ << ", " << j_2__ << ")";
                        stan::lang::rethrow_located(std::runtime_error(std::string("Error initializing variable glmat: ") + msg__.str()), current_statement_begin__, prog_reader__());
                    }
                }
            }
            // model body
            current_statement_begin__ = 79;
            lp_accum__.add(uniform_log(alpha, 0.0, drbound));
            current_statement_begin__ = 80;
            lp_accum__.add(log_sum_exp(glmat));
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
        lp_accum__.add(lp__);
        return lp_accum__.sum();
    } // log_prob()
    template <bool propto, bool jacobian, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto,jacobian,T_>(vec_params_r, vec_params_i, pstream);
    }
    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("alpha");
        names__.push_back("glmat");
    }
    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(5);
        dims__.push_back(5);
        dimss__.push_back(dims__);
    }
    template <typename RNG>
    void write_array(RNG& base_rng__,
                     std::vector<double>& params_r__,
                     std::vector<int>& params_i__,
                     std::vector<double>& vars__,
                     bool include_tparams__ = true,
                     bool include_gqs__ = true,
                     std::ostream* pstream__ = 0) const {
        typedef double local_scalar_t__;
        vars__.resize(0);
        stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
        static const char* function__ = "model_f1_g4_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        double alpha = in__.scalar_lub_constrain(0.0, drbound);
        vars__.push_back(alpha);
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        if (!include_tparams__ && !include_gqs__) return;
        try {
            // declare and define transformed parameters
            current_statement_begin__ = 66;
            validate_non_negative_index("glmat", "5", 5);
            validate_non_negative_index("glmat", "5", 5);
            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> glmat(5, 5);
            stan::math::initialize(glmat, DUMMY_VAR__);
            stan::math::fill(glmat, DUMMY_VAR__);
            // do transformed parameters statements
            current_statement_begin__ = 67;
            for (int i = 1; i <= 5; ++i) {
                {
                current_statement_begin__ = 68;
                validate_non_negative_index("p1", "3", 3);
                std::vector<local_scalar_t__  > p1(3, local_scalar_t__(DUMMY_VAR__));
                stan::math::initialize(p1, DUMMY_VAR__);
                stan::math::fill(p1, DUMMY_VAR__);
                stan::math::assign(p1,segfreq4(alpha, i, pstream__));
                current_statement_begin__ = 69;
                for (int j = 1; j <= 5; ++j) {
                    {
                    current_statement_begin__ = 70;
                    validate_non_negative_index("p2", "3", 3);
                    std::vector<local_scalar_t__  > p2(3, local_scalar_t__(DUMMY_VAR__));
                    stan::math::initialize(p2, DUMMY_VAR__);
                    stan::math::fill(p2, DUMMY_VAR__);
                    stan::math::assign(p2,segfreq4(alpha, i, pstream__));
                    current_statement_begin__ = 71;
                    validate_non_negative_index("q", "5", 5);
                    Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> q(5);
                    stan::math::initialize(q, DUMMY_VAR__);
                    stan::math::fill(q, DUMMY_VAR__);
                    stan::math::assign(q,convolve(p1, p2, 4, 3, pstream__));
                    current_statement_begin__ = 72;
                    stan::model::assign(glmat, 
                                stan::model::cons_list(stan::model::index_uni(i), stan::model::cons_list(stan::model::index_uni(j), stan::model::nil_index_list())), 
                                ((multinomial_log(x, q) + get_base1(p1_gl, i, "p1_gl", 1)) + get_base1(p2_gl, j, "p2_gl", 1)), 
                                "assigning variable glmat");
                    }
                }
                }
            }
            if (!include_gqs__ && !include_tparams__) return;
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            // write transformed parameters
            if (include_tparams__) {
                size_t glmat_j_2_max__ = 5;
                size_t glmat_j_1_max__ = 5;
                for (size_t j_2__ = 0; j_2__ < glmat_j_2_max__; ++j_2__) {
                    for (size_t j_1__ = 0; j_1__ < glmat_j_1_max__; ++j_1__) {
                        vars__.push_back(glmat(j_1__, j_2__));
                    }
                }
            }
            if (!include_gqs__) return;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    template <typename RNG>
    void write_array(RNG& base_rng,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool include_tparams = true,
                     bool include_gqs = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng, params_r_vec, params_i_vec, vars_vec, include_tparams, include_gqs, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }
    std::string model_name() const {
        return "model_f1_g4";
    }
    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "alpha";
        param_names__.push_back(param_name_stream__.str());
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            size_t glmat_j_2_max__ = 5;
            size_t glmat_j_1_max__ = 5;
            for (size_t j_2__ = 0; j_2__ < glmat_j_2_max__; ++j_2__) {
                for (size_t j_1__ = 0; j_1__ < glmat_j_1_max__; ++j_1__) {
                    param_name_stream__.str(std::string());
                    param_name_stream__ << "glmat" << '.' << j_1__ + 1 << '.' << j_2__ + 1;
                    param_names__.push_back(param_name_stream__.str());
                }
            }
        }
        if (!include_gqs__) return;
    }
    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "alpha";
        param_names__.push_back(param_name_stream__.str());
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            size_t glmat_j_2_max__ = 5;
            size_t glmat_j_1_max__ = 5;
            for (size_t j_2__ = 0; j_2__ < glmat_j_2_max__; ++j_2__) {
                for (size_t j_1__ = 0; j_1__ < glmat_j_1_max__; ++j_1__) {
                    param_name_stream__.str(std::string());
                    param_name_stream__ << "glmat" << '.' << j_1__ + 1 << '.' << j_2__ + 1;
                    param_names__.push_back(param_name_stream__.str());
                }
            }
        }
        if (!include_gqs__) return;
    }
}; // model
}  // namespace
typedef model_f1_g4_namespace::model_f1_g4 stan_model;
#ifndef USING_R
stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}
#endif
#endif
