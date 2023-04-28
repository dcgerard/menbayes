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
#ifndef USE_STANC3
#define USE_STANC3
#endif
#include <rstan/rstaninc.hpp>
// Code generated by %%NAME%% %%VERSION%%
#include <stan/model/model_header.hpp>
namespace model_marg_dr_pp_g4_namespace {
inline void validate_positive_index(const char* var_name, const char* expr,
                                    int val) {
  if (val < 1) {
    std::stringstream msg;
    msg << "Found dimension size less than one in simplex declaration"
        << "; variable=" << var_name << "; dimension size expression=" << expr
        << "; expression value=" << val;
    std::string msg_str(msg.str());
    throw std::invalid_argument(msg_str.c_str());
  }
}
inline void validate_unit_vector_index(const char* var_name, const char* expr,
                                       int val) {
  if (val <= 1) {
    std::stringstream msg;
    if (val == 1) {
      msg << "Found dimension size one in unit vector declaration."
          << " One-dimensional unit vector is discrete"
          << " but the target distribution must be continuous."
          << " variable=" << var_name << "; dimension size expression=" << expr;
    } else {
      msg << "Found dimension size less than one in unit vector declaration"
          << "; variable=" << var_name << "; dimension size expression=" << expr
          << "; expression value=" << val;
    }
    std::string msg_str(msg.str());
    throw std::invalid_argument(msg_str.c_str());
  }
}
using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using std::pow;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::model_base_crtp;
using stan::model::rvalue;
using stan::model::cons_list;
using stan::model::index_uni;
using stan::model::index_max;
using stan::model::index_min;
using stan::model::index_min_max;
using stan::model::index_multi;
using stan::model::index_omni;
using stan::model::nil_index_list;
using namespace stan::math;
using stan::math::pow; 
stan::math::profile_map profiles__;
static int current_statement__= 0;
static const std::vector<string> locations_array__ = {" (found before start of program)",
                                                      " (in 'marg_dr_pp_g4', line 61, column 2 to column 36)",
                                                      " (in 'marg_dr_pp_g4', line 62, column 2 to column 27)",
                                                      " (in 'marg_dr_pp_g4', line 65, column 2 to column 15)",
                                                      " (in 'marg_dr_pp_g4', line 66, column 2 to column 15)",
                                                      " (in 'marg_dr_pp_g4', line 67, column 2 to column 14)",
                                                      " (in 'marg_dr_pp_g4', line 68, column 2 to column 31)",
                                                      " (in 'marg_dr_pp_g4', line 69, column 2 to column 31)",
                                                      " (in 'marg_dr_pp_g4', line 70, column 2 to column 29)",
                                                      " (in 'marg_dr_pp_g4', line 71, column 2 to column 47)",
                                                      " (in 'marg_dr_pp_g4', line 72, column 2 to column 37)",
                                                      " (in 'marg_dr_pp_g4', line 73, column 2 to column 36)",
                                                      " (in 'marg_dr_pp_g4', line 55, column 2 to column 20)",
                                                      " (in 'marg_dr_pp_g4', line 56, column 2 to column 36)",
                                                      " (in 'marg_dr_pp_g4', line 57, column 2 to column 26)",
                                                      " (in 'marg_dr_pp_g4', line 58, column 2 to column 26)",
                                                      " (in 'marg_dr_pp_g4', line 11, column 4 to column 16)",
                                                      " (in 'marg_dr_pp_g4', line 32, column 11 to line 34, column 5)",
                                                      " (in 'marg_dr_pp_g4', line 29, column 6 to column 17)",
                                                      " (in 'marg_dr_pp_g4', line 30, column 6 to column 17)",
                                                      " (in 'marg_dr_pp_g4', line 31, column 6 to column 17)",
                                                      " (in 'marg_dr_pp_g4', line 28, column 23 to line 32, column 5)",
                                                      " (in 'marg_dr_pp_g4', line 28, column 11 to line 34, column 5)",
                                                      " (in 'marg_dr_pp_g4', line 25, column 6 to column 25)",
                                                      " (in 'marg_dr_pp_g4', line 26, column 6 to column 33)",
                                                      " (in 'marg_dr_pp_g4', line 27, column 6 to column 34)",
                                                      " (in 'marg_dr_pp_g4', line 24, column 23 to line 28, column 5)",
                                                      " (in 'marg_dr_pp_g4', line 24, column 11 to line 34, column 5)",
                                                      " (in 'marg_dr_pp_g4', line 21, column 6 to column 61)",
                                                      " (in 'marg_dr_pp_g4', line 22, column 6 to column 42)",
                                                      " (in 'marg_dr_pp_g4', line 23, column 6 to column 61)",
                                                      " (in 'marg_dr_pp_g4', line 20, column 23 to line 24, column 5)",
                                                      " (in 'marg_dr_pp_g4', line 20, column 11 to line 34, column 5)",
                                                      " (in 'marg_dr_pp_g4', line 17, column 6 to column 34)",
                                                      " (in 'marg_dr_pp_g4', line 18, column 6 to column 33)",
                                                      " (in 'marg_dr_pp_g4', line 19, column 6 to column 25)",
                                                      " (in 'marg_dr_pp_g4', line 16, column 23 to line 20, column 5)",
                                                      " (in 'marg_dr_pp_g4', line 16, column 11 to line 34, column 5)",
                                                      " (in 'marg_dr_pp_g4', line 13, column 6 to column 17)",
                                                      " (in 'marg_dr_pp_g4', line 14, column 6 to column 17)",
                                                      " (in 'marg_dr_pp_g4', line 15, column 6 to column 17)",
                                                      " (in 'marg_dr_pp_g4', line 12, column 16 to line 16, column 5)",
                                                      " (in 'marg_dr_pp_g4', line 12, column 4 to line 34, column 5)",
                                                      " (in 'marg_dr_pp_g4', line 35, column 4 to column 13)",
                                                      " (in 'marg_dr_pp_g4', line 10, column 46 to line 36, column 3)",
                                                      " (in 'marg_dr_pp_g4', line 42, column 11 to column 14)",
                                                      " (in 'marg_dr_pp_g4', line 42, column 4 to column 18)",
                                                      " (in 'marg_dr_pp_g4', line 44, column 6 to column 38)",
                                                      " (in 'marg_dr_pp_g4', line 45, column 6 to column 34)",
                                                      " (in 'marg_dr_pp_g4', line 46, column 6 to column 17)",
                                                      " (in 'marg_dr_pp_g4', line 48, column 8 to column 38)",
                                                      " (in 'marg_dr_pp_g4', line 47, column 25 to line 49, column 7)",
                                                      " (in 'marg_dr_pp_g4', line 47, column 6 to line 49, column 7)",
                                                      " (in 'marg_dr_pp_g4', line 43, column 23 to line 50, column 5)",
                                                      " (in 'marg_dr_pp_g4', line 43, column 4 to line 50, column 5)",
                                                      " (in 'marg_dr_pp_g4', line 51, column 4 to column 13)",
                                                      " (in 'marg_dr_pp_g4', line 41, column 58 to line 52, column 3)"};
template <typename T0__, typename T1__>
Eigen::Matrix<stan::promote_args_t<T0__,
T1__>, -1, 1>
segfreq4(const T0__& alpha, const T1__& xi, const int& g,
         std::ostream* pstream__) {
  using local_scalar_t__ = stan::promote_args_t<T0__, T1__>;
  const static bool propto__ = true;
  (void) propto__;
  local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
  (void) DUMMY_VAR__;  // suppress unused var warning
  
  try {
    Eigen::Matrix<local_scalar_t__, -1, 1> p;
    p = Eigen::Matrix<local_scalar_t__, -1, 1>(3);
    stan::math::fill(p, DUMMY_VAR__);
    
    current_statement__ = 42;
    if (logical_eq(g, 0)) {
      current_statement__ = 38;
      assign(p, cons_list(index_uni(1), nil_index_list()), 1.0,
        "assigning variable p");
      current_statement__ = 39;
      assign(p, cons_list(index_uni(2), nil_index_list()), 0.0,
        "assigning variable p");
      current_statement__ = 40;
      assign(p, cons_list(index_uni(3), nil_index_list()), 0.0,
        "assigning variable p");
    } else {
      current_statement__ = 37;
      if (logical_eq(g, 1)) {
        current_statement__ = 33;
        assign(p, cons_list(index_uni(1), nil_index_list()),
          (0.5 + (0.25 * alpha)), "assigning variable p");
        current_statement__ = 34;
        assign(p, cons_list(index_uni(2), nil_index_list()),
          (0.5 - (0.5 * alpha)), "assigning variable p");
        current_statement__ = 35;
        assign(p, cons_list(index_uni(3), nil_index_list()), (alpha / 4.0),
          "assigning variable p");
      } else {
        current_statement__ = 32;
        if (logical_eq(g, 2)) {
          current_statement__ = 28;
          assign(p, cons_list(index_uni(1), nil_index_list()),
            ((0.5 * alpha) + ((0.25 * (1 - alpha)) * (1 - xi))),
            "assigning variable p");
          current_statement__ = 29;
          assign(p, cons_list(index_uni(2), nil_index_list()),
            ((0.5 * (1 - alpha)) * (1 + xi)), "assigning variable p");
          current_statement__ = 30;
          assign(p, cons_list(index_uni(3), nil_index_list()),
            ((0.5 * alpha) + ((0.25 * (1 - alpha)) * (1 - xi))),
            "assigning variable p");
        } else {
          current_statement__ = 27;
          if (logical_eq(g, 3)) {
            current_statement__ = 23;
            assign(p, cons_list(index_uni(1), nil_index_list()),
              (alpha / 4.0), "assigning variable p");
            current_statement__ = 24;
            assign(p, cons_list(index_uni(2), nil_index_list()),
              (0.5 - (0.5 * alpha)), "assigning variable p");
            current_statement__ = 25;
            assign(p, cons_list(index_uni(3), nil_index_list()),
              (0.5 + (0.25 * alpha)), "assigning variable p");
          } else {
            current_statement__ = 22;
            if (logical_eq(g, 4)) {
              current_statement__ = 18;
              assign(p, cons_list(index_uni(1), nil_index_list()), 0.0,
                "assigning variable p");
              current_statement__ = 19;
              assign(p, cons_list(index_uni(2), nil_index_list()), 0.0,
                "assigning variable p");
              current_statement__ = 20;
              assign(p, cons_list(index_uni(3), nil_index_list()), 1.0,
                "assigning variable p");
            } else {
              
            }
          }
        }
      }
    }
    current_statement__ = 43;
    return p;
  } catch (const std::exception& e) {
    stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
  }
  
}
struct segfreq4_functor__ {
template <typename T0__, typename T1__>
Eigen::Matrix<stan::promote_args_t<T0__,
T1__>, -1, 1>
operator()(const T0__& alpha, const T1__& xi, const int& g,
           std::ostream* pstream__)  const 
{
return segfreq4(alpha, xi, g, pstream__);
}
};
template <typename T0__, typename T1__>
Eigen::Matrix<stan::promote_args_t<stan::value_type_t<T0__>,
stan::value_type_t<T1__>>, -1, 1>
convolve(const T0__& p1_arg__, const T1__& p2_arg__, const int& K,
         const int& khalf, std::ostream* pstream__) {
  using local_scalar_t__ = stan::promote_args_t<stan::value_type_t<T0__>,
          stan::value_type_t<T1__>>;
  const auto& p1 = to_ref(p1_arg__);
  const auto& p2 = to_ref(p2_arg__);
  const static bool propto__ = true;
  (void) propto__;
  local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
  (void) DUMMY_VAR__;  // suppress unused var warning
  
  try {
    current_statement__ = 45;
    validate_non_negative_index("q", "K + 1", (K + 1));
    Eigen::Matrix<local_scalar_t__, -1, 1> q;
    q = Eigen::Matrix<local_scalar_t__, -1, 1>((K + 1));
    stan::math::fill(q, DUMMY_VAR__);
    
    current_statement__ = 54;
    for (int k = 1; k <= (K + 1); ++k) {
      int iup;
      iup = std::numeric_limits<int>::min();
      
      current_statement__ = 47;
      iup = std::min((k - 1), (khalf - 1));
      int ilo;
      ilo = std::numeric_limits<int>::min();
      
      current_statement__ = 48;
      ilo = std::max(0, (k - khalf));
      current_statement__ = 49;
      assign(q, cons_list(index_uni(k), nil_index_list()), 0.0,
        "assigning variable q");
      current_statement__ = 52;
      for (int i = ilo; i <= iup; ++i) {
        current_statement__ = 50;
        assign(q, cons_list(index_uni(k), nil_index_list()),
          (q[(k - 1)] + (p1[((i + 1) - 1)] * p2[((k - i) - 1)])),
          "assigning variable q");}}
    current_statement__ = 55;
    return q;
  } catch (const std::exception& e) {
    stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
  }
  
}
struct convolve_functor__ {
template <typename T0__, typename T1__>
Eigen::Matrix<stan::promote_args_t<stan::value_type_t<T0__>,
stan::value_type_t<T1__>>, -1, 1>
operator()(const T0__& p1, const T1__& p2, const int& K, const int& khalf,
           std::ostream* pstream__)  const 
{
return convolve(p1, p2, K, khalf, pstream__);
}
};
#include <stan_meta_header.hpp>
class model_marg_dr_pp_g4 final : public model_base_crtp<model_marg_dr_pp_g4> {
private:
  std::vector<int> x;
  double drbound;
  int g1;
  int g2;
 
public:
  ~model_marg_dr_pp_g4() { }
  
  inline std::string model_name() const final { return "model_marg_dr_pp_g4"; }
  inline std::vector<std::string> model_compile_info() const noexcept {
    return std::vector<std::string>{"stanc_version = %%NAME%%3 %%VERSION%%", "stancflags = "};
  }
  
  
  model_marg_dr_pp_g4(stan::io::var_context& context__,
                      unsigned int random_seed__ = 0,
                      std::ostream* pstream__ = nullptr) : model_base_crtp(0) {
    using local_scalar_t__ = double ;
    boost::ecuyer1988 base_rng__ = 
        stan::services::util::create_rng(random_seed__, 0);
    (void) base_rng__;  // suppress unused var warning
    static const char* function__ = "model_marg_dr_pp_g4_namespace::model_marg_dr_pp_g4";
    (void) function__;  // suppress unused var warning
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning
    
    try {
      int pos__;
      pos__ = std::numeric_limits<int>::min();
      
      pos__ = 1;
      current_statement__ = 12;
      context__.validate_dims("data initialization","x","int",
          context__.to_vec(5));
      x = std::vector<int>(5, std::numeric_limits<int>::min());
      
      current_statement__ = 12;
      assign(x, nil_index_list(), context__.vals_i("x"),
        "assigning variable x");
      current_statement__ = 12;
      for (int sym1__ = 1; sym1__ <= 5; ++sym1__) {
        current_statement__ = 12;
        current_statement__ = 12;
        check_greater_or_equal(function__, "x[sym1__]", x[(sym1__ - 1)], 0);}
      current_statement__ = 13;
      context__.validate_dims("data initialization","drbound","double",
          context__.to_vec());
      drbound = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 13;
      drbound = context__.vals_r("drbound")[(1 - 1)];
      current_statement__ = 13;
      current_statement__ = 13;
      check_greater_or_equal(function__, "drbound", drbound, 0.0);
      current_statement__ = 13;
      current_statement__ = 13;
      check_less_or_equal(function__, "drbound", drbound, 1.0);
      current_statement__ = 14;
      context__.validate_dims("data initialization","g1","int",
          context__.to_vec());
      g1 = std::numeric_limits<int>::min();
      
      current_statement__ = 14;
      g1 = context__.vals_i("g1")[(1 - 1)];
      current_statement__ = 14;
      current_statement__ = 14;
      check_greater_or_equal(function__, "g1", g1, 0);
      current_statement__ = 14;
      current_statement__ = 14;
      check_less_or_equal(function__, "g1", g1, 4);
      current_statement__ = 15;
      context__.validate_dims("data initialization","g2","int",
          context__.to_vec());
      g2 = std::numeric_limits<int>::min();
      
      current_statement__ = 15;
      g2 = context__.vals_i("g2")[(1 - 1)];
      current_statement__ = 15;
      current_statement__ = 15;
      check_greater_or_equal(function__, "g2", g2, 0);
      current_statement__ = 15;
      current_statement__ = 15;
      check_less_or_equal(function__, "g2", g2, 4);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    num_params_r__ = 0U;
    
    try {
      num_params_r__ += 1;
      num_params_r__ += 1;
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
  }
  template <bool propto__, bool jacobian__, typename VecR, typename VecI, stan::require_vector_like_t<VecR>* = nullptr, stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr>
  inline stan::scalar_type_t<VecR> log_prob_impl(VecR& params_r__,
                                                 VecI& params_i__,
                                                 std::ostream* pstream__ = nullptr) const {
    using T__ = stan::scalar_type_t<VecR>;
    using local_scalar_t__ = T__;
    T__ lp__(0.0);
    stan::math::accumulator<T__> lp_accum__;
    static const char* function__ = "model_marg_dr_pp_g4_namespace::log_prob";
(void) function__;  // suppress unused var warning
    stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning
    
    try {
      local_scalar_t__ alpha;
      alpha = DUMMY_VAR__;
      
      current_statement__ = 1;
      alpha = in__.scalar();
      current_statement__ = 1;
      if (jacobian__) {
        current_statement__ = 1;
        alpha = stan::math::lub_constrain(alpha, 0, drbound, lp__);
      } else {
        current_statement__ = 1;
        alpha = stan::math::lub_constrain(alpha, 0, drbound);
      }
      local_scalar_t__ xi;
      xi = DUMMY_VAR__;
      
      current_statement__ = 2;
      xi = in__.scalar();
      current_statement__ = 2;
      if (jacobian__) {
        current_statement__ = 2;
        xi = stan::math::lub_constrain(xi, 0, 1, lp__);
      } else {
        current_statement__ = 2;
        xi = stan::math::lub_constrain(xi, 0, 1);
      }
      {
        Eigen::Matrix<local_scalar_t__, -1, 1> p1;
        p1 = Eigen::Matrix<local_scalar_t__, -1, 1>(3);
        stan::math::fill(p1, DUMMY_VAR__);
        
        Eigen::Matrix<local_scalar_t__, -1, 1> p2;
        p2 = Eigen::Matrix<local_scalar_t__, -1, 1>(3);
        stan::math::fill(p2, DUMMY_VAR__);
        
        Eigen::Matrix<local_scalar_t__, -1, 1> q;
        q = Eigen::Matrix<local_scalar_t__, -1, 1>(5);
        stan::math::fill(q, DUMMY_VAR__);
        
        current_statement__ = 6;
        assign(p1, nil_index_list(), segfreq4(alpha, xi, g1, pstream__),
          "assigning variable p1");
        current_statement__ = 7;
        assign(p2, nil_index_list(), segfreq4(alpha, xi, g2, pstream__),
          "assigning variable p2");
        current_statement__ = 8;
        assign(q, nil_index_list(), convolve(p1, p2, 4, 3, pstream__),
          "assigning variable q");
        current_statement__ = 9;
        lp_accum__.add(uniform_lpdf<false>(alpha, 0.0, drbound));
        current_statement__ = 10;
        lp_accum__.add(beta_lpdf<false>(xi, 1.0, 2.0));
        current_statement__ = 11;
        lp_accum__.add(multinomial_lpmf<false>(x, q));
      }
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    lp_accum__.add(lp__);
    return lp_accum__.sum();
    } // log_prob_impl() 
    
  template <typename RNG, typename VecR, typename VecI, typename VecVar, stan::require_vector_like_vt<std::is_floating_point, VecR>* = nullptr, stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr, stan::require_std_vector_vt<std::is_floating_point, VecVar>* = nullptr>
  inline void write_array_impl(RNG& base_rng__, VecR& params_r__,
                               VecI& params_i__, VecVar& vars__,
                               const bool emit_transformed_parameters__ = true,
                               const bool emit_generated_quantities__ = true,
                               std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    vars__.resize(0);
    stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
    static const char* function__ = "model_marg_dr_pp_g4_namespace::write_array";
(void) function__;  // suppress unused var warning
    (void) function__;  // suppress unused var warning
    double lp__ = 0.0;
    (void) lp__;  // dummy to suppress unused var warning
    stan::math::accumulator<double> lp_accum__;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning
    
    try {
      double alpha;
      alpha = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 1;
      alpha = in__.scalar();
      current_statement__ = 1;
      alpha = stan::math::lub_constrain(alpha, 0, drbound);
      double xi;
      xi = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 2;
      xi = in__.scalar();
      current_statement__ = 2;
      xi = stan::math::lub_constrain(xi, 0, 1);
      vars__.emplace_back(alpha);
      vars__.emplace_back(xi);
      if (logical_negation((primitive_value(emit_transformed_parameters__) ||
            primitive_value(emit_generated_quantities__)))) {
        return ;
      } 
      if (logical_negation(emit_generated_quantities__)) {
        return ;
      } 
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    } // write_array_impl() 
    
  template <typename VecVar, typename VecI, stan::require_std_vector_t<VecVar>* = nullptr, stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr>
  inline void transform_inits_impl(const stan::io::var_context& context__,
                                   VecI& params_i__, VecVar& vars__,
                                   std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    vars__.clear();
    vars__.reserve(num_params_r__);
    
    try {
      int pos__;
      pos__ = std::numeric_limits<int>::min();
      
      pos__ = 1;
      double alpha;
      alpha = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 1;
      alpha = context__.vals_r("alpha")[(1 - 1)];
      double alpha_free__;
      alpha_free__ = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 1;
      alpha_free__ = stan::math::lub_free(alpha, 0, drbound);
      double xi;
      xi = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 2;
      xi = context__.vals_r("xi")[(1 - 1)];
      double xi_free__;
      xi_free__ = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 2;
      xi_free__ = stan::math::lub_free(xi, 0, 1);
      vars__.emplace_back(alpha_free__);
      vars__.emplace_back(xi_free__);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    } // transform_inits_impl() 
    
  inline void get_param_names(std::vector<std::string>& names__) const {
    
    names__.clear();
    names__.emplace_back("alpha");
    names__.emplace_back("xi");
    } // get_param_names() 
    
  inline void get_dims(std::vector<std::vector<size_t>>& dimss__) const {
    dimss__.clear();
    dimss__.emplace_back(std::vector<size_t>{});
    
    dimss__.emplace_back(std::vector<size_t>{});
    
    } // get_dims() 
    
  inline void constrained_param_names(
                                      std::vector<std::string>& param_names__,
                                      bool emit_transformed_parameters__ = true,
                                      bool emit_generated_quantities__ = true) const
    final {
    
    param_names__.emplace_back(std::string() + "alpha");
    param_names__.emplace_back(std::string() + "xi");
    if (emit_transformed_parameters__) {
      
    }
    
    if (emit_generated_quantities__) {
      
    }
    
    } // constrained_param_names() 
    
  inline void unconstrained_param_names(
                                        std::vector<std::string>& param_names__,
                                        bool emit_transformed_parameters__ = true,
                                        bool emit_generated_quantities__ = true) const
    final {
    
    param_names__.emplace_back(std::string() + "alpha");
    param_names__.emplace_back(std::string() + "xi");
    if (emit_transformed_parameters__) {
      
    }
    
    if (emit_generated_quantities__) {
      
    }
    
    } // unconstrained_param_names() 
    
  inline std::string get_constrained_sizedtypes() const {
    stringstream s__;
    s__ << "[{\"name\":\"alpha\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"xi\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"}]";
    return s__.str();
    } // get_constrained_sizedtypes() 
    
  inline std::string get_unconstrained_sizedtypes() const {
    stringstream s__;
    s__ << "[{\"name\":\"alpha\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"xi\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"}]";
    return s__.str();
    } // get_unconstrained_sizedtypes() 
    
  
    // Begin method overload boilerplate
    template <typename RNG>
    inline void write_array(RNG& base_rng,
                            Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                            Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                            const bool emit_transformed_parameters = true,
                            const bool emit_generated_quantities = true,
                            std::ostream* pstream = nullptr) const {
      std::vector<double> vars_vec(vars.size());
      std::vector<int> params_i;
      write_array_impl(base_rng, params_r, params_i, vars_vec,
          emit_transformed_parameters, emit_generated_quantities, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i) {
        vars.coeffRef(i) = vars_vec[i];
      }
    }
    template <typename RNG>
    inline void write_array(RNG& base_rng, std::vector<double>& params_r,
                            std::vector<int>& params_i,
                            std::vector<double>& vars,
                            bool emit_transformed_parameters = true,
                            bool emit_generated_quantities = true,
                            std::ostream* pstream = nullptr) const {
      write_array_impl(base_rng, params_r, params_i, vars, emit_transformed_parameters, emit_generated_quantities, pstream);
    }
    template <bool propto__, bool jacobian__, typename T_>
    inline T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
                       std::ostream* pstream = nullptr) const {
      Eigen::Matrix<int, -1, 1> params_i;
      return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
    }
    template <bool propto__, bool jacobian__, typename T__>
    inline T__ log_prob(std::vector<T__>& params_r,
                        std::vector<int>& params_i,
                        std::ostream* pstream = nullptr) const {
      return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
    }
  
    inline void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream = nullptr) const final {
      std::vector<double> params_r_vec(params_r.size());
      std::vector<int> params_i;
      transform_inits_impl(context, params_i, params_r_vec, pstream);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i) {
        params_r.coeffRef(i) = params_r_vec[i];
      }
    }
    inline void transform_inits(const stan::io::var_context& context,
                                std::vector<int>& params_i,
                                std::vector<double>& vars,
                                std::ostream* pstream = nullptr) const final {
      transform_inits_impl(context, params_i, vars, pstream);
    }        
};
}
using stan_model = model_marg_dr_pp_g4_namespace::model_marg_dr_pp_g4;
#ifndef USING_R
// Boilerplate
stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}
stan::math::profile_map& get_stan_profile_data() {
  return model_marg_dr_pp_g4_namespace::profiles__;
}
#endif
#endif
