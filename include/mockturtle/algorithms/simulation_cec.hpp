/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2021  EPFL
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/*!
  \file simulation_cec.hpp
  \brief Simulation-based CEC

  EPFL CS-472 2021 Final Project Option 2
*/

#pragma once

#include <kitty/constructors.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/operations.hpp>

#include "../utils/node_map.hpp"
#include "miter.hpp"
#include "simulation.hpp"

namespace mockturtle
{

/* Statistics to be reported */
struct simulation_cec_stats
{
  /*! \brief Split variable (simulation size). */
  uint32_t split_var{ 0 };

  /*! \brief Number of simulation rounds. */
  uint32_t rounds{ 0 };
};

namespace detail
{

template<class Ntk>
class simulation_cec_impl
{
public:
  using pattern_t = unordered_node_map<kitty::dynamic_truth_table, Ntk>;
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;

public:
  explicit simulation_cec_impl( Ntk& ntk, simulation_cec_stats& st )
      : _ntk( ntk ),
        _st( st )
  {
  }
 
  bool run()
  {
    _st.split_var = compute_splitting_var(_ntk);
    //std::cout<< "st_split_var est" << _st.split_var << std::endl;
    _st.rounds = compute_rounds(_ntk.num_pis(),_st.split_var);
    //std::cout<< "_st.rounds est" << _st.rounds << std::endl;

    pattern_t patterns(_ntk);
    init_patterns(_ntk, _st.split_var, patterns);

    /*simulate(_ntk,patterns,_st.split_var);*/

    default_simulator<kitty::dynamic_truth_table> sim(_st.split_var);

    simulate_nodes(_ntk,patterns,sim);


    /*equivalent(_ntk,patterns);*/

    if(!equivalent(_ntk,patterns)){
      return(false);
    }

    //default_simulator<kitty::dynamic_truth_table> sim(_st.split_var);

    for (uint32_t num_r = 1; num_r <_st.rounds; num_r++){
      free(patterns);
      update_pattern(patterns,num_r);
      simulate_nodes(_ntk,patterns, sim);
      if(!equivalent(_ntk,patterns)){
      return(false);
    }
    }

    return true;
  }

private:

uint32_t compute_splitting_var ( Ntk& _n)
{ 
  uint32_t n ,v;

  n= _n.num_pis();
  v= _n._storage->nodes.size();
  if(n<=6)
  {
    return n;
  }
  else
  {
    uint32_t z = log2((1<<29)/v - 32) + 3;
    uint32_t split_var = std::min(n,z);
    return split_var;
  }
}
      
uint32_t compute_rounds(uint32_t n,uint32_t m)
{
  uint32_t rounds = 1<<(n-m);
  return rounds;
}

void init_patterns(Ntk& _net, uint32_t n, pattern_t& patterns)
{
  /*pattern_t  patterns(_net);*/
  
  _net.foreach_pi([&]( auto const& m, auto k){
  kitty::dynamic_truth_table tt (n);
  if (k<n)
  {
    kitty::create_nth_var(tt,k);
  }
  patterns[m]=tt;

});

}
   void free( pattern_t& patterns ){
    _ntk.foreach_gate( [&]( auto const& n )
    {
       patterns.erase(n);
    } );
  }
/*void simulate(Ntk& _net, pattern_t patterns,uint32_t n)
{
  default_simulator<kitty::dynamic_truth_table> sim(n);
  simulate_nodes(_net,patterns,sim);
}
*/

bool equivalent(Ntk& _net,pattern_t& patterns){
  bool check = true;
    _net.foreach_po([&](auto const& m){
    if (_net.is_complemented(m))
    {
     if(!is_const0(~patterns[m])){
       check = false;
     }
    }
    else
    {
      if(!is_const0(patterns[m])){
        check = false;
      }
    }
  });
  return check;
}

void update_pattern( pattern_t& patterns, uint32_t& num_r ){
    uint32_t nb_rounds= num_r;
      _ntk.foreach_pi( [&]( auto const& m, auto k ) 
      {
        if (k >= _st.split_var ){
          if (nb_rounds % 2 == 1) {
             if ( is_const0(patterns[m]) ) patterns[m] = ~patterns[m];
          }
        else {
           if ( !is_const0(patterns[m]) ) patterns[m] = ~patterns[m];
        }
        nb_rounds /= 2;
        }
      } );
  }


private:
  Ntk& _ntk;
  simulation_cec_stats& _st;

  /* you can add other attributes here */
};
 
} // namespace detail

/* Entry point for users to call */

/*! \brief Simulation-based CEC.
 *
 * This function implements a simulation-based combinational equivalence checker.
 * The implementation creates a miter network and run several rounds of simulation
 * to verify the functional equivalence. For memory and speed reasons this approach
 * is limited up to 40 input networks. It returns an optional which is nullopt,
 * if the network has more than 40 inputs.
 */
template<class Ntk>
std::optional<bool> simulation_cec( Ntk const& ntk1, Ntk const& ntk2, simulation_cec_stats* pst = nullptr )
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( has_num_pis_v<Ntk>, "Ntk does not implement the num_pis method" );
  static_assert( has_size_v<Ntk>, "Ntk does not implement the size method" );
  static_assert( has_get_node_v<Ntk>, "Ntk does not implement the get_node method" );
  static_assert( has_foreach_pi_v<Ntk>, "Ntk does not implement the foreach_pi method" );
  static_assert( has_foreach_po_v<Ntk>, "Ntk does not implement the foreach_po method" );
  static_assert( has_foreach_node_v<Ntk>, "Ntk does not implement the foreach_node method" );
  static_assert( has_is_complemented_v<Ntk>, "Ntk does not implement the is_complemented method" );

  simulation_cec_stats st;

  bool result = false;

  if ( ntk1.num_pis() > 40 )
    return std::nullopt;

  auto ntk_miter = miter<Ntk>( ntk1, ntk2 );

  if ( ntk_miter.has_value() )
  {
    detail::simulation_cec_impl p( *ntk_miter, st );
    result = p.run();
  }

  if ( pst )
    *pst = st;

  return result;
}

} // namespace mockturtle