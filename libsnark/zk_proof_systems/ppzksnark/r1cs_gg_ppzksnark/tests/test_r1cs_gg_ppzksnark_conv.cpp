/** @file
 *****************************************************************************
 Test program that exercises the ppzkSNARK (first generator, then
 prover, then verifier) on a synthetic R1CS instance.

 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#include <cassert>
#include <cstdio>

#include <libff/common/profiling.hpp>
#include <libff/common/utils.hpp>

#include <libsnark/common/default_types/r1cs_gg_ppzksnark_pp.hpp>
#include <libsnark/relations/constraint_satisfaction_problems/r1cs/examples/r1cs_examples.hpp>
#include <libsnark/zk_proof_systems/ppzksnark/r1cs_gg_ppzksnark/examples/run_r1cs_gg_ppzksnark.hpp>

#include <libsnark/gadgetlib1/examples/simple_example.hpp>

using namespace libsnark;

#ifndef NDEBUG

// test with original examples (generated from field/binary input)
template<typename ppT>
void test_r1cs_gg_ppzksnark(size_t num_constraints, size_t input_size)
// void test_r1cs_gg_ppzksnark(size_t i_h, size_t i_w, size_t k_h, size_t k_w)
{
    libff::print_header("(enter) Test R1CS GG-ppzkSNARK");

    const bool test_serialization = true;

    // r1cs_example<libff::Fr<ppT> > example = generate_r1cs_example_with_binary_input<libff::Fr<ppT> >(num_constraints, input_size);
    r1cs_example<libff::Fr<ppT> > example = generate_r1cs_example_with_field_input<libff::Fr<ppT> >(num_constraints, input_size);
    
    const bool bit = run_r1cs_gg_ppzksnark<ppT>(example, test_serialization);
    
    assert(bit);

    libff::print_header("(leave) Test R1CS GG-ppzkSNARK");
}

// test with examples (generated from matrix computation )
template<typename ppT>
void test_r1cs_gg_ppzksnark(size_t i_h, size_t i_w, size_t k_h, size_t k_w)
{
    libff::print_header("(enter) Test R1CS GG-ppzkSNARK");

    const bool test_serialization = false;

    // for mat_mul
    // r1cs_example<libff::Fr<ppT> > example = generate_r1cs_example_with_matrix<libff::Fr<ppT> >(i_h,i_w,k_h,k_w);
    
    // ours (with crpc only)
    // r1cs_example<libff::Fr<ppT> > example = generate_r1cs_example_with_matrix_opt<libff::Fr<ppT> >(i_h,i_w,k_h,k_w);
    
    // vcnn
    // r1cs_example<libff::Fr<ppT> > example = generate_r1cs_example_with_matrix_opt_2<libff::Fr<ppT> >(i_h,i_w,k_h,k_w);

    /*don't use 3*/
    // r1cs_example<libff::Fr<ppT> > example = generate_r1cs_example_with_matrix_opt_3<libff::Fr<ppT> >(i_h,i_w,k_h,k_w);

    // with opt + prefix sum
    // r1cs_example<libff::Fr<ppT> > example = generate_r1cs_example_with_matrix_opt_4<libff::Fr<ppT> >(i_h,i_w,k_h,k_w);

    // with only prefix sum
    r1cs_example<libff::Fr<ppT> > example = generate_r1cs_example_with_matrix_opt_5<libff::Fr<ppT> >(i_h,i_w,k_h,k_w);


    // for conv
    // r1cs_example<libff::Fr<ppT> > example = generate_r1cs_example_with_conv_2_opt<libff::Fr<ppT> >(i_h,i_w,k_h,k_w);
    // r1cs_example<libff::Fr<ppT> > example = generate_r1cs_example_with_conv_2<libff::Fr<ppT> >(i_h,i_w,k_h,k_w);

    // for pooling
    // r1cs_example<libff::Fr<ppT> > example = generate_r1cs_example_from_pooling_matrix<libff::Fr<ppT> >(i_h,i_w,k_h,k_w);

    const bool bit = run_r1cs_gg_ppzksnark<ppT>(example, test_serialization);
    assert(bit);

    libff::print_header("(leave) Test R1CS GG-ppzkSNARK");
}

// test with examples (generated from matrix computation )
template<typename ppT>
void test_r1cs_gg_ppzksnark(size_t n)
{
    libff::print_header("(enter) Test R1CS GG-ppzkSNARK");

    // const bool test_serialization = true;
    const bool test_serialization = false;

    // r1cs_example<libff::Fr<ppT> > example = gen_r1cs_example_from_relu<libff::Fr<ppT> >(n);
    // r1cs_example<libff::Fr<ppT> > example = gen_r1cs_example_from_protoboard<libff::Fr<ppT> >(n);
    r1cs_example<libff::Fr<ppT> > example = gen_r1cs_example_from_relu_vector<libff::Fr<ppT> >(n);

    // with only prefix sum
    // r1cs_example<libff::Fr<ppT> > example = gen_r1cs_example_from_relu_vector_opt_5<libff::Fr<ppT> >(n);

    const bool bit = run_r1cs_gg_ppzksnark<ppT>(example, test_serialization);
    assert(bit);

    libff::print_header("(leave) Test R1CS GG-ppzkSNARK");
}

int main()
{
    default_r1cs_gg_ppzksnark_pp::init_public_params();
    libff::start_profiling();
    // row == output dimension, com == input dimension, col == #tokens
    size_t row=128, com=64, col=49;
    // size_t row=2, com=3, col=2;
    test_r1cs_gg_ppzksnark<default_r1cs_gg_ppzksnark_pp>(row, com, com, col);
}
#else // NDEBUG
int main()
{
    printf("All tests here depend on assert() which is disabled by -DNDEBUG. Please recompile and run again.\n");
}
#endif // NDEBUG
