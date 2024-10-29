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

// #include <libsnark/models/r1cs_examples/generate_from_model.hpp>
#include <libsnark/models/r1cs_examples/generate_from_mixer.hpp>

using namespace libsnark;

#ifndef NDEBUG

// test with original examples (generated from field/binary input)
template<typename ppT>
void test_r1cs_gg_ppzksnark()
// void test_r1cs_gg_ppzksnark(size_t i_h, size_t i_w, size_t k_h, size_t k_w)
{
    libff::print_header("(enter) Test R1CS GG-ppzkSNARK");

    const bool test_serialization = false;
    // const bool test_serialization = true;

    // r1cs_example<libff::Fr<ppT> > example = generate_r1cs_example_from_model<libff::Fr<ppT> >();
    // r1cs_example<libff::Fr<ppT> > example = generate_r1cs_example_from_model<libff::Fr<ppT> >();
    // r1cs_example<libff::Fr<ppT> > example = generate_r1cs_example_from_fc<libff::Fr<ppT> >();
    // r1cs_example<libff::Fr<ppT> > example = generate_r1cs_example_from_component<libff::Fr<ppT> >();
    // r1cs_example<libff::Fr<ppT> > example = generate_r1cs_example_from_mlp<libff::Fr<ppT> >();
    // r1cs_example<libff::Fr<ppT> > example = generate_r1cs_example_from_cnn<libff::Fr<ppT> >();
    // r1cs_example<libff::Fr<ppT> > example = generate_r1cs_example_from_layer<libff::Fr<ppT> >();

    r1cs_example<libff::Fr<ppT> > example = generate_r1cs_example_from_pool_mixer<libff::Fr<ppT> >();


    // one example is a circuts for one computation, testing generator, prover and verifier altogether here
    const bool bit = run_r1cs_gg_ppzksnark<ppT>(example, test_serialization);
    
    assert(bit);

    libff::print_header("(leave) Test R1CS GG-ppzkSNARK");
}

int main()
{
    default_r1cs_gg_ppzksnark_pp::init_public_params();
    libff::start_profiling();
    test_r1cs_gg_ppzksnark<default_r1cs_gg_ppzksnark_pp>();
}
#else // NDEBUG
int main()
{
    printf("All tests here depend on assert() which is disabled by -DNDEBUG. Please recompile and run again.\n");
}
#endif // NDEBUG
