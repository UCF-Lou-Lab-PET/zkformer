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

using namespace libsnark;

#ifndef NDEBUG
template<typename ppT>
void test_r1cs_gg_ppzksnark(size_t num_constraints,
                         size_t input_size)
{
    libff::print_header("(enter) Test R1CS GG-ppzkSNARK");

    const bool test_serialization = true;
    // r1cs_example<libff::Fr<ppT> > example = generate_r1cs_example_with_binary_input<libff::Fr<ppT> >(num_constraints, input_size);
    r1cs_example<libff::Fr<ppT> > example = generate_r1cs_example_with_field_input<libff::Fr<ppT> >(num_constraints, input_size);
    const bool bit = run_r1cs_gg_ppzksnark<ppT>(example, test_serialization);
    assert(bit);

    libff::print_header("(leave) Test R1CS GG-ppzkSNARK");
}

int main()
{
    default_r1cs_gg_ppzksnark_pp::init_public_params();
    libff::start_profiling();

    test_r1cs_gg_ppzksnark<default_r1cs_gg_ppzksnark_pp>(1000, 100);
}
#else // NDEBUG
int main()
{
    printf("All tests here depend on assert() which is disabled by -DNDEBUG. Please recompile and run again.\n");
}
#endif // NDEBUG
