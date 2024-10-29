/** @file
 *****************************************************************************

 Implementation of functions to sample R1CS examples with prescribed parameters
 (according to some distribution).

 See r1cs_examples.hpp .

 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include <cassert>

#include <libff/common/utils.hpp>
#include <libsnark/gadgetlib1/gadgets/basic_gadgets.hpp>
#include <libsnark/relations/constraint_satisfaction_problems/r1cs/r1cs.hpp>
#include <libsnark/relations/constraint_satisfaction_problems/r1cs/examples/r1cs_examples.hpp>

#include <libsnark/models/transformer/components.hpp>

namespace libsnark {

// input[128, 28, 28]
// a token mixer outputs [128, 28, 28]
// a template for test
template<typename FieldT>
r1cs_example<FieldT> generate_r1cs_example_from_pool_mixer(){

    libff::enter_block("Call to generate_r1cs_example_from_pool_mixer");
    // working with pb rather than cs itself
    protoboard<FieldT> pb;
    pb_tensor<FieldT> model_output(1, 0);
    model_output.allocate(pb, model_output.size(), "model output");
    std::cout << "==>model output: " << model_output << std::endl;

    size_t input_c=1, input_h=320, input_w=49;
    pb_tensor<FieldT> model_input(input_c*input_h*input_w, 0);
    model_input.shape = {input_c, input_h, input_w};
    model_input.allocate(pb, model_input.size(), " model input ");
    std::cout << "==>model input: " << model_input << std::endl;


    libff::enter_block("NET: token mixer gadget");
    //token mixer: : [64, 56, 56] --> [64, 56, 56]
    pb_tensor<FieldT> mixed_feature(model_input.size(), 0);
    mixed_feature.shape = model_input.shape;
    mixed_feature.allocate(pb, mixed_feature.size(), "output of token mixer");
    pool_token_mixer_gadget<FieldT> mixer(pb, model_input, mixed_feature, "token mixer");
    // attention_gadget<FieldT> mixer(pb, model_input, mixed_feature, input_c/32, "attention mixer");
    // linear_attention_gadget<FieldT> mixer(pb, model_input, mixed_feature, input_c/32, "attention mixer");
    mixer.generate_r1cs_constraints();


    libff::leave_block("NET: token mixer gadget");
    std::cout << "==>token mixer output: " << mixed_feature << std::endl;

    libff::enter_block("Generate constraints");
    libff::leave_block("Generate constraints");

    libff::enter_block("Generate witness");
    for(size_t i=0; i<model_input.shape[0]; i++){
        for( size_t j=0; j<model_input.shape[1] * model_input.shape[2]; j++){
            pb.val(model_input[i*model_input.shape[1] * model_input.shape[2]+j]) = FieldT((j+1)%64);
        }
    }
    mixer.generate_r1cs_witness();


    libff::leave_block("Generate witness");



    libff::leave_block("Call to generate_r1cs_example_from_pool_mixer");
    pb.set_input_sizes(model_output.size() + model_input.size());
    r1cs_constraint_system<FieldT> cs(pb.get_constraint_system());
    return r1cs_example<FieldT>(std::move(cs), std::move(pb.primary_input()), std::move(pb.auxiliary_input()));

}



// // a template for test
// template<typename FieldT>
// r1cs_example<FieldT> generate_r1cs_example_from_fc(){

//     libff::enter_block("Call to generate_r1cs_example");
//     // working with pb rather than cs itself
//        protoboard<FieldT> pb;
    //     pb_tensor<FieldT> model_output(1, 0);
    // model_output.allocate(pb, model_output.size(), "model output");
    // std::cout << "output " << model_output << " of size " << model_output.size() << std::endl;

//     libff::enter_block("NET: patch embedding gadget");
   
//     libff::leave_block("NET: patch embedding gadget");


//     libff::enter_block("Generate constraints");



//     libff::leave_block("Generate constraints");




//     libff::enter_block("Generate witness");



//     libff::leave_block("Generate witness");



//     libff::leave_block("Call to generate_r1cs_example_from_model");
//     pb.set_input_sizes(model_output.size() + model_input.size());
//     r1cs_constraint_system<FieldT> cs(pb.get_constraint_system());
//     return r1cs_example<FieldT>(std::move(cs), std::move(pb.primary_input()), std::move(pb.auxiliary_input()));

// }



}