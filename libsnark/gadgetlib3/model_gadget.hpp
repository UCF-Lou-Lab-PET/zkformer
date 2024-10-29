#include <libsnark/gadgetlib3/matrix_gadget.hpp>
#include <libsnark/gadgetlib3/convolution_gadget.hpp>
#include <libsnark/gadgetlib3/relu_gadget.hpp>
#include <libsnark/gadgetlib3/pool_gadget.hpp>
/**
 * generally speaking,
 * all these gadgets following the following form
 * pb_tensor is a data structure to hold both the value and size information
 * for a computation with learnable parameters (weights)
 *      gadget(pb, input, (input_size), weight, (weight_size), output, "")
 * for a computation with out weights (such as relu)
 *      gadget(pb, input, output, (additional information))
 * 
 * in all hidden layers, the output can be assigned by gen_witness
 * while the model's output is given in the primary part
*/