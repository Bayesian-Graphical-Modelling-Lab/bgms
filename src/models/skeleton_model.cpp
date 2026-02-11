// /**
//  * ==================================================
//  *  Header dependencies for SkeletonVariables
//  * ==================================================
//  *
//  * This file defines a template ("skeleton") for variable-type classes
//  * in the bgms codebase. The included headers constitute the minimal
//  * set of dependencies required by any variable model, independent of
//  * its specific statistical formulation.
//  *
//  * Each include serves a specific purpose:
//  *  - <memory>                 : ownership and cloning via std::unique_ptr
//  *  - base_model.h             : abstract interface for variable models
//  *  - adaptiveMetropolis.h     : Metropolis–Hastings proposal mechanism
//  *  - rng_utils.h              : reproducible random number generation
//  */

// #include <memory>

// #include "models/base_model.h"
// #include "models/adaptive_metropolis.h"
// #include "rng/rng_utils.h"


// /**
//  * ==================================================
//  *  SkeletonVariables class
//  * ==================================================
//  *
//  * A class in C++ represents a concrete type that bundles together:
//  *  - internal state, including observed data, model parameters,
//  *    and auxiliary objects that maintain sampling state,
//  *  - and functions (methods) that act on this state to evaluate
//  *    the log posterior, its gradient, and to perform sampling updates.
//  *
//  * In the bgms codebase, each statistical variable model is implemented
//  * as a C++ class that stores the model state and provides the methods
//  * required for inference.
//  *
//  * SkeletonVariables defines a template for such implementation classes.
//  * It specifies the structure and interface that concrete implementations
//  * of variable models must follow, without imposing a particular
//  * statistical formulation.
//  *
//  * SkeletonVariables inherits from BaseModel, which defines the common
//  * interface for all variable-model implementations in bgms:
//  *
//  *   class SkeletonVariables : public BaseModel
//  *
//  * Inheriting from BaseModel means that SkeletonVariables must provide
//  * a fixed set of functions (such as log_posterior and do_one_mh_step)
//  * that the rest of the codebase relies on. As a result, code elsewhere
//  * in bgms can interact with SkeletonVariables through the BaseModel
//  * interface, without needing to know which specific variable model
//  * implementation is being used.
//  */
// class SkeletonVariables : public BaseModel {

//   /*
//    * The 'public:' label below marks the beginning of the part of the class
//    * that is accessible from outside the class.
//    *
//    * Functions and constructors declared under 'public' are intended to be
//    * called by other components of the bgms codebase, such as samplers,
//    * model-selection routines, and result containers.
//    *
//    * Together, these public members define how a variable-model
//    * implementation can be created, queried, and updated by external code.
//    */
//   public:

//   /*
//    * Constructors are responsible for establishing the complete internal
//    * state of the object.
//    *
//    * A SkeletonVariables object represents a fully specified variable-model
//    * implementation at a given point in an inference procedure. Therefore,
//    * all information required to evaluate the log posterior and to perform
//    * sampling updates must be stored within the object itself.
//    *
//    * After construction, the object is expected to be immediately usable:
//    * no additional initialization steps are required before it can be
//    * queried, updated, or copied.
//    *
//    * The constructor below uses a constructor initializer list to define
//    * how base classes and data members are constructed. The initializer
//    * list is evaluated before the constructor body runs and is used to:
//    *  - construct the BaseModel subobject,
//    *  - initialize data members directly from constructor arguments.
//    *
//    * This ensures that all base-class and member invariants are established
//    * before any additional derived quantities are computed in the
//    * constructor body.
//    */
//   SkeletonVariables(
//     const arma::mat& observations,
//     const arma::mat& inclusion_probability,
//     const arma::imat& initial_edge_indicators,
//     const bool edge_selection = true
//   )
//     : BaseModel(edge_selection),
//       observations_(observations),
//       inclusion_probability_(inclusion_probability),
//       edge_indicators_(initial_edge_indicators)
//   {
//     /*
//      * The constructor body initializes derived state that depends on
//      * already-constructed members, such as dimensions, parameter vectors,
//      * and sampling-related objects.
//      */

//     n_ = observations_.n_rows;
//     p_ = observations_.n_cols;

//     // Dimension of the parameter vector (model-specific).
//     // For the skeleton, we assume one parameter per variable.
//     dim_ = p_;

//     // Initialize parameter vector
//     parameters_.zeros(dim_);

//     // Initialize adaptive Metropolis–Hastings sampler
//     adaptive_mh_ = AdaptiveMetropolis(dim_);
//   }

//   /*
//    * Copy constructor.
//    *
//    * The copy constructor creates a new SkeletonVariables object that is an
//    * exact copy of an existing one. This includes not only the observed data
//    * and model parameters, but also all internal state required for inference,
//    * such as sampler state and random number generator state.
//    *
//    * Copying is required in bgms because variable-model objects are duplicated
//    * during inference, for example when running multiple
//    * chains, or storing and restoring model states.
//    *
//    * The copy is performed using a constructor initializer list to ensure
//    * that the BaseModel subobject and all data members are constructed
//    * directly from their counterparts in the source object.
//    *
//    * After construction, the new object is independent from the original
//    * but represents the same model state.
//    */
//   SkeletonVariables(const SkeletonVariables& other)
//     : BaseModel(other),
//       observations_(other.observations_),
//       inclusion_probability_(other.inclusion_probability_),
//       edge_indicators_(other.edge_indicators_),
//       n_(other.n_),
//       p_(other.p_),
//       dim_(other.dim_),
//       parameters_(other.parameters_),
//       adaptive_mh_(other.adaptive_mh_),
//       rng_(other.rng_)
//   {}

//   // --------------------------------------------------
//   // Polymorphic copy
//   // --------------------------------------------------
//   std::unique_ptr<BaseModel> clone() const override {
//     return std::make_unique<SkeletonVariables>(*this);
//   }

//   // --------------------------------------------------
//   // Capabilities
//   // --------------------------------------------------
//   bool has_log_posterior() const override { return true; }
//   bool has_gradient() const override { return true; }
//   bool has_adaptive_mh() const override { return true; }

//   // --------------------------------------------------
//   // Log posterior (LIKELIHOOD + PRIOR)
//   // --------------------------------------------------
//   double log_posterior(const arma::vec& parameters) override {
//     // Skeleton: flat log-density
//     // Real models will:
//     //  - unpack parameters
//     //  - compute likelihood
//     //  - add priors
//     return 0.0;
//   }

//   // --------------------------------------------------
//   // Gradient of LOG (LIKELIHOOD * PRIOR)
//   // --------------------------------------------------
//   void gradient(const arma::vec& parameters) override {
//     // Skeleton:
//   }

//   // --------------------------------------------------
//   // One Metropolis–Hastings step
//   // --------------------------------------------------
//   void do_one_mh_step(arma::vec& parameters) override {
//     // Skeleton:
//   }

//   // --------------------------------------------------
//   // Required interface
//   // --------------------------------------------------
//   size_t parameter_dimension() const override {
//     return dim_;
//   }

//   void set_seed(int seed) override {
//     rng_ = SafeRNG(seed);
//   }

// protected:
//   // --------------------------------------------------
//   // Data
//   // --------------------------------------------------
//   arma::mat  observations_;
//   arma::mat  inclusion_probability_;
//   arma::imat edge_indicators_;

//   // --------------------------------------------------
//   // Dimensions
//   // --------------------------------------------------
//   size_t n_   = 0;   // number of observations
//   size_t p_   = 0;   // number of variables
//   size_t dim_ = 0;   // dimension of parameter vector

//   // --------------------------------------------------
//   // Parameters & MCMC machinery
//   // --------------------------------------------------
//   arma::vec parameters_;
//   AdaptiveMetropolis adaptive_mh_;
//   SafeRNG rng_;
// };
