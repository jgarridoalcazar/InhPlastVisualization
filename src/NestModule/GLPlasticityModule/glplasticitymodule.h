/*
 *  glplasticitymodule.h
 *
 *  This file is based on the example module distributed with NEST.
 *  
 *  Modified by: Jesœs Garrido (jgarridoalcazar at gmail.com) in 2014.
 */

#ifndef GLPLASTICITYMODULE_H
#define GLPLASTICITYMODULE_H

#include "dynmodule.h"
#include "slifunction.h"

namespace nest
{
  class Network;
}

// Put your stuff into your own namespace.
namespace mynest {

/**
 * Class defining your model.
 * @note For each model, you must define one such class, with a unique name.
 */
class GLPlasticityModule : public DynModule
{
public:

  // Interface functions ------------------------------------------

  /**
   * @note The constructor registers the module with the dynamic loader.
   *       Initialization proper is performed by the init() method.
   */
  GLPlasticityModule();

  /**
   * @note The destructor does not do much in modules. Proper "downrigging"
   *       is the responsibility of the unregister() method.
   */
  ~GLPlasticityModule();

  /**
   * Initialize module by registering models with the network.
   * @param SLIInterpreter* SLI interpreter
   * @param nest::Network*  Network with which to register models
   * @note  Parameter Network is needed for historical compatibility
   *        only.
   */
  void init(SLIInterpreter*, nest::Network*);

  /**
   * Return the name of your model.
   */
  const std::string name(void) const;

  /**
   * Return the name of a sli file to execute when mymodule is loaded.
   * This mechanism can be used to define SLI commands associated with your
   * module, in particular, set up type tries for functions you have defined.
   */
  const std::string commandstring(void) const;

public:

  // Classes implementing your functions -----------------------------

  };
} // namespace mynest

#endif
