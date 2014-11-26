/*
 *  mymodule.cpp
 *
 *  This file is part of NEST.
 *
 *  Copyright (C) 2004 The NEST Initiative
 *
 *  NEST is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  NEST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

// include necessary NEST headers
#include "config.h"
#include "network.h"
#include "model.h"
#include "dynamicloader.h"
#include "genericmodel.h"
#include "generic_connector.h"
#include "booldatum.h"
#include "integerdatum.h"
#include "tokenarray.h"
#include "exceptions.h"
#include "sliexceptions.h"
#include "nestmodule.h"

// include headers with your own stuff
#include "glplasticitymodule.h"
#include "iaf_cond_exp_ip.h"
#include "iaf_cond_exp_sto_ip.h"
#include "iaf_cond_exp_sym.h"
#include "stdp_sym_connection_hom.h"

// -- Interface to dynamic module loader ---------------------------------------

/*
 * The dynamic module loader must be able to find your module.
 * You make the module known to the loader by defining an instance of your
 * module class in global scope. This instance must have the name
 *
 * <modulename>_LTX_mod
 *
 * The dynamicloader can then load modulename and search for symbol "mod" in it.
 */

mynest::GLPlasticityModule glplasticitymodule_LTX_mod;

// -- DynModule functions ------------------------------------------------------

mynest::GLPlasticityModule::GLPlasticityModule()
  {
#ifdef LINKED_MODULE
     // register this module at the dynamic loader
     // this is needed to allow for linking in this module at compile time
     // all registered modules will be initialized by the main app's dynamic loader
     nest::DynamicLoaderModule::registerLinkedModule(this);
#endif
   }

mynest::GLPlasticityModule::~GLPlasticityModule()
   {}

   const std::string mynest::GLPlasticityModule::name(void) const
   {
     return std::string("Granular-layer plasticity NEST module"); // Return name of the module
   }

   const std::string mynest::GLPlasticityModule::commandstring(void) const
   {
     // Instruct the interpreter to load mymodule-init.sli
     return std::string("(glplasticitymodule-init) run");
   }

  //-------------------------------------------------------------------------------------

  void mynest::GLPlasticityModule::init(SLIInterpreter *i, nest::Network*)
  {
	/* Register a neuron or device model.
	   Give node type as template argument and the name as second argument.
	   The first argument is always a reference to the network.
	   Return value is a handle for later unregistration.
	*/
	nest::register_model<iaf_cond_exp_sto_ip>(nest::NestModule::get_network(),
	                                          "iaf_cond_exp_sto_ip");

	/* Register a neuron or device model.
       Give node type as template argument and the name as second argument.
       The first argument is always a reference to the network.
       Return value is a handle for later unregistration.
    */
    nest::register_model<iaf_cond_exp_ip>(nest::NestModule::get_network(),
                                        "iaf_cond_exp_ip");

    /* Register a neuron or device model.
	   Give node type as template argument and the name as second argument.
	   The first argument is always a reference to the network.
	   Return value is a handle for later unregistration.
	*/
	nest::register_model<iaf_cond_exp_sym>(nest::NestModule::get_network(),
											  "iaf_cond_exp_sym");


    /* Register a synapse type.
       Give synapse type as template argument and the name as second argument.
       The first argument is always a reference to the network.
    */
   nest::register_prototype_connection_commonproperties < STDPSymConnectionHom,STDPSymHomCommonProperties>
   	   	   	   	   	   	   	   	   	   	   	   (nest::NestModule::get_network(), "stdp_sym_synapse_hom");

  }  // GLPlasticityModule::init()


