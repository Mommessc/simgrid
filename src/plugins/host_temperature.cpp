
/*

https://en.wikipedia.org/wiki/Function_pointer


This is the C version

void HostTemperature::update()
{
    this->my_update_function();
}

my_update_function est de type void (*my_update_function)();

à l'init :

if (host->get_propery("temperature_mode") == "A")
    this->my_update_function = *update_mode_A;
else if (host->get_property("temperature_mode") == "B")
    this-my_update_function = *update_mode_B;

(maybe this is *update_mode_A or without the *)


C++ version:
#include <functional>

use std::function<double(double)> &f (this is an argument of another function in the given example)


*/

/* Copyright (c) 2010-2020. The SimGrid Team. All rights reserved.          */

/* This program is free software; you can redistribute it and/or modify it
 * under the terms of the license (GNU LGPL) which comes with this package. */

#include "simgrid/Exception.hpp"
#include "simgrid/plugins/temperature.h"
#include "simgrid/plugins/energy.h"
#include "simgrid/s4u/Engine.hpp"
#include "src/include/surf/surf.hpp"
#include "src/plugins/vm/VirtualMachineImpl.hpp"
#include "src/surf/cpu_interface.hpp"

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>

#include <vector>
#include <string>
#include <cmath>

SIMGRID_REGISTER_PLUGIN(host_temperature, "host temperature.", &sg_host_temperature_plugin_init)

/** @defgroup plugin_host_temperature plugin_host_temperature Plugin Host Temperature

  @beginrst

This is the temperature plugin, enabling to compute the temperature of a host.
To activate this plugin, first call :cpp:func:`sg_host_temperature_plugin_init()` *before* loading your platform,
and call :cpp:func:`sg_host_get_temperature()` to retrieve the temperature of a given host.

Note that this plugin relies on the host energy plugin to provide the power and energy consumption of hosts.
As each host consumes power, some energy is dissipated as heat.
This plugin computes the temperature of a host and the energy dissipated as heat during a given period of time
from the previous temperature of the host and the temperature of the surrounding air.



The master/workers mechanism can be used to aggregate the energy consumption of multiple hosts to compute a single temperature,
when for example multiple hosts represents multiple cores of a single machine and we are only interseted in the temperature of the machine.
The function :cpp:func:`sg_host_init_master_workers()` must be called *after* loading the platform.

  @endrst
*/

XBT_LOG_NEW_DEFAULT_SUBCATEGORY(surf_temperature, surf, "Logging specific to the SURF temperature plugin");

namespace simgrid {
namespace plugin {

class HostTemperature {
public:
  static simgrid::xbt::Extension<simgrid::s4u::Host, HostTemperature> EXTENSION_ID;

  explicit HostTemperature(simgrid::s4u::Host *ptr);
  ~HostTemperature();

  double getAggregatedConsumedEnergy();
  double getHostTemperature();
  //double getLastDissipatedEnergy();
  double getSurroundingTemperature();
  void setSurroundingTemperature(double temperature);
  void update();
  void initMasterWorkers();

private:
  simgrid::s4u::Host * host = nullptr;
  std::vector<simgrid::s4u::Host *> worker_hosts; /*< The list of worker hosts to whom I'm the master. Empty if i'm not a master. >*/
  simgrid::s4u::Host * master_host = nullptr; /*< The master host if I'm a worker, nullptr if I'm a master. >*/

  double air_temperature = 20.0;/*< Temperature of the surrounding air of the host (in Celsius) >*/
  double host_temperature = air_temperature; /*< Current temperature of the host (in Celsius) >*/
  
  double host_conductive_coeff = 1.6; /*< host thermal conductivity * area / thickness (in Watt/Kelvin), with area and thickness of the surface of heat exchange b/w the host and the air.>*/
  double host_specific_heat = 0.39; /*< Specific heat of the host in joules/(gram*kelvin)>*/
  double host_mass = 150; /*< Host mass in gram >*/

public:
  double total_energy;           /*< Amount of energy used so far >*/
  //double last_dissipated_energy; /*< Amount of energy dissipated during the previous timestep */ // TODO maybe it does not make sense (regarding the update calls and stuff)
  double last_updated;           /*< Timestamp of the last update event >*/
  bool is_master = true;         /*< Whether I'm a master host in the sense of the temperature plugin >*/
};

simgrid::xbt::Extension<simgrid::s4u::Host, HostTemperature> HostTemperature::EXTENSION_ID;

/* Computes the temperature. Called lazily on need. */
void HostTemperature::update()
{
    if (dynamic_cast<simgrid::s4u::VirtualMachine*>(host)) // Ignore virtual machines
        return;

    if (!is_master) // Ignore calls on worker hosts
        return;

    double start_time  = this->last_updated;
    double finish_time = surf_get_clock();
    double elapsed_time = finish_time - start_time;

    if (elapsed_time < 1){
        /* we do not update if the time period is too small */
        return;
    }

    // Get the new energy of me and all my workers, if any
    double new_energy = sg_host_get_consumed_energy(host);
    for (simgrid::s4u::Host * worker : worker_hosts) {
        new_energy += sg_host_get_consumed_energy(worker);
    }

    double consumed_energy_at_period = new_energy - this->total_energy; // Energy consumed in this time period by me and my workers

    XBT_DEBUG("[Update temperature of %s] period=[%.2f-%.2f] Last energy was %f Joules, consumed %f Joules during this period.",
              host->get_cname(), start_time, finish_time, this->total_energy, consumed_energy_at_period);


    /*
        TODO call to the lumped model to compute the temperature
        Need also to keep track of the energy dissipated for this timestep (?)
     */
    double new_temperature = 20;

    this->host_temperature = new_temperature;
    this->last_updated = finish_time;
    this->total_energy = new_energy;
}

void HostTemperature::initMasterWorkers()
{
    if (is_master) {
        // Let's retrieve my list of worker hosts, if any!
        const char* worker_hosts_str = host->get_property("temperature_worker_hosts");
        if (worker_hosts_str != nullptr && strcmp(worker_hosts_str, "")) {
            std::vector<std::string> worker_names;
            boost::split(worker_names, worker_hosts_str, boost::is_any_of(","));
            worker_hosts.reserve(worker_names.size());
            for (const std::string host_name : worker_names)
            {
                worker_hosts.push_back(simgrid::s4u::Host::by_name(host_name)); // Throws a error if the host does not exist
            }
        }
    } else {
        // If I'm a worker, retrieve the Host of my master
        const char* master_host_str = host->get_property("temperature_master_host");
        xbt_assert(master_host_str != nullptr,"Host %s has no master host for the temperature plugin", host->get_cname());
        this->master_host = simgrid::s4u::Host::by_name(master_host_str);
    }
}

HostTemperature::HostTemperature(simgrid::s4u::Host* ptr) : host(ptr), last_updated(surf_get_clock()), worker_hosts()
{
    // First check if this host is a "master" in the sense of the temperature plugin
    const char* temperature_role_str = host->get_property("temperature_role");

    if ((temperature_role_str != nullptr) && (!strcmp(temperature_role_str, "worker")) ) {
        // I'm just a worker, do nothing
        this->is_master = false;
        return;
    }
    // Else I'm a master (by default if there is no "temperature_role" property for this host)
    this->is_master = true;
    this->total_energy = sg_host_get_consumed_energy(host);

    // Then read some properties from the XML platform file
    const char* air_temp_str = host->get_property("temperature_surroundings");
    if (air_temp_str != nullptr) {
        try {
          this->air_temperature = std::stof(std::string(air_temp_str));
        } catch (std::invalid_argument& ia) {
          throw std::invalid_argument(std::string("Invalid value for property 'temperature_surroundings' of host ") + host->get_cname() +
                                      ": " + air_temp_str);
        }
    } else {
        XBT_DEBUG("Host '%s': Missing value for property 'temperature_air', using default value (%.2f).", host->get_cname(), this->air_temperature);
    }

    const char* host_temp_str = host->get_property("temperature_host");
    if (host_temp_str != nullptr) {
        try {
            this->host_temperature = std::stof(std::string(host_temp_str));
        } catch (std::invalid_argument& ia) {
            throw std::invalid_argument(std::string("Invalid value for property 'temperature_host' of host ") + host->get_cname() +
                                        ": " + host_temp_str);
        }
    } else{
        XBT_DEBUG("Host '%s': Missing value for property 'temperature_host', using value of 'temperature_air' (%.2f).", host->get_cname(),
                 this->air_temperature);
        this->host_temperature = this->air_temperature;
    }

    /* Array of values for the conductive coefficient, specific heat and mass of host.
       Example: <prop ip="temperature_coeffs" value="1.6, 0.39, 150"/> */
    const char* coefficients_str = host->get_property("temperature_coeffs");
    if (coefficients_str != nullptr) {
        std::vector<std::string> all_coeffs;
        boost::split(all_coeffs, coefficients_str, boost::is_any_of(","));
        xbt_assert(all_coeffs.size() == 3,
                  "Incorrect size for 'temperature_coeffs' property of host '%s'."
                  "This property should contain 3 values separated by commas.", host->get_cname());
        this->host_conductive_coeff = std::stof(std::string(all_coeffs[0]));
        this->host_specific_heat = std::stof(std::string(all_coeffs[1]));
        this->host_mass = std::stof(std::string(all_coeffs[2]));
    } else {
        XBT_WARN("Host '%s': Missing value for property 'temperature_coeffs', using default values.", host->get_cname());
    }
}

HostTemperature::~HostTemperature() = default;

double HostTemperature::getAggregatedConsumedEnergy()
{
    if (this->is_master) {
        if (last_updated < surf_get_clock()) // We need to simcall this as it modifies the environment
            simgrid::kernel::actor::simcall(std::bind(&HostTemperature::update, this));

        return this->total_energy;
    } else {
        return master_host->extension<HostTemperature>()->getAggregatedConsumedEnergy();
    }
}


double HostTemperature::getHostTemperature()
{
    if (this->is_master) {
        if (last_updated < surf_get_clock()) // We need to simcall this as it modifies the environment
          simgrid::kernel::actor::simcall(std::bind(&HostTemperature::update, this));

        return host_temperature;
    } else {
        return master_host->extension<HostTemperature>()->getHostTemperature();
    }
}

double HostTemperature::getSurroundingTemperature()
{
    if (this->is_master) {
        if (last_updated < surf_get_clock()) // We need to simcall this as it modifies the environment
          simgrid::kernel::actor::simcall(std::bind(&HostTemperature::update, this));

        return air_temperature;
    } else {
        return master_host->extension<HostTemperature>()->getSurroundingTemperature();
    }
}


void HostTemperature::setSurroundingTemperature(double temperature)
{
    if (this->is_master) {
        if (last_updated < surf_get_clock()) // We need to simcall this as it modifies the environment
        simgrid::kernel::actor::simcall(std::bind(&HostTemperature::update, this));

        this->air_temperature = temperature;
    } else {
        XBT_WARN("Trying to set the surrounding temperature of Host %s which is a worker host. Doing nothing", host->get_cname());
    }
}

} // namespace plugin
} // namespace simgrid

using simgrid::plugin::HostTemperature;

/* **************************** events  callback *************************** */
static void on_creation(simgrid::s4u::Host& host)
{
    if (dynamic_cast<simgrid::s4u::VirtualMachine*>(&host)) // Ignore virtual machines
      return;

  host.extension_set(new HostTemperature(&host));
}

static void on_action_state_change(simgrid::kernel::resource::CpuAction const& action,
                                   simgrid::kernel::resource::Action::State /*previous*/)
{
    for (simgrid::kernel::resource::Cpu* const& cpu : action.cpus()) {
        simgrid::s4u::Host* host = cpu->get_host();
        if (host != nullptr) {
            // If it's a VM, take the corresponding PM
            simgrid::s4u::VirtualMachine* vm = dynamic_cast<simgrid::s4u::VirtualMachine*>(host);
            if (vm) // If it's a VM, take the corresponding PM
                host = vm->get_pm();

            // Get the host_temperature extension for the relevant host
            HostTemperature* host_temperature = host->extension<HostTemperature>();

            if (host_temperature->last_updated < surf_get_clock())
                host_temperature->update();
        }
    }
}

/* This callback is fired either when the host changes its state (on/off) ("onStateChange") or its speed
 * (because the user changed the pstate, or because of external trace events) ("onSpeedChange") */
static void on_host_change(simgrid::s4u::Host const& host)
{
    if (dynamic_cast<simgrid::s4u::VirtualMachine const*>(&host)) // Ignore virtual machines
        return;

    HostTemperature* host_temperature = host.extension<HostTemperature>();

    host_temperature->update();
}

static void on_host_destruction(simgrid::s4u::Host const& host)
{
    if (dynamic_cast<simgrid::s4u::VirtualMachine const*>(&host)) // Ignore virtual machines
        return;

    HostTemperature* host_temperature = host.extension<HostTemperature>();

    if(!host_temperature->is_master)
        return;
  
    host_temperature->update();
    double tot_energy = host_temperature->total_energy;
    if (tot_energy > 1e6){
        XBT_INFO("[Host %s] Host temperature: %.2f °C, Surrounding temperature: %.2f °C, Total energy consumed: %.0f kJ",
                host.get_cname(), host_temperature->getHostTemperature(), host_temperature->getSurroundingTemperature(), (tot_energy / 1000.0));
    }
    else{
        XBT_INFO("[Host %s] Host temperature: %.2f °C, Surroudning temperature: %.2f °C, Total energy consumed %.0f J",
                host.get_cname(), host_temperature->getHostTemperature(), host_temperature->getSurroundingTemperature(), tot_energy);
    }
}

static void on_simulation_end()
{
    //Do nothing.
    return;
}

/* **************************** Public interface *************************** */

/** \ingroup plugin_temperature
 * \brief Enable host temperature plugin
 * \details Enable temperature plugin to get the temperature of each cpu. Call this function before #MSG_init().
 */
void sg_host_temperature_plugin_init()
{
    if (HostTemperature::EXTENSION_ID.valid())
        return;

    if (not sg_host_energy_is_inited()) // TODO CHECK THAT
        sg_host_energy_plugin_init();
    // Else the energy_plugin was already inited.

    HostTemperature::EXTENSION_ID = simgrid::s4u::Host::extension_create<HostTemperature>();

    simgrid::s4u::Host::on_creation.connect(&on_creation);
    simgrid::s4u::Host::on_state_change.connect(&on_host_change);
    simgrid::s4u::Host::on_destruction.connect(&on_host_destruction);
    simgrid::kernel::resource::CpuAction::on_state_change.connect(&on_action_state_change);
}

static void ensure_plugin_inited()
{
    if (not HostTemperature::EXTENSION_ID.valid())
        throw simgrid::xbt::InitializationError("The Temperature plugin is not active. Please call sg_host_temperature_plugin_init()"
                                                "before calling any function related to that plugin.");
}


/** @ingroup plugin_temperature
 *  @brief Returns the current temperature of the (group of) hosts (in Celcius)
 *
 *  Please note that since the temperature is lazily updated, it may require a simcall to update it.
 *  The result is that the actor requesting this value will be interrupted,
 *  the value will be updated in kernel mode before returning the control to the requesting actor.
 */
double sg_host_get_host_temperature(const_sg_host_t host)
{
    ensure_plugin_inited();
    return host->extension<HostTemperature>()->getHostTemperature();
}

/** @ingroup plugin_temperature
 *  @brief Returns the surrounding temperature of the (group of) hosts (in Celcius)
 *
 *  Please note that since the temperature is lazily updated, it may require a simcall to update it.
 *  The result is that the actor requesting this value will be interrupted,
 *  the value will be updated in kernel mode before returning the control to the requesting actor.
 */
double sg_host_get_surrounding_temperature(const_sg_host_t host)
{
    ensure_plugin_inited();
    return host->extension<HostTemperature>()->getSurroundingTemperature();
}

/** @ingroup plugin_temperature
 *  @brief Returns the total consumed energy of the (group of) hosts (in Joules)
 *
 *  Please note that since the temperature is lazily updated, it may require a simcall to update it.
 *  The result is that the actor requesting this value will be interrupted,
 *  the value will be updated in kernel mode before returning the control to the requesting actor.
 */
double sg_host_get_aggregated_consumed_energy(const_sg_host_t host)
{
    ensure_plugin_inited();
    return host->extension<HostTemperature>()->getAggregatedConsumedEnergy();
}

/** @ingroup plugin_temperature
 *  @brief Sets the surrounding temperature to the temperature given in argument.
 *  Please note that since the temperature is lazily updated, it may require a simcall to update it.
 *  The result is that the actor requesting this value will be interrupted,
 *  the value will be updated in kernel mode before returning the control to the requesting actor.
 */
void sg_host_set_surrounding_temperature(const_sg_host_t host, double temperature)
{
  ensure_plugin_inited();
  host->extension<HostTemperature>()->setSurroundingTemperature(temperature);
}

/** @ingroup plugin_temperature
 *  @brief updates the temperatures of the host
 */
void sg_host_update_temperatures(const_sg_host_t host)
{
  ensure_plugin_inited();
  host->extension<HostTemperature>()->update();
}

/** @ingroup plugin_temperature
 *  @brief init the temperature master/workers information on this host
 */
void sg_host_init_master_workers(const_sg_host_t host)
{
  ensure_plugin_inited();
  host->extension<HostTemperature>()->initMasterWorkers();
}