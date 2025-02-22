/* Copyright (c) 2006-2019. The SimGrid Team. All rights reserved.          */

/* This program is free software; you can redistribute it and/or modify it
 * under the terms of the license (GNU LGPL) which comes with this package. */

#include "simgrid/kernel/routing/NetPoint.hpp"
#include "simgrid/s4u/Actor.hpp"
#include "simgrid/s4u/Engine.hpp"
#include "simgrid/s4u/Exec.hpp"
#include "simgrid/s4u/VirtualMachine.hpp"
#include "src/plugins/vm/VirtualMachineImpl.hpp"
#include "src/simix/smx_private.hpp"
#include "src/surf/HostImpl.hpp"

#include <algorithm>
#include <string>

XBT_LOG_NEW_DEFAULT_SUBCATEGORY(s4u_host, s4u, "Logging specific to the S4U hosts");
XBT_LOG_EXTERNAL_CATEGORY(surf_route);

namespace simgrid {
namespace xbt {
template class Extendable<s4u::Host>;
} // namespace xbt

namespace s4u {

xbt::signal<void(Host&)> Host::on_creation;
xbt::signal<void(Host const&)> Host::on_destruction;
xbt::signal<void(Host const&)> Host::on_state_change;
xbt::signal<void(Host const&)> Host::on_speed_change;

Host::Host(const std::string& name) : name_(name)
{
  xbt_assert(Host::by_name_or_null(name_) == nullptr, "Refusing to create a second host named '%s'.", get_cname());
  Engine::get_instance()->host_register(name_, this);
  new surf::HostImpl(this);
}

Host::~Host()
{
  xbt_assert(currently_destroying_, "Please call h->destroy() instead of manually deleting it.");

  delete pimpl_;
  if (pimpl_netpoint != nullptr) // not removed yet by a children class
    Engine::get_instance()->netpoint_unregister(pimpl_netpoint);
  delete pimpl_cpu;
  delete mounts_;
}

/** @brief Fire the required callbacks and destroy the object
 *
 * Don't delete directly a host, call h->destroy() instead.
 *
 * This is cumbersome but this is the simplest solution to ensure that the onDestruction() callback receives a valid
 * object (because of the destructor order in a class hierarchy).
 */
void Host::destroy()
{
  if (not currently_destroying_) {
    currently_destroying_ = true;
    on_destruction(*this);
    Engine::get_instance()->host_unregister(std::string(name_));
    delete this;
  }
}

Host* Host::by_name(const std::string& name)
{
  return Engine::get_instance()->host_by_name(name);
}
Host* Host::by_name_or_null(const std::string& name)
{
  return Engine::get_instance()->host_by_name_or_null(name);
}

Host* Host::current()
{
  kernel::actor::ActorImpl* self = SIMIX_process_self();
  if (self == nullptr)
    xbt_die("Cannot call Host::current() from the maestro context");
  return self->get_host();
}

void Host::turn_on()
{
  if (not is_on()) {
    kernel::actor::simcall([this] {
      this->pimpl_cpu->turn_on();
      this->pimpl_->turn_on();
      on_state_change(*this);
    });
  }
}

/** @brief Stop the host if it is on */
void Host::turn_off()
{
  if (is_on()) {
    kernel::actor::simcall([this] {
      for (VirtualMachine* const& vm : vm::VirtualMachineImpl::allVms_)
        if (vm->get_pm() == this) {
          vm->shutdown();
          vm->turn_off();
        }
      this->pimpl_cpu->turn_off();
      this->pimpl_->turn_off();

      on_state_change(*this);
    });
  }
}

bool Host::is_on() const
{
  return this->pimpl_cpu->is_on();
}

int Host::get_pstate_count() const
{
  return this->pimpl_cpu->get_pstate_count();
}

/**
 * @brief Return a copy of the list of actors that are executing on this host.
 *
 * Daemons and regular actors are all mixed in this list.
 */
std::vector<ActorPtr> Host::get_all_actors()
{
  return pimpl_->get_all_actors();
}

/** @brief Returns how many actors (daemonized or not) have been launched on this host */
int Host::get_actor_count()
{
  return pimpl_->get_actor_count();
}

/**
 * @brief Find a route toward another host
 *
 * @param dest [IN] where to
 * @param links [OUT] where to store the list of links (must exist, cannot be nullptr).
 * @param latency [OUT] where to store the latency experienced on the path (or nullptr if not interested)
 *                It is the caller responsibility to initialize latency to 0 (we add to provided route)
 * @pre links!=nullptr
 *
 * walk through the routing components tree and find a route between hosts
 * by calling each "get_route" function in each routing component.
 */
void Host::route_to(Host* dest, std::vector<Link*>& links, double* latency)
{
  std::vector<kernel::resource::LinkImpl*> linkImpls;
  this->route_to(dest, linkImpls, latency);
  for (kernel::resource::LinkImpl* const& l : linkImpls)
    links.push_back(&l->piface_);
}

/** @brief Just like Host::routeTo, but filling an array of link implementations */
void Host::route_to(Host* dest, std::vector<kernel::resource::LinkImpl*>& links, double* latency)
{
  kernel::routing::NetZoneImpl::get_global_route(pimpl_netpoint, dest->pimpl_netpoint, links, latency);
  if (XBT_LOG_ISENABLED(surf_route, xbt_log_priority_debug)) {
    XBT_CDEBUG(surf_route, "Route from '%s' to '%s' (latency: %f):", get_cname(), dest->get_cname(),
               (latency == nullptr ? -1 : *latency));
    for (auto const& link : links)
      XBT_CDEBUG(surf_route, "Link %s", link->get_cname());
  }
}

/** @brief Returns the networking zone englobing that host */
NetZone* Host::get_englobing_zone()
{
  return pimpl_netpoint->get_englobing_zone()->get_iface();
}

void Host::send_to(Host* dest, double byte_amount)
{
  std::vector<Host*> m_host_list   = {this, dest};
  std::vector<double> flops_amount = {0, 0};
  std::vector<double> bytes_amount = {0, byte_amount, 0, 0};
  this_actor::parallel_execute(m_host_list, flops_amount, bytes_amount);
}

/** Get the properties assigned to a host */
const std::unordered_map<std::string, std::string>* Host::get_properties() const
{
  return this->pimpl_->get_properties();
}

/** Retrieve the property value (or nullptr if not set) */
const char* Host::get_property(const std::string& key) const
{
  return this->pimpl_->get_property(key);
}

void Host::set_property(const std::string& key, const std::string& value)
{
  kernel::actor::simcall([this, &key, &value] { this->pimpl_->set_property(key, value); });
}

void Host::set_properties(const std::map<std::string, std::string>& properties)
{
  kernel::actor::simcall([this, &properties] { this->pimpl_->set_properties(properties); });
}

/** Specify a profile turning the host on and off according to a exhaustive list or a stochastic law.
 * The profile must contain boolean values. */
void Host::set_state_profile(kernel::profile::Profile* p)
{
  return kernel::actor::simcall([this, p] { pimpl_cpu->set_state_profile(p); });
}
/** Specify a profile modeling the external load according to a exhaustive list or a stochastic law.
 *
 * Each event of the profile represent a peak speed change that is due to external load. The values are given as a rate
 * of the initial value. This means that the actual value is obtained by multiplying the initial value (the peek speed
 * at this pstate level) by the rate coming from the profile.
 */
void Host::set_speed_profile(kernel::profile::Profile* p)
{
  return kernel::actor::simcall([this, p] { pimpl_cpu->set_speed_profile(p); });
}

/** @brief Get the peak processor speed (in flops/s), at the specified pstate  */
double Host::get_pstate_speed(int pstate_index) const
{
  return this->pimpl_cpu->get_pstate_peak_speed(pstate_index);
}

/** @brief Get the peak computing speed in flops/s at the current pstate, NOT taking the external load into account.
 *
 *  The amount of flops per second available for computing depends on several things:
 *    - The current pstate determines the maximal peak computing speed (use @ref get_pstate_speed() to retrieve the
 *      computing speed you would get at another pstate)
 *    - If you declared an external load (with @ref simgrid::surf::Cpu::set_speed_profile()), you must multiply the
 * result of get_speed() by get_available_speed() to retrieve what a new computation would get.
 *
 *  The remaining speed is then shared between the executions located on this host.
 *  You can retrieve the amount of tasks currently running on this host with @ref get_load().
 *
 *  The host may have multiple cores, and your executions may be able to use more than a single core.
 *
 *  Finally, executions of priority 2 get twice the amount of flops than executions of priority 1.
 */
double Host::get_speed() const
{
  return this->pimpl_cpu->get_speed(1.0);
}
/** @brief Returns the current computation load (in flops per second)
 *
 * The external load (coming from an availability trace) is not taken in account.
 * You may also be interested in the load plugin.
 */
double Host::get_load() const
{
  return this->pimpl_cpu->get_load();
}
/** @brief Get the available speed ratio, between 0 and 1.
 *
 * This accounts for external load (see @ref simgrid::surf::Cpu::set_speed_profile()).
 */
double Host::get_available_speed() const
{
  return this->pimpl_cpu->get_speed_ratio();
}

/** @brief Returns the number of core of the processor. */
int Host::get_core_count() const
{
  return this->pimpl_cpu->get_core_count();
}

/** @brief Set the pstate at which the host should run */
void Host::set_pstate(int pstate_index)
{
  kernel::actor::simcall([this, pstate_index] { this->pimpl_cpu->set_pstate(pstate_index); });
}
/** @brief Retrieve the pstate at which the host is currently running */
int Host::get_pstate() const
{
  return this->pimpl_cpu->get_pstate();
}

std::vector<Disk*> Host::get_disks() const
{
  return kernel::actor::simcall([this] { return this->pimpl_->get_disks(); });
}

void Host::add_disk(Disk* disk)
{
  kernel::actor::simcall([this, disk] { this->pimpl_->add_disk(disk); });
}

void Host::remove_disk(const std::string& disk_name)
{
  kernel::actor::simcall([this, disk_name] { this->pimpl_->remove_disk(disk_name); });
}
/**
 * @ingroup simix_storage_management
 * @brief Returns the list of storages attached to a host.
 * @return a vector containing all storages attached to the host
 */
std::vector<const char*> Host::get_attached_storages() const
{
  return kernel::actor::simcall([this] { return this->pimpl_->get_attached_storages(); });
}

std::unordered_map<std::string, Storage*> const& Host::get_mounted_storages()
{
  if (mounts_ == nullptr) {
    mounts_ = new std::unordered_map<std::string, Storage*>();
    for (auto const& m : this->pimpl_->storage_) {
      mounts_->insert({m.first, &m.second->piface_});
    }
  }
  return *mounts_;
}

ExecPtr Host::exec_async(double flops)
{
  return this_actor::exec_init(flops);
}

void Host::execute(double flops)
{
  execute(flops, 1.0 /* priority */);
}

void Host::execute(double flops, double priority)
{
  this_actor::exec_init(flops)->set_priority(1 / priority)->start()->wait();
}

} // namespace s4u
} // namespace simgrid

/* **************************** Public C interface *************************** */
size_t sg_host_count()
{
  return simgrid::s4u::Engine::get_instance()->get_host_count();
}
/** @brief Returns the host list
 *
 * Uses sg_host_count() to know the array size.
 *
 * @return an array of @ref sg_host_t containing all the hosts in the platform.
 * @remark The host order in this array is generally different from the
 * creation/declaration order in the XML platform (we use a hash table
 * internally).
 * @see sg_host_count()
 */
sg_host_t* sg_host_list()
{
  xbt_assert(sg_host_count() > 0, "There is no host!");
  std::vector<simgrid::s4u::Host*> hosts = simgrid::s4u::Engine::get_instance()->get_all_hosts();

  sg_host_t* res = xbt_new(sg_host_t, hosts.size());
  memcpy(res, hosts.data(), sizeof(sg_host_t) * hosts.size());

  return res;
}

const char* sg_host_get_name(sg_host_t host)
{
  return host->get_cname();
}

void* sg_host_extension_get(sg_host_t host, size_t ext)
{
  return host->extension(ext);
}

size_t sg_host_extension_create(void (*deleter)(void*))
{
  return simgrid::s4u::Host::extension_create(deleter);
}

sg_host_t sg_host_by_name(const char* name)
{
  return simgrid::s4u::Host::by_name_or_null(name);
}

xbt_dynar_t sg_hosts_as_dynar()
{
  std::vector<simgrid::s4u::Host*> list = simgrid::s4u::Engine::get_instance()->get_all_hosts();

  auto last = std::remove_if(begin(list), end(list), [](const simgrid::s4u::Host* host) {
    return not host || not host->pimpl_netpoint || not host->pimpl_netpoint->is_host();
  });
  std::sort(begin(list), last,
            [](const simgrid::s4u::Host* a, const simgrid::s4u::Host* b) { return a->get_name() < b->get_name(); });

  xbt_dynar_t res = xbt_dynar_new(sizeof(sg_host_t), nullptr);
  std::for_each(begin(list), last, [res](sg_host_t host) { xbt_dynar_push_as(res, sg_host_t, host); });
  return res;
}

// ========= Layering madness ==============*

// ========== User data Layer ==========
void* sg_host_data(sg_host_t host)
{
  return host->get_data();
}
void sg_host_data_set(sg_host_t host, void* userdata)
{
  host->set_data(userdata);
}
void* sg_host_user(sg_host_t host) // deprecated
{
  return host->get_data();
}
void sg_host_user_set(sg_host_t host, void* userdata) // deprecated
{
  host->set_data(userdata);
}
void sg_host_user_destroy(sg_host_t host) // deprecated
{
  host->set_data(nullptr);
}

// ========= storage related functions ============
xbt_dict_t sg_host_get_mounted_storage_list(sg_host_t host)
{
  xbt_assert((host != nullptr), "Invalid parameters");
  xbt_dict_t res = xbt_dict_new_homogeneous(nullptr);
  for (auto const& elm : host->get_mounted_storages()) {
    const char* mount_name = elm.first.c_str();
    sg_storage_t storage   = elm.second;
    xbt_dict_set(res, mount_name, (void*)storage->get_cname(), nullptr);
  }

  return res;
}

xbt_dynar_t sg_host_get_attached_storage_list(sg_host_t host)
{
  xbt_dynar_t storage_dynar               = xbt_dynar_new(sizeof(const char*), nullptr);
  std::vector<const char*> storage_vector = host->get_attached_storages();
  for (auto const& name : storage_vector)
    xbt_dynar_push(storage_dynar, &name);
  return storage_dynar;
}

// =========== user-level functions ===============
// ================================================
/** @brief Returns the total speed of a host */
double sg_host_speed(sg_host_t host)
{
  return host->get_speed();
}

/** @brief Return the speed of the processor (in flop/s) at a given pstate. See also @ref plugin_energy.
 *
 * @param  host host to test
 * @param pstate_index pstate to test
 * @return Returns the processor speed associated with pstate_index
 */
double sg_host_get_pstate_speed(sg_host_t host, int pstate_index)
{
  return host->get_pstate_speed(pstate_index);
}

/** @ingroup m_host_management
 * @brief Return the number of cores.
 *
 * @param host a host
 * @return the number of cores
 */
int sg_host_core_count(sg_host_t host)
{
  return host->get_core_count();
}

double sg_host_get_available_speed(sg_host_t host)
{
  return host->get_available_speed();
}

/** @brief Returns the number of power states for a host.
 *
 *  See also @ref plugin_energy.
 */
int sg_host_get_nb_pstates(sg_host_t host)
{
  return host->get_pstate_count();
}

/** @brief Gets the pstate at which that host currently runs.
 *
 *  See also @ref plugin_energy.
 */
int sg_host_get_pstate(sg_host_t host)
{
  return host->get_pstate();
}
/** @brief Sets the pstate at which that host should run.
 *
 *  See also @ref plugin_energy.
 */
void sg_host_set_pstate(sg_host_t host, int pstate)
{
  host->set_pstate(pstate);
}

/** @ingroup m_host_management
 *
 * @brief Start the host if it is off
 *
 * See also #sg_host_is_on() to test the current state of the host and @ref plugin_energy
 * for more info on DVFS.
 */
void sg_host_turn_on(sg_host_t host)
{
  host->turn_on();
}

/** @ingroup m_host_management
 *
 * @brief Stop the host if it is on
 *
 * See also #MSG_host_is_on() to test the current state of the host and @ref plugin_energy
 * for more info on DVFS.
 */
void sg_host_turn_off(sg_host_t host)
{
  host->turn_off();
}

/** @ingroup m_host_management
 * @brief Determine if a host is up and running.
 *
 * See also #sg_host_turn_on() and #sg_host_turn_off() to switch the host ON and OFF and @ref plugin_energy for more
 * info on DVFS.
 *
 * @param host host to test
 * @return Returns true if the host is up and running, and false if it's currently down
 */
int sg_host_is_on(sg_host_t host)
{
  return host->is_on();
}

/** @deprecated */
int sg_host_is_off(sg_host_t host)
{
  return not host->is_on();
}

/** @brief Get the properties of a host */
xbt_dict_t sg_host_get_properties(sg_host_t host)
{
  xbt_dict_t as_dict = xbt_dict_new_homogeneous(xbt_free_f);
  const std::unordered_map<std::string, std::string>* props = host->get_properties();
  if (props == nullptr)
    return nullptr;
  for (auto const& elm : *props) {
    xbt_dict_set(as_dict, elm.first.c_str(), xbt_strdup(elm.second.c_str()), nullptr);
  }
  return as_dict;
}

/** @ingroup m_host_management
 * @brief Returns the value of a given host property
 *
 * @param host a host
 * @param name a property name
 * @return value of a property (or nullptr if property not set)
 */
const char* sg_host_get_property_value(sg_host_t host, const char* name)
{
  return host->get_property(name);
}

void sg_host_set_property_value(sg_host_t host, const char* name, const char* value)
{
  host->set_property(name, value);
}

/**
 * @brief Find a route between two hosts
 *
 * @param from where from
 * @param to where to
 * @param links [OUT] where to store the list of links (must exist, cannot be nullptr).
 */
void sg_host_route(sg_host_t from, sg_host_t to, xbt_dynar_t links)
{
  std::vector<simgrid::s4u::Link*> vlinks;
  from->route_to(to, vlinks, nullptr);
  for (auto const& link : vlinks)
    xbt_dynar_push(links, &link);
}
/**
 * @brief Find the latency of the route between two hosts
 *
 * @param from where from
 * @param to where to
 */
double sg_host_route_latency(sg_host_t from, sg_host_t to)
{
  std::vector<simgrid::s4u::Link*> vlinks;
  double res = 0;
  from->route_to(to, vlinks, &res);
  return res;
}
/**
 * @brief Find the bandwitdh of the route between two hosts
 *
 * @param from where from
 * @param to where to
 */
double sg_host_route_bandwidth(sg_host_t from, sg_host_t to)
{
  double min_bandwidth = -1.0;

  std::vector<simgrid::s4u::Link*> vlinks;
  from->route_to(to, vlinks, nullptr);
  for (auto const& link : vlinks) {
    double bandwidth = link->get_bandwidth();
    if (bandwidth < min_bandwidth || min_bandwidth < 0.0)
      min_bandwidth = bandwidth;
  }
  return min_bandwidth;
}

void sg_host_send_to(sg_host_t from, sg_host_t to, double byte_amount)
{
  from->send_to(to, byte_amount);
}

/** @brief Displays debugging information about a host */
void sg_host_dump(sg_host_t host)
{
  XBT_INFO("Displaying host %s", host->get_cname());
  XBT_INFO("  - speed: %.0f", host->get_speed());
  XBT_INFO("  - available speed: %.2f", sg_host_get_available_speed(host));
  const std::unordered_map<std::string, std::string>* props = host->get_properties();

  if (not props->empty()) {
    XBT_INFO("  - properties:");
    for (auto const& elm : *props) {
      XBT_INFO("    %s->%s", elm.first.c_str(), elm.second.c_str());
    }
  }
}

/** @brief Return the list of actors attached to a host.
 *
 * @param host a host
 * @param whereto a dynar in which we should push actors living on that host
 */
void sg_host_get_actor_list(sg_host_t host, xbt_dynar_t whereto)
{
  auto actors = host->get_all_actors();
  for (auto& actor : actors)
    xbt_dynar_push(whereto, &actor);
}

sg_host_t sg_host_self()
{
  smx_actor_t process = SIMIX_process_self();
  return (process == nullptr) ? nullptr : process->get_host();
}

/* needs to be public and without simcall for exceptions and logging events */
const char* sg_host_self_get_name()
{
  sg_host_t host = sg_host_self();
  if (host == nullptr || SIMIX_process_self() == simix_global->maestro_process)
    return "";

  return host->get_cname();
}

double sg_host_load(sg_host_t host)
{
  return host->get_load();
}
